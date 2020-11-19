#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANEUPLOIDY_TEST
Builds a dictionary that lists genomic windows that contain at least two reads
and gives the associated log-likelihood BPH/SPH ratio (LLR).
BPH (Both Parental Homologs) correspond to the presence of three unmatched
haplotypes, while SPH (Single Parental Homolog) correspond to chromosome gains
involving identical homologs.
Daniel Ariad (daniel@ariad.org)
Nov 2, 2020
"""

import collections, time, pickle, argparse, re, sys, random, inspect, os
from MAKE_OBS_TAB import read_impute2
from LLR_CALCULATOR import wrapper_func_of_create_LLR as get_LLR

from itertools import product
from functools import reduce
from operator import not_, and_
from statistics import mean, variance, pstdev

def mean_and_var(sample):
    """ Calculates the mean and the sample standard deviation. """
    m = mean(sample)
    var = variance(sample, xbar=m)
    return m, var 

def bools2int(x):
        """ Transforms a tuple/list of bools to a int. """
        return int(''.join('%d'%i for i in x),2)
    
def build_reads_dict(obs_tab,leg_tab):
    """ Returns a dictionary that lists read IDs of reads that overlap with
        SNPs and gives the alleles in each read. """

    reads = collections.defaultdict(list)

    for (pos, ind, read_id, base) in obs_tab:
        if base in leg_tab[ind][2:]:
            reads[read_id].append((pos,base))

    return reads

def build_score_dict(reads_dict,obs_tab,leg_tab,hap_tab,min_HF):
    """ Returns a dicitonary lists read_IDs and gives their score. The scoring
    algorithm scores each read according to the number of differet haplotypes
    that the reference panel supports at the chromosomal region that overlaps
    with the read. Only bialleic SNP with a minor allele frequeancy above 
    0.01 are considered for the calculation, since they are unlikely to affect
    the score. In addition, only haplotypes with a frequnecy between min_HF
    and 1-min_HF add to the score of a read. """

    N = len(hap_tab[0])

    hap_dict = dict()
    for (pos, ind, read_id, base) in obs_tab:
        if pos not in hap_dict and (0.01 <= hap_tab[ind].count(1)/N <= 0.99): #Include only biallelic SNPs with MAF of at least 0.01. 
            hap_dict[pos] = (bools2int(hap_tab[ind]), bools2int(map(not_,hap_tab[ind])))
 
    score_dict = dict()
    for read_id in reads_dict:
        haplotypes = (hap_dict[pos] for pos,base in reads_dict[read_id] if pos in hap_dict)
        score_dict[read_id] = sum(min_HF <= bin(reduce(and_,hap)).count('1')/N <= (1-min_HF)
                                  for hap in product(*haplotypes) if len(hap)!=0)

    return score_dict

def pick_reads(reads_dict,score_dict,read_IDs,min_reads,max_reads):
    """ Draws up to max_reads reads from a given genomic window. In addition,
        if the number of reads in a given genomic window is less than the
        minimal requirment then the genomic window would not be considered."""

    if len(read_IDs)<max(3,min_reads):
        haplotypes = None
    else:
        drawn_read_IDs = random.sample(read_IDs, min(len(read_IDs)-1,max_reads))
        haplotypes = tuple(reads_dict[read_ID] for read_ID in drawn_read_IDs)
    
    return haplotypes
    
def iter_windows(obs_tab,leg_tab,score_dict,window_size,offset,max_reads,minimal_score):
    """ Returns an iterator over the genomic windows together with read IDs of
        the reads that overlap with SNPs in the genomic window. Only reads with
        a score larger than one are considered. """

    adaptive, window_size = (False, int(window_size)) if window_size else (True, 50000)
    
    offset = int(offset)
    
    aux_dict = collections.defaultdict(list) ### aux_dict is a dictionary that lists chromosome positions of SNPs and gives a list of read IDs for all the reads that overlap with the SNP.  
    for (pos, ind, read_id, base) in obs_tab:
        if base in leg_tab[ind][2:]:
            aux_dict[pos].append(read_id)
            
    first, last = obs_tab[0][0]+offset, obs_tab[-1][0]+window_size
    a, b, readIDs_in_window = first, first+window_size, set()
    
    for pos in aux_dict:
        if pos<first: continue   
        while b<last:
            if a<=pos<b:
                readIDs_in_window.update(read_ID for read_ID in aux_dict[pos] if minimal_score<=score_dict[read_ID])
                break
            elif adaptive and 0<len(readIDs_in_window)<2*max_reads and b-a<350000:
                b += 10000
            else:
                yield ((a,b-1), readIDs_in_window)
                a, b, readIDs_in_window = b, b+window_size, set() 
                
def aneuploidy_test(obs_filename,leg_filename,hap_filename,window_size,subsamples,offset,min_reads,max_reads,minimal_score,min_HF,output_filename,**kwargs):
    """ Returns a dictionary that lists the boundaries of approximately
    independent genomic windows. For each genomic window it gives the
    associated log-likelihood BPH/SPH ratio (LLR)."""

    a = time.time()
    random.seed(a=None, version=2) #I should set a=None after finishing to debug the code.

    with open(obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        info = pickle.load(f)
        
    hap_tab = read_impute2(hap_filename, filetype='hap')
    leg_tab = read_impute2(leg_filename, filetype='leg')

    reads_dict = build_reads_dict(obs_tab,leg_tab)
    score_dict = build_score_dict(reads_dict,obs_tab,leg_tab,hap_tab,min_HF)
    windows_dict = {window: read_IDs for window,read_IDs in iter_windows(obs_tab,leg_tab,score_dict,window_size,offset,max_reads,minimal_score)}       
    
    
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    path = os.path.dirname(os.path.abspath(filename))
    model = kwargs.get('model', path + '/MODELS/' + ('MODELS18D.p' if max_reads>16 else 'MODELS16D.p'))
    LLR = get_LLR(obs_tab, leg_tab, hap_tab, model)
    LLR_dict = collections.defaultdict(list)

    for k in range(subsamples):
        sys.stdout.write('\r')
        sys.stdout.write(f"[{'=' * (33*(k+1)//subsamples):{33}s}] {int(100*(k+1)/subsamples)}% ")
        sys.stdout.flush()
        
        windows_dict_picked = {window: pick_reads(reads_dict,score_dict,read_IDs,min_reads,max_reads)
                           for window,read_IDs in windows_dict.items()}

        for window,haplotypes in windows_dict_picked.items():
            LLR_dict[window].append(LLR(*haplotypes) if haplotypes!=None else None)
 
    ############################### STATISTICS ################################
    reads_per_window_dict = {window:len(read_IDs) for window,read_IDs in windows_dict.items()}
    reads_per_window_filtered = [l for l in reads_per_window_dict.values() if l>1]
    reads_mean = mean(reads_per_window_filtered)
    reads_std = pstdev(reads_per_window_filtered, mu=reads_mean)
    
    LLR_stat = {window: mean_and_var(LLRs) for window,LLRs in LLR_dict.items() if None not in LLRs}
    
    M, V = zip(*LLR_stat.values())
    mean_LLR, std_of_mean_LLR = sum(M)/len(M), sum(V)**.5/len(V) #The standard deviation is calculated according to the Bienaymé formula.

    num_of_windows = len(M)
    fraction_of_negative_LLRs = sum([1 for i in M if i<0])/len(M)
    
    print('\nFilename: %s' % obs_filename)
    print('Depth: %.2f, Mean and standard error of meaningful reads per genomic window: %.1f, %.1f.' % (info['depth'], reads_mean, reads_std))
    print('Number of genomic windows: %d, Fraction of genomic windows with a negative LLR: %.3f' % (num_of_windows,fraction_of_negative_LLRs))
    print('Mean LLR: %.3f, Standard error of the mean LLR: %.3f' % ( mean_LLR,  std_of_mean_LLR))
    ###########################################################################
    
    info.update({'window_size': window_size,
                 'subsamples': subsamples,
                 'offset': offset,
                 'min_reads': min_reads,
                 'max_reads': max_reads,
                 'minimal_score': minimal_score,
                 'min_HF': min_HF,
                 'runtime': time.time()-a})

    info['statistics'] = {'mean': mean_LLR, 'std': std_of_mean_LLR,
                          'num_of_windows': num_of_windows,
                          'fraction_of_negative_LLRs': fraction_of_negative_LLRs,
                          'reads_mean': reads_mean, 'reads_std': reads_std,
                          'reads_per_window_dict': reads_per_window_dict}

    if output_filename!=None:
        output_dir = kwargs.get('output_dir', '')
        default_filename = re.sub('(.*)obs','\\1LLR', obs_filename.split('/')[-1],1)
        output_filename = default_filename if output_filename=='' else output_filename
        with open( output_dir + output_filename, "wb") as f:
            pickle.dump(LLR_dict, f, protocol=4)
            pickle.dump(info, f, protocol=4)

    b = time.time()
    print('Done calculating LLRs for all the genomic windows in %.3f sec.' % ((b-a)))
    return LLR_dict, info

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
    description='Builds a dictionary that lists genomic windows that contain'
                'at least two reads and gives the associated log-likelihood '
                'BPH/SPH ratio (LLR). BPH (Both Parental Homologs) correspond'
                'to the presence of three unmatched haplotypes, while SPH'
                '(Single Parental Homolog) correspond to chromosome gains'
                'involving identical homologs.')
    parser.add_argument('obs_filename', metavar='OBS_FILENAME', type=str,
                        help='A pickle file created by MAKE_OBS_TAB, containing base observations at known SNP positions.')
    parser.add_argument('leg_filename', metavar='LEG_FILENAME', type=str,
                        help='IMPUTE2 legend file')
    parser.add_argument('hap_filename', metavar='HAP_FILENAME', type=str,
                        help='IMPUTE2 haplotype file')
    parser.add_argument('-b', '--window-size', type=int,
                        metavar='INT', default='100000',
                        help='Specifies the size of the genomic window. The default value is 100 kbp. When given a zero-size genomic window, it adjusts the size of the window according to the local depth coverage.')
    parser.add_argument('-s', '--subsamples', type=int,
                        metavar='INT', default='32',
                        help='Sets the number of subsamples per genomic window. The default value is 32.')
    parser.add_argument('-o', '--offset', type=int,
                        metavar='INT', default=0,
                        help='Shifts all the genomic windows by the requested base pairs. The default value is 0.')
    parser.add_argument('-m', '--min-reads', type=int, metavar='INT', default=3,
                        help='Takes into account only genomic windows with at least INT reads, admitting non-zero score. The default value is 3.')
    parser.add_argument('-M', '--max-reads', type=int, metavar='INT', default=16,
                        help='Selects up to INT reads from each genomic windows. The default value is 16.')
    parser.add_argument('-l', '--min-HF', type=int, metavar='FLOAT', default=0.05,
                        help='Only haplotypes with a frequnecy between FLOAT and 1-FLOAT add to the score of a read. The default value is 0.05.')
    parser.add_argument('-c', '--min-score', type=int, metavar='INT', default=16,
                        help='Consider only reads that reach the minimal score. The default value is 2.')
    parser.add_argument('-O', '--output-filename', type=str, metavar='output_filename',  default='',
                        help='The output filename. The default is the input filename with the extension \".obs.p\" replaced by \".LLR.p\".')
    args = parser.parse_args()


    LLR_dict, info = aneuploidy_test(**vars(args))

    sys.exit(0)
else:
    print("The module ANEUPLOIDY_TEST was imported.")

### END OF FILE ###


"""
if __name__ == "__main__":
    print("Executed when invoked directly")
    a = time.time()
    random.seed(a=None, version=2) #I should set a=None after finishing to debug the code.

    args = dict(obs_filename = 'results_EUR/mixed3haploids.X0.05.HG00096A.HG00096B.HG00097A.chr21.recomb.%.2f.obs.p' % 0,
                hap_filename = '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.hap',
                leg_filename = '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.legend',
                window_size = 5e4,
                subsamples = 100,
                offset = 0,
                min_reads = 2,
                max_reads = 16,
                output_filename = None)

    with open(args['obs_filename'], 'rb') as f:
        obs_tab = pickle.load(f)
        info = pickle.load(f)
        
    hap_tab = read_impute2(args['hap_filename'], filetype='hap')
    leg_tab = read_impute2(args['leg_filename'], filetype='leg')

    reads_dict = build_reads_dict(obs_tab,leg_tab)
    score_dict = build_score_dict(reads_dict,obs_tab,leg_tab,hap_tab)
    windows_dict = {window: read_IDs for window,read_IDs in iter_windows(obs_tab,leg_tab,score_dict,args['window_size'],args['offset'],args['max_reads'])}
    
    reads_per_window_dict = {window:len(read_IDs) for window,read_IDs in windows_dict.items()}
    reads_per_window_filtered = [l for l in reads_per_window_dict.values() if l>1]
    reads_mean = mean(reads_per_window_filtered)
    reads_std = pstdev(reads_per_window_filtered, mu=reads_mean)
    
    sys.exit(0)


    LLR = get_LLR(obs_tab, leg_tab, hap_tab, 'MODELS/MODELS18D.p' if args['max_reads']>16 else 'MODELS/MODELS16D.p')
    LLR_dict = collections.defaultdict(list)

    for k in range(args['subsamples']):
        sys.stdout.write('\r')
        sys.stdout.write(f"[{'=' * int(k+1):{args['subsamples']}s}] {int(100*(k+1)/args['subsamples'])}% ")
        sys.stdout.flush()
        
        windows_dict_picked = {window: pick_reads(reads_dict,score_dict,read_IDs,args['min_reads'],args['max_reads'])
                           for window,read_IDs in windows_dict.items()}

        for window,haplotypes in windows_dict_picked.items():
            LLR_dict[window].append(LLR(*haplotypes) if haplotypes!=None else None)
 
    LLR_stat = {window: mean_and_var(LLRs) for window,LLRs in LLR_dict.items() if None not in LLRs}
    
    M, V = zip(*LLR_stat.values())
    mean, std = sum(M)/len(M), sum(V)**.5/len(V) #The standard deviation is calculated according to the Bienaymé formula.

    num_of_windows = len(M)
    fraction_of_negative_LLRs = sum([1 for i in M if i<0])/len(M)
    
    
    
    print('\nFilename: %s' % args['obs_filename'])
    print('Depth: %.2f, Number of genomic windows: %d, Fraction of genomic windows with a negative LLR: %.3f' % (info['depth'], num_of_windows,fraction_of_negative_LLRs))
    print('Mean LLR: %.3f, Standard error of the mean LLR: %.3f' % ( mean, std))
    
    info.update({'window_size': args['window_size'],
                 'offset': args['offset'],
                 'min_reads': args['min_reads'],
                 'max_reads': args['max_reads'],
                 'runtime': time.time()-a})

    info['statistics'] = {'mean': mean, 'std': std,
                          'num_of_genomic_windows': num_of_windows,
                          'fraction_of_negative_LLRs': fraction_of_negative_LLRs}

    if args['output_filename']!=None:
        default_filename = re.sub('(.*)obs','\\1LLR', args['obs_filename'].split('/')[-1],1)
        output_filename = default_filename if args['output_filename']=='' else args['output_filename']
        with open( output_filename, "wb") as f:
            pickle.dump(LLR_dict, f, protocol=4)
            pickle.dump(info, f, protocol=4)
    
    
    #filename = 'results_EUR/mixed3haploids.X0.10.HG00096A.HG00096B.HG00097A.chr21.recomb.%.2f.LLR.p' % (i * 0.1)
    #with open(filename, 'rb') as f:
    #    LLR_dict = pickle.load(f)
    #    info = pickle.load(f)
    #info['statistics']['reads_per_window_dict'] = reads_per_window_dict
    #info['statistics']['reads_mean'] = reads_mean
    #info['statistics']['reads_var'] = reads_var
    #with open( filename, "wb") as f:
    #    pickle.dump(LLR_dict, f, protocol=4)
    #    pickle.dump(info, f, protocol=4)
    #print(i)
    
    b = time.time()
    print('Done calculating LLRs for all the genomic window in %.3f sec.' % ((b-a)))
"""
