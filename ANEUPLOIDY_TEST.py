#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANEUPLOIDY_TEST

Builds a dictionary that lists genomic windows that contain at least three
reads and gives the likelihoods to observese these reads under various
aneuploidy landscapes --- i.e., monosomy, disomy, SPH and BPH.

BPH (Both Parental Homologs) correspond to the presence of three unmatched
haplotypes, while SPH (Single Parental Homolog) correspond to chromosome gains
involving identical homologs.
Daniel Ariad (daniel@ariad.org)
Dec 22, 2020
"""

import collections, time, pickle, argparse, re, sys, random, os, bz2, gzip
from MAKE_OBS_TAB import read_impute2
from LIKELIHOODS_CALCULATOR_HOMOGENOUES import examine_homogeneous
from LIKELIHOODS_CALCULATOR_ADMIXED import examine_admixed


from itertools import product, starmap
from functools import reduce
from operator import and_
from statistics import mean, variance, pstdev
from math import log

try:
    from gmpy2 import popcount
except ModuleNotFoundError:
    print('caution: the module gmpy2 is missing.')
    def popcount(x):
        """ Counts non-zero bits in positive integer. """
        return bin(x).count('1')

try:
    from math import comb
except ImportError:
    print('caution: cound not import comb from the math module.')
    def comb(n, k):
        """ Return the number of ways to choose k items from n items without repetition and without order. """
        if not 0 <= k <= n:
            return 0
        b = 1
        for t in range(min(k, n-k)):
            b *= n
            b //= t+1
            n -= 1    
        return b

class examine:
    """ Chooses the set of the statistical models to calculate the LLR based on
        the reference panel and the number of haplotypes to subsample. """
   
    def __init__(self, obs_tab, leg_tab, hap_tab, sam_tab, models_dict, number_of_haplotypes):
        g = examine_homogeneous if all(row[2] == sam_tab[0][2] for row in sam_tab) else examine_admixed
        self.E = g(obs_tab, leg_tab, hap_tab, sam_tab, models_dict, number_of_haplotypes)
        
    def get_likelihoods(self, *x):
        F = {2: self.E.likelihoods2, 3: self.E.likelihoods3, 4: self.E.likelihoods4 }.get(len(x), self.E.likelihoods)
        return F(*x)

def mean_and_var(data):
    """ Calculates the mean and variance. """
    m = mean(data)
    var = variance(data, xbar=m)
    return m, var 

def mean_and_std(data):
    """ Calculates the mean and population standard deviation. """
    m = mean(data)
    std = pstdev(data, mu=m)
    return m, std 

def summarize(M,V):
    """ Calculates chromosome-wide statistics of the LLRs """
    result =  {'mean': mean(M),
               'std_of_mean': sum(V)**.5/len(V),  #The standard deviation is calculated according to the Bienaymé formula.
               'fraction_of_negative_LLRs': [i<0 for i in M].count(1)/len(M)}
    return result

def LLR(y,x):
    """ Calculates the logarithm of y over x and deals with edge cases. """
    if x and y:
        result = log(y/x)
    elif x and not y:
        result = -1.23456789 
    elif not x and y:
        result = +1.23456789 
    elif not x and not y:
        result = 0 
    else:
        result = None    
    return result

def invert(x,n):
    """ Inverts the bits of a positive integer. """
    return  x ^ ((1 << n) - 1)

def mismatches(obs_tab,leg_tab):
    """ Calculates the fraction of observed alleles at known SNP positions that
    mismatched the legend. """ 
    
    mismatches = 0    
    if not len(obs_tab): raise Exception('error: obs_tab is empty.')
    for (pos, ind, read_id, base) in obs_tab:
        chr_id, pos2, ref, alt = leg_tab[ind]
        if pos!=pos2:
            raise Exception('error: the line numbers in obs_tab refer to the wrong chromosome positions in leg_tab.')
        elif base!=alt and base!=ref:
            mismatches += 1
    result = mismatches/len(obs_tab)
    
    return result
        
def build_reads_dict(obs_tab,leg_tab):
    """ Returns a dictionary that lists read IDs of reads that overlap with
        SNPs and gives the alleles in each read. """

    reads = collections.defaultdict(list)
    
    for (pos, ind, read_id, base) in obs_tab:
        if base in leg_tab[ind][2:]:
            reads[read_id].append((pos,base))
            
    return reads

def build_score_dict(reads_dict,obs_tab,hap_tab,number_of_haplotypes,min_HF):
    """ Returns a dicitonary lists read_IDs and gives their score. The scoring
    algorithm scores each read according to the number of differet haplotypes
    that the reference panel supports at the chromosomal region that overlaps
    with the read. Only bialleic SNP with a minor allele frequeancy above 
    0.01 are considered for the calculation, since they are unlikely to affect
    the score. In addition, only haplotypes with a frequnecy between min_HF
    and 1-min_HF add to the score of a read. """

    N = number_of_haplotypes
    b = (1 << number_of_haplotypes) - 1 #### equivalent to int('1'*number_of_haplotypes,2)


    hap_dict = dict()
    for (pos, ind, read_id, base) in obs_tab:
        if pos not in hap_dict and (0.01 <= popcount(hap_tab[ind])/N <= 0.99): #Include only biallelic SNPs with MAF of at least 0.01. 
            hap_dict[pos] = (hap_tab[ind], hap_tab[ind] ^ b) ### ^b flips all bits of the binary number, hap_tab[ind] using bitwise xor operator. 

    score_dict = dict()
    for read_id in reads_dict:
        haplotypes = (hap_dict[pos] for pos,base in reads_dict[read_id] if pos in hap_dict)
        score_dict[read_id] = sum(min_HF <= popcount(reduce(and_,hap))/N <= (1-min_HF)
                                  for hap in product(*haplotypes) if len(hap)!=0)

    return score_dict

def iter_windows(obs_tab,leg_tab,score_dict,window_size,offset,min_reads,
                 max_reads,minimal_score):
    """ Returns an iterator over the genomic windows together with read IDs of
        the reads that overlap with SNPs in the genomic window. Only reads with
        a score larger than one are considered. """

    max_dist = 100000 #maximal distance between consecutive observed alleles.
    max_win_size = 350000 #maximal genomic window size
    initial_win_size = 50000 #initial genomic window size
    
    adaptive, window_size = (False, int(window_size)) if window_size else (True, initial_win_size)
    
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
            elif adaptive and 0<len(readIDs_in_window)<min_reads and b-pos<=max_dist and b-a<=max_win_size:
                b += 10000
            else:
                yield ((a,b-1), readIDs_in_window)
                a, b, readIDs_in_window = b, b+window_size, set() 

def pick_reads(reads_dict,score_dict,read_IDs,max_reads):
    """ Draws up to max_reads reads from a given genomic window. """

    drawn_read_IDs = random.sample(read_IDs, min(len(read_IDs)-1,max_reads))
    haplotypes = tuple(reads_dict[read_ID] for read_ID in drawn_read_IDs)
    
    return haplotypes

def effective_number_of_subsamples(num_of_reads,min_reads,max_reads,subsamples):
    """ Ensures that the number of requested subsamples is not larger than the
    number of unique subsamples. """ 
    
    if  min_reads <= num_of_reads > max_reads:
        eff_subsamples = min(comb(num_of_reads,max_reads),subsamples)
    elif min_reads <= num_of_reads <= max_reads:
        eff_subsamples = min(num_of_reads,subsamples)
    else:
        eff_subsamples = 0
        
    return eff_subsamples
        
def bootstrap(obs_tab, leg_tab, hap_tab, sam_tab, number_of_haplotypes,
              models_dict, window_size, subsamples, offset, min_reads,
              max_reads, minimal_score, min_HF):
    """ Applies a bootstrap approach in which: (i) the resample size is smaller
    than the sample size and (ii) resampling is done without replacement. """
    
    random.seed(a=None, version=2) #I should set a=None after finishing to debug the code.
    min_reads == max(min_reads,3) # Due to the bootstrap approach, min_reads must be at least 3.
    max_reads == max(max_reads,2) # Our statistical models require at least 2 reads.    
    

    reads_dict = build_reads_dict(obs_tab, leg_tab)
    score_dict = build_score_dict(reads_dict, obs_tab, hap_tab, number_of_haplotypes, min_HF)
    windows_dict = dict(iter_windows(obs_tab, leg_tab, score_dict, window_size, offset, min_reads, max_reads, minimal_score))       
    
    E = examine(obs_tab, leg_tab, hap_tab, sam_tab, models_dict, number_of_haplotypes)
    
    likelihoods = {}
    
    for k,(window,read_IDs) in enumerate(windows_dict.items()):    
        sys.stdout.write(f"\r[{'=' * (33*(k+1)//len(windows_dict)):{33}s}] {int(100*(k+1)/len(windows_dict))}%"); sys.stdout.flush()
        
        effN = effective_number_of_subsamples(len(read_IDs),min_reads,max_reads,subsamples)
        if effN>0:
            likelihoods[window] = tuple(E.get_likelihoods(*pick_reads(reads_dict,score_dict,read_IDs,max_reads)) for _ in range(effN))
    
    return likelihoods, windows_dict
        
def statistics(likelihoods,windows_dict,mismatched_alleles):
    """ Compares likelihoods of different aneuploidy scenarios and extracts
    useful information about the genmoic windows. """
    
    if likelihoods:
        window_size_mean, window_size_std = mean_and_std([j-i for (i,j) in likelihoods])    
        reads_mean, reads_std = mean_and_std([len(read_IDs) for window,read_IDs in windows_dict.items() if window in likelihoods])
        num_of_windows = len(likelihoods)
        
        
        pairs = (('BPH','SPH'), ('BPH','disomy'), ('disomy','SPH'), ('SPH','monosomy')); _ = {};
        LLRs_per_genomic_window = {(i,j): {window:  mean_and_var([*starmap(LLR, ((_[i], _[j]) for _['monosomy'], _['disomy'], _['SPH'], _['BPH'] in L))])
                           for window,L in likelihoods.items()} for i,j in pairs}
        
        LLRs_per_chromosome = {pair: summarize(*zip(*stat.values())) for pair,stat in LLRs_per_genomic_window.items()}
        
        result = {'num_of_windows': num_of_windows,
                  'reads_mean': reads_mean, 
                  'reads_std': reads_std,
                  'window_size_mean': window_size_mean,
                  'window_size_std': window_size_std,
                  'LLRs_per_genomic_window': LLRs_per_genomic_window,
                  'LLRs_per_chromosome': LLRs_per_chromosome,
                  'mismatched_alleles': mismatched_alleles}
    else:
        result = {'num_of_windows': 0,
                  'mismatched_alleles': mismatched_alleles}
    return result

def print_summary(obs_filename,info):
    S = info['statistics']
    print('\nFilename: %s' % obs_filename)
    print('Depth: %.2f, Chromosome ID: %s, Mean and standard error of meaningful reads per genomic window: %.1f, %.1f.' % (info['depth'], info['chr_id'], S.get('reads_mean',0), S.get('reads_std',0)))
    print('Number of genomic windows: %d, Mean and standard error of genomic window size: %d, %d.' % (S.get('num_of_windows',0),S.get('window_size_mean',0),S.get('window_size_std',0)))
    print('Ancestry: %s.' % info['ancestry'])
    if S.get('LLRs_per_chromosome',None):
        for (i,j), L in S['LLRs_per_chromosome'].items():
            print(f"--- LLR between {i:s} and {j:s} ----")        
            print(f"Mean LLR: {L['mean']:.3f}, Standard error of the mean LLR: {L['std_of_mean']:.3f}")
            print(f"Fraction of genomic windows with a negative LLR: {L['fraction_of_negative_LLRs']:.3f}")

def save_results(likelihoods,info,compress,obs_filename,output_filename,output_dir):
    """ Saves the likelihoods together with information about the chromosome
        number, depth of coverage, ancestry, statistics of the genomic windows 
        and flags that were used. Also, data compression is supported in gzip
        and bzip2 formats. """
        
    Open = {'bz2': bz2.open, 'gz': gzip.open}.get(compress, open)
    ext = ('.'+compress) * (compress in ('bz2','gz'))
    obs_filename_stripped =  obs_filename.rsplit('/', 1).pop()
    default_output_filename = re.sub('(.*)obs.p(.*)',f'\\1LLR.p{ext:s}', obs_filename_stripped, 1)
    output_filename = default_output_filename * (output_filename=='')
    output_dir = re.sub('/$','',output_dir)+'/' #undocumented option
    if output_dir!='' and not os.path.exists(output_dir): os.makedirs(output_dir)
    with Open(output_dir + output_filename, "wb") as f:
        pickle.dump(likelihoods, f, protocol=4)
        pickle.dump(info, f, protocol=4)       
    return output_dir + output_filename
            
def aneuploidy_test(obs_filename,leg_filename,hap_filename,sam_filename,
                    window_size,subsamples,offset,min_reads,max_reads,
                    minimal_score,min_HF,output_filename,compress,**kwargs):
    """ Returns a dictionary that lists the boundaries of approximately
    independent genomic windows. For each genomic window it gives the
    likelihood of four scenarios, namely, monosomy, disomy, SPH and BPH.
    It also returns a dictionary with various run parameters. """

    time0 = time.time()
    
    path = os.path.realpath(__file__).rsplit('/', 1)[0] + '/MODELS/'
    models_filename = kwargs.get('model', path + ('MODELS18.p' if max_reads>16 else ('MODELS16.p' if max_reads>12 else 'MODELS12.p')))

    Open = {'bz2': bz2.open, 'gzip': gzip.open}.get(obs_filename.rpartition('.')[-1], open)    
    with Open(obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        info = pickle.load(f)
        
    leg_tab = read_impute2(leg_filename, filetype='leg')
    hap_tab, number_of_haplotypes = read_impute2(hap_filename, filetype='hap')
    sam_tab = read_impute2(sam_filename, filetype='sam')
    
    ancestry = {row[2] for row in sam_tab}
    if len(ancestry)>2: print('warning: individuals in the sample file are associated with more than two populations.')
    
    load_model = bz2.BZ2File if models_filename[-4:]=='pbz2' else open
    with load_model(models_filename, 'rb') as f:
        models_dict = pickle.load(f)
    
    mismatched_alleles = mismatches(obs_tab,leg_tab)

    likelihoods, windows_dict = bootstrap(obs_tab, leg_tab, hap_tab, sam_tab, number_of_haplotypes, models_dict, window_size, subsamples, offset, min_reads, max_reads, minimal_score, min_HF)
     
    info.update({'ancestry': ancestry,
                 'window_size': window_size,
                 'subsamples': subsamples,
                 'offset': offset,
                 'min_reads': min_reads,
                 'max_reads': max_reads,
                 'minimal_score': minimal_score,
                 'min_HF': min_HF,
                 'statistics': statistics(likelihoods,windows_dict,mismatched_alleles),
                 'runtime': time.time()-time0})
    
    if output_filename!=None:
        save_results(likelihoods,info,compress,obs_filename,output_filename,kwargs.get('output_dir', 'results'))
    print_summary(obs_filename,info)
    
    time1 = time.time()
    print('Done calculating LLRs for all the genomic windows in %.3f sec.' % ((time1-time0)))
    
    return likelihoods, info

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
    parser.add_argument('sam_filename', metavar='SAM_FILENAME', type=str,
                        help='IMPUTE2 samples file')
    parser.add_argument('-w', '--window-size', type=int,
                        metavar='INT', default='100000',
                        help='Specifies the size of the genomic window. The default value is 100 kbp. When given a zero-size genomic window, it adjusts the size of the window to include min-reads reads.')
    parser.add_argument('-s', '--subsamples', type=int,
                        metavar='INT', default='32',
                        help='Sets the number of subsamples per genomic window. The default value is 32.')
    parser.add_argument('-o', '--offset', type=int,
                        metavar='INT', default=0,
                        help='Shifts all the genomic windows by the requested base pairs. The default value is 0.')
    parser.add_argument('-m', '--min-reads', type=int, metavar='INT', default=6,
                        help='Takes into account only genomic windows with at least INT reads, admitting non-zero score. The minimal value is 3, while the default is 6.')
    parser.add_argument('-M', '--max-reads', type=int, metavar='INT', default=4,
                        help='Selects up to INT reads from each genomic windows in each bootstrap sampling. The minimal value is 2, while the default value is 4.')
    parser.add_argument('-F', '--min-HF', type=int, metavar='FLOAT', default=0.05,
                        help='Only haplotypes with a frequnecy between FLOAT and 1-FLOAT add to the score of a read. The default value is 0.05.')
    parser.add_argument('-S', '--min-score', type=int, metavar='INT', default=16,
                        help='Consider only reads that reach the minimal score. The default value is 2.')
    parser.add_argument('-O', '--output-filename', type=str, metavar='output_filename',  default='',
                        help='The output filename. The default is the input filename with the extension \".obs.p\" replaced by \".LLR.p\".')
    parser.add_argument('-C', '--compress', metavar='gz/bz2/unc', type=str, default='unc',
                        help='Output compressed via gzip, bzip2 or uncompressed. Default is uncompressed.')
    args = parser.parse_args()


    LLR_dict, info = aneuploidy_test(**vars(args))

    sys.exit(0)
else:
    print("The module ANEUPLOIDY_TEST was imported.")

### END OF FILE ###
