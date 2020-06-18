#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

ANEUPLOIDY_TEST

Builds a dictionary that lists linkage disequilibrium (LD) blocks that contain
at least two reads and gives the associated log-likelihood BPH/SPH ratio (LLR).
BPH (Both Parental Homologs) correspond to the presence of three unmatched
haplotypes, while SPH (Single Parental Homolog) correspond to chromosome gains
involving identical homologs.

Daniel Ariad (daniel@ariad.org)
May 11st, 2020
"""

import collections, time, pickle, statistics, argparse, re, sys, operator, random

from MAKE_OBS_TAB import read_impute2
from LLR_CALCULATOR import wrapper_func_of_create_LLR as get_LLR

def mean_and_std(sample):
    """ Calculates the mean and standard deviation of normally distributed
        random variables. """
    mean = statistics.mean(sample)
    std = statistics.pstdev(sample, mu=mean)/len(sample)**.5
    return mean, std

def jackknifing(sample,weights):
    """ Given sample elements and the weight of each element, the jackknife
    estimator, the jackknife standard error and the Jackknife bias (RB)
    are calculated. More information about delete-m jackknife for unequal m
    can be found in F.M.Busing et al. (1999), [DOI:10.1023/A:1008800423698]. """
    
    N = len(sample)
    t0 = sum(sample) / N
    H = [1/w for w in weights]
    ###H = [N for w in weights]
    T = [sum(sample[:i]+sample[i+1:])/(N-1) for i in range(N)] 
    pseudo_values = [h*t0-(h-1)*t for t,h in zip(T,H)] 
    jackknife_estimator = sum((p/h for p,h in zip(pseudo_values,H)))
    jackknife_variance = sum(((p-jackknife_estimator)**2/(h-1) for p,h in zip(pseudo_values,H)))/N
    jackknife_bias = (N-1)*t0-sum((t*(h-1)/h for t,h in zip(T,H)))
    return jackknife_estimator, jackknife_variance**.5 , jackknife_bias

def build_aux_dict(obs_tab,leg_tab):
    """ Returns a dictionary that lists chromosome positions of SNPs and gives
        a dictionary that lists the observed bases at these postions. The
        nested dictionary gives all the read IDs that contain the observed base
        at the specific chromosome position. """
    
    nested = lambda: collections.defaultdict(list)
    aux_dict = collections.defaultdict(nested)

    for (pos, ind, read_id, base) in obs_tab:
        if base in leg_tab[ind][2:]:
            aux_dict[pos][base].append(read_id)  
       
    return aux_dict 

def build_reads_dict(obs_tab,leg_tab):
    """ Returns a dictionary that lists read IDs of reads that overlap with
        SNPs and gives the alleles in each read. """
    
    reads = collections.defaultdict(list)

    for (pos, ind, read_id, base) in obs_tab:
        if base in leg_tab[ind][2:]:
            reads[read_id].append((pos,base))  
       
    return reads

def pick_reads(reads_dict,read_IDs,min_reads,max_reads):
    """ Picks randomly up to N reads out of all the reads in a given LD block. 
        In addition, if the number of reads in a given LD block is less than
        the minimal requirment then the block would not be considered."""
    
    if len(read_IDs) < max(2,min_reads): return None
    reads = tuple(reads_dict[read_ID] for read_ID in read_IDs)
    haplotypes = reads if len(reads)<=max_reads else random.sample(reads,max_reads)    
    
    return haplotypes
    
def build_blocks_dict(aux_dict,block_size,offset):
    """ Returns a dictionary that lists LD blocks and gives the read IDs of 
        reads that overlap with SNPs in the block."""
    
    first, *_, last = iter(aux_dict)
    boundaries = tuple(range(int(first+offset), int(last), int(block_size)))
    blocks = [(i,j-1) for i,j in zip(boundaries,boundaries[1:])]
    
    blocks_dict = collections.defaultdict(set)
    
    blocks_iterator = iter(blocks)
    block = next(blocks_iterator, None)
    for pos in aux_dict:
        while block:
            if block[0]<=pos<=block[1]:
                for read_IDs in aux_dict[pos].values(): blocks_dict[block].update(read_IDs)
                break
            block = next(blocks_iterator, None)
    return blocks_dict   

def aneuploidy_test(obs_filename,leg_filename,hap_filename,block_size,offset,min_reads,max_reads,output_filename):
    """ Returns a dictionary that lists the boundaries of approximately
    independent blocks of linkage disequilibrium (LD). For each LD block it
    gives the associated log-likelihood BPH/SPH ratio (LLR)."""
    
    a = time.time()

    hap_tab = read_impute2(hap_filename, filetype='hap')
    leg_tab = read_impute2(leg_filename, filetype='leg')
    
    with open(obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        info = pickle.load(f)
        
    aux_dict = build_aux_dict(obs_tab,leg_tab)
    blocks_dict = build_blocks_dict(aux_dict,block_size,offset)
    
    random.seed(a=0, version=2) #I should set a=None after finishing to debug the code.
    reads_dict = build_reads_dict(obs_tab,leg_tab)
    blocks_dict_picked = {block: pick_reads(reads_dict,read_IDs,min_reads,max_reads) 
                               for block,read_IDs in blocks_dict.items()}
    
    LLR = get_LLR(obs_tab, leg_tab, hap_tab, 'MODELS/MODELS16A.pbz2')
    LLR_dict = {block: LLR(*haplotypes) if haplotypes!=None else None for block,haplotypes in blocks_dict_picked.items()}
    
    population = tuple(value for value in LLR_dict.values() if value!=None)    
    mean, std = mean_and_std(population)
    
    w0 = [len(blocks_dict[key]) for key,value in LLR_dict.items() if value!=None]
    W = sum(w0)
    weights = [i/W for i in w0]
    jk_mean, jk_std, jk_bias = jackknifing(population,weights)
    
    num_of_LD_blocks = len(population)
    fraction_of_negative_LLRs = sum([1 for i  in population if i<0])/len(population)
    
    print('Depth: %.2f, Mean: %.3f, STD: %.3f' % (info['depth'], mean, std))
    print('Jackknife estimator: %.3f, Jackknife standard error: %.3f, Jackknife bias: %.3f' % (jk_mean, jk_std, jk_bias))
    print('Number of LD blocks: %d, Fraction of LD blocks with a negative LLR: %.3f' % (num_of_LD_blocks,fraction_of_negative_LLRs))
    
    info.update({'block_size': block_size,
                 'offset': offset,
                 'min_reads': min_reads,
                 'max_reads': max_reads,
                 'runtime': time.time()-a})
    
    info['statistics'] = {'mean': mean, 'std': std,
                          'jk_mean': jk_mean, 'jk_std': jk_std, 'jk_bias': jk_bias,
                          'num_of_LD_blocks': num_of_LD_blocks,
                          'fraction_of_negative_LLRs': fraction_of_negative_LLRs} 
    
    if output_filename!=None:
        default_filename = re.sub('(.*)obs','\\1LLR', obs_filename.split('/')[-1],1)    
        output_filename = default_filename if output_filename=='' else output_filename 
        with open( output_filename, "wb") as f:
            pickle.dump(obs_tab, f, protocol=4)
            pickle.dump(info, f, protocol=4)    
    
    b = time.time()
    print('Done calculating LLRs for all the LD block in %.3f sec.' % ((b-a)))
    return LLR_dict, info

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
    description='Builds a dictionary that lists linkage disequilibrium (LD) '
                'blocks that contain at least two reads and gives the '
                'associated log-likelihood BPH/SPH ratio (LLR). BPH (Both '
                'Parental Homologs) correspond to the presence of three '
                'unmatched haplotypes, while SPH (Single Parental Homolog) '
                'correspond to chromosome gains involving identical homologs.')
    parser.add_argument('obs_filename', metavar='OBS_FILENAME', type=str, 
                        help='A pickle file created by MAKE_OBS_TAB, containing base observations at known SNP positions.')
    parser.add_argument('leg_filename', metavar='LEG_FILENAME', type=str, 
                        help='IMPUTE2 legend file')
    parser.add_argument('hap_filename', metavar='HAP_FILENAME', type=str, 
                        help='IMPUTE2 haplotype file')
    parser.add_argument('-s', '--block-size', type=int, 
                        metavar='INT', default='100000', 
                        help='Specifies the typical size of the LD block. The default value is 10^5.')
    parser.add_argument('-o', '--offset', type=int, 
                        metavar='INT', default=0,
                        help='Shifts all the LD blocks by the requested base pairs. The default value is 0.')
    parser.add_argument('-m', '--min-reads', type=int, metavar='INT', default=2,
                        help='Takes into account only LD blocks with at least INT alleles. The default value is 2.')
    parser.add_argument('-M', '--max-reads', type=int, metavar='INT', default=2,
                        help='Randomly choose up to INT reads from each LD blocks. The default value is 16.'
                            'Filters a,b and c are applied sequentially, one after the other.')
    parser.add_argument('-O', '--output-filename', type=str, metavar='output_filename',  default='',
                        help='The output filename. The default is the input filename with the extension \".obs.p\" replaced by \".LLR#.p\".')
    args = parser.parse_args()
    
    
    LLR_dict, info = aneuploidy_test(**vars(args))
            
    sys.exit(0)
else: 
    print("The module ANEUPLOIDY_TEST was imported.")





"""  
if __name__ == "__main__": 
    print("Executed when invoked directly")
    a = time.time()
    
    args = dict(obs_filename = 'results_HapMix_EXT/mixed2haploids.X0.5.SRR10393062.SRR151495.0-2.hg38.obs.p',
                hap_filename = '../build_reference_panel/ref_panel.HapMix_EXT.hg38.BCFtools/chr21_HapMix_EXT_panel.hap',
                leg_filename = '../build_reference_panel/ref_panel.HapMix_EXT.hg38.BCFtools/chr21_HapMix_EXT_panel.legend',
                block_size = 2e5,
                offset = 0,
                min_reads = 16,
                max_reads = 16,
                output_filename = None)
     
    hap_tab = read_impute2(args['hap_filename'], filetype='hap')
    leg_tab = read_impute2(args['leg_filename'], filetype='leg')
    
    with open(args['obs_filename'], 'rb') as f:
        obs_tab = pickle.load(f)
        info = pickle.load(f)
        
    aux_dict = build_aux_dict(obs_tab,leg_tab)
    blocks_dict = build_blocks_dict(aux_dict,args['block_size'],args['offset'])
    
    random.seed(a=0, version=2) #I should set a=None after finishing to debug the code.
    reads_dict = build_reads_dict(obs_tab,leg_tab)
    blocks_dict_picked = {block: pick_reads(reads_dict,read_IDs,args['min_reads'],args['max_reads']) 
                               for block,read_IDs in blocks_dict.items()}
    
    #sys.exit(0)
    LLR = get_LLR(obs_tab, leg_tab, hap_tab, 'MODELS/MODELS16A.pbz2')
    LLR_dict = {block: LLR(*haplotypes) if haplotypes!=None else None for block,haplotypes in blocks_dict_picked.items()}
    
    population = tuple(value for value in LLR_dict.values() if value!=None)    
    mean, std = mean_and_std(population)
    
    w0 = [len(blocks_dict[key]) for key,value in LLR_dict.items() if value!=None]
    W = sum(w0)
    weights = [i/W for i in w0]
    jk_mean, jk_std, jk_bias = jackknifing(population,weights)
    
    num_of_LD_blocks = len(population)
    fraction_of_negative_LLRs = sum([1 for i  in population if i<0])/len(population)
    
    print('Depth: %.2f, Mean: %.3f, STD: %.3f' % (info['depth'], mean, std))
    print('Jackknife estimator: %.3f, Jackknife standard error: %.3f, Jackknife bias: %.3f' % (jk_mean, jk_std, jk_bias))
    print('Number of LD blocks: %d, Fraction of LD blocks with a negative LLR: %.3f' % (num_of_LD_blocks,fraction_of_negative_LLRs))
    
    info.update({'block_size': args['block_size'],
                 'offset': args['offset'],
                 'min_reads': args['min_reads'],
                 'max_reads': args['max_reads']})
    
    info['statistics'] = {'mean': mean, 'std': std,
                          'jk_mean': jk_mean, 'jk_std': jk_std, 'jk_bias': jk_bias,
                          'num_of_LD_blocks': num_of_LD_blocks,
                          'fraction_of_negative_LLRs': fraction_of_negative_LLRs} 
    
    output_filename, obs_filename = args['output_filename'], args['obs_filename']
    
    if output_filename!=None:
        default_filename = re.sub('(.*)obs','\\1LLR', obs_filename.split('/')[-1],1)    
        output_filename = default_filename if output_filename=='' else output_filename 
        with open( output_filename, "wb") as f:
            pickle.dump(obs_tab, f, protocol=4)
            pickle.dump(info, f, protocol=4)    
    
    b = time.time()
    print('Done calculating LLRs for all the LD block in %.3f sec.' % ((b-a)))
"""
