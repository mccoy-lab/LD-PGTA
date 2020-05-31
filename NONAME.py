#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 16:25:49 2020

@author: ariad
"""
import collections,time,pickle,statistics

from LLR_CALCULATOR import read_impute2
from LLR_CALCULATOR import wrapper_func_of_create_LLR_2 as get_LLR

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
        observed bases and their corresponding read name."""
        
    aux_dict = collections.defaultdict(list)
        
    for (pos, ind, read_id, allele) in obs_tab:
        if allele in leg_tab[ind][2:]:
            aux_dict[pos].append((read_id,allele))  
       
    return aux_dict 

def group_alleles(aux_dict,positions,**thresholds):
    """ For each chromosome position within the tuple positions, all the
    observed alleles together with their associted read id are extracted from
    the dictionary aux_dict, which was created by build_aux_dict. 
    Then, all the alleles are grouped into haplotypes (essetially tuples),
    according to the reads they origined from. All the haplotypes together form
    a tuple, named HAPLOTYPES, which is sorted in descending order according to
    the number of alleles in each haplotype. Moreover, only the 16 longest
    haplotypes are kept in HAPLOTYPES.
    In addition, overlapping haplotypes that have common alleles are identified.
    For each group of overlapping haplotypes a set of indices is formed; The
    haplotypes are indexed by their position in the tuple HAPLOTYPES. All these
    sets are kept within a tuple named OVERLAPS. Finally, the tuples HAPLOTYPES
    and OVERLAPS are returned. """
    
    reads = collections.defaultdict(list)
    for pos in positions:
        for (read_id,allele) in aux_dict[pos]:
            reads[read_id].append((pos,allele))
    
    ######################## APPLY THRESHOLDS AND FILTERS #####################    
    if thresholds.get('min_alleles_per_win',None):
        total_num_of_alleles = sum(len(alleles) for alleles in reads.values())
        if total_num_of_alleles<thresholds['min_alleles_per_win']: return None
    
    if thresholds.get('min_alleles_per_win',None):    
        for read_id,alleles in tuple(reads.items()):
            if len(alleles)<thresholds['min_alleles_per_read']: del reads[read_id]
    
    if thresholds.get('min_alleles_per_win',2):        
        if len(reads)<thresholds['min_reads_per_win']: return None
    ###########################################################################
    
    
    read_IDs, HAPLOTYPES = zip(*sorted(reads.items(), key=lambda x: len(x[1]), reverse=True)[:16])
    
    X = [{read_IDs.index(read_id) for (read_id,_) in aux_dict[pos] 
          if read_id in read_IDs}
                for pos in positions]

    for i in range(len(X)-1,-1,-1):
        for j in range(i-1,-1,-1):
            if not X[i].isdisjoint(X[j]):
                X[j].update(X.pop(i))
                break
    
    OVERLAPS = tuple(x for x in X if len(x)>1)
            
    return HAPLOTYPES, OVERLAPS
    
def build_windows_dict(positions,window_size,offset):
    """ Returns a dictionary that lists windows and gives all the SNP positions
        within them."""
    
    a = positions[0]-(window_size-1)+offset
    b = positions[-1]+(window_size-1)
    boundaries = tuple(range(int(a), int(b), int(window_size)))
    windows = [(i,j-1) for i,j in zip(boundaries,boundaries[1:])]
    
    windows_dict = collections.defaultdict(list)
    
    windows_iterator = iter(windows)
    window = next(windows_iterator, None)
    for p in positions:
        while window:
            if window[0]<=p<=window[1]:
                windows_dict[window].append(p)
                break
            window = next(windows_iterator, None)
    return windows_dict   

def main(obs_filename,leg_filename,hap_filename,window_size,offset,min_reads_per_window,min_alleles_per_read):
    """ TBA """
    a = time.time()
    
    hap_tab = read_impute2(hap_filename, filetype='hap')
    leg_tab = read_impute2(leg_filename, filetype='legend')
    
    with open(obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        info = pickle.load(f)
      
    aux_dict = build_aux_dict(obs_tab, leg_tab)
    
    args = dict(positions = tuple(aux_dict.keys()),
                window_size =  window_size,
                offset =  offset) 
    
    thresholds = dict(min_alleles_per_win = 1,
                      min_alleles_per_read = 1,
                      min_reads_per_win = 2)
    
    windows_dict = build_windows_dict(**args)
    windows_dict_grouped = {window: group_alleles(aux_dict,positions,**thresholds) for window,positions in windows_dict.items()}
    LLR = get_LLR(obs_tab, leg_tab, hap_tab)
    LLR_dict = {window: LLR(*arg) if arg!=None else None for window,arg in windows_dict_grouped.items()}
    
    b = time.time()
    print('Done in %.3f sec.' % ((b-a)))
    return LLR_dict
     
if __name__ == "__main__": 
    print("Executed when invoked directly")
    a = time.time()
    
    obs_filename = 'results/mixed2haploids.X0.1.SRR10393062.SRR151495.0-2.hg38.OBS.p'
    hap_filename = '../build_reference_panel/ref_panel.HapMix.hg38.BCFtools/chr21_HapMix_panel.hap'
    leg_filename = '../build_reference_panel/ref_panel.HapMix.hg38.BCFtools/chr21_HapMix_panel.legend'
    #models_filename = 'MODELS.p'
    
    
    hap_tab = read_impute2(hap_filename, filetype='hap')
    leg_tab = read_impute2(leg_filename, filetype='legend')
    
    with open(obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        info = pickle.load(f)
      
    aux_dict = build_aux_dict(obs_tab, leg_tab)
    
    args = {'positions': tuple(aux_dict.keys()),
            'window_size': 1e5,
            'offset': 0 } 
    
    windows_dict = build_windows_dict(**args)
    
    thresholds = dict(min_alleles_per_win = 1,
                      min_alleles_per_read = 1,
                      min_reads_per_win = 4)
    
    windows_dict_grouped = {window: group_alleles(aux_dict,positions,**thresholds) for window,positions in windows_dict.items()}
    LLR = get_LLR(obs_tab, leg_tab, hap_tab)
    LLR_dict = {window: LLR(*arg) for window,arg in windows_dict_grouped.items() if arg!=None}
    LLR_dict_without_nones = {key: value for key,value in LLR_dict.items() if value!=None} 
    
    w0 = [len(windows_dict[key]) for key in LLR_dict_without_nones]; W = sum(w0)
    weights = [i/W for i in w0]
    population = tuple(LLR_dict_without_nones.values())
    print(jackknifing(population,weights))
    print(len(population),sum([1 for i  in population if i<0])/len(population) )
    print(mean_and_std(population))
    #print(sum(LLR_dict_without_nones.values()))
        
    #A = group_alleles(aux_dict,positions[:20])
    
    b = time.time()

    print('Done in %.3f sec.' % ((b-a)))
    
else: 
    print("Executed when imported")