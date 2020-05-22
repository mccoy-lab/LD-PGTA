#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 16:25:49 2020

@author: ariad
"""
import collections,time,pickle

from LLR_CALCULATION import wraps_create_LLR2

def read_impute2(impute2_filename,**kwargs):
    """ Reads an IMPUTE2 file format (SAMPLE/LEGEND/HAPLOTYPE) and builds a list
        of lists, containing the dataset. """
    
    filetype = kwargs.get('filetype', None)
    with open(impute2_filename, 'r') as impute2_in:
        if filetype == 'legend':
            impute2_in.readline()   # Bite off the header
            def parse(x): 
                y=x.strip().split()
                y[0] = 'chr'+y[0].split(':')[0]
                y[1]=int(y[1])
                return y
        elif filetype == 'hap':
            def parse(x):
                return [i=='1' for i in x.strip().split()]
        else:
            def parse(x):
                return x.strip().split() # Trailing whitespaces stripped from the ends of the string. Then it splits the string into a list of words.
       
        impute2_tab = [parse(line) for line in impute2_in]
    return impute2_tab

def build_aux_dict(obs_tab,leg_tab):
    """ Returns a dictionary that lists chromosome positions of SNPs and gives
        observed bases and their corresponding read name."""
        
    aux_dict = collections.defaultdict(list)
        
    for (pos, ind, read_id, read) in obs_tab:
        if read in leg_tab[ind][2:]:
            aux_dict[pos].append((read_id,read))  
       
    return aux_dict 

def group_alleles(aux_dict,positions):
    """ Returns a dictionary that lists reads IDs and gives all the SNPs within
        the corresponding read."""
    
    alleles = collections.defaultdict(list)
    for pos in positions:
        for (read_id,read) in aux_dict[pos]:
            alleles[read_id].append((pos,read))
    grouped_alleles = sorted(alleles.values(), key=len, reverse=True)
    return grouped_alleles

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
     
if __name__ == "__main__": 
    print("Executed when invoked directly")
    a = time.time()
    
    
    obs_filename = 'results/mixed2haploids.X0.05.SRR10393062.SRR151495.0-2.hg38.OBS.p'
    hap_filename = '../make_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.hap'
    leg_filename = '../make_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.legend'
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
    
    windows_dict_grouped = {key: group_alleles(aux_dict,value) for key,value in windows_dict.items()} 
    
    #windows_dict_filtered = {key: value for key,value in windows_dict_grouped.items() if len(value)>1}
    
    #LLR = wraps_create_LLR1(obs_filename,leg_filename,hap_filename,models_filename)
    LLR = wraps_create_LLR2(obs_tab, leg_tab, hap_tab)
    LLR_dict = {key: LLR(*value[:16]) for key,value in windows_dict_grouped.items() if len(value)>1}
    
    K = [i for i in LLR_dict.values() if i!=None]
    
    print(sum(K))
    
    ########### SEND TO LLR ##########
    
    #A = group_alleles(aux_dict,positions[:20])
    
    b = time.time()

    print('Done in %.3f sec.' % ((b-a)))
    
else: 
    print("Executed when imported")