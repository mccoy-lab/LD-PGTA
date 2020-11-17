#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 10:31:24 2020

@author: ariad
"""
import sys, argparse, itertools, collections, operator, pickle, math, time, random, re

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

def build_hap_dict(obs_tab,leg_tab,hap_tab):
    """ Returns a dictionary that lists chromosome positions of SNPs and gives
        their relevent row from haplotypes table. The row is stored as a tuple
        of booleans, where True represents the observed allele. We denote
        the returned dictionary as the reference panel."""
        
    hap_dict = dict()
    mismatches = 0
        
    for (pos, ind, read_id, read) in obs_tab:
        chr_id, pos2, ref, alt = leg_tab[ind]
        if pos!=pos2:
            raise Exception('error: the line numbers in obs_tab refer to the wrong chromosome positions in leg_tab.')         
        if read==alt:
            hap_dict[(pos,read)] = tuple(hap_tab[ind])
        elif read==ref:
            hap_dict[(pos,read)] = tuple(map(operator.not_,hap_tab[ind]))
        else:
            mismatches += 1
    
    print('%%%.2f of the reads matched known alleles.' % (100*(1-mismatches/len(obs_tab))))
    
    return hap_dict

def samples_info():
    with open('igsr_samples.txt', 'r') as f:
        f.readline()
        tab = [line.strip('\n').split('\t') for line in f]
    samples_dict = {row[0]:row[1:] for row in tab}
    return samples_dict

def main(obs_filename,leg_filename,hap_filename, sam_filename):
    """ main function. """
    time0 = time.time()
    
    with open(obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        #info = pickle.load(f)
        
    leg_tab = read_impute2( leg_filename ,filetype='legend')
    hap_tab = read_impute2( hap_filename ,filetype='hap')
    sam_tab = read_impute2( sam_filename ,filetype='samples')
    
    hap_dict = build_hap_dict(obs_tab, leg_tab, hap_tab)
    
    genotypes = (' '.join(genotype) for genotype in itertools.product(next(zip(*sam_tab)),('A','B')))
    hap = zip(*iter(hap_dict.values()))
    hap_indv_dict = dict(zip(genotypes,hap))
    sd = samples_info()
    #positions = tuple(hap_pos_dict.keys())
    
    A1 = [(i,sum(j)) for i,j in hap_indv_dict.items()]
    A2 = sorted(A1, key=operator.itemgetter(1))[-len(A1)//250:]
    A3 = tuple(zip(*A2))[0]
    A4 = [sd[i[:-2]][3] for i in A3]
    A5 = collections.Counter(A4)
    A6 = {a: 100*b/len(A4) for a,b in dict(A5).items()}
    
    time1 = time.time()
    print('Done in %.2f sec.' % (time1-time0))
    
    return A6

def demo(filename):
    mydict = {'obs_filename': 'results_ALL_EXT/'+filename,
              'leg_filename': '../build_reference_panel/ref_panel.ALL_EXT.hg38.BCFtools/chr21_ALL_EXT_panel.legend', 
              'hap_filename': '../build_reference_panel/ref_panel.ALL_EXT.hg38.BCFtools/chr21_ALL_EXT_panel.hap', 
              'sam_filename': '../build_reference_panel/ref_panel.ALL_EXT.hg38.BCFtools/chr21_ALL_EXT_panel.samples'}
    ancestry =  main(**mydict)
    
    return ancestry
    
    
       
    
