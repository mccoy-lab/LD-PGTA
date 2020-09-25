#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 00:32:35 2020

@author: ariad
"""

from math import log
from operator import countOf
from itertools import compress, islice, combinations
from collections import defaultdict
from time import time
from sys import stdout

def read_impute2(impute2_filename,**kwargs):
    """ Reads an IMPUTE2 file format (SAMPLE/LEGEND/HAPLOTYPE) and builds a list
        of lists, containing the dataset. """
    
    filetype = kwargs.get('filetype', None)
    with open(impute2_filename, 'r') as impute2_in:
        if filetype == 'leg':
            impute2_in.readline()   # Bite off the header
            def parse(x): 
                y=x.strip().split()
                y[0] = 'chr'+y[0].split(':')[0]
                y[1]=int(y[1])
                return tuple(y)
        elif filetype == 'hap':
            def parse(x):
                return tuple(i=='1' for i in x.strip().split())
        else:
            def parse(x):
                return tuple(x.strip().split()) # Trailing whitespaces stripped from the ends of the string. Then it splits the string into a list of words.
       
        impute2_tab = tuple(parse(line) for line in impute2_in)
    return impute2_tab

def transpose(dic):
    result = defaultdict(dict)
    for i,t in dic.items():
        for j,k in t.items():
            result[j][i] = k
    return result

def symmetrize(dic):
    result = defaultdict(dict)
    for i,t in dic.items():
        for j,k in t.items():
            result[j][i] = k
            result[i][j] = k
    return result

def remove(dic,pos):
    result = {i:{j: k for j,k in t.items() if j not in pos} for i,t in dic.items() if i not in pos} if len(pos) else dic
    return result

def aux(a,b):
    if a<=0.01 and b<=0.01:
        result = True
    elif (a<=1e-16 and b>=1e-16) or (a>=1e-16 and b<=1e-16):
        result = False
    else:
        result = -0.5<log(a/b)<0.5
    return result
    
def stability(A,B):
    result = {Ak: countOf((aux(a,b) for a,b in zip(Av.values(),Bv.values())),True)/(len(Av) or 1) 
                  for ((Ak, Av), (Bk, Bv)) in zip(A.items(),B.items())}
    return result

def get_common(leg_tab,hap_tab):
    superpop = {'AFR':(0,659),'AMR':(660,1006),'EAS':(1007,1510),'EUR':(1511,2013),'SAS':(2014,2502)}
    common = [all(countOf(r[a:b+1],True)%(b+1-a) for (a,b) in superpop.values()) for r in hap_tab]
    positions = [i[1] for i in leg_tab]
    hap_dict = {sp: {pos: hap[a:b+1] for pos,hap in zip(compress(positions,common),compress(hap_tab,common))} for sp,(a,b) in superpop.items()}
    return hap_dict

def build_LLR2_dict(hap_dict,population,max_dist):
    """ Builds a dictionary that contains the LLR of all SNP pairs (pos1,pos2)
    with pos1<pos2 and min_dist<pos2-pos1<mix_dist. The LLRs are stored in the
    nested dictionary, result[pos1][pos2]. 
    """

    result = defaultdict(dict)
    time0 = time()
    M, N = len(hap_dict[population]), len(hap_dict[population][next(iter(hap_dict[population]))])
    #freq = {pos: countOf(hap,True)/N for pos,hap in hap_dict[population].items()}
    
    hap_dict_pop = {pos: int(''.join('%d'%i for i in hap),2) for (pos,hap) in hap_dict[population].items()}     
    freq = {pos:  bin(hap).count('1')/N for pos,hap in hap_dict_pop.items()}

    for i,(pos1,hap1) in enumerate(hap_dict_pop.items()):
        for (pos2,hap2) in islice(hap_dict_pop.items(),i,None):
            if pos2-pos1 >= max_dist: break
            #joint_freq = countOf(compress(hap1,hap2),True)/N
            joint_freq = bin(hap1&hap2).count('1')/N
            result[pos1][pos2] = joint_freq/(freq[pos1]*freq[pos2])
            
        if (i+1)%1000==0:
            stdout.write('\r')
            stdout.write(f"Proceeded {(i+1)} / {M} rows. Average time per 1000 rows: {'%.2f' % (1000*(time() - time0)/(i+1))} sec")
            stdout.flush()

    return result

if __name__=='__main__':
    
    hap_filename = '../ref_panel/ALL_panel.hg38.BCFtools/chr21_ALL_panel.hap'
    leg_filename = '../ref_panel/ALL_panel.hg38.BCFtools/chr21_ALL_panel.legend'
    #leg_tab,hap_tab,obs_tab = load(obs_filename,leg_filename,hap_filename)
    leg_tab = read_impute2(leg_filename, filetype='leg')
    hap_tab = read_impute2(hap_filename, filetype='hap')  
    
    hap_dict = get_common(leg_tab,hap_tab)
    print('---')
    R = {i: build_LLR2_dict(hap_dict,i,100000) for i in ('AFR','EUR','SAS','EAS','AMR')} 
    
    import numpy as np
    
    R2 = {i: symmetrize(R[i]) for i in ('AFR','EUR','SAS','EAS','AMR')} 
    threshold = step = 1e-2
    while(threshold<1):
        t0 = time()
        q = stability(R2['EUR'],R2['AFR'])
        pos = {*compress(q, (x<threshold for x in q.values()))}
        step *= 2 if len(pos)<(len(R2['EUR'])//1000) else 0.5    
        threshold += step
        R2 = {j: remove(R2[j],pos) for j in ('AFR','EUR')}
        t1 = time()
        print(threshold,step,len(pos),t1-t0)
        
    print(len(R2['EUR']))
    POSITIONS = tuple(R2['EUR'].keys())
    print(R2['EUR'][POSITIONS[500]])
    print(R2['AFR'][POSITIONS[500]])
    
    pos = set(R2['SAS']).difference(R2['EUR'])
    R3 = {j: remove(R2[j],pos) for j in ('AFR','EUR','SAS','EAS','AMR')}
