#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 00:32:35 2020

@author: ariad
"""
from math import log
from operator import countOf
from itertools import compress, islice
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
       
        for line in impute2_in:
            yield parse(line)

def get_common(leg_iter,hap_iter):
    superpop = {'AFR':(0,659),'AMR':(660,1006),'EAS':(1007,1510),'EUR':(1511,2013),'SAS':(2014,2502)}
    hap_dict = defaultdict(dict)
    for leg, hap in zip(leg_iter,hap_iter):
        if all(5<countOf(hap[a:b+1],True)<(b+1-a-5) for (a,b) in superpop.values()):
            for sp,(a,b) in superpop.items():
                hap_dict[sp][leg[1]] = int(''.join('%d'%i for i in hap[a:b+1]),2) 
    N = len(hap)
    return hap_dict, N
            
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

def keep(dic,pos):
    result = {i:{j: k for j,k in t.items() if j in pos} for i,t in dic.items() if i in pos}
    return result

def aux(a,b):
    if a<=1e-16 and b<=1e-16:
        result = True
    elif (a<1e-16 and b>1e-16) or (a>1e-16 and b<1e-16):
        result = False
    else:
        result = -0.7<log(a/b)<0.7
    return result
    
def stability(A,B):
    result = {Ak: countOf(map(aux,Av.values(),Bv.values()),True)/(len(Av) or 1)
              for ((Ak, Av), (Bk, Bv)) in zip(A.items(),B.items())}
    return result
    
def build_LLR2_dict(hap_dict,N,max_dist):
    """ Builds a dictionary that contains the LLR of all SNP pairs (pos1,pos2)
    with pos1<pos2 and min_dist<pos2-pos1<mix_dist. The LLRs are stored in the
    nested dictionary, result[pos1][pos2]. 
    """

    time0 = time()
    result = defaultdict(dict)
    M = len(hap_dict)
    freq = {pos:  bin(hap).count('1')/N for pos,hap in hap_dict.items()}

    for i,(pos1,hap1) in enumerate(hap_dict.items()):
        for (pos2,hap2) in islice(hap_dict.items(),i,None):
            if pos2-pos1 >= max_dist: break
            joint_freq = bin(hap1&hap2).count('1')/N
            result[pos1][pos2] = joint_freq/(freq[pos1]*freq[pos2])
            
        if (i+1)%1000==0:
            stdout.write('\r')
            stdout.write(f"Proceeded {(i+1)} / {M} rows. Average time per 1000 rows: {'%.2f' % (1000*(time() - time0)/(i+1))} sec")
            stdout.flush()
    stdout.write('\n')
    return result

def filtering(A,B):
    A, B = symmetrize(A), symmetrize(B) 
    threshold = step = 1e-4
    while(threshold<1):
        t0 = time()
        q = stability(A, B)
        pos = {*compress(q, (x<threshold for x in q.values()))}
        step *= 2 if len(pos)<(len(A)//2000) else 0.5    
        threshold += step
        A, B = remove(A,pos), remove(B,pos)
        t1 = time()
        print(threshold, step, len(pos), t1-t0)
    return A, B

def filter_impute2(read_filename, write_filename, lines):
    with open(read_filename, 'r') as impute2_in:
        with open(write_filename, 'w') as impute2_out:
                for i,line in enumerate(impute2_in):
                    if i in lines:
                        impute2_out.write(line)
    return 0

if __name__=='__main__':
    time0 = time()

    hap_filename = '../ref_panel/ALL_panel.hg38.BCFtools/chr21_ALL_panel.hap'
    leg_filename = '../ref_panel/ALL_panel.hg38.BCFtools/chr21_ALL_panel.legend'
    #leg_tab,hap_tab,obs_tab = load(obs_filename,leg_filename,hap_filename)
    leg_tab = read_impute2(leg_filename, filetype='leg')
    hap_tab = read_impute2(hap_filename, filetype='hap')      
    hap_dict, N = get_common(leg_tab,hap_tab)
    EUR = build_LLR2_dict(hap_dict['EUR'],N,50000)
    POSITIONS = {*EUR}
    for sp in ('AFR','AMR','EAS','SAS'):
        print('SUPERPOPULATION: %s' % sp)
        EUR, SP = filtering(EUR, keep(build_LLR2_dict(hap_dict[sp], N, 50000), POSITIONS))
        POSITIONS = {*EUR}
        print('Done.',len(POSITIONS))
    lines = [i for i,j in enumerate(leg_tab) if j[1] in POSITIONS]
    filter_impute2('../ref_panel/ALL_panel.hg38.BCFtools/chr21_ALL_panel.hap', 'chr21_COMMON_panel.hap', lines)
    filter_impute2('../ref_panel/ALL_panel.hg38.BCFtools/chr21_ALL_panel.legend', 'chr21_COMMON_panel.legend', lines)
    time1 = time()
    print(f'Done in {(time1-time0):.2f} sec.')