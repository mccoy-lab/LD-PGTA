#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Checks the rate of heterozygosity. 

Created on Thu Jun  4 18:51:19 2020

@author: ariad
"""
import collections, operator, os, pickle, itertools

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

def build_hap_dict(leg_tab,hap_tab):
    N = len(hap_tab[0])
    hap_dict = {pos: tuple(hap_tab[ind]) 
                     for ind,(chr_id, pos, ref, alt) in enumerate(leg_tab) if 0.05<hap_tab[ind].count(True)/N<0.95}
    
    return hap_dict

            
def HIST(x):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    #ax.hist(x,bins=int(len(x)**.5),histtype='step', linewidth=2.2, label='TBA')
    result = ax.hist(x,bins=600, linewidth=0.2, label='TBA', histtype='stepfilled')
    ax.set_xlabel('Positions')
    ax.set_ylabel('Counts')
    ax.set_title('distribution of SNPs along chromosome 21' )
    ax.set_xticks(range(min(x), max(x),5000000)) #(max(x)-min(x))//20))
    #ax.set_xticks(range(min(x), max(x),100000)) #(max(x)-min(x))//20))
    #ax.legend()
    plt.show()
    return result


def check_var(hap_tab):
    import numpy as np
    K = lambda x:  np.sum(0<np.mod(np.sum(x,axis=1),x.shape[1]))
    A = np.array(hap_tab)
    B0 = [np.sum(A[:,i]!=A[:,i+1]) for i in range(0,len(A[0,:]),2)]
    B = [K(A[:,i:i+2]) for i in range(0,len(A[0,:]),2)]
    C = [K(A[:,i:i+4]) for i in range(0,len(A[0,:]),4)]
    b0 = np.sum(B0)/len(B0)/len(A)
    b = np.sum(B)/len(B)/len(A)
    c = np.sum(C)/len(C)/len(A)
    return (b0,b,c)

def check_var2(hap_dict):
    from operator import ne
    N = len(next(iter(hap_dict.values())))//2
    result = {pos: sum(map(ne,hap[1::2],hap[0::2]))/N for pos,hap in hap_dict.items()}
    return result

def check_var3(hap_dict):
    N = len(next(iter(hap_dict.values())))//3
    C = lambda *x: not (all(x) or not any(x))
    result = {pos: sum(map(C,hap[0::3],hap[1::3],hap[2::3]))/N for pos,hap in hap_dict.items()}
    return result
    
if __name__=='__main__':
    leg_filename = '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.legend'
    hap_filename = '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.hap'
    leg_tab = read_impute2(leg_filename,filetype='leg')
    hap_tab = read_impute2(hap_filename,filetype='hap')
    hap_dict = build_hap_dict(leg_tab,hap_tab)
    A = check_var2(hap_dict)
    B = check_var3(hap_dict)
    from matplotlib import pyplot
    pyplot.scatter(A.keys(),[1/i if i else -10 for i in A.values()],s=1)
    pyplot.scatter(B.keys(),[1/i if i else -10 for i in B.values()],s=1)
