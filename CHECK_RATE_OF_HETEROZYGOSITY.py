#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Checks the rate of heterozygosity. 

Created on Thu Jun  4 18:51:19 2020

@author: ariad
"""
import collections, operator, os, pickle, itertools, gzip

def read_impute2(filename,**kwargs):
    """ Reads an IMPUTE2 file format (LEGEND/HAPLOTYPE/SAMPLE) and builds a list
        of lists, containing the dataset. """
    
    filetype = kwargs.get('filetype', None)
    
    def leg_format(line):
        rs_id, pos, ref, alt = line.strip().split()
        return ('chr'+rs_id[:2].rstrip(':'), int(pos), ref, alt)  
    
    def sam_format(line):
          sample_id, group1, group2, sex = line.strip().split(' ')
          return (sample_id, group1, group2, int(sex))  
      
    with (gzip.open(filename,'rt') if filename[-3:]=='.gz' else open(filename, 'r')) as impute2_in:
        if filetype == 'leg': 
            impute2_in.readline()   # Bite off the header
            result = tuple(map(leg_format,impute2_in))
            
        elif filetype == 'hap':
            firstline = impute2_in.readline()   # Get first line
            a0 = int(firstline.replace(' ', ''), 2)
            a1 = (int(line.replace(' ', ''), 2) for line in impute2_in)
            hap_tab = (a0, *a1)
            number_of_haplotypes = len(firstline.strip().split())
            result = hap_tab, number_of_haplotypes
            
        elif filetype == 'sam': 
            impute2_in.readline()   # Bite off the header
            result = tuple(map(sam_format,impute2_in))
        
        else:
            result = tuple(line.strip().split() for line in impute2_in)
    
    return result 


def build_hap_dict(leg_tab,hap_tab,number_of_haplotypes):
    N = number_of_haplotypes
    hap_dict = {pos: hap
                     for (chr_id, pos, ref, alt), hap in zip(leg_tab,hap_tab) if 0.05<bin(hap).count('1')/N<0.95}
    
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

def check_var2(hap_dict,number_of_haplotypes):
    N = number_of_haplotypes/2
    result = {pos: bin(hap ^ (hap >> 1))[2::2].count('1')/N for pos,hap in hap_dict.items()}
    return result

def check_var3(hap_dict,number_of_haplotypes):
    N = number_of_haplotypes/3
    result = {pos: bin( (hap ^ (hap >> 1)) | (hap ^ (hap >> 2)) | ((hap >> 1) ^ (hap >> 2)) )[2::3].count('1')/N for pos,hap in hap_dict.items()}
    return result
    
if __name__=='__main__':
    leg_filename = '../../build_reference_panel/EUR_panel.hg38.BCFtools/chr21_EUR_panel.legend.gz'
    hap_filename = '../../build_reference_panel/EUR_panel.hg38.BCFtools/chr21_EUR_panel.hap.gz'
    leg_tab = read_impute2(leg_filename,filetype='leg')
    hap_tab, number_of_haplotypes = read_impute2(hap_filename,filetype='hap')
    hap_dict = build_hap_dict(leg_tab,hap_tab,number_of_haplotypes)
    A = check_var2(hap_dict,number_of_haplotypes)
    B = check_var3(hap_dict,number_of_haplotypes)
    #from matplotlib import pyplot
    #pyplot.scatter(A.keys(),[1/i if i else -10 for i in A.values()],s=1)
    #pyplot.scatter(B.keys(),[1/i if i else -10 for i in B.values()],s=1)
