#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 23:18:10 2020

@author: ariad
"""

import pickle

from operator import not_, itemgetter, countOf
from itertools import combinations
from math import log
from time import time

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

def build_hap_dict(obs_tab,leg_tab,hap_tab):
    hap_dict = dict()
    for (pos, ind, read_id, base) in obs_tab:
        chr_id, pos2, ref, alt = leg_tab[ind]
        if base==alt:
            hap_dict[(pos,base)] = hap_tab[ind]
        elif base==ref:
            hap_dict[(pos,base)] = tuple(map(not_,hap_tab[ind]))
    return hap_dict

def intrenal_hap_dict(*alleles,hap_dict):
    hap = dict()
    for i, X in enumerate(alleles,start=65):
        if type(X[0])==int:
            hap[chr(i)] = hap_dict[X]
        elif len(X)==1:
            hap[chr(i)] = hap_dict[X[0]]
        else:
            hap[chr(i)] = tuple(map(all, zip(*itemgetter(*X)(hap_dict))))
    return hap

def joint_frequencies_combo(*alleles,hap_dict,norm_const,normalize):
    N = norm_const if normalize else 1
    hap = intrenal_hap_dict(*alleles,hap_dict=hap_dict)
    result = {c:  A.count(True) / N  for c,A in hap.items() }
    for r in range(2,len(alleles)+1):
        for A in combinations(hap, r):
            result[''.join(A)] = countOf(map(all,zip(*itemgetter(*A)(hap))),True) / N
    return result

def LLR2(*alleles,hap_dict,N):
    F = joint_frequencies_combo(*alleles,hap_dict=hap_dict,norm_const=N,normalize=True)
    a, b, ab = F['A'], F['B'], F['AB']
    BPH = (ab+2*a*b)/3 #The probability of three unmatched haplotypes.
    SPH = (5*ab+4*a*b)/9 #The probability of two identical haplotypes out three.
    result = log(BPH/SPH)
    return result

def LLR3(*alleles,hap_dict,N):
    F = joint_frequencies_combo(*alleles,hap_dict=hap_dict,norm_const=N,normalize=True)
    a, b, c, ab, ac, bc, abc = F['A'], F['B'], F['C'], F['AB'], F['AC'], F['BC'], F['ABC']
    BPH = (abc+2*(ab*c+ac*b+bc*a+a*b*c))/9 #The probability of three unmatched haplotypes.
    SPH = abc/3+2*(ab*c+ac*b+bc*a)/9  #The probability of two identical haplotypes out three.
    result = None if SPH<1e-16 else log(BPH/SPH)
    return result

def LLR4(*alleles,hap_dict,N):
    F = joint_frequencies_combo(*alleles,hap_dict=hap_dict,norm_const=N,normalize=True)
    a, b, c, d = F['A'], F['B'], F['C'], F['D'],
    ab, ac, ad, bc, bd, cd = F['AB'], F['AC'], F['AD'], F['BC'], F['BD'], F['CD'],
    abc, abd, acd, bcd = F['ABC'], F['ABD'], F['ACD'], F['BCD']
    abcd = F['ABCD']
    BPH = (abcd+2*(ab*c*d+a*bd*c+a*bc*d+ac*b*d+a*b*cd+ad*b*c+abc*d+a*bcd+acd*b+abd*c+ab*cd+ad*bc+ac*bd))/27  #The probability of three unmatched haplotypes.
    SPH = (17*abcd+10*(abc*d+bcd*a+acd*b+abd*c)+8*(ab*cd+ad*bc+ac*bd))/81  #The probability of two identical haplotypes out three.
    result = None if SPH<1e-16 or BPH<1e-16 else log(BPH/SPH)
    return result

if __name__ != "__main__":
    print("The module was imported.")
else:
    print("Executed when invoked directly")
    a = time()
    obs_filename = 'results_EUR/mixed2haploids.X0.5.HG00096.HG00096.B.hg38.obs.p'
    hap_filename = '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.hap'
    leg_filename = '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.legend'

    hap_tab = read_impute2(hap_filename, filetype='hap')
    leg_tab = read_impute2(leg_filename, filetype='leg')


    with open(obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        #info = pickle.load(f)

    hap_dict = build_hap_dict(obs_tab, leg_tab, hap_tab)
    #aux_dict = build_aux_dict(obs_tab, leg_tab)

    positions = tuple(hap_dict.keys())
    N = len(hap_tab[0])
    print('-----joint_frequencies_combo-----')
    print(joint_frequencies_combo(positions[0],hap_dict=hap_dict,norm_const=N,normalize=False))
    print(joint_frequencies_combo(positions[:4],hap_dict=hap_dict,norm_const=N,normalize=False))
    print('-----LLR4-Haplotypes-----')
    pos = (positions[:4],positions[4:8],positions[8:12],positions[12:16])
    print(LLR4(*pos,hap_dict=hap_dict,N=N))
    print('-----LLR2-----')
    print(LLR2(*positions[:2],hap_dict=hap_dict,N=N))
    print('-----LLR3-----')
    print(LLR3(*positions[:3],hap_dict=hap_dict,N=N))
    print('-----LLR4-----')
    print(LLR4(*positions[:4],hap_dict=hap_dict,N=N))

    b = time()

    print('Done in %.3f sec.' % ((b-a)))
