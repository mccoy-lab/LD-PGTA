#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAKE_STATISTICAL_MODEL

Builds statistical models for aneuploidy cells with BPH and SPH.
BPH (Both Parental Homologs) correspond to the presence of three unmatched
haplotypes, while SPH (Single Parental Homolog) correspond to chromosome gains
involving identical homologs.

Daniel Ariad (daniel@ariad.org)
AUG 31, 2020
"""

from itertools import product
from collections import defaultdict
from math import gcd as greatest_common_divisor
from pickle import dump
from time import time
from bz2 import BZ2File

def ENGINE(number_of_reads,degeneracies):
    """ Generates polysomy statistical models for n-reads, based of a list of
    the degeneracy of each homolog. """
    
    degeneracies_dict = {i:w for i,w in enumerate(degeneracies) if w>0}
    model = defaultdict(int)
    for sequence in product(degeneracies_dict, repeat=number_of_reads):
        haplotypes = defaultdict(list)
        weight = 1
        for read_ind,hap in enumerate(sequence): 
            haplotypes[hap].append(read_ind)
            weight *=  degeneracies_dict[hap]
        key = tuple(tuple(indices) for indices in haplotypes.values())
        model[key] += weight
    return model

def COMPACT(model,number_of_reads,degeneracies):
    """ The partitioning of reads can be encoded efficiently using the
    occupation basis. In this representation all the reads are enumerated.
    Each subset of reads is represented by a binary sequence, where the
    $i^\mathrm{th}$ element is one when the $i^\mathrm{th}$ read is included
    in the subset and zero otherwise. In other words, bits in a binary sequence
    map to whether a read is included in the subset. """
    
    compact = {i+1: defaultdict(list) for i in range(len(degeneracies))}
    T = sum(degeneracies)**number_of_reads
    while(len(model)!=0):
        haplotypes, weight = model.popitem()
        HAPLOTYPES = tuple(sum(1 << x for x in h) for h in sorted(haplotypes, key=len))
        gcd = greatest_common_divisor(weight,T)
        compact[len(HAPLOTYPES)][weight//gcd,T//gcd].append(HAPLOTYPES)
    for k1 in compact:
        compact[k1] = {k2:tuple(v) for k2,v in compact[k1].items()}
    return compact

def COMMON(model):
    """ In order to identify common multiples in the statistical model, we
    group partitions of reads according to their subset with the smallest
    cardinality. """
       
    nested = lambda: defaultdict(list)
    temp = defaultdict(nested)
    for normalized_weight, triplets in model[3].items():
        for haplotypes in triplets:
            temp[normalized_weight][haplotypes[0]].append((haplotypes[1],haplotypes[2]))
    for k1 in temp:
        model[3][k1] = {k2:tuple(v) for k2,v in temp[k1].items()} 
    return model

def SPH(number_of_reads):
    """ Builds a statistical model for n-reads under the SPH scenario. """

    degeneracies = (2, 1)
    model = ENGINE(number_of_reads,degeneracies)
    return COMPACT(model,number_of_reads,degeneracies)

def BPH(number_of_reads):
    """ Builds a statistical model for n-reads under the BPH scenario. """
    
    degeneracies = (1, 1, 1)
    model = ENGINE(number_of_reads,degeneracies)
    compact = COMPACT(model,number_of_reads,degeneracies)
    common = COMMON(compact)
    return(common)    
    
def BUILD(x):
    """ Build and store a dictionary with BPH and SPH statistical models for
    2-16 reads. """
    
    models = dict()
    for i in range(2,x+1):
        print('Building the statistical model for %d reads.' % i)
        a = time()
        models[i]= {'BPH': BPH(i), 'SPH': SPH(i)}
        b = time()
        print('Done building in %.3f sec.' % ((b-a)))
    with open( 'MODELS.p', 'wb') as f:
        dump(models, f, protocol=4)
    with BZ2File( 'MODELS.pbz2', 'wb') as f:
        dump(models, f, protocol=4)
    return models

if __name__ == "__main__":
    print('The module MAKE_STATISTICAL_MODEL was invoked directly.')
    models = BUILD(18)
else:
    print('The module MAKE_STATISTICAL_MODEL was imported.')