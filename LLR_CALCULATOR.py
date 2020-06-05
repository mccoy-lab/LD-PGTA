#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
LLR_CALCULATOR

This script creates a function that calculates log-likelihood BPH/SPH ratio
(LLR) for a tuple of SNPs. BPH (Both Parental Homologs) correspond to the
presence of three unmatched haplotypes, while SPH (Single Parental Homolog)
correspond to chromosome gains involving identical homologs.  

Daniel Ariad (daniel@ariad.org)
May 4th, 2020
"""

import operator, itertools, pickle, math, os, sys

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

def create_frequencies(hap_dict):
    """ Returns the functions joint_frequencies, which extracts from the
        dictionary hap_dict the joint frequencies of observed alleles. """
    
    N = len(next(iter(hap_dict.values())))

    def intrenal_hap_dict(*alleles):
        """ This function allows treatment of alleles and haplotypes on an
        equal footing. This is done in three steps: (1) All the alleles and
        haplotypes are enumerated. (2) For each given haplotype, tuples in the
        reference panels, corresponding to the haplotype's alleles, are
        intersected. (3) A dictionary that lists all the alleles and haplotypes
        by their index is returned. The dictionary gives for each allele and
        haplotype their associated tuple and intersected tuple, respectively. """ 
        hap = dict()
            
        for i, X in enumerate(alleles):
            if type(X[0])==tuple: #Checks if X is a tuple/list of alleles.
                n = len(X)
                if n==1: 
                    hap[i] = hap_dict[X[0]] 
                elif n==2:
                    hap[i] = tuple(map(operator.and_,hap_dict[X[0]],hap_dict[X[1]]))
                else:
                    hap[i] = tuple(map(all, zip(*(hap_dict[allele] for allele in X))))
            elif type(X[0])==int: #Checks if X is a single allele.
                hap[i] = hap_dict[X]
            else:
                raise Exception('error: joint_frequencies only accepts alleles and tuple/list of alleles.')
       
        return hap                

    def joint_frequencies(*alleles):
        """ Based on the reference panel, it calculates joint frequencies of
            observed alleles. The function arguments are alleles, that is,
            tuples of position and base, e.g., (100,'T'), (123, 'A') and
            (386, 'C'). Each allele is enumerated according to the order it was
            received by the function. The function returns a dictionary that
            lists all the possible subgroups of the given alleles. Each key in
            the dictionary is a tuple of intergers that are lexicographically
            sorted. Moreover, each integer within the keys corresponds to an
            enumerated allele. For each subgroup of alleles the dictionary
            gives the joint frequencies in a given population. The function
            arguments can also include haplotypes, that is, tuples of alleles;
            Haplotypes are treated in the same manner as alleles. """
           
        hap = intrenal_hap_dict(*alleles)   
           
        result = {(i,): A.count(True) / N  for i,A in hap.items() }
        
        for i in itertools.combinations(hap.items(), 2):
            a,b = zip(*i)
            result[a] = sum(itertools.compress(*b)) / N 
        
        for i in itertools.combinations(hap.items(), 3):
            a,b = zip(*i)
            result[a] = sum(itertools.compress(b[0], map(operator.and_,b[1],b[2]))) / N
        
        for r in range(4,len(alleles)):
            for i in itertools.combinations(hap.items(), r):
                a,b = zip(*i)
                result[a] = sum(itertools.compress(b[0], map(all,zip(*b[1:])))) / N
        
        if len(alleles)>=4:
            a,b = zip(*hap.items())
            result[a] = sum(itertools.compress(b[0], map(all,zip(*b[1:])))) / N
        
        return result

    return joint_frequencies


def create_LLR(models_dict,joint_frequencies):
    """ This function receives the dictionary models_dict with the
    statisitcal models and the function frequncies, which calculates
    joint frequncies. Based on these arguments it creates the function
    LLR, which calculates the log-likelihood BPH/SPH ratio."""
    
    def LLR(*alleles):
        """ Calculates the log-likelihood BPH/SPH ratio for a tuple that contains
        alleles and haplotypes. """
                
        l = len(alleles)
        freq = joint_frequencies(*alleles)
        BPH = sum(A[0]/A[1] * sum(math.prod(freq[b] for b in B) for B in C)
                   for A,C in models_dict[l]['BPH'].items())
        SPH = sum(A[0]/A[1] * sum(math.prod(freq[b] for b in B) for B in C) 
                   for A,C in models_dict[l]['SPH'].items())
        
        result = None if SPH<1e-16 else math.log(BPH/SPH)
                
        return result
    
    return LLR

def wrapper_func_of_create_LLR(obs_tab,leg_tab,hap_tab):
    """ Wraps the fuction create_LLR. It receives an observations array, legend
        array and haplotypes array. Based on the given data it creates and
        returns the function LLR."""
        
    if not os.path.isfile('MODELS.p'): raise Exception('Error: MODELS file does not exist.')
    with open('MODELS.p', 'rb') as models:
        models_dict = pickle.load(models)
    
    LLR = create_LLR(models_dict,create_frequencies(build_hap_dict(obs_tab, leg_tab, hap_tab)))
    return LLR

def wrapper_func_of_create_LLR_for_debugging(obs_filename,leg_filename,hap_filename,models_filename):
    """ Wraps the function create_LLR. It receives an observations file, IMPUTE2
        legend file, IMPUTE2 haplotypes file, and the statistical model. Based
        on the given data it creates and returns the LLR function."""
    
    from MAKE_OBS_TAB import read_impute2
    
    if not os.path.isfile(obs_filename): raise Exception('Error: OBS file does not exist.')
    if not os.path.isfile(leg_filename): raise Exception('Error: LEGEND file does not exist.')
    if not os.path.isfile(hap_filename): raise Exception('Error: HAP file does not exist.')
    if not os.path.isfile(models_filename): raise Exception('Error: MODELS file does not exist.')

    leg_tab = read_impute2(leg_filename, filetype='leg')
    hap_tab = read_impute2(hap_filename, filetype='hap')
    with open(obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        #info = pickle.load(f)
    with open(models_filename, 'rb') as f:
        models_dict = pickle.load(f)
    
    hap_dict = build_hap_dict(obs_tab, leg_tab, hap_tab)
    frequencies = create_frequencies(hap_dict)
    LLR = create_LLR(models_dict,frequencies)     
    
    ###LLR = create_LLR(models_dict,create_frequencies(build_hap_dict(obs_tab, leg_tab, hap_tab))) #This line replaces the three lines above.
    return LLR
        
###############################################################################

if __name__ != "__main__": 
    print("The module LLR_CALCULATOR was imported.")   
else:
    print("Executed when the module LLR_CALCULATOR is invoked directly")
    sys.exit(0)
   