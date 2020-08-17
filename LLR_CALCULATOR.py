#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
LLR_CALCULATOR

This script creates a function that calculates log-likelihood BPH/SPH ratio
(LLR) for a given tuple of reads that originated form the same LD block.
BPH (Both Parental Homologs) correspond to the presence of three unmatched
haplotypes, while SPH (Single Parental Homolog) correspond to chromosome gains
involving identical homologs.

Daniel Ariad (daniel@ariad.org)
Aug 2th, 2020
"""

import pickle, os, sys, bz2

from functools import reduce
from operator import not_, and_, itemgetter
from itertools import combinations
from math import log

try:
    from gmpy2 import popcount
except:
    print('caution: cound not import the gmpy2 module.')
    def popcount(x):
        """ Counts non-zero bits in positive integer. """
        return bin(x).count('1')

def build_hap_dict(obs_tab,leg_tab,hap_tab):
    """ Returns a dictionary that lists chromosome positions of SNPs and gives
        their relevent row from haplotypes table. The row is stored as a tuple
        of booleans, where True represents the observed allele. We denote
        the returned dictionary as the reference panel."""

    hap_dict = dict()
    mismatches = 0

    bools2int = lambda x: int(''.join(chr(48+i) for i in x),2)

    for (pos, ind, read_id, base) in obs_tab:
        chr_id, pos2, ref, alt = leg_tab[ind]
        if pos!=pos2:
            raise Exception('error: the line numbers in obs_tab refer to the wrong chromosome positions in leg_tab.')
        if base==alt:
            hap_dict[(pos,base)] = bools2int(hap_tab[ind])
        elif base==ref:
            hap_dict[(pos,base)] = bools2int(map(not_,hap_tab[ind]))
        else:
            mismatches += 1

    print('%%%.2f of the reads matched known alleles.' % (100*(1-mismatches/len(obs_tab))))

    return hap_dict

def create_frequencies(hap_dict):
    """ Returns the function combo_joint_frequencies, which extracts from the dictionary
        hap_dict the joint frequencies of observed alleles. """

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
                    hap[1 << i] = hap_dict[X[0]]
                elif n==2:
                    hap[1 << i] = hap_dict[X[0]] & hap_dict[X[1]]
                else:
                    hap[1 << i] = reduce(and_,itemgetter(*X)(hap_dict))

            elif type(X[0])==int: #Checks if X is a single allele.
                hap[1 << i] = hap_dict[X]
            else:
                raise Exception('error: joint_frequencies only accepts alleles and tuple/list of alleles.')

        return hap

    def joint_frequencies_combo(*alleles):
        """ Based on the reference panel, it calculates unnormalized joint
            frequencies of observed alleles. The function arguments are alleles,
            that is, tuples of position and base, e.g., (100,'T'), (123, 'A')
            and (386, 'C'). Each allele is enumerated according to the order it
            was received by the function. The function returns a dictionary that
            lists all the possible subgroups of the given alleles. Each key in
            the dictionary is a tuple of intergers that are lexicographically
            sorted. Moreover, each integer within the keys corresponds to an
            enumerated allele. For each subgroup of alleles the dictionary
            gives the joint frequencies in a given population. The function
            arguments can also include haplotypes, that is, tuples of alleles;
            Haplotypes are treated in the same manner as alleles. """

        hap = intrenal_hap_dict(*alleles)

        result = {c: popcount(A) for c,A in hap.items() }

        for C in combinations(hap, 2):
            result[C[0]|C[1]] = popcount(hap[C[0]]&hap[C[1]])

        for C in combinations(hap, 3):
            result[C[0]|C[1]|C[2]] = popcount(hap[C[0]]&hap[C[1]]&hap[C[2]])

        for r in range(4,len(alleles)):
            for C in combinations(hap, r):
                result[sum(C)] = popcount(reduce(and_,itemgetter(*C)(hap)))

        if len(alleles)>=4:
            result[sum(hap.keys())] = popcount(reduce(and_,hap.values()))

        return result

    return joint_frequencies_combo

def create_LLR(models_dict,joint_frequencies_combo,D):
    """ This function receives the dictionary models_dict with the
    statisitcal models and the function frequncies, which calculates
    joint frequncies. Based on these arguments it creates the function
    LLR, which calculates the log-likelihood BPH/SPH ratio."""

    def LLR(*alleles):
        """ Calculates the log-likelihood BPH/SPH ratio for a given tuple of
        alleles and haplotypes. """
        
        F = joint_frequencies_combo(*alleles)
        N = len(alleles)

        ### BPH ###
        (((A0, A1),((B0,),)),) = models_dict[N]['BPH'][1].items()
        BPH = A0/A1 * F[B0] / D

        BPH += sum(A0/A1 * sum(F[B0] * F[B1] for (B0, B1) in C)
                   for (A0, A1), C in models_dict[N]['BPH'][2].items()) / D**2

        if N>2:
            BPH += sum(A0/A1 * sum(F[B0] * sum(F[B1] * F[B2] for (B1, B2) in C[B0]) for B0 in C)
                       for (A0, A1), C in models_dict[N]['BPH'][3].items()) / D**3

        ### SPH ###
        (((A0, A1),((B0,),)),) = models_dict[N]['SPH'][1].items()
        SPH = A0/A1 * F[B0] / D

        SPH += sum(A0/A1 * sum(F[B0] * F[B1] for (B0, B1) in C)
                   for (A0, A1), C in models_dict[N]['SPH'][2].items()) / D**2

        result = 0.123456789 if SPH<1e-16 else log(BPH/SPH)
        return result

    return LLR

def wrapper_func_of_create_LLR(obs_tab,leg_tab,hap_tab,models_filename):
    """ Wraps the fuction create_LLR. It receives an observations array, legend
        array and haplotypes array. Based on the given data it creates and
        returns the function LLR."""

    if not os.path.isfile(models_filename): raise Exception('Error: MODELS file does not exist.')
    ###with open(models_filename, 'rb') as f:
    ###with bz2.BZ2File(models_filename, 'rb') as f:
    load_model = bz2.BZ2File if models_filename.split('.')[-1]=='pbz2' else open
    with load_model(models_filename, 'rb') as f:
        models_dict = pickle.load(f)

    LLR = create_LLR(models_dict,create_frequencies(build_hap_dict(obs_tab, leg_tab, hap_tab)),len(hap_tab[0]))
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

    ###with open(models_filename, 'rb') as f:
    ###with bz2.BZ2File(models_filename, 'rb') as f:
    load_model = bz2.BZ2File if models_filename.split('.')[-1]=='pbz2' else open
    with load_model(models_filename, 'rb') as f:
        models_dict = pickle.load(f)

    hap_dict = build_hap_dict(obs_tab, leg_tab, hap_tab)
    joint_frequencies_combo = create_frequencies(hap_dict)
    LLR = create_LLR(models_dict,joint_frequencies_combo,len(hap_tab[0]))

    ###LLR = create_LLR(models_dict,create_frequencies(build_hap_dict(obs_tab, leg_tab, hap_tab))) #This line replaces the three lines above.
    return LLR


"""
if __name__ != "__main__":
    print('The module LLR_CALCULATOR was imported.')
else:
    print('The module LLR_CALCULATOR was invoked directly')
    sys.exit(0)

###############################   END OF FILE   ###############################

"""

if __name__ != "__main__":
    print("The module LLR_CALCULATOR was imported.")
else:
    print("Executed when invoked directly")
    #sys.exit(0)
    import time
    from MAKE_OBS_TAB import read_impute2
    a = time.time()
    obs_filename = 'results_HapMix_EXT/mixed2haploids.X0.01.SRR10393062.SRR151495.0-2.hg38.obs.p'
    hap_filename = '../build_reference_panel/ref_panel.HapMix_EXT.hg38.BCFtools/chr21_HapMix_EXT_panel.hap'
    leg_filename = '../build_reference_panel/ref_panel.HapMix_EXT.hg38.BCFtools/chr21_HapMix_EXT_panel.legend'

    hap_tab = read_impute2(hap_filename, filetype='hap')
    leg_tab = read_impute2(leg_filename, filetype='leg')


    with open(obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        #info = pickle.load(f)

    #with bz2.BZ2File('MODELS/MODELS16D.pbz2', 'rb') as f:
    with open('MODELS/MODELS16D.p', 'rb') as f:
        models_dict = pickle.load(f)

    hap_dict = build_hap_dict(obs_tab, leg_tab, hap_tab)
    #aux_dict = build_aux_dict(obs_tab, leg_tab)

    positions = tuple(hap_dict.keys())

    #frequencies, frequency = create_frequencies(hap_dict)
    N = len(hap_tab[0])
    frequencies = create_frequencies(hap_dict,N)
    def frequencies2(*x):
        return {bin(a)[2:]:b for a,b in frequencies(*x).items()}

    LLR = create_LLR(models_dict,frequencies)

    pos = (positions[:4],positions[4:8],positions[8:12],positions[12:16])

    print(frequencies2(positions[0]))
    print(frequencies2(positions[:4]))
    print('-----')
    print(pos)
    print(frequencies2(*pos))
    print(LLR(*pos))
    print('-----')
    print(positions[:2])
    print(frequencies2(*positions[:2]))
    print(LLR(*positions[:2]))
    print('-----')
    print(positions[:3])
    print(frequencies2(*positions[:3]))
    print(LLR(*positions[:3]))
    print('-----')
    print(positions[:4])
    print(frequencies2(*positions[:4]))
    print(LLR(*positions[:4]))

    b = time.time()

    print('Done in %.3f sec.' % ((b-a)))
