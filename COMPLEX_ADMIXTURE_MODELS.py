#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
LIKELIHOODS_CALCULATOR_GENERALIZED

Given reads that originated form the same genomic window and a reference panel
of two populations, the likelihood of observed reads under four scenarios,
namely, monosomy, disomy, SPH and BPH is calculated. This module is for complex
admixtures, where each descendant haplotype has a certain probability to
originate from one of two ancestral populations.

BPH (Both Parental Homologs) correspond to the presence of three unmatched
haplotypes, while SPH (Single Parental Homolog) correspond to chromosome gains
involving identical homologs.

Daniel Ariad (daniel@ariad.org)
Aug 10, 2021
"""

import pickle, os, sys, bz2, collections

from functools import reduce
from operator import and_, itemgetter
from itertools import combinations

leg_tuple = collections.namedtuple('leg_tuple', ('chr_id', 'pos', 'ref', 'alt')) #Encodes the rows of the legend table
sam_tuple = collections.namedtuple('sam_tuple', ('sample_id', 'group1', 'group2', 'sex')) #Encodes the rows of the samples table
obs_tuple = collections.namedtuple('obs_tuple', ('pos', 'read_id', 'base')) #Encodes the rows of the observations table
admixture_tuple = collections.namedtuple('admixture_tuple', ('group2', 'proportion')) #Encodes the assumed admixture

try:
    from gmpy2 import popcount
except ModuleNotFoundError:
    print('caution: the module gmpy2 is missing.')
    def popcount(x):
        """ Counts non-zero bits in positive integer. """
        return bin(x).count('1')

class examine_generalized:
    """ Based on two IMPUTE2 arrays, which contain the legend and haplotypes,
    and a dictionary with statisitcal models (models_dict), it allows to
    calculate the likelihoods of observed alleles under various statistical
    models (monosomy, disomy, SPH and BPH). """
   
    
    def __init__(self, obs_tab, leg_tab, hap_tab, sam_tab, models_dict, total_number_of_haplotypes, admixture):
        """ Initialize the attributes of the class. """
        
        if len(leg_tab)!=len(hap_tab): 
            raise Exception('Error: the number of SNPs in the LEGEND file differ from the number of SNPs in the HAP file.')
            
        if total_number_of_haplotypes!=2*len(sam_tab): 
            raise Exception('Error: the number of diploid samples in the SAMPLE file differ from the number of haplotypes in the HAP file.')

        self.total_number_of_haplotypes_in_reference_panel = total_number_of_haplotypes
        self.models_dict = models_dict
        self.hap_dict, self.fraction_of_matches = self.build_hap_dict(obs_tab, leg_tab, hap_tab)
        self.flags, self.num_of_hap_in_ref_subpanel, self.name2id = self.subpanels(sam_tab)
        self.admixture = admixture #The first element in the tuple is the group2 name and second is the haplotype proportion.      
        self.proportions = {self.name2id[admixture.group2]: admixture.proportion, 1-self.name2id[admixture.group2]: 1-admixture.proportion}
        ##########print(str({ancestry: admixture.proportion if ancestry==admixture.group2 else 1-admixture.proportion for ancestry in self.name2id})[1:-1])
        
    def subpanels(self, sam_tab):
        """ Differentiates between the two groups that compose the reference
        panel. Then, all the haplotypes that are associated with each group are
        flagged using a binary representation marks and counted. """
        
        differentiate = [row.group2 == sam_tab[0].group2 for row in sam_tab for i in (1,2)]
        flag0 = sum(v<<i for i, v in enumerate(differentiate[::-1]))
        flag1 = flag0 ^ ((1 << self.total_number_of_haplotypes_in_reference_panel) - 1)
        flags = (flag0, flag1)       
        name2id = {sam_tab[0].group2:0,
                   sam_tab[differentiate[::2].index(False,1)].group2: 1}
        N0 = differentiate.count(True)
        N1 = self.total_number_of_haplotypes_in_reference_panel - N0
        number_of_haplotypes_in_reference_subpanel = (N0,N1)
        return flags, number_of_haplotypes_in_reference_subpanel, name2id
    
    def build_hap_dict(self, obs_tab, leg_tab, hap_tab):
        """ Returns a dictionary that lists SNP alleles and gives their
        relevent row from haplotypes table. The row is stored as bits, where
        True means that the haplotype contains the allele. We denote the
        returned dictionary as the reference panel. """
    
        hap_dict = dict()
        mismatches = 0
        combined = {pos: (ref,alt,hap) for (chr_id,pos,ref,alt),hap in zip(leg_tab, hap_tab)}
        missing = 3*(None,)

        b = (1 << self.total_number_of_haplotypes_in_reference_panel) - 1 #### equivalent to int('1'*number_of_haplotypes,2)
        
        for (pos, read_id, base) in obs_tab:
            ref, alt, hap = combined.get(pos, missing)
            if base==alt:
                hap_dict[(pos,base)] = hap
            elif base==ref:
                hap_dict[(pos,base)] = hap ^ b ### ^b flips all bits of the binary number, hap_tab[ind] using bitwise xor operator. 
            else:
                mismatches += 1
    
        fraction_of_matches = 1-mismatches/len(obs_tab)
        
        print('Generalized algorithm: %.2f%% of the observed alleles matched the reference panel.' % (100*fraction_of_matches))
    
        return hap_dict, fraction_of_matches

    def intrenal_hap_dict(self, *alleles):
        """ This function allows treatment of alleles and haplotypes on an
        equal footing. This is done in three steps: (1) All the alleles and
        haplotypes are enumerated. (2) For each given haplotype, tuples in the
        reference panel, associated with the haplotype's alleles, are
        intersected. (3) A dictionary that lists all the alleles and haplotypes
        by their index is returned. The dictionary gives for each allele and
        haplotype their associated tuple and intersected tuple, respectively. """
        hap = dict()

        for i, X in enumerate(alleles):
            if type(X[0])==tuple: #Checks if X is a tuple/list of alleles.
                n = len(X)
                if n==1:
                    hap[1 << i] = self.hap_dict[X[0]]
                elif n==2:
                    hap[1 << i] = self.hap_dict[X[0]] & self.hap_dict[X[1]]
                else:
                    hap[1 << i] = reduce(and_,itemgetter(*X)(self.hap_dict))

            elif type(X[0])==int: #Checks if X is a single allele.
                hap[1 << i] = self.hap_dict[X]
            else:
                raise Exception('error: joint_frequencies only accepts alleles and tuple/list of alleles.')

        return hap

    def joint_frequencies_combo(self, *alleles, group2_id, normalize):
        """ Based on the reference panel, it calculates joint frequencies of
            observed alleles. The function arguments are alleles, that is, 
            tuples of position and base, e.g., (100,'T'), (123, 'A') and 
            (386, 'C'). Each allele is enumerated according to the order it
            was received by the function. The function returns a dictionary that
            lists all the possible subgroups of the given alleles. Each key in
            the dictionary is a tuple of intergers that are lexicographically
            sorted. Moreover, each integer within the keys corresponds to an
            enumerated allele. For each subgroup of alleles the dictionary
            gives the joint frequencies in a given population. The function
            arguments can also include haplotypes, that is, tuples of alleles;
            Haplotypes are treated in the same manner as alleles. """

        flag = self.flags[group2_id]
        hap = self.intrenal_hap_dict(*alleles)

        result = {c: popcount(A&flag) for c,A in hap.items()}

        for C in combinations(hap, 2):
            result[C[0]|C[1]] = popcount(hap[C[0]]&hap[C[1]]&flag)

        for C in combinations(hap, 3):
            result[C[0]|C[1]|C[2]] = popcount(hap[C[0]]&hap[C[1]]&hap[C[2]]&flag)

        for r in range(4,len(alleles)):
            for C in combinations(hap, r):
                result[sum(C)] = popcount(reduce(and_,itemgetter(*C)(hap))&flag)

        if len(alleles)>=4:
            result[sum(hap.keys())] = popcount(reduce(and_,hap.values())&flag)

        if normalize:
            N = self.num_of_hap_in_ref_subpanel[group2_id]
            result = {k: v/N for k,v in result.items()}
        
        return result
    
    def likelihoods(self, *alleles):
        """ Calculates the likelihood to observe a set with alleles
        and haplotypes under four scenarios, namely, monosomy, disomy, SPH
        and BPH. """
        
        models = self.models_dict[len(alleles)]
        F = self.joint_frequencies_combo(*alleles, group2_id=0, normalize=True) 
        G = self.joint_frequencies_combo(*alleles, group2_id=1, normalize=True) 

        m = self.proportions[0] #The probability to draw a haplotype that is associated with group2_id=0.
        n = self.proportions[1] #The probability to draw a haplotype that is associated with group2_id=1.
        
        ### BPH ###
        (((A0, A1),((B0,),)),) = models['BPH'][1].items()
        
        BPH = (A0 / A1)  * (m * F[B0] + n * G[B0])  


        BPH += sum( sum((m * F[B0] + n * G[B0]) * (m * F[B1] + n * G[B1]) for (B0, B1) in C) * A0 / A1
                   for (A0, A1), C in models['BPH'][2].items()) 
      
        if len(alleles)>2:
            BPH += sum( sum((m * F[B0] + n * G[B0]) * sum( (m * F[B1] + n * G[B1]) * (m * F[B2] + n * G[B2])
                  for (B1, B2) in C[B0]) for B0 in C) * A0 / A1
                     for (A0, A1), C in models['BPH'][3].items()) 

        ### SPH ###
        (((A0, A1),((B0,),)),) = models['SPH'][1].items()
        SPH = (A0 / A1)  * (m * F[B0] + n * G[B0]) 

        SPH += sum( sum((m * F[B0] + n * G[B0]) * (m * F[B1] + n * G[B1]) for (B0, B1) in C) * A0 / A1
                   for (A0, A1), C in models['SPH'][2].items()) 

        
        ### DIPLOIDY ###
        (((A0, A1),((B0,),)),) = models['DISOMY'][1].items()
        DISOMY = (A0 / A1)  * (m * F[B0] + n * G[B0]) 

        DISOMY += sum( sum( (m * F[B0] + n * G[B0]) * (m * F[B1] + n * G[B1]) for (B0, B1) in C) * A0 / A1
                   for (A0, A1), C in models['DISOMY'][2].items()) 

        ### MONOSOMY ###
        ((B0,),) = models['MONOSOMY'][1][(1,1)]
        MONOSOMY = m * F[B0] + n * G[B0]
        
        result = (MONOSOMY, DISOMY, SPH, BPH)
        return result
    
    def likelihoods2(self, *alleles):
        """ Calculates the likelihood to observe two alleles/haplotypes
        under four scenarios, namely, monosomy, disomy, SPH and BPH. """
        
        m = self.proportions[0] #The probability to draw a haplotype that is associated with group2_id=0.
        n = self.proportions[1] #The probability to draw a haplotype that is associated with group2_id=1.
        
        F = self.joint_frequencies_combo(*alleles, group2_id=0, normalize=True) 
        G = self.joint_frequencies_combo(*alleles, group2_id=1, normalize=True) 
        a, b, ab = m*F[1]+n*G[1] , m*F[2]+n*G[2], m*F[3]+n*G[3] #The probability of a haplotype to be associated with group2_id=0 is m.
        
        BPH = (ab+2*a*b)/3 #The likelihood of three unmatched haplotypes.
        SPH = (5*ab+4*a*b)/9 #The likelihood of two identical haplotypes out three.
        DISOMY = (ab+a*b)/2 #The likelihood of diploidy.
        MONOSOMY = ab #The likelihood of monosomy.
        
        return MONOSOMY, DISOMY, SPH, BPH

    def likelihoods3(self, *alleles):
        """ Calculates the likelihood to observe three alleles/haplotypes
        under four scenarios, namely, monosomy, disomy, SPH and BPH. """
        
        m = self.proportions[0] #The probability to draw a haplotype that is associated with group2_id=0.
        n = self.proportions[1] #The probability to draw a haplotype that is associated with group2_id=1.
        
        F = self.joint_frequencies_combo(*alleles, group2_id=0, normalize=True) 
        G = self.joint_frequencies_combo(*alleles, group2_id=1, normalize=True) 
        a, b, ab, c, ac, bc, abc = (m*F[i]+n*G[i] for i in range(1,8))        
        
        BPH = (abc+2*(ab*c+ac*b+bc*a+a*b*c))/9 #The likelihood of three unmatched haplotypes.
        SPH = abc/3+2*(ab*c+ac*b+bc*a)/9  #The likelihood of two identical haplotypes out three.
        DISOMY = (abc+ab*c+ac*b+bc*a)/4 #The likelihood of diploidy.
        MONOSOMY = abc #The likelihood of monosomy.

        return MONOSOMY, DISOMY, SPH, BPH
    
    def likelihoods4(self, *alleles):
        """ Calculates the likelihood to observe four alleles/haplotypes
        under four scenarios, namely, monosomy, disomy, SPH and BPH. """
        
        m = self.proportions[0] #The probability to draw a haplotype that is associated with group2_id=0.
        n = self.proportions[1] #The probability to draw a haplotype that is associated with group2_id=1.
        
        F = self.joint_frequencies_combo(*alleles, group2_id=0, normalize=True) 
        G = self.joint_frequencies_combo(*alleles, group2_id=1, normalize=True) 
        a, b, ab, c, ac, bc, abc, d, ad, bd, abd, cd, acd, bcd, abcd = (m*F[i]+n*G[i] for i in range(1,16))
        
        BPH = (abcd+2*(ab*c*d+a*bd*c+a*bc*d+ac*b*d+a*b*cd+ad*b*c+abc*d+a*bcd+acd*b+abd*c+ab*cd+ad*bc+ac*bd))/27  #The likelihood of three unmatched haplotypes.
        SPH = (17*abcd+10*(abc*d+bcd*a+acd*b+abd*c)+8*(ab*cd+ad*bc+ac*bd))/81  #The likelihood of two identical haplotypes out three.
        DISOMY = (abcd+abc*d+bcd*a+acd*b+abd*c+ab*cd+ad*bc+ac*bd)/8 #The likelihood of diploidy.
        MONOSOMY = abcd #The likelihood of monosomy.

        return MONOSOMY, DISOMY, SPH, BPH
    
    def get_likelihoods(self, *x):
        """ Uses the optimal function to calculate the likelihoods. 
        In general, self.likelihoods can get less than five alleles but the
        dedicated functions are optimized to a certain number of alleles. """
        
        l = len(x)
        if l==2:
            result = self.likelihoods2(*x)
        elif l==3:
            result = self.likelihoods3(*x)
        elif l==4:
            result = self.likelihoods4(*x)
        else: 
            result = self.likelihoods(*x)
        return result
        
def wrapper_of_examine_for_debugging(obs_filename,leg_filename,hap_filename,sample_filename,models_filename,admixture):
    """ Wrapper function of the class examine_generalized. It receives an observations
    file, IMPUTE2 legend file, IMPUTE2 haplotypes file, IMPUTE2 samples file,
    and a file with four statistical models. Based on the given data it creates
    and returns an instance of the class. """

    from MAKE_OBS_TAB import read_impute2

    if not os.path.isfile(obs_filename): raise Exception('Error: OBS file does not exist.')
    if not os.path.isfile(leg_filename): raise Exception('Error: LEGEND file does not exist.')
    if not os.path.isfile(hap_filename): raise Exception('Error: HAP file does not exist.')
    if not os.path.isfile(sample_filename): raise Exception('Error: SAMPLE file does not exist.')
    if not os.path.isfile(models_filename): raise Exception('Error: MODELS file does not exist.')

    leg_tab = read_impute2(leg_filename, filetype='leg')
    hap_tab, total_number_of_haplotypes = read_impute2(hap_filename, filetype='hap')
    sam_tab  = read_impute2(sample_filename, filetype='sam')
    
    load_obs = bz2.BZ2File if obs_filename[-6:]=='.p.bz2' else open
    with load_obs(obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        #info = pickle.load(f)

    load_model = bz2.BZ2File if models_filename[-6:]=='.p.bz2' else open
    with load_model(models_filename, 'rb') as f:
        models_dict = pickle.load(f)

    return examine_generalized(obs_tab, leg_tab, hap_tab, sam_tab, models_dict, total_number_of_haplotypes, admixture)

if __name__ != "__main__":
    print('The module LIKELIHOODS_CALCULATOR_GENERALIZED was imported.')
else:
    print('The module LIKELIHOODS_CALCULATOR_GENERALIZED was invoked directly')
    sys.exit(0)

###############################   END OF FILE   ###############################
"""


if __name__ != "__main__":
    print("The module LIKELIHOODS_CALCULATOR was imported.")
else:
    print("Executed when invoked directly")
    #sys.exit(0)
    import time, random
    t0 = time.time()
    obs_filename = 'results/SWI-L-10-27-May-2020_S38.chr6.obs.p.bz2'
    hap_filename = '../build_reference_panel/EAS_EUR_panel.hg38.BCFtools/chr6_EAS_EUR_panel.hap.gz'
    leg_filename = '../build_reference_panel/EAS_EUR_panel.hg38.BCFtools/chr6_EAS_EUR_panel.legend.gz'
    sam_filename = '../build_reference_panel/samples_per_panel/EAS_EUR_panel.samples'
    models_filename = 'MODELS/MODELS16.p'
    admixture = admixture_tuple('EUR',0.8)
    
    A = wrapper_of_examine_for_debugging(obs_filename,leg_filename,hap_filename,sam_filename,models_filename,admixture)

    alleles = tuple(A.hap_dict.keys())

    frequencies0 = lambda *x: {bin(a)[2:]:b for a,b in A.joint_frequencies_combo(*x,group2_id=0,normalize=True).items()}
    frequencies1 = lambda *x: {bin(a)[2:]:b for a,b in A.joint_frequencies_combo(*x,group2_id=1,normalize=True).items()}

    likelihoods = A.likelihoods
    likelihoods2 = A.likelihoods2
    likelihoods3 = A.likelihoods3
    likelihoods4 = A.likelihoods4

    random.seed(a=2021, version=2)
    x = random.randrange(len(alleles)-16)
    haplotypes = (alleles[x:x+4],alleles[x+4:x+8],alleles[x+8:x+12],alleles[x+12:x+16])

    print('-----joint_frequencies_combo-----')    
    print(frequencies0(alleles[x+0]))
    print(frequencies0(*alleles[x:x+4]))
    print(frequencies1(alleles[x+0]))
    print(frequencies1(*alleles[x:x+4]))
    print('-----likelihoods4-haplotypes-----')
    print(haplotypes)
    print(frequencies0(*haplotypes))
    print(frequencies1(*haplotypes))
    print(likelihoods(*haplotypes))
    print(likelihoods4(*haplotypes))
    print('-----likelihoods2-----')
    print(alleles[x:x+2])
    print(frequencies0(*alleles[x:x+2]))
    print(frequencies1(*alleles[x:x+2]))
    print(likelihoods(*alleles[x:x+2]))
    print(likelihoods2(*alleles[x:x+2]))
    print('-----likelihoods3-----')
    print(alleles[x:x+3])
    print(frequencies0(*alleles[x:x+3]))
    print(frequencies1(*alleles[x:x+3]))
    print(likelihoods(*alleles[x:x+3]))
    print(likelihoods3(*alleles[x:x+3]))
    print('-----likelihoods4-----')
    print(alleles[x:x+4])
    print(frequencies0(*alleles[x:x+4]))
    print(frequencies1(*alleles[x:x+4]))
    print(likelihoods(*alleles[x:x+4]))
    print(likelihoods4(*alleles[x:x+4]))

    t1 = time.time()

    print('Done in %.3f sec.' % ((t1-t0)))
"""
