#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DISTANT_ADMIXTURE_MODELS

Given reads that originated form the same genomic window and a reference panel
of two populations, the likelihood of observed reads under four scenarios,
namely, monosomy, disomy, SPH and BPH is calculated. This module is for distant
admixtures, where each descendant haplotype has a certain probability to
originate from one of two ancestral populations.

BPH (Both Parental Homologs) correspond to the presence of three unmatched
haplotypes, while SPH (Single Parental Homolog) correspond to chromosome gains
involving identical homologs.

Daniel Ariad (daniel@ariad.org)
Aug 10, 2021
"""

import pickle, os, sys, bz2, collections, gzip

from functools import reduce
from operator import and_, itemgetter
from itertools import combinations

leg_tuple = collections.namedtuple('leg_tuple', ('chr_id', 'pos', 'ref', 'alt')) #Encodes the rows of the legend table
sam_tuple = collections.namedtuple('sam_tuple', ('sample_id', 'group1', 'group2', 'sex')) #Encodes the rows of the samples table
obs_tuple = collections.namedtuple('obs_tuple', ('pos', 'read_id', 'base')) #Encodes the rows of the observations table
spanel_tuple = collections.namedtuple('spanel_tuple', ('group2', 'flag', 'hap_number', 'proportion')) #Encodes the rows of the observations table
try:
    from gmpy2 import popcount
except ModuleNotFoundError:
    print('caution: the module gmpy2 is missing.')
    def popcount(x):
        """ Counts non-zero bits in positive integer. """
        return bin(x).count('1')

class distant_admixture:
    """ Based on the statisitcal models (models_dict) and the reference panel
    (leg_tab, hap_tab and sam_tab), it allows to calculate the likelihoods of
    observed alleles under various statistical models (monosomy, disomy, SPH
    and BPH). """

    def __init__(self, obs_tab, leg_tab, hap_tab, sam_tab, models_dict, total_number_of_haplotypes, ancestral_makeup):
        """ Initialize the attributes of the class. """

        if len(leg_tab)!=len(hap_tab):
            raise Exception('Error: the number of SNPs in the LEGEND file differ from the number of SNPs in the HAP file.')

        if total_number_of_haplotypes!=2*len(sam_tab):
            raise Exception('Error: the number of diploid samples in the SAMPLE file differ from the number of haplotypes in the HAP file.')
        
        if sum(ancestral_makeup.values())<0.999:
            raise Exception('Error: ancestry proportions must sum to one.')
            
        self.total_number_of_haplotypes_in_reference_panel = total_number_of_haplotypes
        self.sub_panels = self.extract_subpanels(sam_tab)
        self.ancestral_makeup = ancestral_makeup #Keys correspond to group2 names and values to haplotype proportions.
        
        for group2 in ancestral_makeup:
            if group2 not in self.sub_panels: 
                raise Exception('Error: %s is missing from the reference panel.' % group2)
        
        self.proportions = tuple(spanel_tuple(group2, *self.sub_panels[group2], proportion) for group2, proportion in ancestral_makeup.items())
        self.models_dict = models_dict
        self.hap_dict, self.fraction_of_matches = self.build_hap_dict(obs_tab, leg_tab, hap_tab)

    def extract_subpanels(self, sam_tab):
        """ Differentiates between the two groups that compose the reference
        panel. Then, all the haplotypes that are associated with each group are
        flagged using a binary representation marks and counted. """

        groups = list({j.group2 for j in sam_tab})
        flags = [] #occupation representation of samples that are asscoiated with a group.
        N = [] #number of haplotypes in reference subpanel.
        for group2 in groups[:-1]:
            differentiate = [row.group2 == group2 for row in sam_tab for i in (1,2)]
            N.append(differentiate.count(True))
            flag = sum(v<<i for i, v in enumerate(reversed(differentiate)))
            flags.append(flag)
        
        last_N = self.total_number_of_haplotypes_in_reference_panel - sum(N)
        N.append(last_N)
        last_flag = sum(flags) ^ ((1 << self.total_number_of_haplotypes_in_reference_panel) - 1)
        flags.append(last_flag)
        
        result = {group2: (flag, hap_number) for group2, flag, hap_number in zip(groups,flags,N)}
        return result

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

        print('Algorithm for distant admixtures: %.2f%% of the observed alleles matched the reference panel.' % (100*fraction_of_matches))

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

    def effective_joint_frequency(self, hap):
        """ The probability to draw a haplotype from one of the groups2. """
        return sum(i.proportion * popcount(hap & i.flag) / i.hap_number for i in self.proportions)

    def joint_frequencies_combo(self, *alleles):
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

        hap = self.intrenal_hap_dict(*alleles)

        result = {c: self.effective_joint_frequency(A) for c,A in hap.items()}

        for C in combinations(hap, 2):
            result[C[0]|C[1]] = self.effective_joint_frequency(hap[C[0]]&hap[C[1]])

        for C in combinations(hap, 3):
            result[C[0]|C[1]|C[2]] = self.effective_joint_frequency(hap[C[0]]&hap[C[1]]&hap[C[2]])

        for r in range(4,len(alleles)):
            for C in combinations(hap, r):
                result[sum(C)] = self.effective_joint_frequency(reduce(and_,itemgetter(*C)(hap)))

        if len(alleles)>=4:
            result[sum(hap.keys())] = self.effective_joint_frequency(reduce(and_,hap.values()))
        return result

    def likelihoods(self, *alleles):
        """ Calculates the likelihood to observe a set with alleles
        and haplotypes under four scenarios, namely, monosomy, disomy, SPH
        and BPH. """

        models = self.models_dict[len(alleles)]
        
        F = self.joint_frequencies_combo(*alleles)

        ### BPH ###
        (((A0, A1),((B0,),)),) = models['BPH'][1].items()

        BPH = (A0 / A1)  * F[B0]


        BPH += sum( sum(F[B0] * F[B1] for (B0, B1) in C) * A0 / A1
                   for (A0, A1), C in models['BPH'][2].items())

        if len(alleles)>2:
            BPH += sum( sum(F[B0] * sum( F[B1] * F[B2] for (B1, B2) in C[B0]) for B0 in C) * A0 / A1
                     for (A0, A1), C in models['BPH'][3].items())

        ### SPH ###
        (((A0, A1),((B0,),)),) = models['SPH'][1].items()
        SPH = (A0 / A1) * F[B0]

        SPH += sum( sum(F[B0] * F[B1] for (B0, B1) in C) * A0 / A1
                   for (A0, A1), C in models['SPH'][2].items())

        ### DIPLOIDY ###
        (((A0, A1),((B0,),)),) = models['DISOMY'][1].items()
        DISOMY = (A0 / A1)  * F[B0]

        DISOMY += sum( sum( F[B0] * F[B1] for (B0, B1) in C) * A0 / A1
                   for (A0, A1), C in models['DISOMY'][2].items())

        ### MONOSOMY ###
        ((B0,),) = models['MONOSOMY'][1][(1,1)]
        MONOSOMY = F[B0]

        result = (MONOSOMY, DISOMY, SPH, BPH)
        return result

    def likelihoods2(self, *alleles):
        """ Calculates the likelihood to observe two alleles/haplotypes
        under four scenarios, namely, monosomy, disomy, SPH and BPH. """

        F = self.joint_frequencies_combo(*alleles)
        a, b, ab = F[1] , F[2], F[3] #The probability of a haplotype to be associated with group2_id=0 is m.

        BPH = (ab+2*a*b)/3 #The likelihood of three unmatched haplotypes.
        SPH = (5*ab+4*a*b)/9 #The likelihood of two identical haplotypes out three.
        DISOMY = (ab+a*b)/2 #The likelihood of diploidy.
        MONOSOMY = ab #The likelihood of monosomy.

        return MONOSOMY, DISOMY, SPH, BPH

    def likelihoods3(self, *alleles):
        """ Calculates the likelihood to observe three alleles/haplotypes
        under four scenarios, namely, monosomy, disomy, SPH and BPH. """

        F = self.joint_frequencies_combo(*alleles)
        a, b, ab, c, ac, bc, abc = (F[i] for i in range(1,8))

        BPH = (abc+2*(ab*c+ac*b+bc*a+a*b*c))/9 #The likelihood of three unmatched haplotypes.
        SPH = abc/3+2*(ab*c+ac*b+bc*a)/9  #The likelihood of two identical haplotypes out three.
        DISOMY = (abc+ab*c+ac*b+bc*a)/4 #The likelihood of diploidy.
        MONOSOMY = abc #The likelihood of monosomy.

        return MONOSOMY, DISOMY, SPH, BPH

    def likelihoods4(self, *alleles):
        """ Calculates the likelihood to observe four alleles/haplotypes
        under four scenarios, namely, monosomy, disomy, SPH and BPH. """

        F = self.joint_frequencies_combo(*alleles)
        a, b, ab, c, ac, bc, abc, d, ad, bd, abd, cd, acd, bcd, abcd = (F[i] for i in range(1,16))

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

def wrapper_of_distant_admixture_for_debugging(obs_filename,leg_filename,hap_filename,sample_filename,models_filename,admixture):
    """ Wrapper function of the class 'distant_admixture'. It receives an
    observations file, legend file, haplotypes file, samples file and a file
    with the statistical models. Based on the given data it creates and returns
    an instance of the class. """

    if not os.path.isfile(obs_filename): raise Exception('Error: OBS file does not exist.')
    if not os.path.isfile(leg_filename): raise Exception('Error: LEGEND file does not exist.')
    if not os.path.isfile(hap_filename): raise Exception('Error: HAP file does not exist.')
    if not os.path.isfile(sample_filename): raise Exception('Error: SAMPLE file does not exist.')
    if not os.path.isfile(models_filename): raise Exception('Error: MODELS file does not exist.')

    load = lambda filename: {'bz2': bz2.open, 'gz': gzip.open}.get(filename.rsplit('.',1)[1], open)  #Adjusts the opening method according to the file extension.

    open_hap = load(hap_filename)
    with open_hap(hap_filename,'rb') as hap_in:
        hap_tab, total_number_of_haplotypes = pickle.load(hap_in)

    open_leg = load(leg_filename)
    with open_leg(leg_filename,'rb') as leg_in:
        leg_tab = pickle.load(leg_in)

    open_samp = load(sample_filename)
    with open_samp(sample_filename,'rb') as samp_in:
        sam_tab = pickle.load(samp_in)

    open_obs = load(obs_filename)
    with open_obs(obs_filename, 'rb') as obs_in:
        obs_tab = pickle.load(obs_in)
        #info = pickle.load(f)

    open_model = load(models_filename)
    with open_model(models_filename, 'rb') as model_in:
        models_dict = pickle.load(model_in)

    ###return obs_tab, leg_tab, hap_tab, sam_tab, models_dict, total_number_of_haplotypes, admixture
    return distant_admixture(obs_tab, leg_tab, hap_tab, sam_tab, models_dict, total_number_of_haplotypes, admixture)

if __name__ != "__main__":
    print('The module DISTANT_ADMIXTURE_MODELS was imported.')
else:
    print('The module DISTANT_ADMIXTURE_MODELS was invoked directly')
    sys.exit(0)

###############################   END OF FILE   ###############################
"""


if __name__ != "__main__":
    print("The module DISTANT_ADMIXTURE_MODELS was imported.")
else:
    print('The module DISTANT_ADMIXTURE_MODELS was invoked directly')
    #sys.exit(0)
    import time, random
    t0 = time.time()
    obs_filename = 'test/test.obs.p'
    hap_filename = 'test/test.hap.p'
    leg_filename = 'test/test.leg.p'
    sam_filename = 'test/test.sam.p'
    models_filename = 'MODELS/MODELS16.p'
    admixture = {'group0': 1}

    A = wrapper_of_distant_admixture_for_debugging(obs_filename,leg_filename,hap_filename,sam_filename,models_filename,admixture)

    alleles = tuple(A.hap_dict.keys())

    frequencies = lambda *x: {bin(a)[2:]:b for a,b in A.joint_frequencies_combo(*x).items()}

    likelihoods = A.likelihoods
    likelihoods2 = A.likelihoods2
    likelihoods3 = A.likelihoods3
    likelihoods4 = A.likelihoods4

    random.seed(a=2021, version=2)
    x = random.randrange(len(alleles)-16)
    haplotypes = (alleles[x:x+4],alleles[x+4:x+8],alleles[x+8:x+12],alleles[x+12:x+16])

    print('-----joint_frequencies_combo-----')
    print(frequencies(alleles[x+0]))
    print(frequencies(*alleles[x:x+4]))
 
    print('-----likelihoods4-haplotypes-----')
    print(haplotypes)
    print(frequencies(*haplotypes))
    print(likelihoods(*haplotypes))
    print(likelihoods4(*haplotypes))
    print('-----likelihoods2-----')
    print(alleles[x:x+2])
    print(frequencies(*alleles[x:x+2]))
    print(likelihoods(*alleles[x:x+2]))
    print(likelihoods2(*alleles[x:x+2]))
    print('-----likelihoods3-----')
    print(alleles[x:x+3])
    print(frequencies(*alleles[x:x+3]))
    print(likelihoods(*alleles[x:x+3]))
    print(likelihoods3(*alleles[x:x+3]))
    print('-----likelihoods4-----')
    print(alleles[x:x+4])
    print(frequencies(*alleles[x:x+4]))
    print(likelihoods(*alleles[x:x+4]))
    print(likelihoods4(*alleles[x:x+4]))

    t1 = time.time()

    print('Done in %.3f sec.' % ((t1-t0)))
"""
