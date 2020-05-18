#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
LLR_CALCULATION

This script creates a function that calculates log-likelihood BPH/SPH ratio
(LLR) for a tuple of SNPs. BPH (Both Parental Homologs) correspond to the
presence of three unmatched haplotypes, while SPH (Single Parental Homolog)
correspond to chromosome gains involving identical homologs.  

Daniel Ariad (daniel@ariad.org)
May 4th, 2020
"""

import operator, itertools, pickle, math, functools, time

def prod(iterable):
    """ Calculates the product of all the elements in the input iterable. """
    return functools.reduce(operator.mul, iterable, 1)

def build_hap_dict(obs_tab,leg_tab,hap_tab):
    """ Returns a dictionary that lists chromosome positions of SNPs and gives
        their relevent row from haplotypes table. The row is stored as a tuple
        of booleans, where True represents the observed allele."""
        
    hap_dict = dict()
        
    for (pos, ind, read_id, read) in obs_tab:
        chr_id, pos2, ref, alt = leg_tab[ind]
        if pos!=pos2:
            raise Exception('error: the line numbers in obs_tab refer to the wrong chromosome positions in leg_tab.')         
        if read==alt:
            hap_dict[pos] = tuple(hap_tab[ind]) 
        elif read==ref:
            hap_dict[pos] = tuple(map(operator.not_,hap_tab[ind]))        
    
    print('%%%.2f of the reads matched known alleles.' % (100*len(hap_dict)/len(obs_tab)))
    
    return hap_dict

def read_impute2(impute2_filename,**kwargs):
    """ Reads an IMPUTE2 file format (SAMPLE/LEGEND/HAPLOTYPE) and builds a list
        of lists, containing the dataset. """
    
    filetype = kwargs.get('filetype', None)
    with open(impute2_filename, 'r') as impute2_in:
        if filetype == 'legend':
            impute2_in.readline()   # Bite off the header
            def parse(x): 
                y=x.strip().split()
                y[0] = 'chr'+y[0].split(':')[0]
                y[1]=int(y[1])
                return y
        elif filetype == 'hap':
            def parse(x):
                return [i=='1' for i in x.strip().split()]
        else:
            def parse(x):
                return x.strip().split() # Trailing whitespaces stripped from the ends of the string. Then it splits the string into a list of words.
       
        impute2_tab = [parse(line) for line in impute2_in]
    return impute2_tab


def create_frequencies(hap_dict):
    """ Returns the function frequencies, which extracts from the dictionary
        hap_dict the joint frequencies of observed alleles. """
    
    N = len(next(iter(hap_dict.values())))

    def frequencies(*positions):
        """ Based on the reference panel, it calculates joint frequencies of
            observed alleles. In addition, it caches the joint frequencies in a
            dictionary for future requests. """
                    
        hap = dict()
        
        for i, pos in enumerate(sorted(positions)):
            if type(pos)==tuple:
                n = len(pos)
                if n==1: 
                    hap[i] = hap_dict[pos[0]] 
                elif n==2:
                    hap[i] = tuple(map(operator.and_,pos[0],pos[1]))
                else:
                    hap[i] = tuple(map(all,zip(*pos)))
            elif type(pos)==int:
                hap[i] = hap_dict[pos]
            else:
                raise Exception('error: the function, frequencies only receives integers and tuples of integers.')
                   
        result = {(i,): A.count(True) / N  for i,A in hap.items() }
        
        for i in itertools.combinations(hap.items(), 2):
            a,b = zip(*i)
            result[a] = sum(itertools.compress(*b)) / N 
        
        for i in itertools.combinations(hap.items(), 3):
            a,b = zip(*i)
            result[a] = sum(itertools.compress(b[0], map(operator.and_,b[1],b[2]))) / N
        
        for r in range(4,len(positions)+1):
            for i in itertools.combinations(hap.items(), r):
                a,b = zip(*i)
                result[a] = sum(itertools.compress(b[0], map(all,zip(*b[1:])))) / N
        
        return result
    
    return frequencies

def create_LLR(models_dict,frequencies):
    """ This function receives the dictionary models_dict with the
    statisitcal models and the function frequncies, which calculates
    joint frequncies. Based on these arguments it creates the function
    LLR, which calculates the log-likelihood BPH/SPH ratio."""
    
    def LLR(*positions):
        """ Calculates the log-likelihood BPH/SPH ratio for a tuple of SNPs.
        BPH (Both Parental Homologs) denotes three unmatched haplotypes,
        while SPH (Single Parental Homolog) denotes two matched haplotypes
        out of three. """
        l = len(positions)
        freq = frequencies(*positions)
        BPH = (A[0]/A[1]*prod((freq[b] for b in B)) for B,A in models_dict[l]['BPH'] )
        SPH = (A[0]/A[1]*prod((freq[b] for b in B)) for B,A in models_dict[l]['SPH'] )
        return math.log(sum(BPH)/sum(SPH))
    
    return LLR

def main(obs_filename,leg_filename,hap_filename,models_filename):
    """ This function receives a table of observation, IMPUTE2 legend file,
        IMPUTE2 haplotypes file, and the statistical model. Based on the given
        data it creates and returns the LLR function."""
        
    hap_tab = read_impute2(leg_filename, filetype='hap')
    leg_tab = read_impute2(hap_filename, filetype='legend')
    with open(obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        #info = pickle.load(f)
    with open(models_filename, 'rb') as f:
        models_dict = pickle.load(f)
    
    #hap_dict = build_hap_dict(obs_tab, leg_tab, hap_tab)
    #frequencies = create_frequencies(hap_dict)
    #LLR = create_LLR(models_dict,frequencies)     
    
    LLR = create_LLR(models_dict,create_frequencies(build_hap_dict(obs_tab, leg_tab, hap_tab))) #This line replaces the three lines above.
    return LLR
    
if __name__ == "__main__": 
    print("Executed when invoked directly")
    a = time.time()

    hap_tab = read_impute2('../make_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.hap', filetype='hap')
    leg_tab = read_impute2('../make_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.legend', filetype='legend')
    

    with open('results/mixed2haploids.X0.1.SRR10393062.SRR151495.0-2.hg38.OBS.p', 'rb') as f:
        obs_tab = pickle.load(f)
        #info = pickle.load(f)
    
    with open('MODELS.p', 'rb') as f:
        models_dict = pickle.load(f)
    
    hap_dict = build_hap_dict(obs_tab, leg_tab, hap_tab)
    
    positions = tuple(hap_dict.keys())
    
    frequencies = create_frequencies(hap_dict)
    
    LLR = create_LLR(models_dict,frequencies) 
    
    #frequencies(*positions[:16])
    #print(positions[:2])
    #print(frequencies(*positions[:2]))
    #print(LLR(*positions[:2]))
    #print(positions[:3])
    #print(frequencies(*positions[:3]))
    #print(LLR(*positions[:3]))
    #print(positions[:4])
    #print(frequencies(*positions[:4]))
    #print(LLR(*positions[:4]))
    b = time.time()

    print('Done in %.3f sec.' % ((b-a)))
    
else: 
    print("Executed when imported")