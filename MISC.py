#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 18:51:19 2020

@author: ariad
"""
import collections, operator, os, pickle, itertools
#def find_overlaps(obs_filename, leg_filename):
#    """ Given an observation table, the fraction of overlapping reads
#        is calculated; Altough the reads overlap in position, they do not
#        necessary have common alleles. """
#    
#    leg_tab = read_impute2(leg_filename, filetype='leg')
#    obs_tab = pickle.load(open(obs_filename, 'rb'))
#    
#    aux_dict = collections.defaultdict(list)
#    for (pos, ind, read_id, base) in obs_tab:
#        if base in leg_tab[ind][2:]:
#            aux_dict[pos].append(read_id)   
#            
#    X = [x for x in aux_dict.values() if len(x)>1]
#
#    for i in range(len(X)-1,-1,-1):
#       for j in range(i-1,-1,-1):
#            if not X[i].isdisjoint(X[j]):
#                X[j].update(X.pop(i))
#                break
#    
#    OVERLAPS = sum(len(x) for x in X if len(x)>1)
#    TOTAL = len(set(itertools.chain.from_iterable(aux_dict.values())))
#    return OVERLAPS,TOTAL
#"""    
def load_llr(filename):
    with open('results_HapMix_EXT/'+filename, 'rb') as f:
        llr = pickle.load(f)
        info = pickle.load(f)
    print('%f\t%f\t%f\t%f\t%f\t%f\t\t\t%f\t%f' % (info['depth'],info['statistics']['jk_mean'],info['statistics']['jk_std'],info['statistics']['jk_bias'],info['statistics']['mean'],info['statistics']['std'],info['statistics']['num_of_LD_blocks'],info['statistics']['fraction_of_negative_LLRs']))
    return 0

def load(obs_filename,leg_filename,hap_filename):
    """ Wraps the function create_LLR. It receives an observations file, IMPUTE2
        legend file, IMPUTE2 haplotypes file, and the statistical model. Based
        on the given data it creates and returns the LLR function."""
    
    from MAKE_OBS_TAB import read_impute2
    
    if not os.path.isfile(obs_filename): raise Exception('Error: OBS file does not exist.')
    if not os.path.isfile(leg_filename): raise Exception('Error: LEGEND file does not exist.')
    if not os.path.isfile(hap_filename): raise Exception('Error: HAP file does not exist.')

    leg_tab = read_impute2(leg_filename, filetype='leg')
    hap_tab = read_impute2(hap_filename, filetype='hap')
    with open(obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        #info = pickle.load(f)
   
    return leg_tab,hap_tab,obs_tab

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

def build_aux_dict(obs_tab,leg_tab):
    """ Returns a dictionary that lists chromosome positions of SNPs and gives
        a dictionary that lists the observed bases at these postions. The
        nested dictionary gives all the read IDs that contain the observed base
        at the specific chromosome position. """
    
    nested = lambda: collections.defaultdict(list)
    aux_dict = collections.defaultdict(nested)

    for (pos, ind, read_id, base) in obs_tab:
        if base in leg_tab[ind][2:]:
            aux_dict[pos][base].append(read_id)  
       
    return aux_dict 

def build_reads_dict(obs_tab,leg_tab):
    """ Returns a dictionary that lists chromosome positions of SNPs and gives
        a dictionary that lists the observed bases at these postions. The
        nested dictionary gives all the read IDs that contain the observed base
        at the specific chromosome position. """
    
    reads = collections.defaultdict(list)

    for (pos, ind, read_id, base) in obs_tab:
        if base in leg_tab[ind][2:]:
            reads[read_id].append((pos,base))  
       
    return reads

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

def create_frequency(hap_dict):
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

    def joint_frequency(*alleles):
        """ Based on the reference panel, it calculates joint frequencies of
            observed alleles. The function arguments are alleles, that is,
            tuples of position and base, e.g., (100,'T'), (123, 'A') and
            (386, 'C'). The function arguments can also include haplotypes,
            that is, tuples of alleles; Haplotypes are treated in the same
            manner as alleles. """
           
        hap = intrenal_hap_dict(*alleles)   
        
        if len(hap)==1:
            result = hap[0].count(True) / N  
        
        elif len(hap)==2:
            result = sum(itertools.compress(*hap)) / N 
        
        elif len(hap)==3:
            result = sum(itertools.compress(hap[0], map(operator.and_,hap[1],hap[2]))) / N
        elif len(hap)>=4:
            result = sum(itertools.compress(hap[0], map(all,zip(*hap[1:])))) / N
        
        return result

    return joint_frequency

def read_statistics(reads_dict):
    return collections.Counter(len(i) for i in reads_dict.values())

def overlapping_reads(aux_dict,n):
    A = tuple(set(a) for pos in aux_dict for base in aux_dict[pos] if len(a:=aux_dict[pos][base])>=n)
    B = set(tuple(sorted(i)) for i in A) #overlaps atleast once
    C = set(tuple(sorted(a)) for i,j in zip(A[:-1],A[1:]) if len(a:=i.intersection(j))>1)  #overlaps atleast twice
    D = set(tuple(sorted(a)) for i,j,k in zip(A[:-2],A[1:-1],A[2:]) if len(a:=i.intersection(j,k))>1)  #overlaps atleast three times
    E = set(tuple(sorted(a)) for i,j,k,l in zip(A[:-3],A[1:-2],A[2:-1],A[3:]) if len(a:=i.intersection(j,k,l))>1)  #overlaps atleast four times
    return B,C,D,E

def build_blocks_dict(positions,block_size,offset):
    """ Returns a dictionary that lists blocks and gives all the SNP positions
        within them."""
    
    a = positions[0]-(block_size-1)+offset
    b = positions[-1]+(block_size-1)
    boundaries = tuple(range(int(a), int(b), int(block_size)))
    blocks = [(i,j-1) for i,j in zip(boundaries,boundaries[1:])]
    
    blocks_dict = collections.defaultdict(list)
    
    blocks_iterator = iter(blocks)
    block = next(blocks_iterator, None)
    for p in positions:
        while block:
            if block[0]<=p<=block[1]:
                blocks_dict[block].append(p)
                break
            block = next(blocks_iterator, None)
    return blocks_dict   

def HIST(x):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    #ax.hist(x,bins=int(len(x)**.5),histtype='step', linewidth=2.2, label='TBA')
    ax.hist(x,bins=100000, linewidth=0.2, label='TBA', histtype='stepfilled')
    ax.set_xlabel('Positions')
    ax.set_ylabel('Counts')
    ax.set_title('distribution of SNPs along chromosome 21' )
    ax.set_xticks(range(min(x), max(x),5000000)) #(max(x)-min(x))//20))
    #ax.set_xticks(range(min(x), max(x),100000)) #(max(x)-min(x))//20))
    #ax.legend()
    plt.show()
    
if __name__=='__main__':
    obs_filename = 'results_HapMix_EXT/mixed2haploids.X0.5.SRR10393062.SRR151495.0-2.hg38.obs.p'
    hap_filename = '../build_reference_panel/ref_panel.HapMix_EXT.hg38.BCFtools/chr21_HapMix_EXT_panel.hap'
    leg_filename = '../build_reference_panel/ref_panel.HapMix_EXT.hg38.BCFtools/chr21_HapMix_EXT_panel.legend'
    leg_tab,hap_tab,obs_tab = load(obs_filename,leg_filename,hap_filename)
    hap_dict = build_hap_dict(obs_tab,leg_tab,hap_tab)
    reads = build_reads_dict(obs_tab,leg_tab)
    aux_dict = build_aux_dict(obs_tab,leg_tab)
    frequency = create_frequency(hap_dict)
    block_dict = build_blocks_dict(tuple(aux_dict),100000,0)
