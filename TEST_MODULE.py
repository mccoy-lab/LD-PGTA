#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TEST_MODULE

Generates a random reference panel as well as a random observations table
for a single LD block. 

Daniel Ariad (daniel@ariad.org)
Sep 30, 2021
"""

import collections, random, pickle, os

random.seed(a=2021,version=2)

leg_tuple = collections.namedtuple('leg_tuple', ('chr_id', 'pos', 'ref', 'alt')) #Encodes the rows of the legend table
sam_tuple = collections.namedtuple('sam_tuple', ('sample_id', 'group1', 'group2', 'sex')) #Encodes the rows of the samples table
obs_tuple = collections.namedtuple('obs_tuple', ('pos', 'read_id', 'base')) #Encodes the rows of the observations table

def generate_legned(SNPs=24,depth=0.1):
    chr_id = 'chrTEST'
    POS = sorted(random.sample(range(int(SNPs/depth)), SNPs))
    REF, ALT = zip(*(random.sample('ATCG',k=2) for i in range(SNPs)))
    legend = tuple(leg_tuple(chr_id,pos,ref,alt) for pos,ref,alt in zip(POS,REF,ALT))
    return legend

def generate_haplotype(SNPs=24,number_of_haplotypes=10):
    max_int = 2**number_of_haplotypes
    hap_tab = tuple(random.randrange(1,max_int+1) for i in range(SNPs))     
    return hap_tab, number_of_haplotypes

def generate_obs(legend,alleles_per_read=4):
    reads_generator = (i for i in range(len(legend)//alleles_per_read) for j in range(alleles_per_read))
    obs_tab = tuple(obs_tuple(pos,'read'+str(read_id),random.choice(alleles)) 
     for read_id,(chr_id,pos,*alleles) in zip(reads_generator,legend))
    return obs_tab

def generate_sam(groups=1,number_of_haplotypes=10):
    sample_ids = ('sample'+str(i) for i in range(number_of_haplotypes//2))
    groups2 = ('group'+str(i) for i in random.choices(range(groups),k=number_of_haplotypes//2))
    sex = random.choices(range(2),k=number_of_haplotypes//2)
    sam_tab = tuple(sam_tuple(i,'N/A',j,k) for i,j,k in zip(sample_ids,groups2,sex))
    return sam_tab

def print_haplotypes(hap_tab,number_of_haplotypes):
    return tuple(tuple(bin(h)[2:].zfill(number_of_haplotypes)) for h in hap_tab)

leg_tab = generate_legned(SNPs=1000)
hap_tab, number_of_haplotypes = generate_haplotype(SNPs=1000,number_of_haplotypes=250)  ### number of haplotypes must be even, because diploids are simulated.
obs_tab = generate_obs(leg_tab,alleles_per_read=1)
sam_tab = generate_sam(groups=2,number_of_haplotypes=250)

if not os.path.exists('test'): os.makedirs('test')

info = {'redo-BAQ': False,
        'handle-multiple-observations' : 'all',
        'min-bq': 30,
        'min-mq' :  30,
        'max-depth' :  0,
        'chr_id': 'chrTEST',
        'depth': 0.1}

with open('test/test.obs.p','wb') as f:
    pickle.dump(obs_tab,f)
    pickle.dump(info,f)
    
with open('test/test.leg.p','wb') as f:
    pickle.dump(leg_tab,f)

with open('test/test.hap.p','wb') as f:
    pickle.dump((hap_tab,number_of_haplotypes),f)
    
with open('test/test.sam.p','wb') as f:
    pickle.dump(sam_tab,f)

"""
operator, itertools
 
def get_joint_freq(hap_tab,number_of_haplotypes,invert,*include):
    nt = lambda p: 1-p
    new_hap_tab = [ ([*map(int,i)] if not j else [*map(nt,map(int,i))])
             for i,j in zip(print_haplotypes(hap_tab,number_of_haplotypes),invert)]
    #print(z)
    ig = operator.itemgetter(*include)
    t = [all(i) for i in zip(*ig(new_hap_tab))]
    return sum(t)/len(t)

def print_freq(leg_tab,obs_tab,hap_tab,number_of_haplotypes):
    inv = [base==ref for (pos, read_id, base),(chr_id, pos, ref, alt) in zip(obs_tab,leg_tab)]
    for i in range(1,2**4):
        X = []
        for j in bin(i)[2:]:
            X.extend([j=='1']*4)
        #print(bin(i)[2:])
        #print(''.join(map(str,map(int,X))).zfill(16))
        print(bin(i)[2:], get_joint_freq(hap_tab,number_of_haplotypes,inv,*[k for k,l in enumerate(X[::-1]) if l]))
    
        
print_freq(leg_tab,obs_tab,hap_tab,number_of_haplotypes)
"""