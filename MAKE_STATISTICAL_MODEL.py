#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAKE_STATISTICAL_MODEL

Builds statistical models for aneuploidy cells with BPH and SPH.
BPH (Both Parental Homologs) correspond to the presence of three unmatched
haplotypes, while SPH (Single Parental Homolog) correspond to chromosome gains
involving identical homologs.

Daniel Ariad (daniel@ariad.org)
June 5th, 2020
"""
import itertools, collections, math, pickle, time, bz2

def engine(N,degeneracies):
    model = collections.defaultdict(int)
    for C in itertools.product(''.join(i for i in degeneracies.keys()),repeat=N):
        groups = collections.defaultdict(list)
        weight = 1
        for ind,letter in enumerate(C):
            groups[letter].append(ind)
            weight *= degeneracies[letter]
        haplotypes = tuple(tuple(i) for i in groups.values())
        key = sorted(haplotypes,key=lambda x: (len(x), x[0] if len(x)>0 else 0))
        model[tuple(key)]+=weight
    
    compact = collections.defaultdict(list)
    T = sum(degeneracies.values())**N
    while(len(model)!=0):
        haplotypes, weight = model.popitem()
        HAPLOTYPES = tuple(''.join(chr(x+65) for x in h) for h in haplotypes)
        gcd = math.gcd(weight,T)
        compact[weight//gcd,T//gcd].append(HAPLOTYPES)
    for key in compact: compact[key] = tuple(compact[key])
    return compact

def SPH(N):
    return engine(N,{'a': 1, 'b': 2})

def BPH(N):
    return engine(N,{'a': 1, 'b': 1, 'c': 1})
    
def BUILD(x):
    models = dict()
    for i in range(2,x+1):
        print('Building the statistical model for %d reads.' % i)
        a = time.time()
        models[i]= {'BPH': BPH(i), 'SPH': SPH(i)}
        b = time.time()
        print('Done building in %.3f sec.' % ((b-a)))
    with open( 'MODELS.p', 'wb') as f:
        pickle.dump(models, f, protocol=4)
    with bz2.BZ2File( 'MODELS.pbz2', 'wb') as f:
        pickle.dump(models, f, protocol=4)
    return models

if __name__ == "__main__":
    print("Executed when invoked directly")
    models = BUILD(17)
else:
    print("Executed when imported")
   
    