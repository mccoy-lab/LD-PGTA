#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 13:05:06 2020

@author: ariad
"""
import itertools, operator, collections, math, pickle, time

def SPH_V1(N):
    bank = collections.defaultdict(int)
    for C in itertools.product('ab',repeat=N):
        group1, group2, weight = list(), list(), 1
        for ind,letter in enumerate(C):
            if letter=='a':
                group1.append(ind)
                weight *= 1
            elif letter=='b':
                group2.append(ind)
                weight *= 2
                
        if len(group1)<len(group2):
             key = (tuple(group1),tuple(group2))
        elif len(group1)>len(group2):
             key = (tuple(group2),tuple(group1))
        else:
            if group1[0]<group2[0]:
                key = (tuple(group1),tuple(group2))
            else:
                key = (tuple(group2),tuple(group1))
                
        bank[key]+=weight
    return bank

def SPH_V2(N):
    bank = collections.defaultdict(int)
    w = {'a': 1, 'b': 2}
    for C in itertools.product('ab',repeat=N):
        groups = collections.defaultdict(list)
        weight = 1
        for ind,letter in enumerate(C):
            groups[letter].append(ind)
            weight *= w[letter]
        groups2 = tuple(tuple(i) for i in groups.values())
        key = sorted(groups2,key=lambda x: (len(x), x[0] if len(x)>0 else 0))
        bank[tuple(key)]+=weight
    return bank            

def BPH_V2(N):
    bank = collections.defaultdict(int)
    w = {'a': 1, 'b': 1, 'c': 1}
    for C in itertools.product('abc',repeat=N):
        groups = collections.defaultdict(list)
        weight = 1
        for ind,letter in enumerate(C):
            groups[letter].append(ind)
            weight *= w[letter]
        groups2 = tuple(tuple(i) for i in groups.values())
        key = sorted(groups2,key=lambda x: (len(x), x[0] if len(x)>0 else 0))
        bank[tuple(key)]+=weight
    return bank    

def BPH_V3(N):
    model = collections.defaultdict(int)
    O = tuple(itertools.permutations(range(3)))+((3,0,0),(0,3,0),(0,0,3),(1,1,1))
    for i in itertools.product(range(len(O)),repeat=N): 
        d = {j: {'a': O[k][0], 'b': O[k][1], 'c': O[k][2]} for j,k in enumerate(i)}
        for C in itertools.product('abc',repeat=N):
            groups = collections.defaultdict(list)
            weight = 1
            for ind,letter in enumerate(C):
                groups[letter].append(ind)
                weight *= d[ind][letter]
            groups2 = tuple(tuple(i) for i in groups.values())
            key = sorted(groups2,key=lambda x: (len(x), x[0] if len(x)>0 else 0))
            model[tuple(key)]+=weight
    #return model  
    compact = collections.defaultdict(list)
    T = sum(model.values())
    for haplotypes,weight in model.items():
        compact[weight//(gcd:=math.gcd(weight,T)),T//gcd].append(haplotypes)
    return {key: tuple(values) for key,values in compact.items()} 

def SPH_V3(N):
    model = collections.defaultdict(int)
    for i in itertools.product((0,1),repeat=N): 
        d = {j: {'a': 2+k, 'b': 1-k} for j,k in enumerate(i)}
        for C in itertools.product('ab',repeat=N):
            groups = collections.defaultdict(list)
            weight = 1
            for ind,letter in enumerate(C):
                groups[letter].append(ind)
                weight *= d[ind][letter]
            groups2 = tuple(tuple(i) for i in groups.values())
            key = sorted(groups2,key=lambda x: (len(x), x[0] if len(x)>0 else 0))
            model[tuple(key)]+=weight
    #return model  
    compact = collections.defaultdict(list)
    T = sum(model.values())
    for haplotypes,weight in model.items():
        compact[weight//(gcd:=math.gcd(weight,T)),T//gcd].append(haplotypes)
    return {key: tuple(values) for key,values in compact.items()} 

def SPH_V4(N):
    model = collections.defaultdict(int)
    O = ((1,1),)###tuple(itertools.permutations((0,1,2)))+((0,0,3),(3,0,0),(0,3,0))
    for i in itertools.product(range(len(O)),repeat=N): 
        d = {j: {str(m): O[k][m] for m in range(len(O[0]))} for j,k in enumerate(i)}
        for C in itertools.product(''.join(i for i in d[0].keys()),repeat=N):
            groups = collections.defaultdict(list)
            weight = 1
            for ind,letter in enumerate(C):
                groups[letter].append(ind)
                weight *= d[ind][letter]
            groups2 = tuple(tuple(i) for i in groups.values())
            key = sorted(groups2,key=lambda x: (len(x), x[0] if len(x)>0 else 0))
            model[tuple(key)]+=weight
    #return model  
    compact = collections.defaultdict(list)
    T = sum(model.values())
    for haplotypes,weight in model.items():
        compact[weight//(gcd:=math.gcd(weight,T)),T//gcd].append(haplotypes)
    return {key: tuple(values) for key,values in compact.items()} 

def BPH_V4(N):
    model = collections.defaultdict(int)
    O = ((1,1,1),)###tuple(itertools.permutations((0,1,2)))+((0,0,3),(3,0,0),(0,3,0))
    for i in itertools.product(range(len(O)),repeat=N): 
        d = {j: {str(m): O[k][m] for m in range(len(O[0]))} for j,k in enumerate(i)}
        for C in itertools.product(''.join(i for i in d[0].keys()),repeat=N):
            groups = collections.defaultdict(list)
            weight = 1
            for ind,letter in enumerate(C):
                groups[letter].append(ind)
                weight *= d[ind][letter]
            groups2 = tuple(tuple(i) for i in groups.values())
            key = sorted(groups2,key=lambda x: (len(x), x[0] if len(x)>0 else 0))
            model[tuple(key)]+=weight
    #return model  
    compact = collections.defaultdict(list)
    T = sum(model.values())
    for haplotypes,weight in model.items():
        compact[weight//(gcd:=math.gcd(weight,T)),T//gcd].append(haplotypes)
    return {key: tuple(values) for key,values in compact.items()} 


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
    for haplotypes,weight in model.items():
        compact[weight//(gcd:=math.gcd(weight,T)),T//gcd].append(haplotypes)
    return {key: tuple(values) for key,values in compact.items()} 

def SPH(N):
    return engine(N,{'a': 1, 'b': 2})

def BPH(N):
    return engine(N,{'a': 1, 'b': 1, 'c': 1})
   
        
def BUILD(x):
    models = dict()
    for i in range(2,x+1):
        print('Building the statistical model for %d reads.' % i)
        a = time.time()
        models[i]= {'BPH': BPH(i), 'SPH': SPH_V4(i)}
        b = time.time()
        print('Done building in %.3f sec.' % ((b-a)))
    with open( 'MODELS.p', "wb") as f:
            pickle.dump(models, f, protocol=4)
    return models

if __name__ == "__main__":
    print("Executed when invoked directly")
    #models = BUILD(12)
else:
    print("Executed when imported")