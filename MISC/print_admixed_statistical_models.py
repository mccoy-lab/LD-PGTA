#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 22:26:10 2021

@author: ariad
"""

def BPH(n):
    import itertools
    result = []
    for W in itertools.product(['F','G','H'], repeat=n):
        temp = {'F':[],'G':[],'H':[]}
        for i,w in enumerate(W):
            temp[w].append(chr(ord('A')+i))
        temp2 = tuple((i.replace('H','G'),''.join(j)) for i,j in temp.items() if len(j)!=0)
        #temp3 = tuple((i.replace('G','H').replace('F','G').replace('H','F'),''.join(j)) for i,j in temp.items() if len(j)!=0)
        temp3 = tuple((i.replace('H','F'),''.join(j)) for i,j in temp.items() if len(j)!=0)

        result.append(temp2)
        result.append(temp3)
        
    X = tuple(frozenset(k if j=='F' else k.lower() for j,k in i) for i in result)
        
    #X  = {i: result.count(i) for i in result}
    
    import collections
    Y = dict(collections.Counter(X))
    Z = {i:[] for i in range(50)}
    for i,j in Y.items():
        Z[j].append('*'.join(i))
    import math
    P = '+'.join([f'({i//math.gcd(i,2*3**n):d}/{2*3**n//math.gcd(i,2*3**n):d})*('+'+'.join(j)+')' for i,j in Z.items() if len(j)!=0])    
    P2 =  '('+'+'.join([f'{i:d}*('+'+'.join(j)+')' for i,j in Z.items() if len(j)!=0])+f')/{2*3**n:d}'    

    #print(P)    
    #print('+'.join(['*'.join([str(j),*i]) for (i,j) in Y.items()]))
    return P2

def DIPLOID(n):
    import itertools
    result = []
    for W in itertools.product(['F','G'], repeat=n):
        temp = {'F':[],'G':[]}
        for i,w in enumerate(W):
            temp[w].append(chr(ord('A')+i))
        temp2 = tuple((i,''.join(j)) for i,j in temp.items() if len(j)!=0)
        temp3 = tuple((i.replace('G','H').replace('F','G').replace('H','F'),''.join(j)) for i,j in temp.items() if len(j)!=0)
        result.append(temp2)
        result.append(temp3)
        
    X = tuple(frozenset(k if j=='F' else k.lower() for j,k in i) for i in result)
        
    #X  = {i: result.count(i) for i in result}
    
    import collections
    Y = dict(collections.Counter(X))
    Z = {i:[] for i in range(50)}
    for i,j in Y.items():
        Z[j].append('*'.join(i))
    import math
    P = '+'.join([f'({i//math.gcd(i,2*2**n):d}/{2*2**n//math.gcd(i,2*2**n):d})*('+'+'.join(j)+')' for i,j in Z.items() if len(j)!=0])    
    P2 =  '('+'+'.join([f'{i:d}*('+'+'.join(j)+')' for i,j in Z.items() if len(j)!=0])+f')/{2*2**n:d}'    

    #print(P)    
    #print('+'.join(['*'.join([str(j),*i]) for (i,j) in Y.items()]))
    return P2

def SPH(n):
    import itertools
    result = []
    for W in itertools.product(['F','G','G'], repeat=n):
        temp = {'F':[],'G':[],'H':[]}
        for i,w in enumerate(W):
            temp[w].append(chr(ord('A')+i))
        temp2 = tuple((i,''.join(j)) for i,j in temp.items() if len(j)!=0)
        temp3 = tuple((i.replace('G','H').replace('F','G').replace('H','F'),''.join(j)) for i,j in temp.items() if len(j)!=0)
        result.append(temp2)
        result.append(temp3)
        
    X = tuple(frozenset(k if j=='F' else k.lower() for j,k in i) for i in result)
        
    #X  = {i: result.count(i) for i in result}
    
    import collections
    Y = dict(collections.Counter(X))
    Z = {i:[] for i in range(50)}
    for i,j in Y.items():
        Z[j].append('*'.join(i))
    import math
    P = '+'.join([f'({i//math.gcd(i,2*3**n):d}/{2*3**n//math.gcd(i,2*3**n):d})*('+'+'.join(j)+')' for i,j in Z.items() if len(j)!=0])    
    
    P2 = '('+'+'.join([f'{i:d}*('+'+'.join(j)+')' for i,j in Z.items() if len(j)!=0])+f')/{2*3**n:d}' 
    #print(P)    
    #print('+'.join(['*'.join([str(j),*i]) for (i,j) in Y.items()]))
    return P2