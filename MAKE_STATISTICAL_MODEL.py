#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 13:37:57 2020

@author: ariad
"""
import itertools,operator,math,pickle,time

def AUX(X):
    """ Checks if the first non-zero integer in an array is one. """   
    for i in X: 
        if i!=0:
            return True if i==1 else False

def GROUP(X):
    """ divide the indices of the list X according to the elements value."""
    group = {0: [], 1: [], 2: []}
    for i,j in enumerate(X):
        group[j].append(i)
    return tuple([tuple(i) for i in group.values() if len(i)!=0])

def DIV(a,b):
    """ Divides two numbers by their common divisor. """
    c = math.gcd(a,b)
    return (a//c,b//c)

def COMMON(X,t):
    """ Common procedure for both the SPH and BPH functions. """
    A = ((a,len(tuple(b))) for a,b in itertools.groupby(sorted(X)))
    B = ((a,DIV(b,t)) for a,b in A)
    C = ((GROUP(a),b) for a,b in B)
    return tuple(C)
    
def SPH(x):
    """ Returns the SPH probability function """
    A = (i if i[0] else tuple(map(operator.not_,i)) for i in itertools.product([False,True,True],repeat=x))
    result = COMMON(A,3**x)
    return result

def BPH(x):
    """ Returns the BPH probability function """
    A = (tuple(map(lambda x: operator.mod(operator.sub(x,i[0]),3),i)) for i in itertools.product([-1,0,1],repeat=x))
    B = (tuple(map(lambda x: x if AUX(i) else (x if x==0 else (2 if x==1 else 1)) ,i)) for i in A)
    result = COMMON(B,3**x)
    return result

def BUILD(x):
    models = dict()
    for i in range(2,x+1):
        print('Building the statistical model for %d reads.' % i)
        a = time.time()
        models[i]= {'BPH': BPH(i), 'SPH': SPH(i)}
        b = time.time()
        print('Done building in %.3f sec.' % ((b-a)))
    with open( 'MODELS.p', "wb") as f:
            pickle.dump(models, f)
    return models

if __name__ == "__main__": 
    print("Executed when invoked directly")
    BUILD(16)
else: 
    print("Executed when imported")
    