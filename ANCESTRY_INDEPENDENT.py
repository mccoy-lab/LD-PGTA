#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANCESTRY_INDEPENDENT

This script finds which SNPs exist across the five superpopulaions. Then, the
linkage disequilibrium (LD) between these common SNPs is calculated for each 
superpopulation. The LD between neighbor SNPs is compared across
the five superpopulation to track those with a ancestry-independent LD.  

Daniel Ariad (daniel@ariad.org)
Sep 31, 2020
"""
from operator import countOf, itemgetter
from itertools import islice, combinations, combinations_with_replacement
from collections import defaultdict
from time import time
from sys import stdout
from heapq import nlargest
from argparse import ArgumentParser
from sys import exit as sys_exit

def read_impute2(impute2_filename,**kwargs):
    """ Iterates over the rows of an IMPUTE2 file format (SAMPLE/LEGEND/HAPLOTYPE),
        parsers each row and yields a tuple with the parsed data. """
    
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
       
        for line in impute2_in:
            yield parse(line)

def build_LD_matrix(hap_dict,N,max_dist):
    """ Builds a sparse triangular matrix that contains the LD between each
    SNP pairs (pos1,pos2) with pos1<pos2 and pos2-pos1<mix_dist. The matrix is
    stored in the nested dictionary, result[pos1][pos2]. """

    time0 = time()
    result = defaultdict(dict)
    M = len(hap_dict)
    freq = {pos:  bin(hap).count('1')/N for pos,hap in hap_dict.items()}

    for i,(pos1,hap1) in enumerate(hap_dict.items()):
        for (pos2,hap2) in islice(hap_dict.items(),i+1,None):
            if pos2-pos1 >= max_dist: break
            joint_freq = bin(hap1&hap2).count('1')/N
            result[pos1][pos2] = joint_freq/(freq[pos1]*freq[pos2])
            
        if (i+1)%1000==0 or i+1==M:
            stdout.write('\r')
            stdout.write(f"Proceeded {(i+1)} / {M} rows. Average time per 1000 rows: {'%.2f' % (1000*(time() - time0)/(i+1))} sec")
            stdout.flush()
    stdout.write('\n')
    return result
            
def transpose(matrix):
    """ Transposes the LD matrix. """
    result = defaultdict(dict)
    for i,t in matrix.items():
        for j,k in t.items():
            result[j][i] = k
    return result

def symmetrize(matrix):
    """ Returns a symmetric matrix, given a triangular matrix. """
    result = defaultdict(dict)
    for i,t in matrix.items():
        for j,k in t.items():
            result[j][i] = k
            result[i][j] = k
    return result

def remove(matrix,omit):
    """ Return a copy of the LD matrix with specific rows and columns omitted.
        The set, omit contains the numbers of the rows and columns to be
        included. """ 
    result = {i:{j: k for j,k in t.items() if j not in omit} for i,t in matrix.items() if (i not in omit) and (not omit >= t.keys())} if len(omit) else matrix
    return result

def keep(matrix,include):
    """ Return a copy of the LD matrix where only specific rows and columns
        are included. The set, include contains the numbers of the rows and
        columns to be included. """   
    result = {i:{j: k for j,k in t.items() if j in include} for i,t in matrix.items() if (i in include) and (not include.isdisjoint(t.keys()))}
    return result

def replicate(matrix, example_matrix):
    """ Return a copy of the LD matrix where only specific elements are included. 
    The non-zero elements of the example matrix determine which elements of the
    LD matrix would be copied. """
    
    result = {i:{j: matrix[i][j] for j in row} for i,row in example_matrix.items()} if example_matrix else matrix 
    return result

def define_aux(lower_limit,upper_limit):
    """ Given lower and upper limits, creates a function that determines 
    whether two numbers are close to each other """
        
    def aux(a,b):
        """ Checks whether floats a and b have the same order of magnitude. """
        if a<=1e-16 and b<=1e-16:
            result = True
        elif (a<1e-16 and b>1e-16) or (a>1e-16 and b<1e-16):
            result = False
        else:
            result = lower_limit<a/b<upper_limit
        return result
    
    return aux
    
def stability(A,B):
    """ Compares each element in matrix A with the corresponding element in
    matrix B. And returns a dictionary that and gives the number of mismatched
    elements for each matrix row. """
    
    result = {}
    for ((Ak, Av), (Bk, Bv)) in zip(A.items(),B.items()):
            if Ak!=Bk: print('fatal error: dictionaries are out of sync.',Ak,Bk)
            r = countOf(map(aux,Av.values(),Bv.values()),False) / (len(Av) or 1)
            if r!=0:
                result[Ak] = r 
    return result

def get_common(leg_iter,hap_iter):
    """ Creates a reference panel for each superpopulation. The reference panel
    contains only SNPs that are common to all the five superpopulations. """
    
    superpop = {'AFR':(0,659),'AMR':(660,1006),'EAS':(1007,1510),'EUR':(1511,2013),'SAS':(2014,2502)}
    hap_ref = defaultdict(dict)
    hap_alt = defaultdict(dict)

    for leg, hap in zip(leg_iter,hap_iter):
        if all(countOf(hap[a:b+1],True)%(b+1-a) for (a,b) in superpop.values()):
            for sp,(a,b) in superpop.items():
                hap_ref[sp][leg[1]] = int(''.join(f'{not i:d}' for i in hap[a:b+1]),2)
                hap_alt[sp][leg[1]] = int(''.join(f'{i:d}' for i in hap[a:b+1]),2) 
    N = len(hap)
    return hap_ref, hap_alt, N

def filtering_multisteps(A,B,step):
    """ Removes SNPs from LD matrices A and B that do not share the same order
    of LD magnitude with all their neigbour SNPs. In each iteration removes 
    multiple SNPs."""
    
    q = stability(A, B)
    while(len(q)):
        t0 = time()
        pos, values = zip(*nlargest(min(step,len(q)), q.items(), key=itemgetter(1)))
        POSITIONS = {*pos}
        A = remove(A,POSITIONS) 
        B = remove(B,POSITIONS)
        q = stability(A, B)
        t1 = time()
        print(len(A),values[0],values[-1],t1-t0)
    return A, B

def filtering_singlestep(A,B,step):
    """ Removes SNPs from LD matrices A and B that do not share the same order
    of LD magnitude with all their neigbour SNPs. In each iteration removes 
    a single SNP. """
    
    q = stability(A, B)
    while(len(q)):
        t0 = time()
        pos, fraction = max(q.items(), key=itemgetter(1))
        POSITIONS = {pos}
        A = remove(A,POSITIONS) 
        B = remove(B,POSITIONS)
        q = stability(A, B)
        t1 = time()
        print(len(A),fraction,t1-t0)
    return A, B


def write_impute2(read_filename, write_filename, lines, filetype):
    """ Return a partial copy of the IMPUTE2 hap/legend file where only the
    line numbers within the set, lines are copied. """   
    
    with open(read_filename, 'r') as impute2_in:
        with open(write_filename, 'w') as impute2_out:
            if filetype=='legend': impute2_out.write(impute2_in.readline())
            for i,line in enumerate(impute2_in):
                if i in lines:
                    impute2_out.write(line)
    return 0

def main(hap_filename,leg_filename,output_impute2_filename,max_dist,step,lower_limit,upper_limit):
    """ Mainly compares the LD matrices of all the superpopulation to LD matrix
    of the european LD matrix to create a list of SNPs that their LD with their
    neighbor SNPs is stable across superpopulations. """
    global aux
    filtering = filtering_singlestep if step==1 else filtering_multisteps
    time0 = time()
    aux = define_aux(lower_limit,upper_limit)
    leg_iter = read_impute2(leg_filename, filetype='leg')
    hap_iter = read_impute2(hap_filename, filetype='hap')     
    hap_ref, hap_alt, N = get_common(leg_iter,hap_iter)
    hap_dict = {False: hap_ref, True: hap_alt}
    SP0 = False
    for sp0,sp1 in combinations(('EUR','AFR','AMR','EAS','SAS'),2):
        print('Comparing the LD between two superpopulations, %s and %s.' % (sp0,sp1))
        for inv0, inv1 in combinations_with_replacement((True,False),2):
            SP1 = replicate(symmetrize(build_LD_matrix(hap_dict[inv1][sp1], N, max_dist=max_dist)), SP0)
            SP0 = replicate(symmetrize(build_LD_matrix(hap_dict[inv0][sp0], N, max_dist=max_dist)), SP0)
            SP0, SP1 = filtering(SP0, SP1, step=step)
        print('Done.',len(SP0))
    POSITIONS = {*SP0}
    leg_iter = read_impute2(leg_filename, filetype='leg')
    lines = [i for i,j in enumerate(leg_iter) if j[1] in POSITIONS]
    write_impute2(hap_filename, output_impute2_filename+'.hap', lines, 'hap')
    write_impute2(leg_filename, output_impute2_filename+'.legend', lines, 'legend')
    time1 = time()
    print(f'Done in {(time1-time0):.2f} sec.')
    return 0
"""  
if __name__=='__main__':
    parser = ArgumentParser(
        description='Creates a multi-ethnic reference panel.')
    parser.add_argument('leg_filename', metavar='LEG_FILENAME', type=str,
                        help='IMPUTE2 legend file')
    parser.add_argument('hap_filename', metavar='HAP_FILENAME', type=str,
                        help='IMPUTE2 haplotype file')
    parser.add_argument('output_impute2_filename', metavar='FILENAME', type=str,
                        help='IMPUTE2 filename with no extension')
    parser.add_argument('-d', '--max-dist', type=int, 
                        metavar='INT', default=50000, 
                        help='Consider only the LD between SNPs in a range of max-dist. Default value 50,000')
    parser.add_argument('-s', '--step', type=int, 
                        metavar='INT', default=10,
                        help='The number of SNPs to be removed in each iteration. Default value 10.')
    parser.add_argument('-l', '--lower-limit', type=float,
                        metavar='FLOAT', default=0.1,
                        help='The minimal ratio between LD in two different superpopulations for considering them as close to each other. Default value 0.1.')
    parser.add_argument('-u', '--upper-limit', type=float, 
                        metavar='FLOAT', default='10', 
                        help='The maximal ratio between LD in two different superpopulations for considering them as close to each other. Default value 10.')

    x = main(**vars(parser.parse_args()))
    sys_exit(x)
"""    
if __name__=='__main__':
    hap_filename = '../build_reference_panel/ref_panel.ALL.hg38.BCFtools/chr21_ALL_panel.hap'
    leg_filename = '../build_reference_panel/ref_panel.ALL.hg38.BCFtools/chr21_ALL_panel.legend'
    output_impute2_filename = 'chr21_COMMON_panel_3'
    max_dist = 50000
    step = 1
    lower_limit = 0.1
    upper_limit = 10.00
    x = main(hap_filename,leg_filename,output_impute2_filename,max_dist,step,lower_limit,upper_limit)
    sys_exit(x)
