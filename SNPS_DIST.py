#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script plots a histogram of the distance between two nearest-neighbor SNPs. 
The plot reveals that SNPs are not distributed evenly along the chromosome and
tend to bunch together, according to a power law of distance.
(for more details see DOI:10.1111/gtc.12344)


Created on Sun May  3 08:06:57 2020

@author: ariad
"""

def read_impute2(impute2_filename,**kwargs):
    """ Reads an IMPUTE2 file format (SAMPLE/LEGEND/HAPLOTYPE) and builds a list
        of lists, containing the dataset. """
    
    filetype = kwargs.get('filetype', None)
    impute2_in = open(impute2_filename, 'r')
    
    if filetype == 'legend':
        impute2_in.readline()   # Bite off the header
        def parse(x): y=x.strip().split()[1:]; y[0]=int(y[0]); return y
    elif filetype == 'hap':
        def parse(x): return [i=='1' for i in x.strip().split()]
    else:
        def parse(x): return x.strip().split() # Trailing whitespaces stripped from the ends of the string. Then it splits the string into a list of words.
   
    impute2_tab = [parse(line) for line in impute2_in]
    impute2_in.close()
    return impute2_tab


leg_tab = read_impute2('../build_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.legend', filetype='legend')
a,b,c = tuple(zip(*leg_tab))
X = [a[i]-a[i-1] for i in range(1,len(a))]
import matplotlib as mpl
from matplotlib import pyplot as plt
weights = [1/len(a) for i in range(len(a)-1)]
plt.hist(X, bins=range(0,4000,60), facecolor='blue', weights=weights)
plt.xlabel('Distance between nearest neighbor SNPs')
plt.ylabel('Fraction')
plt.title(r'EUR ref panel of CHR21')
plt.show()