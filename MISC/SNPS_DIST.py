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
import matplotlib as mpl
from matplotlib import pyplot as plt
    
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

def plot_dist(SP, chr_id):
    #filename = f'../build_reference_panel/ref_panel.{SP:s}.hg38.BCFtools/{chr_id:s}_{SP:s}_panel.legend'
    filename = f'/Users/ariad/Downloads/{SP:s}_panel.hg38.BCFtools/{chr_id:s}_{SP:s}_panel.legend'
    leg_tab = read_impute2(filename, filetype='legend')
    a,b,c = tuple(zip(*leg_tab))
    X = [a[i]-a[i-1] for i in range(1,len(a))]
    ############
    font_size=18
    fig, ax = plt.subplots(1, 1, figsize=(16, 9))  # setup the plot
    fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.07) 
    weights = [1/len(X)] * len(X)
    result = plt.hist(X, bins=range(0,350,35), facecolor='darkblue', weights=weights)
    plt.xlabel('Distance between nearest neighbor SNPs',fontsize=font_size)
    plt.ylabel('Density of SNPs',fontsize=font_size)
    plt.title(r'Distribution of human SNPs in chromosome 21',fontsize=font_size)
    ax.set_xticks([i for i in range(0,350,35)])
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    plt.show()
    return result

"""
for SP in ('EUR','AFR','AMR','EAS','SAS'):
    for j in (6,8,9,12,13,15,18,21,22):
        print(SP,'chr'+str(j))
        chr_id = 'chr'+str(j)
        A = plot_dist(SP, chr_id)
        plt.close()
        print(A[0][0])
"""
