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

def chr_length(chr_id):
    """ Return the chromosome length for a given chromosome, based on the reference genome hg38.""" 
    #The data of chromosome length was taken from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38
    length_dict = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
                  'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422,
                  'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
                  'chr16': 90338345, 'chr17':  83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
                  'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}
    return length_dict[chr_id]

def plot_dens():
    chr_id = 'chr21'
    num_of_bins = int(chr_length(chr_id)/3e5)
    leg_tab = read_impute2('../build_reference_panel/ref_panel.EUR.hg38.BCFtools/%s_EUR_panel.legend' % chr_id, filetype='legend')
    A,*_= tuple(zip(*leg_tab))
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    font_size=18
    X = [a/chr_length(chr_id) for a in A]
    fig, ax = plt.subplots(1, 1, figsize=(16, 9))  # setup the plot
    fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.07) 
    plt.hist(X, bins=[i/num_of_bins for i in range(num_of_bins)], facecolor='darkblue', weights=[1/len(X)]*len(X))
    plt.xlabel('Normalized chromosomal position',fontsize=font_size)
    plt.ylabel('Density of SNPs',fontsize=font_size)
    plt.title(r'Distribution of human SNPs along chromosome 21',fontsize=font_size)
    plt.show()