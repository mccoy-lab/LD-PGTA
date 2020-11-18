#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PLOT_STREDA
Daniel Ariad (daniel@ariad.org)
Aug 31, 2020
"""

import pickle, statistics

def mean_and_var(sample):
    """ Calculates the mean and the sample standard deviation. """
    mean = statistics.mean(sample)
    var = statistics.variance(sample, xbar=mean)
    return mean, var

def std_of_mean(variances):
    """ Standard error of the mean of uncorrelated variables, based on the
        Bienaymé formula. """
    return sum(variances)**.5/len(variances)

def jackknife_std(sample,weights):
    """ Given sample elements and the weight of each element, the jackknife
    standard deviation is calculated. 
    
    *** More information about delete-m jackknife for unequal m can be found in
    F.M.Busing et al. (1999), [DOI:10.1023/A:1008800423698]. """

    N = len(sample)
    t0 = sum(sample) / N
    S = sum(weights)
    H = [S/w for w in weights]
    T = [sum(sample[:i]+sample[i+1:])/(N-1) for i in range(N)]
    pseudo_values = [h*t0-(h-1)*t for t,h in zip(T,H)]
    jackknife_estimator = sum((p/h for p,h in zip(pseudo_values,H)))
    jackknife_variance = sum(((p-jackknife_estimator)**2/(h-1) for p,h in zip(pseudo_values,H)))/N
    return jackknife_variance**.5
    
def chr_length(chr_id):
    """ Return the chromosome length for a given chromosome, based on the reference genome hg38.""" 
    #The data of chromosome length was taken from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38
    length_dict = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
                  'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422,
                  'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
                  'chr16': 90338345, 'chr17':  83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
                  'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}
    return length_dict[chr_id]

def load_llr(filename):
    """ Loads from a file a dictionary that lists genomic windows that contain
    at least two reads and gives the bootstrap distribution of the 
    log-likelihood BPH/SPH ratios (LLRs). """
    
    with open(filename, 'rb') as f:
        LLR_dict = pickle.load(f)
        info = pickle.load(f)
    
    print('\nFilename: %s' % filename)
    print('Depth: %.2f, Number of genomic windows: %d, Fraction of genomic windows with a negative LLR: %.3f' % (info['depth'], info['statistics']['num_of_windows'],info['statistics']['fraction_of_negative_LLRs']))
    print('Mean and standard error of meaningful reads per genomic window: %.1f, %.1f.' % (info['statistics']['reads_mean'], info['statistics'].get('reads_std',info['statistics'].get('reads_var'))))
    print('Mean LLR: %.3f, Standard error of the mean LLR: %.3f' % ( info['statistics']['mean'], info['statistics']['std']))
    print('Calculation was done in %.3f sec.' % info['runtime'])

    return LLR_dict, info


def plot_streda(LLR_dict,info,N,**kwargs):
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    save = kwargs.get('save', '') 

    font_size=16
       
    M_dict = {block: sum(LLRs)/len(LLRs) for block,LLRs in LLR_dict.items() if None not in LLRs}

    K,V = zip(*M_dict.items())
    B, C = {}, {}
    l = chr_length(info['chr_id'])
    
    I = lambda j,N: j*(len(V)//N)
    F = lambda j,N: (j+1)*(len(V)//N) if j!=N-1 else len(K)-1
    for i in range(1,N+1):
        B[i] = tuple(((i,K[I(j,i)][0]/l),(i,K[F(j,i)][0]/l)) for j in range(i))
        C[i] = tuple(sum(v<0 for v in V[I(j,i):F(j,i)+1])/len(V[I(j,i):F(j,i)]) for j in range(i))

    segs = np.array([j for i in B.values() for j in i ])
    colors = np.array([j for i in C.values() for j in i ])
    fig, ax = plt.subplots(1, 1, figsize=(16, 9))  # setup the plot
    fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.07) 
    ax.set_xlim(segs[0][1][0]-1, segs[-1][1][0]+1)
    ax.set_ylim(segs[0][0][1], segs[0][1][1])
    ax.set_xlabel('Iterations', fontsize=font_size)
    ax.set_ylabel('Normalized chromosomal position', fontsize=font_size, labelpad=5) 

    
    # define the bins and normalize
    bounds = np.linspace(0,1,21,endpoint=True)
    norm = mpl.colors.BoundaryNorm(bounds, 21)
    cmaplistBWR = plt.cm.seismic(bounds)
    cmap = mpl.colors.LinearSegmentedColormap.from_list('BWR', cmaplistBWR, 20)
    
    line_segments = LineCollection(segs, array=colors, linewidth=10, linestyle='solid', cmap=cmap, norm=norm)
    ax.add_collection(line_segments)
    ax.set_title('.'.join(save.split('.')[:-2]),fontsize=font_size) 
    tick = np.linspace(0,1,11,endpoint=True)
    axcb = fig.colorbar(line_segments,ticks=tick, format='%.2f')
    axcb.set_label('Fraction of LD blocks with a negative LLR', fontsize=font_size)
    
    if save!='':
        print('Saving plot...')
        plt.savefig(save, dpi=150, facecolor='w',
                    edgecolor='w', orientation='landscape',
                    format='png', transparent=False, bbox_inches=None,
                    pad_inches=0.3, metadata=None)
        plt.close(fig)
    else:
        plt.show()
    
    return 1
