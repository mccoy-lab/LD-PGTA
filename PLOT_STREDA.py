#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 01:12:11 2020

@author: ariad
"""
import pickle

def load_llr(filename):
    with open('results_EUR/'+filename, 'rb') as f:
        llr = pickle.load(f)
        info = pickle.load(f)
    print('Depth: %.2f, Number of LD blocks: %d, Fraction of LD blocks with a negative LLR: %.3f' % (info['depth'], info['statistics']['num_of_LD_blocks'],info['statistics']['fraction_of_negative_LLRs']))
    print('Mean: %.3f, Standard error: %.3f, Jackknife standard error: %.3f' % ( info['statistics']['mean'], info['statistics']['std'], info['statistics']['jk_std']))
    print('Calculation was done in %.3f sec.' % info['runtime'])

    return llr, info

def plot_streda(LLR_dict0,**kwargs):
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    save = kwargs.get('save', False) 

    font_size=16
    
    LLR_dict = {block: sum(LLRs)/len(LLRs) for block,LLRs in LLR_dict0.items() if None not in LLRs}

    K,V = zip(*((k,v) for k,v in LLR_dict.items() if v!=None))
    B, C = {}, {}
    for i in range(12):
        a = len(V)//(i+1)
        B[i] = tuple(((i,K[j*a][0]),(i,K[min((j+1)*a,len(K)-1)][-1])) for j in range(i+1))
        C[i] = tuple(sum(k<0 for k in V[j*a:(j+1)*a])/len(V[j*a:(j+1)*a]) for j in range(i+1))

    segs = np.array([j for i in B.values() for j in i ])
    colors = np.array([j for i in C.values() for j in i ])
    fig, ax = plt.subplots(1, 1, figsize=(16, 9))  # setup the plot
    fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.07) 
    ax.set_xlim(segs[0][1][0]-1, segs[-1][1][0]+1)
    ax.set_ylim(segs[0][0][1], segs[0][1][1])
    ax.set_xlabel('Iterations', fontsize=font_size)
    ax.set_ylabel('Chromosomal position (bp)', fontsize=font_size, labelpad=5) 

    
    # define the bins and normalize
    bounds = np.linspace(0,1,21,endpoint=True)
    norm = mpl.colors.BoundaryNorm(bounds, 21)
    cmaplistBWR = plt.cm.seismic(bounds)
    cmap = mpl.colors.LinearSegmentedColormap.from_list('BWR', cmaplistBWR, 20)
    
    line_segments = LineCollection(segs, array=colors, linewidth=40, linestyle='solid', cmap=cmap, norm=norm)
    ax.add_collection(line_segments)
    ax.set_title('Tracking down the recombination spot', fontsize=font_size)
    tick = np.linspace(0,1,11,endpoint=True)
    axcb = fig.colorbar(line_segments,ticks=tick, format='%.2f')
    axcb.set_label('Fraction of LD blocks with a negative LLR', fontsize=font_size)
    
    if save:
        print('Saving plot...')
        plt.savefig('conductivity.png', dpi=150, facecolor='w',
                    edgecolor='w', orientation='landscape', papertype='letter',
                    format='png', transparent=False, bbox_inches=None,
                    pad_inches=0.3, frameon=None, metadata=None)
        plt.close(fig)
    else:
        plt.show()
    
    return 1

def analyze(LLR_dict0):
    LLR_dict = {block: sum(LLRs)/len(LLRs) for block,LLRs in LLR_dict0.items() if None not in LLRs}
    K,V = zip(*((k,v) for k,v in LLR_dict.items() if v!=None))
    A = dict()
    B = dict()
    for i in range(12):
        #print(i+1)
        a = len(V)//(i+1)
        A[i] = tuple(sum(k<0 for k in V[j*a:(j+1)*a])/len(V[j*a:(j+1)*a]) for j in range(i+1))
        print('A%d:'%i,A[i])
        B[i] = tuple(((i,K[j*a][0]),(i,K[min((j+1)*a,len(K)-1)][-1])) for j in range(i+1))
        print('B%d:'%i,B[i])
    return tuple(B.values()),tuple(A.values())

def LDblockHIST(LLR_dict0):    
    import matplotlib.pyplot as plt
    LLR_dict = {block: sum(LLRs)/len(LLRs) for block,LLRs in LLR_dict0.items() if None not in LLRs}
    fig, ax = plt.subplots()
    x = [i for i in LLR_dict.values() if i!=None]
    ax.hist(x,bins=int(len(x)**.5),histtype='step', linewidth=2.2, label='LLR distribution accross LD blocks')
    ax.set_xlabel('Aggregated log-likelihood ratio')
    ax.set_ylabel('Counts')
    ax.set_title('LLR distribution accross LD blocks' )
    #ax.legend()
    plt.show()