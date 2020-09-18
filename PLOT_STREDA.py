#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PLOT_STREDA
Daniel Ariad (daniel@ariad.org)
Aug 31, 2020
"""

import pickle, statistics, itertools

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
    with open('results_PRJNA384616/'+filename, 'rb') as f:
        LLR_dict = pickle.load(f)
        info = pickle.load(f)
    
    print('\nFilename: %s' % filename)
    print('Depth: %.2f, Number of LD blocks: %d, Fraction of LD blocks with a negative LLR: %.3f' % (info['depth'], info['statistics']['num_of_LD_blocks'],info['statistics']['fraction_of_negative_LLRs']))
    print('Mean and standard error of meaningful reads per LD block: %.1f, %.1f.' % (info['statistics']['reads_mean'], info['statistics'].get('reads_std',info['statistics'].get('reads_var'))))
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
    F = lambda j: (j+1)*a if j!=i-1 else len(K)-1
    for i in range(1,N+1):
        a = len(V)//i
        B[i] = tuple(((i,K[j*a][0]/l),(i,K[F(j)][0]/l)) for j in range(i))
        C[i] = tuple(sum(v<0 for v in V[j*a:F(j)])/len(V[j*a:F(j)]) for j in range(i))

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

def plot_test(LLR_dict,info,N,**kwargs):
    import numpy as np
    import matplotlib.pyplot as plt
    save = kwargs.get('save', '')
    
    # Create lists for the plot

    TEST = [num_of_reads>info['max_reads']  
            for num_of_reads in info['statistics']['reads_per_LDblock_dict'].values() if num_of_reads>=info['min_reads'] ]    
    
    WEIGHTS = [min(num_of_reads,info['max_reads'])   
               for num_of_reads in info['statistics']['reads_per_LDblock_dict'].values() if num_of_reads>=info['min_reads']]    
    
    LLR_stat = {block: mean_and_var(LLRs)  
                for block,LLRs in LLR_dict.items() if None not in LLRs}
    
    K,M,V = tuple(LLR_stat.keys()), *zip(*LLR_stat.values())
            
    i = lambda j: j*(len(V)//N)
    f = lambda j: (j+1)*(len(V)//N) if j!=N-1 else len(K)-1
    x = lambda p,q: .5*(K[p][0] + K[q][-1])
    y = lambda p,q: statistics.mean(M[p:q])
    e = lambda p,q: std_of_mean(V[p:q]) if all(TEST[p:q]) else jackknife_std(M[p:q],WEIGHTS[p:q])
    
    X,Y,E = ([func(i(j),f(j)) for j in range(N)] for func in (x,y,e))

    widths = tuple((K[f(j)][0]-K[i(j)][0]) for j in range(N))    
    #X_boundaries = tuple(k for j in range(N+1) for k in (K[j*a][0], K[min((j+1)*a,len(K)-1)][-1]))
    X_ticks = [K[i(j)][0] for j in range(N)]+[K[-1][-1]]  
    X_labels = [('%.2f' % (j/chr_length(info['chr_id']))) for j in X_ticks] 
 
    # Build the plot
    fig, ax = plt.subplots(1, 1, figsize=(16, 9))  # setup the plot
    fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.07) 
    ax.bar(X, Y, yerr=E, align='center', alpha=0.5, ecolor='black', capsize=10, width=widths, color=[np.random.rand(3,) for _ in range(N+1)])
    ax.set_ylabel('log-likelihood BPH/SPH ratio')
    ax.set_xlabel('Normalized chromosome position')
    ax.set_xticks(X_ticks)
    ax.set_xticklabels(X_labels)
    ax.set_title('.'.join(save.split('.')[:-2]))
    #ax.xaxis.grid(True)
    #ax.yaxis.grid(True)

    
    for j in range(N):
        plt.text(X[j], .5*(Y[j]-Y[j]/abs(Y[j])*E[j]), '%.2f\u00B1%.2f'% (Y[j],E[j]), ha='center', va='center',color='black',fontsize=8)
    
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
   
