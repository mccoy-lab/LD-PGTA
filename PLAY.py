#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 14:18:44 2020

@author: ariad
"""

import pickle, statistics, re
from statistics import mean, variance
from itertools import starmap
from math import log
from collections import Counter

HOME = '/home' # '/Users'

def chr_length(chr_id):
    """ Return the chromosome length for a given chromosome, based on the reference genome hg38.""" 
    #The data of chromosome length was taken from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38
    length_dict = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
                  'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422,
                  'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
                  'chr16': 90338345, 'chr17':  83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
                  'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}
    return length_dict[chr_id]

def std_of_mean(variances):
    """ Standard error of the mean of uncorrelated variables, based on the
        Bienaymé formula. """
    return sum(variances)**.5/len(variances)

def mean_and_var(data):
    """ Calculates the mean and variance. """
    m = mean(data)
    var = variance(data, xbar=m)
    return m, var 

def LLR(y,x):
    """ Calculates the logarithm of y over x and deals with edge cases. """
    if x and y:
        result = log(y/x)
    elif x and not y:
        result = -1.23456789 
    elif not x and y:
        result = +1.23456789 
    elif not x and not y:
        result = 0 
    else:
        result = None    
    return result

def load_likelihoods(filename):
    """ Loads from a file a dictionary that lists genomic windows that contain
    at least two reads and gives the bootstrap distribution of the 
    log-likelihood BPH/SPH ratios (LLRs). """
    
    with open(filename, 'rb') as f:
        likelihoods = pickle.load(f)
        info = pickle.load(f)
    return likelihoods, info

def show_info(filename,info,pair):
    S = info['statistics']
    X = info['statistics']['LLRs_per_chromosome'][pair]
    print('\nFilename: %s' % filename)
    print('Depth: %.2f, Chromosome ID: %s, Mean and standard error of meaningful reads per genomic window: %.1f, %.1f.' % (info['depth'], info['chr_id'], S['reads_mean'], S['reads_std']))
    print('Number of genomic windows: %d, Mean and standard error of genomic window size:  %d, %d.' % (S['num_of_windows'],S['window_size_mean'],S['window_size_std']))
    print('Mean LLR: %.3f, Standard error of the mean LLR: %.3f' % ( X['mean'],  X['std_of_mean']))
    print('Fraction of genomic windows with a negative LLR: %.3f' % (X['fraction_of_negative_LLRs']))
    print('The calculation was done in %.3f sec.' % info['runtime'])

def confidence(LLR_stat,N,z):
    """ Binning is applied by aggregating the mean LLR of a window across N
        consecutive windows. The boundaries of the bins as well as the mean LLR
        and the standard-error per bin are returned. """
             
    K,M,V = tuple(LLR_stat.keys()), *zip(*LLR_stat.values())
            
    i = lambda j: j*(len(V)//N)
    f = lambda j: (j+1)*(len(V)//N) if j!=N-1 else len(K)
    
    x = lambda p,q: (K[p][0],K[q-1][-1])
    y = lambda p,q: statistics.mean(M[p:q])
    e = lambda p,q: z * std_of_mean(V[p:q]) 
    c = lambda p,q: [u<0 for u in M[p:q]].count(True) / (q-p)
    X,Y,E,C = ([func(i(j),f(j)) for j in range(N) if f(j)-i(j)>0] for func in (x,y,e,c))

    return X,Y,E,C

def bar_plot(info,pair,N,**kwargs):
    """ Plots the mean LLR vs. chromosomal position """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    save = kwargs.get('save', '')
    font_size=16

    # Create lists for the plot    
    X,Y,E,C = confidence(info,pair,N,z=1)
    T = [(x[1]+x[0])/2 for x in X]
    widths = [x[1]-x[0] for x in X]
    #X_boundaries = tuple(k for j in range(N+1) for k in (K[j*a][0], K[min((j+1)*a,len(K)-1)][-1]))
    X_ticks = [x[0] for x in X] + [X[-1][1]]
    X_labels = [('%.2f' % (j/chr_length(info['chr_id']))) for j in X_ticks] 
 
    # Define colormap
    bounds = [0.05*i for i in range(0,21)]  
    norm = mpl.colors.BoundaryNorm(bounds, 21)
    cmaplistBWR = plt.cm.seismic(bounds)
    my_cmap = mpl.colors.LinearSegmentedColormap.from_list('BWR', cmaplistBWR, 20)   
 
    # Build the plot
    fig, ax = plt.subplots(1, 1, figsize=(16, 9))  # setup the plot
    fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.07) 
    rects = ax.bar(T, Y, yerr=E, align='center', alpha=1, edgecolor='white', linewidth=3, capsize=10, width=widths, color=[*map(my_cmap,C)])
    ax.set_ylabel(f'log-likelihood ratio of {pair[0]:s}/{pair[1]:s}')
    ax.set_xlabel('Normalized chromosome position')
    ax.set_xticks(X_ticks)
    ax.set_xticklabels(X_labels)
    ax.set_title('.'.join(save.split('.')[:-2]))
    
    # Add colorbar
    sm = ScalarMappable(cmap=my_cmap, norm=norm)
    sm.set_array([])
    tick = [0.1*i for i in range(0,11)] 
    cbar = plt.colorbar(sm,ticks=tick, format='%.1f')
    cbar.set_label('Fraction of genomic windows with a negative LLR', fontsize=font_size, rotation=270,labelpad=25)
    
    # Add labels to bars
    height = lambda y,e: min(-0.01,y-1.2*y/abs(y)*e) if y>0 else max(0.01,y-1.2*y/abs(y)*e)
    for j in range(N):
        plt.text(T[j], height(Y[j],E[j]), '%.2f\u00B1%.2f'% (Y[j],E[j]), ha='center', va='center',color='black',fontsize=10)
    
    #ax.xaxis.grid(True)
    #ax.yaxis.grid(True)

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

def triple_plot(likelihoods,info,N,**kwargs):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    #from scipy.interpolate import interp1d
    
    save = kwargs.get('save', '')
    pairs = {('BPH','DISOMY'):'saddlebrown', ('DISOMY','SPH'):'darkgreen', ('SPH','MONOSOMY'):'lightskyblue'}

    _ = {};
    LLRs_per_genomic_window = {(i,j): {window:  mean_and_var([*starmap(LLR, ((_[i], _[j]) for _['MONOSOMY'], _['DISOMY'], _['SPH'], _['BPH'] in L))])
                       for window,L in likelihoods.items()} for i,j in pairs}
    
    fig,(ax1)=plt.subplots(1,1)
    fig.set_size_inches(9, 6, forward=True)    
    
    for p,LLR_stat in LLRs_per_genomic_window.items():
        X,Y,E,C = confidence(LLR_stat,N,z=1)
        T = [(x[1]+x[0])/2 for x in X]                
        ax1.plot([X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]],[i for i in Y for j in (1,2)], label=f'LLR of {p[0]:s} to {p[1]:s}',color=pairs[p],linewidth=2)
        ax1.errorbar(T, Y, yerr = E, ecolor='black',marker=None, ls='none',alpha=0.1) 

    ax1.set_ylabel('Log-likelihood ratio')
    ax1.set_xlabel('Normalized chromosome position')
    ax1.set_title(info['chr_id'])
    #Replace ticks along the x-axis 
    X_ticks = [X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]]
    X_labels = [('%.2f' % (j/chr_length(info['chr_id']))) for j in X_ticks] 
    ax1.set_xticks(X_ticks)
    ax1.set_xticklabels(X_labels)
    
    ax1.plot((lambda _: (min(_),max(_)))([j for i in X for j in i]),[0,0],color='red', ls='dotted',alpha=0.5)
    
    # get handles
    handles, labels = ax1.get_legend_handles_labels()
    # remove the errorbars
    #handles = [h[0] for h in handles]
    # use them in the legend
    #ax1.legend(handles, labels, loc='upper left',numpoints=1)
    
    ax1.legend(labels, loc='upper right',numpoints=1)
    
    if save!='':
        print('Saving plot...')
        ax1.set_title(save.rpartition('/')[-1].removesuffix('.png'))
        plt.savefig(save, dpi=150, facecolor='w',
                    edgecolor='w', orientation='landscape',
                    format='png', transparent=False, bbox_inches='tight',
                    pad_inches=0.1, metadata=None)
        plt.close(fig)
    else:
       plt.tight_layout()
       plt.show()
 
def panel_plot(DATA,N,**kwargs):
    import matplotlib as mpl
    mpl.rcParams.update({'figure.max_open_warning': 0})
    import matplotlib.pyplot as plt
    #from scipy.interpolate import interp1d
    
    save = kwargs.get('save', '')
    #pairs = {('BPH','DISOMY'):'saddlebrown', ('DISOMY','SPH'):'darkgreen', ('SPH','MONOSOMY'):'lightskyblue'}
    pairs = {('BPH','DISOMY'):'lightskyblue', ('BPH','SPH'):'darkorange'}

    #pairs = {('DISOMY','MONOSOMY'):'lightskyblue',}

    fig,axs =plt.subplots(4,6)
    fig.set_size_inches(36, 24, forward=True)    
    fig.suptitle(kwargs.get('title', ''), fontsize=16)
    
    AX = [i for j in axs for i in j]
    
    for ax1,(likelihoods,info) in zip(AX,DATA):
        _ = {};
        LLRs_per_genomic_window = {(i,j): {window:  mean_and_var([*starmap(LLR, ((_[i], _[j]) for _['MONOSOMY'], _['DISOMY'], _['SPH'], _['BPH'] in L))])
                       for window,L in likelihoods.items()} for i,j in pairs}
    
        for (p,LLR_stat) in LLRs_per_genomic_window.items():
            X,Y,E,C = confidence(LLR_stat,N,z=1)
            T = [(x[1]+x[0])/2 for x in X]                
            ax1.plot([X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]],[i for i in Y for j in (1,2)], label=f'LLR of {p[0]:s} to {p[1]:s}',color=pairs[p], linewidth=2)
            ax1.errorbar(T, Y, yerr = E, ecolor='black',marker=None, ls='none',alpha=0.1) 
    
        ax1.set_ylabel('Log-likelihood ratio')
        ax1.set_xlabel('Normalized chromosome position')
        ax1.set_title(info['chr_id'])
        #Replace ticks along the x-axis 
        X_ticks = [X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]]
        X_labels = [('%.2f' % (j/chr_length(info['chr_id']))) for j in X_ticks] 
        ax1.set_xticks(X_ticks)
        ax1.set_xticklabels(X_labels)
        
        ax1.plot((lambda _: (min(_),max(_)))([j for i in X for j in i]),[0,0],color='red', ls='dotted',alpha=0.5)
        
        # get handles
        handles, labels = ax1.get_legend_handles_labels()
        # remove the errorbars
        #handles = [h[0] for h in handles]
        # use them in the legend
        #ax1.legend(handles, labels, loc='upper left',numpoints=1)
        
        ax1.legend(labels, loc='upper right',numpoints=1)
        
    if save!='':
        print('Saving plot...')
        ax1.set_title(save.rpartition('/')[-1].removesuffix('.png'))
        plt.savefig(save, dpi=150, facecolor='w',
                    edgecolor='w', orientation='landscape',
                    format='png', transparent=False, bbox_inches='tight',
                    pad_inches=0.1, metadata=None)
        plt.close(fig)
    else:
       plt.tight_layout()
       plt.show()    

def detect_haploids_and_triploids(cases):
    ERRORS = []
    HAPLOIDS = []
    TRIPLOIDS = []
    pairs = (('MONOSOMY','DISOMY'),('BPH','DISOMY'))
    _ = {};
    for case in cases:
        if case['sp']=='AFR' or case['sp']=='AMR':
            continue
        LLRs_per_genomic_window = {p:{} for p in pairs}
        bam_filename = case['filename']
        red_flag=False
        #print(bam_filename)
        for chr_num in case['chr_num']:
            if chr_num=='X':
                continue
            try:
                chr_id = 'chr'+ str(chr_num)
                LLR_filename = re.sub('.bam$','',bam_filename.split('/')[-1]) + f'.{chr_id:s}.LLR.p'
                #show_info(LLR_filename,info,('BPH','SPH'))
                likelihoods, info = load_likelihoods(HOME+'/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/' + LLR_filename)
                if info['statistics']['num_of_windows']<20:
                    red_flag=True
                    continue
                for i,j in pairs:
                    LLRs_per_genomic_window[(i,j)][chr_id] = {window:  mean_and_var([*starmap(LLR, ((_[i], _[j]) for _['MONOSOMY'], _['DISOMY'], _['SPH'], _['BPH'] in L))])
                        for window,L in likelihoods.items()} 
                
            except Exception as error: 
                print(f'ERROR: {LLR_filename:s}', error)
                ERRORS.append((bam_filename.strip().rsplit('/',1).pop(),error))
                continue
        
        if not red_flag:            
            M = {pair: mean([m for c in C.values() for (m,v) in c.values()]) for pair,C in LLRs_per_genomic_window.items()}
            SD = {pair: std_of_mean([v for c in C.values() for (m,v) in c.values()]) for pair,C in LLRs_per_genomic_window.items()}
    
            if M[('BPH','DISOMY')]/SD[('BPH','DISOMY')]>1:
                A = {'mean': M[('BPH','DISOMY')], 'SD': SD[('BPH','DISOMY')] }
                print('Trsiomy:',bam_filename, A)
                case.update(A)
                TRIPLOIDS.append(case)
            elif M[('MONOSOMY','DISOMY')]/SD[('MONOSOMY','DISOMY')]>1:
                A = {'mean': M[('MONOSOMY','DISOMY')], 'SD': SD[('MONOSOMY','DISOMY')] }
                print('Monosomy:',bam_filename, A)
                case.update(A)
                HAPLOIDS.append(case)
                
         
    print(Counter(i[0] for i in ERRORS))
    return HAPLOIDS, TRIPLOIDS, ERRORS

if __name__ == "__main__":
    db_HAPLOIDS = [{'filename': '12751FA-AWL31_14.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12459FA-AU3UR_8.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '14338FA-CFHJ8_14.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '14715FA-CTGFC_21.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '14080FA-CBRH3_24.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '14260FA-CFPGD_19.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '14451FA-CFN45_14.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '12751FA-B5Y5C_14.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '11968FA-AP1KV_6.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '13885FA-C6JPJ_10.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12459FA-ARHR1_1.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13402FA-BK37C_22.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '12459FA-AU3UM_12.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '14417FA-CFMY3_11.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12456FA-ARHR1_16.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13474FA-BRCNL_11.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '13335FA-BK7V7_16.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '13198FA-BK2FN_12.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '11143FA-AK9EU_26.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '13327FA-BK7V7_21.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '12150FA-AW823_LB-8-AMP-17.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13213FA-BK2FN_20.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13213FA-BK2FN_2.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13297FA-B6TD2_12.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '11946FA-AR48V_9.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': 'RES-11733-AW727_LB-8N-AMP16.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13578FA-BRD32_12.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '11946FA-AR48V_3.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13597FA-C476J_35.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13659FA-BJKF9_10.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12073FA-ANTJ1_4.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12394FA-AU3TP_2.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13297FA-B6TD2_17.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12459FA-AU3UM_8.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '11257FA-ANDVJ_11.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '11826FA-AR21R_28.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12155FA-AW829_LB-7-AMP17.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12149FA-AW823_LB-16-AMP17.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12150FA-AW823_LB-8-AMP-17.bam', 'sp': 'EUR','chr_num': [*range(1,23)]+['X']}
                   ]
    
    db_TRIPLOIDS = [{'filename': '12431FA-B22P2_14.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '13809FA-C6JBM_9.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '12149FA-AW829_11-AMP16.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '12149FA-AW738_4.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '14149FA-CB9JP_14.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '14277FA-CFKCF_20.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '10405FA-AAAM9_2.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '12610FA-AV2K5_4.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '10667FA-AFFCM_2.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '12682FA-B8FT9_2.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '11288FA-ALNLD_9.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '13700FA-C4T36_9.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '12716FA-AYW5E_12.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '14236FA-CFJKV_1.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '14328FA-CFHJ8_14.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                     #{'filename': '12742FA-B8FHD_5.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '11510FA-APA12_5.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                     #{'filename': '13175FA-BJNV7_3.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '13175FA-B6RJM_3.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '13762FA-C4T43_3.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '10837FA-AFPAJ_1.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '12279FA-ATVDW_18.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '10771FA-AFPNC_3.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                     #{'filename': '13617FA-C22MH_9.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '12279FA-ATVDW_7.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '13845FA-C6JP9_10.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '14149FA-CB9JP_4.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                     #{'filename': '11657FA-AP1KJ_4.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                     #{'filename': '12073FA-AW12E_7.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                     #{'filename': '13441FA-BRMJT_20.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '13885FA-C6JPJ_8.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                     #{'filename': '13062FA-BJNTW_4.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '10882FA-AJ3U4_3.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                     {'filename': '12478FA-ARHL5_5.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                     #{'filename': '10851FA-AFP9Y_28.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}
                     ]
    
    db_Ydisomy = [{'filename': '12663FA-B8DPD_11.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                  {'filename': '13068FA-BK2G5_23.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']},
                  {'filename': '12405FA-AU3TR_3.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}]
    
    db_DIPLOID = [  {'filename': '10523FA-AFFRU_3.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10523FA-AFFRU_4.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10560FA-AFFPH_3.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10675FA-BJNTV_2.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10675FA-BJNTV_3.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10675FA-BJNTV_4.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10675FA-BJNTV_5.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10686FA-AFFPE_9.bam', 'sp': 'SAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10846FA-AFPAB_3.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10871FA-AJ3U4_12.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10951FA-AJ3WW_3.bam', 'sp': 'SAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10969FA-AJ470_2.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '11522FA-AP925_3.bam', 'sp': 'SAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '11578FA-AR3WC_3.bam', 'sp': 'AMR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '11598FA-AP923_9.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '12662FA-B5Y5R_1.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '12662FA-B5Y5R_3.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '12699FA-B8F4K_4.bam', 'sp': 'AMR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '12789FA-AWL1L_12.bam', 'sp': 'AFR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '12962FA-BK2G8_6.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '14529FA-CM2GK_2.bam', 'sp': 'AMR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': 'GP-CWFRM_8.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': 'MZ-AFFC4_1.bam', 'sp': 'AMR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '13068FA-BK2G5_23.bam', 'sp': 'AMR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10675FA-BJNTV_3c.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': 'MZ-AFFC4_2.bam', 'sp': 'SAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10964FA-AJ470_1.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '13121FA-BK23M_23.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '13086FA-BK2G5_6.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10668FA-AFDL2_2.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '11550FA-AP91V_5.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10722FA-AFFCT_3a.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '12055FA-ANTJ1_15.bam', 'sp': 'SAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '12454FA-AW7BB_2.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10967FA-AJ470_F10.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '11946FA-AR452_9.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '11550FA-AP91V_4.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '13744FA-C4RPY_2.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '13086FA-BK2G5_8.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '10658FA-AFDL2_F10.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '14220FA-CFP2Y_1.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '12446FA-AU3UC_29.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '14212FA-CFP2Y_5.bam', 'sp': 'AMR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '11946FA-AR452_1.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '11944FA-AR452_13.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                {'filename': '11511FA-AP91V_8.bam', 'sp': 'AMR', 'chr_num': [11,17,19,20,21,22]}]
    
    db_DIPLOID_noisy =  [{'filename': '10469FA-A9DBT_5.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '10560FA-AFFPH_5.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '11095FA-AK9FG_5.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '11422FA-APA42_4.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '11517FA-AP925_2.bam', 'sp': 'SAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '11594FA-AP923_15.bam', 'sp': 'AMR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '12853FA-AWKB8_1.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '12862FA-AWKB8_1.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '14529FA-CM2GK_3.bam', 'sp': 'AMR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '14529FA-CM2GK_4.bam', 'sp': 'AMR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '14953FA-CWLHB_6.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '14974FA-CWFV7_13.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '14974FA-CWFV7_14.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '14992FA-CWLYV_1.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '15022FA-J2MJB_11.bam', 'sp': 'AMR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '15042FA-J2JJH_4.bam', 'sp': 'AMR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': 'GP-CWFG9_1.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': 'GP-CWGJL_3.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': 'GP-CWGJL_6.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '11946FA-AR452_7.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': 'RES-11892-AU3TL_5T.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '12446FA-AU3UC_32.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '10722FA-AFFCT_5.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '13068FA-BK2G5_27.bam', 'sp': 'AMR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '11944FA-AR452_15.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '11557FA-AP2W6_4.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '12446FA-AU3UC_19.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '11557FA-AP2W6_1.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)}, 
                     {'filename': '11143FA-AL3GC_27.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '11494FA-AP9FW_6.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '14217FA-CFP2Y_3.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '11180FA-AK9G0_8.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '13126FA-B6RLR_12944FA-1.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '13948FA-C7R7Y_23.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '12881FA-AWK94_5.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                     {'filename': '12446FA-AU3UC_28.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)}]
    
    db_TRIPLOIDS_maybe = [{'filename': '13331FA-BK8NW_2.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14096FA-CBRH3_6.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11749FA-AR3WC_3.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14088FA-C7T6M_15.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12431FA-B22P2_14.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12279FA-ATVDW_6.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13809FA-C6JBM_9.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14119FA-CB9JP_2.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12149FA-AW829_11-AMP16.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12149FA-AW738_4.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14149FA-CB9JP_14.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11494FA-AP9FW_2.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11624FA-ANTHL_6.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13718FA-C4T43_8.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14277FA-CFKCF_20.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14767FA-CTFFN_17.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11270-RES-C7RGH_1T1.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13528FA-BRCJ2_7.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10672FA-AFFPR_7.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13642FA-BJKDV_4.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14473FA-CHTK6_25.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10405FA-AAAM9_2.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12610FA-AV2K5_4.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10667FA-AFFCM_2.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12682FA-B8FT9_2.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11288FA-ALNLD_9.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13700FA-C4T36_9.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11627FA-ANTHL_6.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13809FA-C6JBM_10.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13515FA-BRCJ2_7.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13534FA-BRCJ2_7.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11232FA-AK90E_5.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12279FA-AY8J4_6.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12716FA-AYW5E_12.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12103FA-AW60G_9.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11171FA-AK9FM_19.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14236FA-CFJKV_1.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14328FA-CFHJ8_14.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13564FA-BRCN9_11.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13098FA-B6MTL_3.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12742FA-B8FHD_5.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11510FA-APA12_5.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12566FA-AV2NV_1.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13175FA-BJNV7_3.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13175FA-B6RJM_3.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13762FA-C4T43_3.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10864FA-AFPAN_1.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14444FA-CFN45_17.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10837FA-AFPAJ_1.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14731FA-CTGFC_1.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12103FA-AW60G_5.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12279FA-ATVDW_18.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10771FA-AFPNC_3.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11270FA-ALNKA_21.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11947FA-AR45Y_8.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11070FA-AJY20_3.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11573FA-AP9A4_4.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13617FA-C22MH_9.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14570FA-CLPKK_6-10925.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11826FA-AR21R_25.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13288FA-BJNVW_6.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14473FA-CHTK6_26.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12568FA-AV2NV_2.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11021FA-AJK7A_22.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10635FA-AFDK3_13.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12289FA-AW60A_26.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12279FA-ATVDW_7.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10671FA-AFFPE_24.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12258FA-AU3DR_9.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11270FA-ALNKA_11.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13692FA-C4VK7_15.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11261FA-AL2TY_10.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10417FA-AAALW_7.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11913FA-APHAJ_14.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10749FA-AFP9W_14.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10940FA-AJ47P_5.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12831FA-AWK85_10.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11279FA-ALNLD_23.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13845FA-C6JP9_10.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14149FA-CB9JP_4.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11541FA-AP99W_1.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11657FA-AP1KJ_4.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13699FA-C4T36_1.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12073FA-AW12E_7.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12946FA-AV39W_12.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11059FA-AJY20_27.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13441FA-BRMJT_20.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13885FA-C6JPJ_8.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13143FA-BJNV7_32.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13025FA-BJNV6_4.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13062FA-BJNTW_4.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12258FA-AU3DR_7.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12478FA-ARHL5_3.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12896FA-AWK94_1.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10882FA-AJ3U4_3.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12478FA-ARHL5_5.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13510FA-C4BW2_2.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13334FA-BK8HK_3.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10851FA-AFP9Y_28.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11947FA-AP1GW_10.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}]


    #for i in range(1,23):
        #likelihoods, info = load_likelihoods(f'/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/12751FA-AWL31_14.chr{i:d}.LLR.p')
        #bar_plot(info,('BPH','SPH'),N=10)
        #triple_plot(info,N=10,save=f'/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/12751FA-AWL31_14.chr{i:d}.png')
    
    #DATA = [load_likelihoods(f'/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/12751FA-AWL31_14.chr{i:d}.LLR.p') for i in range(1,23)]
    #panel_plot(DATA,N=10,title='12751FA-AWL31_14')
    #likelihoods, info = load_likelihoods('/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/12751FA-AWL31_14.chr7.LLR.p')
    #triple_plot(likelihoods,info,N=10)
    
    #bob = ['10523FA-AFFRU_3', '10523FA-AFFRU_4', '10560FA-AFFPH_3', '10675FA-BJNTV_2', '10675FA-BJNTV_3', '10675FA-BJNTV_4', '10675FA-BJNTV_5', '10686FA-AFFPE_9', '10846FA-AFPAB_3', '10871FA-AJ3U4_12', '10951FA-AJ3WW_3', '10969FA-AJ470_2', '11522FA-AP925_3', '11578FA-AR3WC_3', '11598FA-AP923_9', '12662FA-B5Y5R_1', '12662FA-B5Y5R_3', '12699FA-B8F4K_4', '12789FA-AWL1L_12', '12962FA-BK2G8_6', '14529FA-CM2GK_2', 'GP-CWFRM_8', 'MZ-AFFC4_1', '13068FA-BK2G5_23', '10675FA-BJNTV_3c', 'MZ-AFFC4_2', '10964FA-AJ470_1', '13121FA-BK23M_23', '13086FA-BK2G5_6', '10668FA-AFDL2_2', '11550FA-AP91V_5', '10722FA-AFFCT_3a', '12055FA-ANTJ1_15', '12454FA-AW7BB_2', '10967FA-AJ470_F10', '11946FA-AR452_9', '11550FA-AP91V_4', '13744FA-C4RPY_2', '13086FA-BK2G5_8', '10658FA-AFDL2_F10', '14220FA-CFP2Y_1', '12446FA-AU3UC_29', '14212FA-CFP2Y_5', '11946FA-AR452_1', '11944FA-AR452_13', '11511FA-AP91V_8']
    
    
    with open(HOME+'/ariad/Dropbox/postdoc_JHU/BlueFuse/Play/diploid_females.p', 'rb') as f:
        db_TEST = pickle.load(f)
        #db_TEST = [i for i in db_TEST if 'AFR'!=i['sp']!='AMR']K = [i for i in db_TEST if 'AFR'!=i['sp']!='AMR']
    
    result = detect_haploids_and_triploids(db_TEST)
    
    #for case in db_TRIPLOIDS:
    #    print(case)
    #    DATA = [load_likelihoods(f"/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/{case['filename'][:-4]:s}.chr{i:s}.LLR.p") for i in [*map(str,range(1,23))]+['X']]
    #    panel_plot(DATA,N=10,title=f"{case['filename'][:-4]:s}")