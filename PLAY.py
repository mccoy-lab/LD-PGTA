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
from collections import Counter, defaultdict

HOME = '/home' #  '/Users' 

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
    colors = {frozenset(('BPH','DISOMY')):'saddlebrown',
              frozenset(('DISOMY','SPH')):'darkgreen',
              frozenset(('SPH','MONOSOMY')):'lightskyblue',
              frozenset(('DISOMY','MONOSOMY')):'lightskyblue',
              frozenset(('BPH','SPH')):'darkorange'}

    fig,axs =plt.subplots(4,6)
    fig.set_size_inches(36, 24, forward=True)    
    fig.suptitle(kwargs.get('title', ''), fontsize=16)
    
    AX = [i for j in axs for i in j]
    
    classification = {(True,True,True): (('BPH','SPH'),),
                      (False,True,True): (('DISOMY','BPH'),('DISOMY','MONOSOMY')),
                      (False,False,True): (('SPH','BPH'),('SPH','MONOSOMY')),
                      (False,False,False): (('SPH','MONOSOMY'),)}

    pairs = (('BPH','DISOMY'), ('DISOMY','SPH'), ('SPH','MONOSOMY'), ('BPH','SPH'))

    default_classification = defaultdict(lambda: pairs, classification)
    
    for ax1,(likelihoods,info) in zip(AX,DATA):
        _ = {};
        LLRs_per_genomic_window = {(i,j): {window:  mean_and_var([*starmap(LLR, ((_[i], _[j]) for _['MONOSOMY'], _['DISOMY'], _['SPH'], _['BPH'] in L))])
                       for window,L in likelihoods.items()} for i,j in pairs}
    
        A = {p: confidence(LLR_stat,N,z=1) for (p,LLR_stat) in LLRs_per_genomic_window.items()}
        key = (*(sum(A[p][1])>0 for p in pairs[:3]),) #2*sum(A[p][2])**.5
        for (a,b) in default_classification[key]:
            LLR_stat = {window:  mean_and_var([*starmap(LLR, ((_[a], _[b]) for _['MONOSOMY'], _['DISOMY'], _['SPH'], _['BPH'] in L))])
                       for window,L in likelihoods.items()}
            X,Y,E,C = confidence(LLR_stat,N,z=1)
            T = [(x[1]+x[0])/2 for x in X]                
            ax1.plot([X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]],[i for i in Y for j in (1,2)], label=f'LLR of {a:s} to {b:s}',color=colors[frozenset((a,b))], linewidth=2)
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
    pairs = (('MONOSOMY','DISOMY'),('DISOMY','BPH'))
    _ = {};
    for case in cases:
        #if case['sp']=='AFR' or case['sp']=='AMR':
        #    continue
        LLRs = {}
        
        bam_filename = case['filename']
        #print(bam_filename)
        for chr_num in case['chr_num']:
            if chr_num=='X':
                continue
            try:
                chr_id = 'chr'+ str(chr_num)
                LLR_filename = re.sub('.bam$','',bam_filename.split('/')[-1]) + f'.{chr_id:s}.LLR.p'
                #show_info(LLR_filename,info,('BPH','SPH'))
                likelihoods, info = load_likelihoods(HOME+'/ariad/Dropbox/postdoc_JHU/LD-PGTA_ecosystem/LD-PGTA_V2/results_ZOUVES/' + LLR_filename)
                if 'ancestry' not in info:
                    print('bob!')
                ####if info['depth'] < 0.007: 
                    #print(f"{info['depth']:f}")
                    #print(f'skipped {bam_filename:s} due to low coverage')
                    #continue
                if not info['statistics']['num_of_windows']:
                    #print(f'skipped {bam_filename:s} due to low coverage')
                    continue
                else:
                    LLRs[chr_id] = likelihoods
            except Exception as error: 
                print(f'ERROR: {LLR_filename:s}', error)
                ERRORS.append((bam_filename.strip().rsplit('/',1).pop(),error))
                continue
        
        if len(LLRs)<22:
        #    #print(len(LLRs))
            continue
        
        #if sum(1 for likelihoods in LLRs.values() for L in likelihoods.values())<1000:
        #    continue
        
        LLRs_aggregation = {(i,j): 
                                [mean_and_var([*starmap(LLR, ((_[i], _[j]) for _['MONOSOMY'], _['DISOMY'], _['SPH'], _['BPH'] in L))])
                                     for likelihoods in LLRs.values() for L in likelihoods.values()] 
                                         for i,j in pairs}
        
        if len(LLRs_aggregation[pairs[0]])<100: continue
                
        M = {pair: mean([m for (m,v) in C]) for pair,C in LLRs_aggregation.items()}
        SD = {pair: std_of_mean([v for (m,v) in C]) for pair,C in LLRs_aggregation.items()}

        if M[('DISOMY','BPH')]/SD[('DISOMY','BPH')]<2 and len(LLRs_aggregation[pairs[0]])>=1000:
            A = {'mean': M[('DISOMY','BPH')], 'SD': SD[('DISOMY','BPH')] }
            print('Triploid',bam_filename.replace('.bam',''), M[('DISOMY','BPH')], SD[('DISOMY','BPH')], case['sp'], len(LLRs_aggregation[pairs[0]]))
            case.update(A)
            TRIPLOIDS.append(case)
        elif M[('MONOSOMY','DISOMY')]/SD[('MONOSOMY','DISOMY')]>1:
            A = {'mean': M[('MONOSOMY','DISOMY')], 'SD': SD[('MONOSOMY','DISOMY')] }
            print('Haploid',bam_filename.replace('.bam',''), M[('MONOSOMY','DISOMY')], SD[('MONOSOMY','DISOMY')], case['sp'], len(LLRs_aggregation[pairs[0]]))
            case.update(A)
            HAPLOIDS.append(case)
        
         
    print(Counter(i[0] for i in ERRORS))
    return HAPLOIDS, TRIPLOIDS, ERRORS

def detect_haploids_and_triploids2(cases):
    ERRORS = []
    HAPLOIDS = []
    TRIPLOIDS = []
    
    for case in cases:
        if case['sp']=='AFR' or case['sp']=='AMR':
            continue
        buffer = []
        bam_filename = case['filename']
        for chr_num in case['chr_num']:
            try:
                chr_id = 'chr'+ str(chr_num)
                LLR_filename = re.sub('.bam$','',bam_filename.split('/')[-1]) + f'.{chr_id:s}.LLR.p'
                likelihoods, info = load_likelihoods('/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/' + LLR_filename)
                #if chr_id!='chrX' and info['statistics']['num_of_windows']<20: raise Exception('Number of genomic windows is below 20')
                if info['statistics']['num_of_windows'] and chr_id!='chrX':
                    #show_info(LLR_filename,info,('SPH','MONOSOMY'))
                    S = info['statistics']['LLRs_per_chromosome']
                    a =  S[('SPH','MONOSOMY')]['mean']/S[('SPH','MONOSOMY')]['std_of_mean']<-1 and S[('BPH','DISOMY')]['mean']/S[('BPH','DISOMY')]['std_of_mean']<-1
                    b =  S[('SPH','MONOSOMY')]['mean']/S[('SPH','MONOSOMY')]['std_of_mean']>1 and S[('BPH','DISOMY')]['mean']/S[('BPH','DISOMY')]['std_of_mean']>1
                    buffer.append(b-a)
            except Exception as error: 
                print(f'ERROR: {LLR_filename:s} ---', error)
                ERRORS.append((bam_filename.strip().split('/')[-1],error))
                continue
        if sum(buffer)<-1:
            print('haploid:',buffer,case['sp'])
            HAPLOIDS.append(case)
        if sum(buffer)>+1:            
            print('triploid:',buffer,case['sp'])
            TRIPLOIDS.append(case)
    print(Counter(i[0] for i in ERRORS))
    return HAPLOIDS, TRIPLOIDS, ERRORS

def about(cases):
    for case in cases:
        bam_filename = case['filename']
        for i in range(1,23):
            chr_id = f'chr{i:d}'
            LLR_filename = re.sub('.bam$','',bam_filename.split('/')[-1]) + f'.{chr_id:s}.LLR.p'
            likelihoods, info = load_likelihoods('/home/ariad/Dropbox/postdoc_JHU/LD-PGTA_ecosystem/LD-PGTA_V2/results_ZOUVES/' + LLR_filename)
            show_info(LLR_filename,info,('SPH','MONOSOMY'))

if __name__ == "__main__":
    #for i in range(1,23):
        #likelihoods, info = load_likelihoods(f'/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/12751FA-AWL31_14.chr{i:d}.LLR.p')
        #bar_plot(info,('BPH','SPH'),N=10)
        #triple_plot(info,N=10,save=f'/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/12751FA-AWL31_14.chr{i:d}.png')
    
    #identifier = '13094FA-B6MTL_3'
    #DATA = [load_likelihoods(f'/home/ariad/Dropbox/postdoc_JHU/LD-PGTA_ecosystem/LD-PGTA_V2/results_ZOUVES/{identifier:s}.chr{i:d}.LLR.p') for i in range(1,23)]
    #panel_plot(DATA,N=10,title=f'{identifier:s}')
    
    #likelihoods, info = load_likelihoods('/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/12751FA-AWL31_14.chr7.LLR.p')
    #triple_plot(likelihoods,info,N=10)
    
    #bob = ['10523FA-AFFRU_3', '10523FA-AFFRU_4', '10560FA-AFFPH_3', '10675FA-BJNTV_2', '10675FA-BJNTV_3', '10675FA-BJNTV_4', '10675FA-BJNTV_5', '10686FA-AFFPE_9', '10846FA-AFPAB_3', '10871FA-AJ3U4_12', '10951FA-AJ3WW_3', '10969FA-AJ470_2', '11522FA-AP925_3', '11578FA-AR3WC_3', '11598FA-AP923_9', '12662FA-B5Y5R_1', '12662FA-B5Y5R_3', '12699FA-B8F4K_4', '12789FA-AWL1L_12', '12962FA-BK2G8_6', '14529FA-CM2GK_2', 'GP-CWFRM_8', 'MZ-AFFC4_1', '13068FA-BK2G5_23', '10675FA-BJNTV_3c', 'MZ-AFFC4_2', '10964FA-AJ470_1', '13121FA-BK23M_23', '13086FA-BK2G5_6', '10668FA-AFDL2_2', '11550FA-AP91V_5', '10722FA-AFFCT_3a', '12055FA-ANTJ1_15', '12454FA-AW7BB_2', '10967FA-AJ470_F10', '11946FA-AR452_9', '11550FA-AP91V_4', '13744FA-C4RPY_2', '13086FA-BK2G5_8', '10658FA-AFDL2_F10', '14220FA-CFP2Y_1', '12446FA-AU3UC_29', '14212FA-CFP2Y_5', '11946FA-AR452_1', '11944FA-AR452_13', '11511FA-AP91V_8']
    
    
    with open(HOME+'/ariad/Dropbox/postdoc_JHU/Zouves-BlueFuse/Play/diploids.p', 'rb') as f:
         db_TEST = pickle.load(f)
         #db_TEST = [i for i in db_TEST if 'AFR'!=i['sp']!='AMR'] 
    
    HAPLOIDS, TRIPLOIDS, ERRORS = detect_haploids_and_triploids(db_TEST)
    
    #for case in db_TRIPLOIDS:
    #    print(case)
    #    DATA = [load_likelihoods(f"/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/{case['filename'][:-4]:s}.chr{i:s}.LLR.p") for i in [*map(str,range(1,23))]+['X']]
    #    panel_plot(DATA,N=10,title=f"{case['filename'][:-4]:s}")
    
    
    #####
    #RPL#
    #####
    
    RPL = ['2019-12-23_RIG003-2_S11.bam', '2020-05-22_RPL005-4_S8.bam',
           '2019-12-23_RIG003-3_S4.bam',   '2020-05-22_RPL006-1_S9.bam',
           '2019-12-23_RIG003-4_S8.bam',   '2020-05-22_RPL006-2_S10.bam',
           '2020-05-22_RPL004-1_S1.bam',   '2020-05-22_RPL006-3_S11.bam',
           '2020-05-22_RPL004-2_S2.bam',   '2020-05-22_RPL006-4_S12.bam',
           '2020-05-22_RPL004-3_S3.bam',   '2020-05-22_RPL007-1_S13.bam',
           '2020-05-22_RPL004-4_S4.bam',   '2020-05-22_RPL007-2_S14.bam',
           '2020-05-22_RPL005-1_S5.bam',   '2020-05-22_RPL007-3_S15.bam',
           '2020-05-22_RPL005-2_S6.bam',   '2020-05-22_RPL007-4_S16.bam',
           '2020-05-22_RPL005-3_S7.bam']
    
    #for R in RPL:
    #    DATA = [load_likelihoods(f"{HOME:s}/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_RPL/{R[:-4]:s}.chr{i:s}.LLR.p") for i in [*map(str,range(1,23))]]
    #    panel_plot(DATA,N=10,title=f"{R[:-4]:s}",save=f"{R[:-4]:s}.png")
    #cases = [{'filename': i ,'sp': 'EUR' , 'chr_num': [*map(str,range(1,23)),'X'] } for i in RPL]
    #detect_haploids_and_triploids(cases)    