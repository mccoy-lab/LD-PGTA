#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PLOT_PANEL

Plots log-likelihood ratio vs. chromosomal position from a LLR file.

May 20, 2020
"""

import pickle, bz2, gzip, collections, math
from statistics import mean, variance
from math import log
from operator import attrgetter, itemgetter
from itertools import accumulate
import argparse, sys

likelihoods_tuple = collections.namedtuple('likelihoods_tuple', ('monosomy', 'disomy', 'SPH', 'BPH'))

def chr_length(chr_id):
    """ Return the chromosome length for a given chromosome, based on the reference genome hg38.""" 
    #The data of chromosome length was taken from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38
    length_dict = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
                  'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422,
                  'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
                  'chr16': 90338345, 'chr17':  83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
                  'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}
    return length_dict[chr_id]

def mean_and_var(x):
    """ Calculates the mean and variance. """
    cache = tuple(x)
    m = mean(cache)
    var = variance(cache, xbar=m)
    return m, var

def mean_and_std_of_mean_of_rnd_var(A):
    """ Calculates the mean and sample standard deviation of the mean of random variables.
        Each row of A represents a random variable, with observations in the columns."""
    if type(A)==dict:
        A = tuple(tuple(i) for i in A.values()) 
    
    M, N = len(A), len(A[0]) ### N is the number of random variables, while M is the number of samples.
    mu = sum(sum(likelihoods_in_window) for likelihoods_in_window in A) / N
    arg = ((sum(sampled_likelihoods) - mu)**2 for sampled_likelihoods in zip(*A))
    std = (sum(arg) / (N - 1))**.5 / M
    mean = mu / M
    return mean, std

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
    log-likelihood ratios (LLRs). """
    
    Open = {'bz2': bz2.open, 'gzip': gzip.open}.get(filename.rpartition('.')[-1], open)
    
    with Open(filename, 'rb') as f:
        likelihoods = pickle.load(f)
        info = pickle.load(f)
    return likelihoods, info

def show_info(filename, info, pairs=(('BPH','SPH'),)):
    S = info['statistics']
    ancestral_makeup = ", ".join("{:.1f}% {}".format(100*v, k) for k, v in info['ancestral_makeup'].items()) if type(info['ancestral_makeup'])==dict else ', '.join(info['ancestral_makeup'])
    matched_alleles = ", ".join("{}: {:.1f}%".format(k,100*v) for k, v in info['statistics']['matched_alleles'].items())

    print('\nFilename: %s' % filename)
    print('\nSummary statistics')
    print('------------------')    
    print('Chromosome ID: %s, Depth: %.2f.' % (info['chr_id'],info['depth']))
    print('Number of genomic windows: %d, Mean and standard error of genomic window size: %d, %d.' % (S.get('num_of_windows',0),S.get('window_size_mean',0),S.get('window_size_std',0)))
    print('Mean and standard error of meaningful reads per genomic window: %.1f, %.1f.' % (S.get('reads_mean',0), S.get('reads_std',0)))
    print('Ancestral makeup: %s, Fraction of alleles matched to the reference panel: %s.' % (ancestral_makeup, matched_alleles))

    for pair in pairs:
        if 'LLRs_per_chromosome' in S and tuple(pair) in S['LLRs_per_chromosome']:
            L = S['LLRs_per_chromosome'][tuple(pair)]
            print(f"--- Chromosome-wide LLR between {pair[0]:s} and {pair[1]:s} ----")
            print(f"Mean LLR: {L['mean_of_mean']:.3f}, Standard error of the mean LLR: {L['std_of_mean']:.3f}")
            print(f"Fraction of genomic windows with a negative LLR: {L['fraction_of_negative_LLRs']:.3f}")

def bin_genomic_windows(windows,chr_id,num_of_bins):
    """ Lists the bins and gives the genomic windows that they contain. """
    bin_size = chr_length(chr_id) / num_of_bins
    result = {}
    j = 0
    
    for i in range(num_of_bins): ### All bins before the first the genomic window are filled with Nones.
        if sum(windows[0])/2 < (i+1)*bin_size:
            break
        result[i/num_of_bins,(i+1)/num_of_bins] = None
    
    for k,(a,b) in enumerate(windows):
        if not bin_size*i <= (a+b)/2 < bin_size*(i+1):
            result[i/num_of_bins,(i+1)/num_of_bins] = (j,k)
            j = k
            for i in range(i+1,num_of_bins): #Proceed to the next non-empty bin; Empty bins are filled with Nones.
                if (a+b)/2 < (i+1)*bin_size:
                    break
                result[i/num_of_bins,(i+1)/num_of_bins] = None
    
    for i in range(i,num_of_bins): ### All bins after the last the genomic window are filled with Nones.
        result[i/num_of_bins,(i+1)/num_of_bins] = (j,k) if j != k else None
        j = k 
    return result

def binning(LLRs_per_window,info,num_of_bins):
    """ Genomic windows are distributed into bins. The LLRs in a genomic windows
    are regarded as samples of a random variable. Within each bin, we calculate
    the mean and population standard deviation of the mean of random variables. 
    The boundaries of the bins as well as the mean LLR and the standard-error
    per bin are returned. """
             
    #K,M,V = tuple(LLR_stat.keys()), *zip(*LLR_stat.values())
    list_of_windows = [*LLRs_per_window.keys()]
    bins = bin_genomic_windows(list_of_windows, info['chr_id'], num_of_bins)
    X = [*bins]
    
    LLR_matrix = [*LLRs_per_window.values()]
    Y, E = [], []
    for C in bins.values():
        if C:
            mean, std = mean_and_std_of_mean_of_rnd_var(LLR_matrix[C[0]:C[1]])
        else:
            mean, std = None, None
        
        Y.append(mean)
        E.append(std)
    
    return X,Y,E

def detect_crossovers(genomic_windows, mean_of_LLRs, variance_of_LLRs, z_score=1.96, lookahead=20):
    """ Detecting crossovers by indetifying transitions between BPH and SPH 
        regions. """
        
    crossovers = {}
    #Scan the chromosome in the 5'-to-3' direction to find crossovers.
    x_coord = tuple(0.5*(a+b) for a,b in genomic_windows)
    acc_means = tuple(accumulate(mean_of_LLRs))
    acc_vars = tuple(accumulate(variance_of_LLRs))
    triple = tuple(zip(x_coord, acc_means, acc_vars))
    
    #maxima and minima candidates are temporarily stored in mx and mn, respectively.
    mx, mn, last_ind = None, None, 0
        
    for index, (x, y, v) in enumerate(triple):
        
        if  mx==None or y > mx:
            mx_index,mx_pos,mx,mx_var = index,x,y,v
        
        if  mn==None or y < mn:
            mn_index,mn_pos,mn,mn_var = index,x,y,v

        if mx!=None and 0 < (mx-y)-z_score*(v-mx_var)**.5 and index-mx_index>=lookahead:
            for x2, y2, v2 in triple[max(mx_index-lookahead,last_ind):last_ind:-1]:
                if 0 < (mx-y2)-z_score*(mx_var-v2)**.5:
                    kappa = min((mx-y2)/(mx_var-v2)**.5,(mx-y)/(v-mx_var)**.5)
                    crossovers[mx_pos] = kappa
                    mx, mn, last_ind = None, None, mx_index # set algorithm to find the next minima
                    break
                
        if mn!=None and 0 < (y-mn)-z_score*(v-mn_var)**.5  and index-mn_index>=lookahead:    
            for x2, y2, v2 in triple[max(mn_index-lookahead,last_ind):last_ind:-1]:
                if  0 < (y2-mn)-z_score*(mn_var-v2)**.5:
                    kappa = min((y2-mn)/(mn_var-v2)**.5,(y-mn)/(v-mn_var)**.5)
                    crossovers[mn_pos] = kappa
                    mx, mn, last_ind = None, None, mn_index # set algorithm to find the next maxima
                    break
                
    return crossovers

def detect_crossovers_v2(genomic_windows, mean_of_LLRs, variance_of_LLRs, z_score=1.96, lookahead=20):
    """ Detecting crossovers by indetifying transitions between BPH and SPH 
        regions. """
        
    crossovers = {}
    #Scan the chromosome in the 5'-to-3' direction to find crossovers.
    x_coord = tuple(0.5*(a+b) for a,b in genomic_windows)
    acc_means = tuple(accumulate(mean_of_LLRs))
    acc_vars = tuple(accumulate(variance_of_LLRs))
    triple = tuple(zip(x_coord, acc_means, acc_vars))
    
    #Maxima and minima candidates are temporarily stored in mx and mn, respectively.
    mx, mn, last_ind, last_extremum = None, None, 0, 0
        
    for index, (x, y, v) in enumerate(triple):
        
        if  mx==None or y > mx:
            mx_index,mx_pos,mx,mx_var = index,x,y,v
        
        if  mn==None or y < mn:
            mn_index,mn_pos,mn,mn_var = index,x,y,v

        if mx!=None and 0 < (mx-y)-z_score*(v-mx_var)**.5 and index-mx_index>=lookahead: #maximum point detected
            for x2, y2, v2 in triple[max(mx_index-lookahead,last_ind):last_ind:-1]:
                if 0 < (mx-y2)-z_score*(mx_var-v2)**.5:
                    if last_extremum == +1: #Recovering skipped minimum point.
                        x0, mx0, mx_var0 = triple[last_ind]    
                        x1, mn1, mn_var1 = min(triple[last_ind+1:mx_index],key=itemgetter(1))
                        kappa = min((mx0-mn1)/(mn_var1-mx_var0)**.5 , (mx-mn1)/(mx_var-mn_var1)**.5 )
                        crossovers[x1] = kappa

                    kappa = min((mx-y2)/(mx_var-v2)**.5,(mx-y)/(v-mx_var)**.5)
                    crossovers[mx_pos] = kappa
                    mx, mn, last_ind, last_extremum = None, None, mx_index, +1 # set algorithm to find the next minimum
                    break
                
        if mn!=None and 0 < (y-mn)-z_score*(v-mn_var)**.5  and index-mn_index>=lookahead: #minimum point detected
            for x2, y2, v2 in triple[max(mn_index-lookahead,last_ind):last_ind:-1]:
                if  0 < (y2-mn)-z_score*(mn_var-v2)**.5:
                    if last_extremum == -1: #Recovering skipped maximum point.
                        x0, mn0, mn_var0 = triple[last_ind]    
                        x1, mx1, mx_var1 = max(triple[last_ind+1:mn_index],key=itemgetter(1))
                        kappa = min((mx1-mn0)/(mx_var1-mn_var0)**.5 , (mx1-mn)/(mn_var-mx_var1)**.5 )
                        crossovers[x1] = kappa
                        
                    kappa = min((y2-mn)/(mn_var-v2)**.5,(y-mn)/(v-mn_var)**.5)
                    crossovers[mn_pos] = kappa
                    
                    mx, mn, last_ind, last_extremum = None, None, mn_index, -1 # set algorithm to find the next maxima
                    break
                
    return crossovers

def capitalize(x):
    return x[0].upper() + x[1:]
    
def panel_plot(DATA,**kwargs):
    """ Creates a multi-panel figure. For each numbered chromosome, a figure 
        depicts the log-likelihood ratio vs. chromosomal position for BPH over
        SPH. """
    
    import matplotlib as mpl
    mpl.rcParams.update({'figure.max_open_warning': 0})
    
    
    scale = kwargs.get('scale', 0.5)
    save = kwargs.get('save', '')
    z_score = kwargs.get('z_score', 1.96)
    lookahead = kwargs.get('lookahead', 30)

    bin_size = kwargs.get('bin_size', 4000000)
    pairs = kwargs.get('pairs', (('BPH','disomy'),('disomy','SPH'),('SPH','monosomy')))
    
    fs=28 * scale
    columns = 6
    rows = math.ceil(len(DATA)/columns)
        
    if save!='':
            mpl.use('Agg')
    else:
        #['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']
        mpl.use('Qt5Agg')
    
    
    import matplotlib.pyplot as plt
    num_of_bins = {'chr'+str(i): chr_length('chr'+str(i))//bin_size for i in [*range(1,23)]+['X','Y']}

    colors = {frozenset(('BPH','disomy')):(177/255,122/255,162/255),
              frozenset(('disomy','SPH')):(242/255,142/255,44/255),
              frozenset(('SPH','monosomy')):(239/255,106/255,92/255),
              frozenset(('disomy','monosomy')):(104/255,162/255,183/255),
              frozenset(('BPH','SPH')):(104/255,162/255,104/255)}

    if len(DATA)>columns:
        fig,axs = plt.subplots(rows ,columns, sharex='col', sharey='row', figsize=(6.666 * columns * scale, 5.625 * rows * scale))
        fig.subplots_adjust(left=0.05, bottom=0.1, right=.99, top=(0.92 if kwargs.get('title',None) else 0.96), wspace=None, hspace=None)
    else:
        fig,axs = plt.subplots(rows ,columns, sharex='none', sharey='row', figsize=( 6.666 * columns * scale, 1.25 * 5.625 * rows * scale))
        fig.subplots_adjust(left=0.05, bottom=0.3, right=.99, top=(0.82 if kwargs.get('title',None) else 0.86), wspace=None, hspace=None)
    
    
    AX = [i for j in axs for i in j] if len(DATA)>columns else axs
    
    H = {}
    YMAX = [0]*len(DATA)
    crossovers = {}
    for a,b in pairs:
        for g,(ax1,(likelihoods,info)) in enumerate(zip(AX,DATA.values())):
    

            LLRs = {window: tuple(LLR(attrgetter(a)(l), attrgetter(b)(l)) for l in likelihoods_in_window)
                               for window,likelihoods_in_window in likelihoods.items()} 
                                                               
            X,Y,E = binning(LLRs,info,num_of_bins[info['chr_id']])
            Y = [(y if y else 0) for y in Y]
            E = [(z_score*e if e else 0) for e in E]
            
            T = [(x[1]+x[0])/2 for x in X]                
            steps_x = [X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]]
            steps_y = [i for i in Y for j in (1,2)]
            H[a,b] = ax1.plot(steps_x, steps_y, label=f'{capitalize(a):s} vs. {capitalize(b):s}' ,color=colors[frozenset((a,b))], linewidth=2, zorder=10, scalex=True, scaley=True, alpha=0.8)
            
            P = [(x[1]-x[0])/2 for x in X]                
            ax1.errorbar(T, Y, xerr = P, ecolor=colors[frozenset((a,b))],marker=None, ls='none',alpha=1, zorder=13, linewidth=5*scale) 
            ax1.errorbar(T, Y, yerr = E, ecolor='black',marker=None, ls='none',alpha=0.2, zorder=15, linewidth=4*scale) 
            
            yabsmax = max(map(abs,Y))
            
            if pairs==(('BPH','SPH'),) or pairs==(('SPH','BPH'),):
                genomic_windows = info['statistics']['LLRs_per_genomic_window'][('BPH','SPH')]
                mean_of_LLRs = [*map(itemgetter(0),genomic_windows.values())]
                variance_of_LLRs = [*map(itemgetter(1),genomic_windows.values())]
                unnormalized_crossovers = detect_crossovers_v2(genomic_windows, mean_of_LLRs, variance_of_LLRs, z_score=z_score, lookahead=lookahead)
                l = chr_length(info['chr_id'])
                crossovers[g] = [pos/l for pos in unnormalized_crossovers] #Normalize position according to the chromosome length.
                
            YMAX[g] = yabsmax if YMAX[g]< yabsmax else YMAX[g]

    for g,(ax1,(identifier,(likelihoods,info))) in enumerate(zip(AX,DATA.items())):
        mean_genomic_window_size = info['statistics']['window_size_mean']/chr_length(info['chr_id']) 
        ymax = max(YMAX[columns*(g//columns):columns*(g//columns+1)])
        ax1.errorbar( 0.88-mean_genomic_window_size, -0.76*ymax,marker=None, ls='none', xerr=25*mean_genomic_window_size, linewidth=2*scale, color='k', capsize=4*scale, zorder=20)
        ax1.text(     0.88-mean_genomic_window_size, -0.82*ymax, '25 GW',  horizontalalignment='center', verticalalignment='top',fontsize=2*fs//3, zorder=20)
        ax1.plot([0,1],[0,0],color='black', ls='dotted',alpha=0.7,zorder=0, linewidth=2*scale, scalex=False, scaley=False)
        ax1.set_title(identifier,fontsize=fs)

    for g,ax1 in enumerate(AX[:len(DATA)]):
        ymax = max(YMAX[columns*(g//columns):columns*(g//columns+1)])
        ax1.set_ylim((-1.01*ymax,+1.01*ymax))
        ax1.set_xlim((0,1)) 
        
        #Replace ticks along the x-axis 
        X_ticks = [i/10 for i in range(0,11,2)]
        X_labels = [('%g' % j) for j in X_ticks] 
        ax1.set_xticks(X_ticks)
        ax1.set_xticklabels(X_labels)
        
        ax1.tick_params(axis='x', labelsize=fs) 
        ax1.tick_params(axis='y', labelsize=fs)
        ax1.xaxis.set_tick_params(width=2*scale)
        ax1.yaxis.set_tick_params(width=2*scale)
        ###ax1.grid(color='black', linestyle='-.', linewidth=1,alpha=0.5)
        for axis in ['top','bottom','left','right']:
            ax1.spines[axis].set_linewidth(2*scale)
            
        if pairs==(('BPH','SPH'),) or pairs==(('SPH','BPH'),):
            for i in crossovers[g]:
                ax1.plot([i,i],[-1.01*ymax,1.01*ymax],color='purple', ls='dotted',alpha=0.7,zorder=19, linewidth=2*scale, scalex=False, scaley=False)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    
    plt.xlabel('Chromosomal position', fontsize=fs,labelpad=23*scale)
    plt.ylabel('Log-likelihood ratio', fontsize=fs,labelpad=45*scale)        
    
    for l in range(1,len(AX)-len(DATA)+1):
        AX[-l].tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False, width=0)
        for axis in ['top','bottom','left','right']:
            AX[-l].spines[axis].set_visible(False)
        AX[-l].xaxis.set_tick_params(labelbottom=True)
    if len(H)>1:
        fig.legend(handles=[i[0] for i in H.values()], title='', loc='lower right', ncol=len(H), fancybox=True,fontsize=fs) # bbox_to_anchor=(.5, .45)
    
    if kwargs.get('title',None): 
        fig.suptitle(kwargs['title'], fontsize=int(1.2*fs), color='black', fontweight="bold") #, y=1.01
            
        
    if save!='':
        print('Saving plot...')
        #####plt.tight_layout()
        extension = 'svg'
        plt.savefig('.'.join([save,extension]), format=extension) # bbox_inches='tight'
        plt.close(fig)
    else:
       #plt.tight_layout()
       plt.show() 

def single_plot(likelihoods,info,**kwargs):
    """ Creates a figure  depicts the log-likelihood ratio vs. chromosomal
        position for (a) BPH over disomy, (b) disomy over SPH and (c) SPH over 
        monosomy. """
        
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    scale = kwargs.get('scale', 1)
    z_score = kwargs.get('z_score', 1.96)
    bin_size = kwargs.get('bin_size', 4000000)
    save = kwargs.get('save', '')
    pairs = kwargs.get('pairs', (('BPH','disomy'),('disomy','SPH'),('SPH','monosomy')))
    
    if save!='':
        mpl.use('Agg')
    else:
        #['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']
        mpl.use('Qt5Agg')
    
    
    num_of_bins = {'chr'+str(i): chr_length('chr'+str(i))//bin_size for i in [*range(1,23)]+['X','Y']}

    
    fs = 24 * scale
    
    colors = {frozenset(('BPH','disomy')):(177/255,122/255,162/255),
              frozenset(('disomy','SPH')):(242/255,142/255,44/255),
              frozenset(('SPH','monosomy')):(239/255,106/255,92/255),
              frozenset(('disomy','monosomy')):(104/255,162/255,183/255),
              frozenset(('BPH','SPH')):(104/255,162/255,104/255)}
    
    LLRs = {(i,j): 
            {window: tuple(LLR(attrgetter(i)(l), attrgetter(j)(l)) for l in likelihoods_in_window)
                       for window,likelihoods_in_window in likelihoods.items()} 
                                                        for i,j in pairs}
        
    fig,(ax1)=plt.subplots(1,1, figsize=(16 * scale, 9 * scale))
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    H = {}
    for p,LLRs_per_genomic_window in LLRs.items():
        X,Y,E = binning(LLRs_per_genomic_window,info,num_of_bins[info['chr_id']])
        Y = [(y if y else 0) for y in Y]
        E = [(z_score*e if e else 0) for e in E]        
        T = [(x[1]+x[0])/2 for x in X]            
        
        ###ax1.plot([X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]],[i for i in Y for j in (1,2)], label=f'{p[0]:s} vs. {p[1]:s}',color=colors[p],linewidth=2)
        
        steps_x = [X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]]
        steps_y = [i for i in Y for j in (1,2)]
        H[p] = ax1.plot(steps_x, steps_y, label=f'{capitalize(p[0]):s} vs. {capitalize(p[1]):s}' ,color=colors[frozenset(p)], linewidth=2*scale, zorder=10, scalex=True, scaley=True, alpha=0.8)
        P = [(x[1]-x[0])/2 for x in X]                
        ax1.errorbar(T, Y, xerr = P, color=colors[frozenset(p)],marker=None, ls='none',alpha=1, zorder=13, linewidth=3*scale) 
        
        
        ax1.errorbar(T, Y, yerr = E, ecolor='black',marker=None, ls='none',alpha=0.2, linewidth=scale, zorder=15) 

    
    ax1.tick_params(axis='x', labelsize=fs) 
    ax1.tick_params(axis='y', labelsize=fs)
    ax1.xaxis.set_tick_params(width=scale)
    ax1.yaxis.set_tick_params(width=scale)
    ###ax1.grid(color='black', linestyle='-.', linewidth=1,alpha=0.5)
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(scale)
        
    #ax1.set_title(info['chr_id'].replace('chr', 'Chromosome '),fontsize=fs)
    ax1.set_title(f"{info['chr_id'].replace('chr', 'Chromosome '):s}, {info['depth']:.2f}x", fontsize=int(1.2 * fs), color='black', fontweight="bold")

    ax1.set_ylabel('Log-likelihood ratio', fontsize=fs,labelpad=2*scale)   
    ax1.set_xlabel('Chromosomal position', fontsize=fs,labelpad=2*scale)

    
    #Replace ticks along the x-axis 
    X_ticks = [i/10 for i in range(0,11,2)]
    X_labels = [('%g' % j) for j in X_ticks] 
    ax1.set_xticks(X_ticks)
    ax1.set_xticklabels(X_labels)
    
    #Y_ticks = [i for i in ax1.get_yticks()]
    #ax1.set_yticks(Y_ticks)
    #ax1.set_yticklabels(f'{j:g}' for j in Y_ticks) 
    
    mean_genomic_window_size = info['statistics']['window_size_mean']/chr_length(info['chr_id']) 
    ymin,ymax = ax1.get_ylim()
    ax1.errorbar( 0.9-mean_genomic_window_size, ymin + 0.08*(ymax-ymin),marker=None, ls='none', xerr=25*mean_genomic_window_size, linewidth=2*scale, color='k', capsize=4*scale)
    ax1.text(   0.9-mean_genomic_window_size, ymin + 0.05*(ymax-ymin), '25 GW',  horizontalalignment='center', verticalalignment='top',fontsize=2*fs//3)
    ax1.plot([0,1],[0,0],color='black', ls='dotted',alpha=0.5)
    ax1.set_ylim((ymin,ymax))
    ax1.set_xlim((0,1))    

    
    
    #handles, labels = ax1.get_legend_handles_labels()
    if len(H)>1:
        ax1.legend(handles=[i[0] for i in H.values()], title='', loc='upper right', ncol=len(H), fancybox=True,fontsize=int(0.8*fs))
    
    if save!='':
        print('Saving plot...')
        #ax1.set_title(save.rpartition('/')[-1].removesuffix('.png'))
        extension = 'svg'
        plt.tight_layout()
        plt.savefig('.'.join([save,extension]), format=extension, bbox_inches='tight')
        plt.close(fig)
        
    else:
       plt.tight_layout()
       plt.show()
       
def wrap_panel_plot_for_single_indv(identifier, **kwargs):
    """ Wraps the function panel_plot to show all the chromosomes from a single individual. """
    
    DATA = {}
    for i in [*range(1,23)]+['X']:
        llr_filename = kwargs.get('work_dir','.').rstrip('/') + '/' + f'{identifier:s}.chr{str(i):s}.LLR.p.bz2'
        likelihoods,info = load_likelihoods(llr_filename) 
        DATA[f"{info['chr_id'].replace('chr', 'Chromosome '):s}, {info['depth']:.2f}x"] = (likelihoods,info)
        show_info(llr_filename,info,kwargs['pairs'])
        kwargs['title'] = identifier
    panel_plot(DATA,**kwargs)
    
    return 0

def wrap_panel_plot_many_cases(filenames, **kwargs):
    """ Wraps the function panel_plot to show a panel with many cases. """
    
    DATA = {}
    for llr_filename in filenames:
        likelihoods,info = load_likelihoods(llr_filename)
        if llr_filename[-6:]=='.LLR.p':
            identifer = llr_filename[:-6].rsplit('/',1).pop()
        elif llr_filename[-10:]=='.LLR.p.bz2':
            identifer = llr_filename[:-10].rsplit('/',1).pop()
        elif llr_filename[-9:]=='.LLR.p.gz':
            identifer = llr_filename[:9].rsplit('/',1).pop()
        else:
            identifer = llr_filename.rsplit('/',1).pop()
        DATA[identifer]=(likelihoods,info)
        show_info(llr_filename,info,kwargs['pairs'])
    panel_plot(DATA, **kwargs)
    return 0

def wrap_single_plot(llr_filename, **kwargs):
    """ Wraps the function single_plot. """
    likelihoods,info = load_likelihoods(llr_filename)
    show_info(llr_filename,info,kwargs['pairs'])
    single_plot(likelihoods, info, **kwargs)
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Plots log-likelihood ratios (LLR) vs. chromosomal position from a LLR file.')
    parser.add_argument('llr_filename', metavar='LLR_FILENAME', type=str, nargs='+',
                        help='One or more LLR files created by ANEUPLOIDY_TEST, containing likelihoods to observese reads under various aneuploidy landscapes .')
    parser.add_argument('-p', '--pairs', type=str, nargs='+', metavar='scenario_A,scenario_B', default=['BPH,SPH'],
                        help='Plots the LLR between scenario A and scenario B along the chromosome. The possible pairs are: BPH,disomy; disomy,SPH; SPH,monosomy; disomy,monosomy; BPH,SPH.'
                             'In addition, giving a list of pairs would plot the LLR of each pair in the same figure, e.g. \"BPH,SPH SPH,MONOSOMY\". The default value is BPH,SPH.')
    parser.add_argument('-b', '--bin-size', type=int, metavar='INT', default=4000000,
                        help='The bin size in which the chromosome is divided. The default value is 4,000,000 bp.')
    parser.add_argument('-z', '--z-score', type=int, metavar='INT', default=1.96,
                        help='The z-score value for the confidence intervals. The default value is 1.96, which corresponds to confidence level of 95\%.')

    kwargs = vars(parser.parse_args())
    kwargs['pairs'] = [j.split(',') for j in kwargs.get('pairs','')]
    
    if  len(kwargs['llr_filename'])==1:
        kwargs['llr_filename'] = kwargs['llr_filename'].pop()
        wrap_single_plot(**kwargs)
    else:
        kwargs['filenames'] = kwargs['llr_filename']
        del kwargs['llr_filename']
        wrap_panel_plot_many_cases(**kwargs)
    sys.exit(0)

else:
    print('The module PLOT_PANEL was imported.')

######## END OF FILE ########    
   




"""
####################################################
# Produce panel plots for all cases in the folders #
####################################################

import os
work_dir = 'results2/'
identifiers = {i.split('.')[0] for i in os.listdir(work_dir) if i[-3:]=='bz2'}
for identifier in identifiers:
    try:
        if not os.path.isfile(work_dir.rstrip('/') + '/' + identifier+'.svg'):
            wrap_panel_plot_for_single_indv(identifier, bin_size=4000000, pairs=(('BPH','SPH'),), save=identifier, work_dir=work_dir)
    except Exception as e:
        print(identifier,e)
"""


#wrap_panel_plot_for_single_indv(identifier='robert', bin_size=4000000, pairs=(('BPH','SPH'),), work_dir='/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/single/')
"""
filenames = [
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr1.x0.100.HG00105A.NA12827B.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr2.x0.100.NA12872A.HG00343A.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr3.x0.100.HG01710B.NA20774A.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr4.x0.100.NA20813A.HG01707A.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr5.x0.100.NA12273B.HG00111B.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr6.x0.100.HG01773A.NA20770A.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr7.x0.100.HG02235B.NA12286B.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr8.x0.100.HG01521A.NA20521A.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr9.x0.100.HG00240B.HG01682A.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr10.x0.100.HG00360B.HG00109A.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr11.x0.100.NA12778A.HG01527A.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr12.x0.100.NA20511B.HG01682A.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr13.x0.100.NA12828B.HG00261B.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr14.x0.100.NA07034B.NA12815A.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr15.x0.100.HG00130A.HG00174A.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr16.x0.100.NA11993B.HG00146A.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr17.x0.100.HG00108B.HG01522A.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr18.x0.100.HG00278A.HG01697B.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr19.x0.100.NA12144A.NA12413B.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr20.x0.100.HG00240B.HG01709B.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr21.x0.100.HG00235A.HG01767B.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chr22.x0.100.NA12249A.NA12891B.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/simulations/simulation_single_indv/simulated.SPH.chrX.x0.100.HG00362A.NA20766B.LLR.p.bz2" 
]

filenames = [
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/COR-A-16-SplitA-11-Aug-2020.chr2.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/COR-A-10-SplitA-11-Aug-2020.chr6.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/WIG-A-4-27-Nov-2020.chr2.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/RUF-A-10-SplitA-14-Jul-2020.chr4.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/KNI-J-M-4-SplitA-15-Jul-2020.chr8.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/GAN-D-4-splitA-12-Feb-2021.chr4.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/ZEU-A-1-RETEST-02-Sep-2021.chr2.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/HOL-K-11-22-Dec-2020.chr5.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/DES-N-1-splitA-28-Aug-2021.chr7.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/JIA-L-3-25-Mar-2020.chr2.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/COR-A-6-SplitA-11-Aug-2020.chr8.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/BOT-S-1-10-Oct-2021.chr3.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/AUG-G-6-03-Dec-2020.chr6.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/BHA-J-7-NextSeq.chr9.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/COM-M-11-1-01-Oct-2020.chr8.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/HOL-K-11-22-Dec-2020.chr7.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/SAN-B-2-splitA-26-May-2020.chr9.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/ZEU-A-1-RETEST-02-Sep-2021.chr8.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/JIA-L-3-25-Mar-2020.chr5.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/HOL-K-11-22-Dec-2020.chr9.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/REI-SM-1-Split-A-NextSeq.chr5.LLR.p.bz2",
"/Users/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/CReATe_results/MES-J-8-SplitA-15-Jun-2020.chr7.LLR.p.bz2"
]

filenames = [
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/CHA-S-11-30-May-2020.chr6.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/MAR-A-1-13-Aug-2020.chr2.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/RESEARCH-E7-8-01-Sep-2020.chr8.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/BOU-M-3-11-Oct-2021.chr1.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/LOC-A-7-20-Dec-2020.chr2.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/HEN-T-3-17-Dec-2020.chr1.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/COM-M-5-1-01-Oct-2020.chr1.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/BOU-JB-5-29-Jun-2020.chr2.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/WRI-A-2-06-Oct-2020.chr3.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/RESEARCH-E6-19-01-Sep-2020.chr2.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/LOH-J-2-28-Nov-2021.chr3.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/MES-A-1-02-Jul-2020.chr2.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/YU-SH-3-16-Feb-2021.chr1.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/DHA-D-3-17-Feb-2021.chr3.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/TRO-S-2-03-Jul-2020.chr4.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/SCH-S-4-12-Jun-2021.chr3.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/GOY-D-12-15-Dec-2020.chr5.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/SCH-P-1-04-Jul-2021.chr4.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/NG-M-K-2-16-Nov-2020.chr12.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/LEP-T-4-24-Jun-2021.chr1.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/CAR-N-1-07-Nov-2020.chr6.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/AND-J-5-02-Nov-2020.chr3.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/COM-M-5-1-01-Oct-2020.chr5.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/COM-M-1-1-01-Oct-2020.chr5.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/RIC-L-7-1-08-Sep-2020.chr3.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/AZA-G-8-03-Oct-2021.chr8.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/NG-M-K-1-16-Jun-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/KHA-H-10-04-Nov-2020.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/IAN-B-H-E-16-05-Dec-2020.chr6.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/KAU-C-3-10-Dec-2020.chr4.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/VIG-B-L-I-4-05-Sep-2021.chr4.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/COL-S-1-09-Apr-2021.chr1.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/GRE-O-11-RETEST-14-Aug-2020.chr7.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/MRO-L-2-23-Jul-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/WU-A-5-11-Aug-2021.chr8.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/CON-J-12-17-Dec-2020.chr3.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/GOY-D-3-15-Dec-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/AND-J-12-02-Nov-2020.chr11.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/MOL-C-1-18-Jun-2020.chr6.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/BOT-S-1-10-Oct-2021.chr4.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/GEL-B-3-10-Jan-2021.chr7.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/POL-T-7-retest-13-Mar-2021.chr10.LLR.p.bz2"
]


filenames = ["/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/ROS-L-8-NextSeq.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/DO-T-1-14-Sep-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/GOY-D-3-15-Dec-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/MRO-L-2-23-Jul-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/CAS-S-27-17-Nov-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/HOL-TR-1-16-Aug-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/NG-M-K-1-16-Jun-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/NG-M-K-2-16-Jun-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/SAI-D-10-30-Nov-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/BRU-M-N-2-14-Apr-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/WHE-B-11-1-17-Sep-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/DHA-D-3-17-Feb-2021.chr18.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/LAV-F-4-06-Jul-2021.chr18.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/SHE-M-9-03-Sep-2021.chr18.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/BRA-S-1-splitA-25-Sep-2020.chr18.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/GOS-L-9-splitA-19-Jun-2020.chr18.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/PUG-J-15-splitA-14-Jun-2020.chr18.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/KNI-J-M-11-SplitA-15-Jul-2020.chr18.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/BEK-N-5-21-Aug-2020.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/GON-S-24-17-Sep-2021.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/KHA-H-10-04-Nov-2020.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/TIL-M-11-18-Dec-2021.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/IAN-B-H-E-16-05-Dec-2020.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/LON-M-6-splitA-25-Nov-2021.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/POL-T-7-retest-13-Mar-2021.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/TUN-BN-8-splitA-12-Dec-2020.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/HAG-R-4-2-NIPT-TOP-16-Nov-2021.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/CHE-X-6-06-Jul-2020.chr12.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/GON-J-7-15-Aug-2021.chr12.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/SAH-P-7-04-Dec-2021.chr12.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/DON-D-16-18-Nov-2021.chr12.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/NG-M-K-2-16-Nov-2020.chr12.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/RIC-L-7-1-08-Sep-2020.chr12.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/MON-AL-FT1-18-Oct-2021.chr12.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/TUC-E-1-splitA-04-Dec-2020.chr12.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/HOL-B-37-splitA-25-Sep-2021.chr12.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/AFR-F-5-13-Dec-2020.chr13.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/DUK-A-2-15-Jul-2020.chr13.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/RAD-A-2-24-Oct-2020.chr13.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/RAD-A-3-03-Aug-2020.chr13.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/RAD-M-2-18-Nov-2020.chr13.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/LIU-TXY-6-11-Oct-2021.chr13.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/RIC-L-2-1-08-Sep-2020.chr13.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/RESEARCH-E5-9-01-Sep-2020.chr13.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/BRA-S-1-20-Oct-2020.chr21.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/TUN-BN-8-splitA-12-Dec-2020.chr21.LLR.p.bz2"
]

filenames = ["/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/ROS-L-8-NextSeq.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/DO-T-1-14-Sep-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/CAS-J-3-14-Jan-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/GAN-C-3-23-Nov-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/GIL-K-8-14-Jun-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/GOY-D-3-15-Dec-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/HAS-A-2-16-Jun-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/JIW-S-1-09-Jun-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/KOU-R-5-10-May-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/LIN-Y-3-02-Aug-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/MAC-S-5-04-Aug-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/MOO-K-2-04-Apr-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/MRO-L-2-23-Jul-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/PAH-Y-4-28-Jun-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/SHU-S-3-05-May-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/BOX-AN-7-23-Apr-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/CAS-S-27-17-Nov-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/HOL-TR-1-16-Aug-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/HYD-D-17-04-Jul-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/KLE-M-12-03-Oct-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/MAT-J-19-08-Sep-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/NG-M-K-1-16-Jun-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/NG-M-K-2-16-Jun-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/NG-M-K-3-11-Jul-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/NG-M-K-3-16-Jun-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/SAI-D-10-30-Nov-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/BRU-M-N-2-14-Apr-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/BRU-M-N-3-14-Apr-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/WHE-B-11-1-17-Sep-2020.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/WAN-H-8-RETEST-24-Jan-2021.chr9.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/LD-PGTA_analysis/results/CReATe_segmentals/JOH-J-18-retest-17-Feb-2021.chr9.LLR.p.bz2"]

wrap_panel_plot_many_cases(filenames,pairs=(('BPH','SPH'),),title='BPH vs. SPH')
#wrap_single_plot(filenames[0],pairs=(('BPH','SPH'),))
"""