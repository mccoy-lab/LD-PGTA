#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PLOT_PANEL

Plots log-likelihood ratio vs. chromosomal position from a LLR file.

May 20, 2020
"""

import pickle, bz2, gzip, collections
from statistics import mean, variance
from math import log
from operator import attrgetter
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
    """ Calculates the mean and population standard deviation of the mean of random variables.
        Each row of A represents a random variable, with observations in the columns."""
    if type(A)==dict:
        A = tuple(tuple(i) for i in A.values()) 
    
    M, N = len(A), len(A[0])
    mu = sum(sum(likelihoods_in_window)/N for likelihoods_in_window in A)
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

def show_info(filename, info, pairs):
    S = info['statistics']
    print('\nFilename: %s' % filename)
    print('\nSummary statistics')
    print('------------------')    
    print('Chromosome ID: %s, Depth: %.2f.' % (info['chr_id'],info['depth']))
    print('Number of genomic windows: %d, Mean and standard error of genomic window size: %d, %d.' % (S.get('num_of_windows',0),S.get('window_size_mean',0),S.get('window_size_std',0)))
    print('Mean and standard error of meaningful reads per genomic window: %.1f, %.1f.' % (S.get('reads_mean',0), S.get('reads_std',0)))
    print('Ancestry: %s, Fraction of alleles matched to the reference panel: %.3f.' % (str(info['ancestry']),info['statistics']['matched_alleles']))

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

def detect_transition(X,Y,E):
    """  Detection of meiotic crossovers based on inferred switches between
         tracts of BPH and SPH trisomy. """
    
    A = [(l,j,k) for i,j,k in zip(X,Y,E) for l in i if j!=None and (abs(j)-k)>0 ]
    if len(A)>2:
        x,y,e = zip(*A)
        result = [(x[i+1]+x[i])/2 for i,(a,b) in enumerate(zip(y[:-1],y[1:])) if b/a<0]
    else:
        result = []
        
    return result


def capitalize(x):
    return x[0].upper() + x[1:]
    
def panel_plot(DATA,num_of_bins_in_chr21,pairs,**kwargs):
    """ Creates a multi-panel figure. For each numbered chromosome, a figure 
        depicts the log-likelihood ratio vs. chromosomal position for BPH over
        SPH. """
    
    scale = 0.5
    z_score = 1
    fs=28 * scale
    rows = 4
    columns = 6
    import matplotlib as mpl
    save = kwargs.get('save', '')
    if save!='':
            mpl.use('Agg')
    else:
        #['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']
        mpl.use('Qt5Agg')
    mpl.rcParams.update({'figure.max_open_warning': 0})
    import matplotlib.pyplot as plt
    #from scipy.interpolate import interp1d
    num_of_bins = {'chr'+str(i): num_of_bins_in_chr21*chr_length('chr'+str(i))//chr_length('chr21') for i in [*range(1,23)]+['X','Y']}

    colors = {frozenset(('BPH','disomy')):(177/255,122/255,162/255),
              frozenset(('disomy','SPH')):(242/255,142/255,44/255),
              frozenset(('SPH','monosomy')):(239/255,106/255,92/255),
              frozenset(('disomy','monosomy')):(104/255,162/255,183/255),
              frozenset(('BPH','SPH')):(104/255,162/255,104/255)}


    fig,axs = plt.subplots(rows ,columns, sharex='col', sharey='row', figsize=(6.666 * columns * scale, 5.625 * rows * scale))
    fig.subplots_adjust(left=0.05, bottom=0.06, right=.99, top=.97, wspace=None, hspace=None)
    #fig.suptitle(kwargs.get('title', ''), fontsize=16)
    #fig.text(0.5, 0, 'Chromosomal position (normalized)',fontsize=28, ha='center')
    #fig.text(0, 0.5, 'Log-likelihood ratio per genomic window', fontsize=28, va='center', rotation='vertical')


    AX = [i for j in axs for i in j]
    
    H = {}
    YMAX = [0]*len(DATA)
    transitions = []
    for a,b in pairs:
        for g,(ax1,(likelihoods,info)) in enumerate(zip(AX,DATA)):
    
            if (a,b) in info['statistics']['LLRs_per_genomic_window']:
                LLRs_per_genomic_window = info['statistics']['LLRs_per_genomic_window'][(a,b)]
            else:
                    
                LLRs = {window: tuple(LLR(attrgetter(a)(l), attrgetter(b)(l)) for l in likelihoods_in_window)
                                   for window,likelihoods_in_window in likelihoods.items()} 
                                                                   
                    
                LLRs_per_genomic_window = {window: mean_and_var(LLRs_in_window) for window, LLRs_in_window in LLRs.items()}
                                                                    
            
            X,Y,E = binning(LLRs_per_genomic_window,info,num_of_bins[info['chr_id']])
            Y = [(y if y else 0) for y in Y]
            E = [(z_score*e if e else 0) for e in E]
            
            T = [(x[1]+x[0])/2 for x in X]                
            steps_x = [X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]]
            steps_y = [i for i in Y for j in (1,2)]
            H[a,b] = ax1.plot(steps_x, steps_y, label=capitalize(f'{a:s} vs. {b:s}') ,color=colors[frozenset((a,b))], linewidth=2, zorder=10, scalex=True, scaley=True, alpha=0.8)
            
            P = [(x[1]-x[0])/2 for x in X]                
            ax1.errorbar(T, Y, xerr = P, ecolor=colors[frozenset((a,b))],marker=None, ls='none',alpha=1, zorder=13, linewidth=5*scale) 
            ax1.errorbar(T, Y, yerr = E, ecolor='black',marker=None, ls='none',alpha=0.2, zorder=15, linewidth=4*scale) 
            
            yabsmax = max(map(abs,Y))
            
            if pairs==(('BPH','SPH'),) or pairs==(('SPH','BPH'),):
                transitions.append(detect_transition(X,Y,E))
                    
            YMAX[g] = yabsmax if YMAX[g]< yabsmax else YMAX[g]

    for g,(ax1,(likelihoods,info)) in enumerate(zip(AX,DATA)):
        mean_genomic_window_size = info['statistics']['window_size_mean']/chr_length(info['chr_id']) 
        ymax = max(YMAX[columns*(g//columns):columns*(g//columns+1)])
        ax1.errorbar( 0.88-mean_genomic_window_size, -0.76*ymax,marker=None, ls='none', xerr=25*mean_genomic_window_size, linewidth=2*scale, color='k', capsize=4*scale, zorder=20)
        ax1.text(     0.88-mean_genomic_window_size, -0.82*ymax, '25 GW',  horizontalalignment='center', verticalalignment='top',fontsize=2*fs//3, zorder=20)
        ax1.plot([0,1],[0,0],color='black', ls='dotted',alpha=0.7,zorder=0, linewidth=2*scale, scalex=False, scaley=False)
        #ax1.set_title(info['chr_id'].replace('chr', 'Chromosome '),fontsize=fs)
        ax1.set_title(f"{info['chr_id'].replace('chr', 'Chromosome '):s}, {info['depth']:.2f}x",fontsize=fs)
        
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
            for i in transitions[g]:
                ax1.plot([i,i],[-1.01*ymax,1.01*ymax],color='purple', ls='dotted',alpha=0.7,zorder=19, linewidth=2*scale, scalex=False, scaley=False)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    ###plt.title(f'{RATIO[0]:s} vs. {RATIO[1]:s}\n', fontsize=int(1.2*fs))
    
    plt.xlabel('Chromosomal position (normalized)', fontsize=fs,labelpad=23*scale)
    plt.ylabel('Log-likelihood ratio (normalized)', fontsize=fs,labelpad=45*scale)        
    
    for l in range(1,len(AX)-len(DATA)+1):
        AX[-l].tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False, width=0)
        for axis in ['top','bottom','left','right']:
            AX[-l].spines[axis].set_visible(False)
        AX[-l].xaxis.set_tick_params(labelbottom=True)
    AX[-1].legend(handles=[i[0] for i in H.values()], title='', bbox_to_anchor=(.5, .45), loc='upper center', ncol=1, fancybox=True,fontsize=fs)
    if kwargs.get('title',None): AX[-1].set_title(kwargs.get('title',None), fontsize=int(fs), y=0.55, color='purple')
            
        
    if save!='':
        print('Saving plot...')
        #ax1.set_title(save.rpartition('/')[-1].removesuffix('.png'))
        plt.tight_layout()
        extension = 'svg'
        plt.savefig('.'.join([save,extension]), format=extension, bbox_inches='tight')
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
    #from scipy.interpolate import interp1d
    
    scale = kwargs.get('scale', 1)
    z_score = kwargs.get('z_score', 1)
    num_of_bins = kwargs.get('num_of_bins', 10)
    save = kwargs.get('save', '')
    pairs = kwargs.get('pairs', (('BPH','DISOMY'),('DISOMY','SPH'),('SPH','MONOSOMY')))
    
    if save!='':
        mpl.use('Agg')
    else:
        #['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']
        mpl.use('Qt5Agg')
    
    
    fs = 24 * scale
    
    colors = {('BPH','DISOMY'):(177/255,122/255,162/255),
              ('DISOMY','SPH'):(242/255,142/255,44/255),
              ('SPH','MONOSOMY'):(239/255,106/255,92/255),
              ('DISOMY','MONOSOMY'):(104/255,162/255,183/255),
              ('BPH','SPH'):(104/255,162/255,104/255)}
    
    LLRs = {(i,j): 
            {window: tuple(LLR(attrgetter(i)(l), attrgetter(j)(l)) for l in likelihoods_in_window)
                       for window,likelihoods_in_window in likelihoods.items()} 
                                                        for i,j in pairs}
        
    #LLRs_per_genomic_window = {(i,j): 
    #        {window: mean_and_var(LLRs_in_window) for window, LLRs_in_window in LLRs[(i,j)].items()}
    #                                                    for (i,j) in pairs}
    
    fig,(ax1)=plt.subplots(1,1, figsize=(16 * scale, 9 * scale))
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    H = {}
    for p,LLRs_per_genomic_window in LLRs.items():
        X,Y,E = binning(LLRs_per_genomic_window,info,num_of_bins)
        Y = [(y if y else 0) for y in Y]
        E = [(z_score*e if e else 0) for e in E]        
        T = [(x[1]+x[0])/2 for x in X]            
        
        ###ax1.plot([X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]],[i for i in Y for j in (1,2)], label=f'{p[0]:s} vs. {p[1]:s}',color=colors[p],linewidth=2)
        
        steps_x = [X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]]
        steps_y = [i for i in Y for j in (1,2)]
        H[p] = ax1.plot(steps_x, steps_y, label=f'{p[0]:s} vs. {p[1]:s}',color=colors[p], linewidth=2*scale, zorder=10, scalex=True, scaley=True, alpha=0.8)
        P = [(x[1]-x[0])/2 for x in X]                
        ax1.errorbar(T, Y, xerr = P, color=colors[p],marker=None, ls='none',alpha=1, zorder=13, linewidth=3*scale) 
        
        
        ax1.errorbar(T, Y, yerr = E, ecolor='black',marker=None, ls='none',alpha=0.2, linewidth=scale, zorder=15) 

    
    ax1.tick_params(axis='x', labelsize=fs) 
    ax1.tick_params(axis='y', labelsize=fs)
    ax1.xaxis.set_tick_params(width=scale)
    ax1.yaxis.set_tick_params(width=scale)
    ###ax1.grid(color='black', linestyle='-.', linewidth=1,alpha=0.5)
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(scale)
        
    #ax1.set_title(info['chr_id'].replace('chr', 'Chromosome '),fontsize=fs)
    ax1.set_title(f"{info['chr_id'].replace('chr', 'Chromosome '):s}, {info['depth']:.2f}x",fontsize=fs)

    ax1.set_ylabel('Log-likelihood ratio (normalized)', fontsize=fs,labelpad=2*scale)   
    ax1.set_xlabel('Chromosomal position (normalized)', fontsize=fs,labelpad=2*scale)

    
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
    ax1.legend(handles=[i[0] for i in H.values()], title='', loc='upper center', ncol=len(H), fancybox=True,fontsize=8*fs//10)
    
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
       
def wrap_panel_plot(identifier,pairs=(('BPH','SPH'),), num_of_bins_in_chr21=5, save='', work_dir=''):
    """ Wraps the function panel_plot. """
    #DATA = {filename: load_likelihoods(filename)}      
    DATA = {f'{identifier:s}.chr{str(i):s}': load_likelihoods(work_dir + f'{identifier:s}.chr{str(i):s}.LLR.p.bz2') for i in [*range(1,23)]+['X']}
    
    for f,(likelihoods,info) in DATA.items():
        show_info(f'{f:s}.LLR.p.bz2',info,pairs)
    panel_plot(DATA.values(),num_of_bins_in_chr21=num_of_bins_in_chr21,pairs=pairs,title=f'{identifier:s}',save=save)
    return 0

def wrap_single_plot(llr_filename,pairs,bins):
    """ Wraps the function single_plot. """
    likelihoods,info =  load_likelihoods(llr_filename)
    show_info(llr_filename,info,pairs)
    single_plot(likelihoods,info,pairs=pairs,num_of_bins=bins)
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Plots log-likelihood ratios (LLR) vs. chromosomal position from a LLR file.')
    parser.add_argument('llr_filename', metavar='LLR_FILENAME', type=str,
                        help='A pickle file created by ANEUPLOIDY_TEST, containing likelihoods to observese reads under various aneuploidy landscapes .')
    parser.add_argument('-p', '--pairs', type=str, nargs='+', metavar='scenario_A,scenario_B scenario_C,scenario_D', default=['BPH,SPH'],
                        help='Plots the LLR between scenario A and scenario B along the chromosome. The possible pairs are: BPH,DISOMY; DISOMY,SPH; SPH,MONOSOMY; DISOMY,MONOSOMY; BPH,SPH.'
                             'In addition, giving a list of pairs would plot the LLR of each pair in the same figure, e.g. \"BPH,SPH SPH,MONOSOMY\". The default value is BPH,SPH.')
    parser.add_argument('-b', '--bins', type=int, metavar='INT', default=10,
                        help='The numbers of bins the chromosome is divided into. The default value is 10.')

    kwargs = vars(parser.parse_args())
    kwargs['pairs'] = [j.split(',') for j in kwargs.get('pairs','')]
    wrap_single_plot(**kwargs)
    sys.exit(0)

else:
    print('The module PLOT_PANEL was imported.')
    
#identifier = '13094FA-B6MTL_3'

####################################################
# Produce panel plots for all cases in the folders #
####################################################

import os
work_dir = 'results2/'
identifiers = {i.split('.')[0] for i in os.listdir(work_dir) if i[-3:]=='bz2'}
for identifier in identifiers:
    try:
        if not os.path.isfile(work_dir+identifier+'.svg'):
            wrap_panel_plot(identifier,pairs=(('BPH','SPH'),),save=identifier,work_dir=work_dir, num_of_bins_in_chr21=20)
    except Exception as e:
        print(identifier,e)