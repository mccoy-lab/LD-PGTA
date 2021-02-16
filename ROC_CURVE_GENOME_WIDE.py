#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROC_CURVE_GENOME_WIDE

Daniel Ariad (daniel@ariad.org)
Dec 20, 2020
"""

import pickle, statistics, collections
from itertools import starmap
from operator import or_
from ANEUPLOIDY_TEST import LLR 
    
def mean_and_var(sample):
    """ Calculates the mean and the sample standard deviation. """
    mean = statistics.mean(sample)
    var = statistics.variance(sample, xbar=mean)
    return mean, var

def std_of_mean(variances):
    """ Standard error of the mean of uncorrelated variables, based on the
        Bienaymé formula. """
    return sum(variances)**.5/len(variances)


def load_llr(filename):
    """ Loads from a file a dictionary that lists genomic windows that contain
    at least two reads and gives the bootstrap distribution of the 
    log-likelihood BPH/SPH ratios (LLRs). """
    
    with open(filename, 'rb') as f:
        likelihoods = pickle.load(f)
        info = pickle.load(f)
    return likelihoods, info

def show_info(filename,info):
    S = info['statistics']
    X = info['statistics']['LLRs_per_chromosome'][('BPH','SPH')]
    print('\nFilename: %s' % filename)
    print('Depth: %.2f, Chromosome ID: %s, Mean and standard error of meaningful reads per genomic window: %.1f, %.1f.' % (info['depth'], info['chr_id'], S['reads_mean'], S['reads_std']))
    print('Number of genomic windows: %d, Mean and standard error of genomic window size:  %d, %d.' % (S['num_of_windows'],S['window_size_mean'],S['window_size_std']))
    print('Mean LLR: %.3f, Standard error of the mean LLR: %.3f' % ( X['mean'],  X['std_of_mean']))
    print('Fraction of genomic windows with a negative LLR: %.3f' % (X['fraction_of_negative_LLRs']))
    print('The calculation was done in %.3f sec.' % info['runtime'])


def build_confidence_dict(criterias, num_of_buckets, ratio, work_dir):
    """ Iterates over all the data files in the folder and creates a dictionary
    the list all the files that fit the criterias and gives their analysis 
    (via the confidence function). """
    
    import glob
    effective_ratio = tuple('BPH' if i=='transitions' else i for i in ratio)
    filenames = glob.glob(work_dir + '*.LLR.p')
    result = {r: {} for r in ratio}

    buffer = {r: collections.defaultdict(dict) for r in ratio}
    for filename in filenames:
        likelihoods, info = load_llr(filename)
        subinfo = {x: info.get(x,None) for x in criterias.keys()}
        ### print(subinfo)
        
        scenario = info.get('scenario','')
        if criterias==subinfo and scenario in ratio:
            show_info(filename,info)  
        
            if effective_ratio in info['statistics']['LLRs_per_genomic_window']:
                LLRs_per_genomic_window = [*info['statistics']['LLRs_per_genomic_window'][effective_ratio].values()]
            else:
                i,j = effective_ratio
                _ = {}
                LLRs_per_genomic_window = [mean_and_var([*starmap(LLR, ((_[i], _[j]) for _['monosomy'], _['disomy'], _['SPH'], _['BPH'] in L))])
                               for window,L in likelihoods.items()]
            buffer[scenario][info['chr_id']][filename] = LLRs_per_genomic_window    
    
    count = {scenario: min(len(buffer[scenario][f'chr{i:d}']) for i in range(1,23)) for scenario in ratio}
        
    result = collections.defaultdict(dict)
    for scenario in ratio:
        for i in range(count[scenario]):
            M = []
            V = []
            for j in range(1,23):
                 filename, LLRs = buffer[scenario][f'chr{j:d}'].popitem()
                 mean, std = zip(*LLRs)
                 M.extend(mean)
                 V.extend(std)
            
            result[scenario][i] = {'mean': statistics.mean(M), 'std': std_of_mean(V)}
                                               
    
    return result

def build_ROC_curve(criterias, positive, ratio, thresholds, num_of_buckets, work_dir):
    """ Creates a nested dictionary that lists bins and thresholds and gives 
        the false and true positive rates. """ 
    
    num_of_buckets = 1
    DATA = build_confidence_dict(criterias, num_of_buckets, ratio, work_dir)
    B,A = DATA[ratio[0]], DATA[ratio[1]]
    
    print(len(B),len(A))
    
    result = {}

    for z in thresholds:
        true_B = [file['mean'] > z * file['std'] for file in B.values()]
        false_B = [file['mean'] > z * file['std'] for file in A.values()]
        
        false_A = [file['mean'] < -z * file['std'] for file in B.values()]
        true_A = [file['mean'] < -z * file['std'] for file in A.values()]
        
        if positive == 'both':
            undetermined_B = [*map(or_,true_B,false_A)].count(False) / len(true_B)
            undetermined_A = [*map(or_,true_A,false_B)].count(False) / len(true_A)
            undetermined = 0.5 * (undetermined_A + undetermined_B)
        
            TPR = 0.5 * (true_B.count(1)/len(true_B) + true_A.count(1)/len(true_A))
            FPR = 0.5 * (false_B.count(1)/len(false_B) + false_A.count(1)/len(false_A))
        elif positive == 'A':
            TPR = true_A.count(1)/len(true_A)
            FPR = false_A.count(1)/len(false_A)
            undetermined = None
        elif positive == 'B':
            TPR = true_B.count(1)/len(true_B)
            FPR = false_B.count(1)/len(false_B)
            undetermined = None
        else:
            break
        
        result[z] = (FPR,TPR,undetermined)
            
    return result

def plot_ROC_curve(x):
    """ Plots the ROC curve. """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    for i in range(len(x)):
        ax.scatter(*zip(*x[i].values()), label=f'bin {i:d}', s=len(x)+1-i)        
    ax.legend()
    ax.grid(True)
    plt.show()
    
def plot_ROC_curve_fig7a(x):
    """ Plots the ROC curve. """
    import matplotlib.pyplot as plt
    colors = {'cyan': (104/256,162/256,183/256),
              'peach': (239/256,106/256,92/256),
              'purple': (177/256,122/256,162/256),
              'orange': (242/256,142/256,44/256)}
    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(111)
    #fig, ax1 = plt.subplots()
    fs = 24
    ax1.set_title('ROC curve of a triploid classfier',fontsize=fs)
    ax2 = ax1.twinx()
    ax1.tick_params(axis='y', labelsize=fs )
    ax2.tick_params(axis='y', labelsize=fs )
    ax1.tick_params(axis='x', labelsize=fs )
    ax1.set_xlabel('False Positive Rate',fontsize=fs)
    ax1.set_ylabel('True Positive Rate',fontsize=fs)
    ax2.set_ylabel('$z$-score', color=colors['peach'], fontsize=fs)

    ax1.set_ylim(-0.02,1.02)
    ax1.set_xlim(-0.02,1.02)
    #ax1.plot([0,1],[0.77,0.77],color='gray',linestyle='--')
    for k,v in x.items():
        #ax1.scatter(*zip(*v.values()), label=f'{k:s}', s=10, color='black', marker='s')    
        X, Y, Z = zip(*v.values())
        ax1.plot(X, Y, label=f'{k:s}', color='black',linewidth=6)   
        #ax2.scatter([*v],[*zip(*v.values())][1], label=f'{k:s}', s=10, color='red')        
        X,Y = zip(*((k,-i) for i,(k,l,m) in v.items() if -11.5<i<-2))
        #ax2.scatter(X,Y, label=f'{k:s}', s=1, color=peach)
        ax2.plot(X,Y, label=f'{k:s}', color=colors['peach'], linewidth=2)         
    #ax1.legend()
    #ax1.grid(True)
    ax2.spines['right'].set_color(colors['peach'])
    ax2.tick_params(axis='y', colors=colors['peach'])
    fig.savefig('fig7a.pdf',bbox_inches='tight', format='pdf')
    #plt.tight_layout()
    plt.show()

def plot_ROC_curve_fig7b(x):
    """ Plots the ROC curve. """
    import matplotlib.pyplot as plt
    colors = {'cyan': (104/256,162/256,183/256),
              'peach': (239/256,106/256,92/256),
              'purple': (177/256,122/256,162/256),
              'orange': (242/256,142/256,44/256)}
    
    fig = plt.figure(figsize=(10, 10))
    fs = 24
    ax = fig.add_subplot(111)
    ax.set_title('Balanced ROC curve',fontsize=fs-1)

    #fig, ax1 = plt.subplots()
    ax.set_xlabel('False Positive Rate',fontsize=fs)
    ax.set_ylabel('True Positive Rate',fontsize=fs)
    ax.tick_params(axis='y', labelsize=fs )
    ax.tick_params(axis='x', labelsize=fs )
    ax.set_ylim(-0.02,1.02)
    ax.set_xlim(-0.02,1.02)
    style = {0: {'color': 'black', 'linewidth':6, 'linestyle':'solid'},
             1: {'color': colors['orange'], 'linewidth':6, 'linestyle':'dashed'}}
    for i, (k,v) in enumerate(x.items()):
        #if i==1: continue
        #ax.scatter(*zip(*v.values()), label=f'{k:s}', marker='s', **style[i])
        X, Y, Z = zip(*v.values())
        ax.plot(X,Y, label=f'{k:s}',  **style[i])        
        
    ax.legend(prop={'size': fs-4})
    #ax.grid(True)
    fig.savefig('fig7b.pdf',bbox_inches='tight', format='pdf')
    plt.tight_layout()
    plt.show()

def plot_ROC_curve_zscore(x):
    """ Plots the ROC curve. """
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(111)
    #fig, ax1 = plt.subplots()
    fs = 24
    ax1.set_title('ROC curve of a diploid classfier', fontsize=fs)
    ax2 = ax1.twinx()
    ax1.tick_params(axis='y', labelsize=fs )
    ax2.tick_params(axis='y', labelsize=fs )
    ax1.tick_params(axis='x', labelsize=fs )
    ax1.set_xlabel('False Positive Rate', fontsize=fs)
    ax1.set_ylabel('True Positive Rate', fontsize=fs)
    ax2.set_ylabel('$z$-score', color='red', fontsize=fs)
    ax1.set_ylim(-0.02,1.02)
    ax1.set_xlim(-0.02,1.02)
    ax1.plot([0.23,0.23],[-0.1,1.1],color='gray',linestyle='--')
    for k,v in x.items():
        ax1.scatter(*zip(*v.values()), label=f'{k:s}', s=10, color='black')        
        #ax2.scatter([*v],[*zip(*v.values())][1], label=f'{k:s}', s=10, color='red')        
        X,Y = zip(*((k,i) for i,(k,l) in v.items() if 2.2<i<8.4))
        ax2.scatter(X,Y, label=f'{k:s}', s=1, color='red')        
    #ax1.legend()
    #ax1.grid(True)
    ax2.spines['right'].set_color('red')
    ax2.tick_params(axis='y', colors='red')
    plt.tight_layout()
    plt.show()    
    
def plot_ROC_curve_undetermined(x):
    """ Plots the ROC curve. """
    import matplotlib.pyplot as plt
    colors = {'cyan': (104/256,162/256,183/256),
              'peach': (239/256,106/256,92/256),
              'purple': (177/256,122/256,162/256),
              'orange': (242/256,142/256,44/256)}
    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(111)
    #fig, ax1 = plt.subplots()
    fs = 24
    ax1.set_title('ROC curve of a triploid classfier',fontsize=fs)
    ax2 = ax1.twinx()
    ax1.tick_params(axis='y', labelsize=fs )
    ax2.tick_params(axis='y', labelsize=fs )
    ax1.tick_params(axis='x', labelsize=fs )
    ax1.set_xlabel('False Positive Rate',fontsize=fs)
    ax1.set_ylabel('True Positive Rate',fontsize=fs)
    ax2.set_ylabel('Undetermined Rate', color=colors['peach'], fontsize=fs)
    ax1.set_ylim(-0.02,1.02)
    ax1.set_xlim(-0.02,1.02)
    #ax1.plot([0,1],[0.77,0.77],color='gray',linestyle='--')
    for k,v in x.items():
        #ax1.scatter(*zip(*v.values()), label=f'{k:s}', s=10, color='black', marker='s')        
        X, Y, Z = zip(*v.values())
        ax1.plot(X,Y, label=f'{k:s}', color='black',linewidth=6)   
        #ax2.scatter([*v],[*zip(*v.values())][1], label=f'{k:s}', s=10, color='red')        
        #ax2.scatter(X,Y, label=f'{k:s}', s=1, color=peach)
        ax2.plot(X,Z, label=f'{k:s}', color=colors['peach'], linewidth=2)         
    #ax1.legend()
    #ax1.grid(True)
    ax2.spines['right'].set_color(colors['peach'])
    ax2.tick_params(axis='y', colors=colors['peach'])
    #fig.savefig('fig7a.pdf',bbox_inches='tight', format='pdf')
    #plt.tight_layout()
    plt.show()
    
def configuration(C):
    C0 = {
         'depth': 0.01,
         'read_length': 36,
         'window_size': 0,
         'min_reads': 6,
         'max_reads': 4,
         'minimal_score': 2,
         'min_HF': 0.05}
    
    C1 = {
         'depth': 0.1,
         'read_length': 36,
         'window_size': 0,
         'min_reads': 32,
         'max_reads': 16,
         'minimal_score': 2,
         'min_HF': 0.05}


    return locals()[f'{C:s}']

if __name__ == "__main__":    
    Z = [i/500 for i in range(-10000,10000)]
    A = {}
    #for SP in ('EUR','EAS','SAS','AFR','AMR'):
    #    R = build_ROC_curve(criterias = configuration('C0'), positive = 'both', ratio=('BPH','SPH'), thresholds = Z, num_of_buckets = 1, work_dir = f'results_{SP:s}/')
    #    A[SP] = R
        
    A = {'EUR': build_ROC_curve(criterias = configuration('C0'), positive = 'A', ratio=('transitions','disomy'), thresholds = Z, num_of_buckets = 1, work_dir = 'results_EUR/')}
    B = {'EUR': build_ROC_curve(criterias = configuration('C0'), positive = 'B', ratio=('transitions','disomy'), thresholds = Z, num_of_buckets = 1, work_dir = 'results_EUR/')}
    C = {'Multinomial classifier of haploids and diploids': build_ROC_curve(criterias = configuration('C0'), positive = 'both', ratio=('monosomy','disomy'), thresholds =  [i/50 for i in range(-10000,10000)], num_of_buckets = 1, work_dir = 'results_EUR/'),
         'Multinomial classifier of triploids and diploids': build_ROC_curve(criterias = configuration('C0'), positive = 'both', ratio=('BPH','disomy'), thresholds = Z, num_of_buckets = 1, work_dir = 'results_EUR/')}

    plot_ROC_curve_fig7b(C)
    plot_ROC_curve_fig7a(B)
    #plot_ROC_curve(C)

else:
    print("The module ROC_curve was imported.")

#from pathlib import Path
#LLR_files = [str(j) for j in Path('/home/ariad/Dropbox/postdoc_JHU/Tools/origin_V2/results_ZOUVES/').rglob('*LLR.p')]
#for filename in LLR_files:
#    LLR_dict, info = load_llr(filename); plot_single_case(LLR_dict, info, N=10, save=filename.replace('LLR.p','png'))


### END OF FILE ###
