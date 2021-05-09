#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROC_CURVE_GENOME_WIDE

Daniel Ariad (daniel@ariad.org)
Dec 20, 2020
"""

import pickle, statistics, collections, bz2, gzip
from itertools import starmap
from operator import or_
from ANEUPLOIDY_TEST import LLR 
from multiprocessing import Process

def mean_and_var(sample):
    """ Calculates the mean and the sample standard deviation. """
    mean = statistics.mean(sample)
    var = statistics.variance(sample, xbar=mean)
    return mean, var

def std_of_mean(variances):
    """ Standard error of the mean of uncorrelated variables, based on the
        Bienaymé formula. """
    return sum(variances)**.5/len(variances)


def load_likelihoods(filename):
    """ Loads from a file a dictionary that lists genomic windows that contain
    at least two reads and gives the bootstrap distribution of the 
    log-likelihood ratios (LLRs). """
    
    Open = {'bz2': bz2.open, 'gzip': gzip.open}.get(filename.rpartition('.')[-1], open)
    
    with Open(filename, 'rb') as f:
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


def build_confidence_dict(criterias, ratio, work_dir):
    """ Iterates over all the data files in the folder and creates a dictionary
    the list all the files that fit the criterias and gives their analysis 
    (via the confidence function). """
    
    import glob
    effective_ratio = tuple('BPH' if i=='transitions' else i for i in ratio)
    ###filenames = glob.glob(work_dir + '*.LLR.p')
    ###filenames = glob.glob(work_dir + f"simulated*{criterias['depth']:.3f}*LLR.p")
    filenamesB = glob.glob(work_dir + f"simulated.{ratio[0]:s}*{criterias['depth']:.3f}*LLR.p")
    filenamesA = glob.glob(work_dir + f"simulated.{ratio[1]:s}*{criterias['depth']:.3f}*LLR.p")
    result = {r: {} for r in ratio}

    buffer = {r: collections.defaultdict(dict) for r in ratio}
    for filename in filenamesA+filenamesB:
        likelihoods, info = load_likelihoods(filename)
        if likelihoods==None: continue
        subinfo = {x: info.get(x,None) for x in criterias.keys()}
        ### print(subinfo)
        
        scenario = info.get('scenario','')
        if criterias==subinfo and scenario in ratio:
            #show_info(filename,info)  
        
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
                 MEANs, VARs = zip(*LLRs)
                 M.extend(MEANs)
                 V.extend(VARs)            
            result[scenario][i] = {'mean': statistics.mean(M), 'std': std_of_mean(V)}
                                               
    
    return result

def build_ROC_curve(B, A, positive, thresholds):
    """ Creates a nested dictionary that lists bins and thresholds and gives 
        the false and true positive rates. """ 
    
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
    

def configuration(C):
    
    result = dict(
        
    C0 = {'depth': 0.01,
         'read_length': 36,
         'window_size': 0,
         'min_reads': 6,
         'max_reads': 4,
         'minimal_score': 2,
         'min_HF': 0.05},
    
    C1 = {'depth': 0.1,
         'read_length': 36,
         'window_size': 0,
         'min_reads': 24,
         'max_reads': 12,
         'minimal_score': 2,
         'min_HF': 0.05},
    
    C2 = {'depth': 0.05,
         'read_length': 36,
         'window_size': 0,
         'min_reads': 18,
         'max_reads': 8,
         'minimal_score': 2,
         'min_HF': 0.05}
    )

    return result[C]

def main(RATIO,SP,matched,mixed,C):
    conf = configuration(C)
    DATA = build_confidence_dict(criterias = conf, ratio=RATIO, work_dir = f'/mybox/simulations/results{mixed:s}_{SP:s}/')
    B,A = DATA['transitions'], DATA['disomy']
    R = build_ROC_curve(B, A, positive = 'A', thresholds = Z)
    with open(f"PREDICTED_RATES_FOR_{RATIO[0]:s}_vs_{RATIO[1]:s}_{matched:s}_{SP:s}_{conf['depth']:g}x.p" , "wb") as f:
        pickle.dump(R, f, protocol=4)
        pickle.dump(configuration(C), f, protocol=4)
        pickle.dump({'transitions':len(B), 'disomy':len(A)}, f, protocol=4)
                

if __name__ == "__main__":
    proc = []
    RATIO = ('transitions','disomy')    
    Z = [i/10 for i in range(-100000,300)] + [i/50 for i in range(-300,300)] + [i/10 for i in range(300,100000)]# Z = [i/500 for i in range(-10000,10000)]
    for C in ('C0','C2','C1'):
        for matched,mixed in (('MISMATCHED','_mixed'),('MATCHED','')):
            for SP in ('EUR','EAS','SAS','AFR','AMR'):
                print(RATIO,SP,matched,mixed,C)
                ###main(RATIO,SP,matched,mixed,C,depth)
                try:
                    p = Process(target=main,args=(RATIO,SP,matched,mixed,C))
                    p.start()
                    proc.append(p)
                except:
                    print('caution: a process failed!')
    for p in proc:
        try:
            p.join()
        except:
            None  
                
    #A = {'EUR': build_ROC_curve(criterias = configuration('C0'), positive = 'A', ratio=('transitions','disomy'), thresholds = Z, num_of_buckets = 1, work_dir = 'results_EUR/')}
    #B = {'EUR': build_ROC_curve(criterias = configuration('C0'), positive = 'B', ratio=('transitions','disomy'), thresholds = Z, num_of_buckets = 1, work_dir = 'results_EUR/')}
    #C = {'Multinomial classifier of haploids and diploids': build_ROC_curve(criterias = configuration('C0'), positive = 'both', ratio=('monosomy','disomy'), thresholds =  [i/50 for i in range(-10000,10000)], num_of_buckets = 1, work_dir = 'results_EUR/'),
    #     'Multinomial classifier of triploids and diploids': build_ROC_curve(criterias = configuration('C0'), positive = 'both', ratio=('BPH','disomy'), thresholds = Z, num_of_buckets = 1, work_dir = 'results_EUR/')}

    #plot_ROC_curve_fig7b(C)
    #plot_ROC_curve_fig7a(B)
    #plot_ROC_curve(C)

else:
    print("The module ROC_curve was imported.")

#from pathlib import Path
#LLR_files = [str(j) for j in Path('/home/ariad/Dropbox/postdoc_JHU/Tools/origin_V2/results_ZOUVES/').rglob('*LLR.p')]
#for filename in LLR_files:
#    LLR_dict, info = load_llr(filename); plot_single_case(LLR_dict, info, N=10, save=filename.replace('LLR.p','png'))


### END OF FILE ###
