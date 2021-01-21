#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROC_CURVE_GENOME_WIDE

Daniel Ariad (daniel@ariad.org)
Dec 20, 2020
"""

import pickle, statistics, collections

    
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
    filenames = glob.glob(work_dir + '*.LLR.p')
    result = {r: {} for r in ratio}

    buffer = {r: collections.defaultdict(dict) for r in ratio}
    for filename in filenames:
        likelihoods, info = load_llr(filename)
        subinfo = {x: info.get(x,None) for x in criterias.keys()}
        ### print(subinfo)
        
        S = info.get('scenario','').upper()
        if S in ('SPH','DISOMY','MONOSOMY'):
            scenario = S
        elif S=='BPH' and info.get('recombination_spot',None)==1.00:
            scenario = 'SPH'
        elif S=='BPH' and info.get('recombination_spot',None)==0.00:
            scenario = 'BPH'
        else:
            scenario = None
        
        if criterias==subinfo and scenario in ratio:
            show_info(filename,info)  
            buffer[scenario][info['chr_id']][filename] = [*info['statistics']['LLRs_per_genomic_window'][ratio].values()]
    
     ### buffer ---> result     
    
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
    
    result = {}

    for z in thresholds:
        true_B = [file['mean'] > z * file['std'] for file in B.values()]
        false_B = [file['mean'] > z * file['std'] for file in A.values()]
        
        false_A = [file['mean'] < -z * file['std'] for file in B.values()]
        true_A = [file['mean'] < -z * file['std'] for file in A.values()]
        
        if positive == 'both':
            TPR = 0.5 * (true_B.count(1)/len(true_B) + true_A.count(1)/len(true_A))
            FPR = 0.5 * (false_B.count(1)/len(false_B) + false_A.count(1)/len(false_A))
        elif positive == 'A':
            TPR = true_A.count(1)/len(true_A)
            FPR = false_A.count(1)/len(false_A)
        elif positive == 'B':
            TPR = true_B.count(1)/len(true_B)
            FPR = false_B.count(1)/len(false_B)
        else:
            break
        
        result[z] = (FPR,TPR)
            
    return result

def plot_ROC_curve(x):
    """ Plots the ROC curve. """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    for k,v in x.items():
        ax.scatter(*zip(*v.values()), label=f'{k:s}', s=10)        
    ax.legend()
    ax.grid(True)
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


    return locals()[f'{C:s}']

if __name__ == "__main__":    
    Z = [i/20 for i in range(-1000,1000)]
    A = {}
    for SP in ('EUR','EAS','SAS','AFR','AMR'):
        R = build_ROC_curve(criterias = configuration('C0'), positive = 'both', ratio=('BPH','SPH'), thresholds = Z, num_of_buckets = 1, work_dir = f'results_{SP:s}/')
        A[SP] = R
    plot_ROC_curve(A)

else:
    print("The module ROC_curve was imported.")

#from pathlib import Path
#LLR_files = [str(j) for j in Path('/home/ariad/Dropbox/postdoc_JHU/Tools/origin_V2/results_ZOUVES/').rglob('*LLR.p')]
#for filename in LLR_files:
#    LLR_dict, info = load_llr(filename); plot_single_case(LLR_dict, info, N=10, save=filename.replace('LLR.p','png'))


### END OF FILE ###
