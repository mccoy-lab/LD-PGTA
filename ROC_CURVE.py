#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROC_CURVE

Daniel Ariad (daniel@ariad.org)
Dec 20, 2020
"""

import pickle, statistics

def chr_length(chr_id):
    """ Return the chromosome length for a given chromosome, based on the reference genome hg38.""" 
    #The data of chromosome length was taken from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38
    length_dict = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
                  'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422,
                  'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
                  'chr16': 90338345, 'chr17':  83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
                  'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}
    return length_dict[chr_id]
    
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

def coordinates(windows,chr_id,num_of_buckets):
    bin_size = chr_length(chr_id) / num_of_buckets
    result = {}
    j = 0
    for i in range(num_of_buckets): ### Fill all the bins before the first genomic window with Nones
        if sum(windows[0])/2 < (i+1)*bin_size:
            break
        result[i/num_of_buckets,(i+1)/num_of_buckets] = None
    
    for k,(a,b) in enumerate(windows):
        if not bin_size*i <= (a+b)/2 < bin_size*(i+1):
            result[i/num_of_buckets,(i+1)/num_of_buckets] = (j,k)
            j = k
            i += 1
    for l in range(i,num_of_buckets): 
        result[l/num_of_buckets,(l+1)/num_of_buckets] = (j,k) if j != k else None
        j = k 
    return result

def confidence(info,num_of_buckets,z,ratio):
    """ Binning is applied by aggregating the mean LLR of a window across N
        consecutive windows. The boundaries of the bins as well as the mean LLR
        and the standard-error per bin are returned. """
        
    LLR_stat = info['statistics']['LLRs_per_genomic_window'][ratio]
     
    K,M,V = tuple(LLR_stat.keys()), *zip(*LLR_stat.values())

    buckets = coordinates(windows=K,chr_id=info['chr_id'],num_of_buckets=num_of_buckets)
    
    X = [*buckets]
    Y = [statistics.mean(M[P[0]:P[1]]) if P else None for P in buckets.values()]
    E = [std_of_mean(V[P[0]:P[1]]) if P else None for P in buckets.values()] 
    
    return X,Y,E

def build_confidence_dict(criterias, num_of_buckets, ratio, work_dir):
    """ Iterates over all the data files in the folder and creates a dictionary
    the list all the files that fit the criterias and gives their analysis 
    (via the confidence function). """
    
    import glob
    effective_ratio = tuple('BPH' if i=='transitions' else i for i in ratio)
    filenames = glob.glob(work_dir + f"simulated.*.{criterias['chr_id']:s}.x{criterias['depth']:.3f}*.LLR.p")
    result = {r: {} for r in ratio}
    for filename in filenames:
        likelihoods, info = load_llr(filename)
        subinfo = {x: info.get(x,None) for x in criterias.keys()}
        ### print(subinfo)
        
        scenario = info.get('scenario','') 
        if criterias==subinfo and scenario in ratio:
            show_info(filename,info)
            result[scenario][filename] = {bucket: {'mean': mean, 'std': std} 
                                               for (bucket,mean,std) in zip(*confidence(info,num_of_buckets=num_of_buckets,z=1,ratio=effective_ratio))}
    return result

def build_ROC_curve(criterias, positive, ratio, thresholds, num_of_buckets, work_dir):
    """ Creates a nested dictionary that lists bins and thresholds and gives 
        the false and true positive rates. """ 
    
    DATA = build_confidence_dict(criterias, num_of_buckets, ratio, work_dir)
    B,A = DATA[ratio[0]], DATA[ratio[1]]
    
    buckets = [*next(iter(B.values()))]
    ###print(buckets)
    
    print(len(A),len(B))
    result = {}
    for bucket in buckets:
        result[bucket] = {}
        for z in thresholds:
            true_B = [file[bucket]['mean'] > z * file[bucket]['std'] for file in B.values() if file[bucket]['mean']!=None]
            false_B = [file[bucket]['mean'] > z * file[bucket]['std'] for file in A.values() if file[bucket]['mean']!=None]
            
            false_A = [file[bucket]['mean'] < -z * file[bucket]['std'] for file in B.values() if file[bucket]['mean']!=None]
            true_A = [file[bucket]['mean'] < -z * file[bucket]['std'] for file in A.values() if file[bucket]['mean']!=None]
            
            if true_B==[] or false_B==[]:
                break
            elif positive == 'both':
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
            
            result[bucket][z] = (FPR,TPR)
            
    return result

def plot_ROC_curve(x):
    """ Plots the ROC curve. """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    for i,j in enumerate(x):
        if x[j]!={}:
            ax.scatter(*zip(*x[j].values()), label=f'bin {str(j):s}', s=len(x)+1-i)        
    ax.legend()
    ax.grid(True)
    plt.show()

def configuration(C):
    
    result = dict(
        
    C0 = {'chr_id': 'chr21',
         'depth': 0.01,
         'read_length': 36,
         'window_size': 0,
         'min_reads': 6,
         'max_reads': 4,
         'minimal_score': 2,
         'min_HF': 0.05},
    
    C1 = {'chr_id': 'chr21',
         'depth': 0.1,
         'read_length': 36,
         'window_size': 0,
         'min_reads': 24,
         'max_reads': 12,
         'minimal_score': 2,
         'min_HF': 0.05},
    
    C2 = {'chr_id': 'chr21',
         'depth': 0.05,
         'read_length': 36,
         'window_size': 0,
         'min_reads': 18,
         'max_reads': 8,
         'minimal_score': 2,
         'min_HF': 0.05}
    )


    return result[C]

if __name__ == "__main__":    
    Z = [i/600 for i in range(-3600,3600)]
    R = build_ROC_curve(criterias = configuration('C0'), positive = 'both', ratio=('BPH','SPH'), thresholds = Z, num_of_buckets = 10, work_dir = 'results_EUR/')
    plot_ROC_curve(R)

else:
    print("The module ROC_curve was imported.")

#from pathlib import Path
#LLR_files = [str(j) for j in Path('/home/ariad/Dropbox/postdoc_JHU/Tools/origin_V2/results_ZOUVES/').rglob('*LLR.p')]
#for filename in LLR_files:
#    LLR_dict, info = load_llr(filename); plot_single_case(LLR_dict, info, N=10, save=filename.replace('LLR.p','png'))


### END OF FILE ###
