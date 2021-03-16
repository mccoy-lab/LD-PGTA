#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BUILD_ROC_CURVES

Based on simulated data that was created by MIX_HAPLOIDS and analyzed by
ANEUPLOIDY_TEST, Balanced ROC curves for predicting aneuploidy are created.

The genome is divided into bins. For each simulated data, the mean LLR and the
standard deviation for a bin are calculated. A positive (negative) prediction
is made if the bounds of the confidence interval lie on the positive (negative)
side of the number line. In order to determine the optimal a confidence
level, the z-score is varied and the balanced positive and negative rates are
calculated for each value that is takes. For each z-score the the balanced
positive and negative rates are averaged across the bins. Then, averaged rates
are used for creating a balanced ROC curve.



Daniel Ariad (daniel@ariad.org)
March 20, 2021
"""

import pickle, statistics
from math import log
from itertools import starmap
from operator import itemgetter
from multiprocessing import Process

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

def mean_and_std(data):
    """ Calculates the mean and population standard deviation. """
    m = statistics.mean(data)
    std = statistics.pstdev(data, mu=m)
    return m, std 


def load_llr(filename):
    """ Loads from a file a dictionary that lists genomic windows that contain
    at least two reads and gives the bootstrap distribution of the 
    log-likelihood BPH/SPH ratios (LLRs). """
    try:
        with open(filename, 'rb') as f:
            likelihoods = pickle.load(f)
            info = pickle.load(f)
    except:
        likelihoods = info = None
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
    """ Lists the buckets and gives the genomic windows that they contain. """
    bin_size = chr_length(chr_id) / num_of_buckets
    result = {}
    j = 0
    
    for i in range(num_of_buckets): ### All buckets before the first the genomic window are filled with Nones.
        if sum(windows[0])/2 < (i+1)*bin_size:
            break
        result[i/num_of_buckets,(i+1)/num_of_buckets] = None
    
    for k,(a,b) in enumerate(windows):
        if not bin_size*i <= (a+b)/2 < bin_size*(i+1):
            result[i/num_of_buckets,(i+1)/num_of_buckets] = (j,k)
            j = k
            for i in range(i+1,num_of_buckets): #Proceed to the next non-empty bucket; Empty buckets are filled with Nones.
                if (a+b)/2 < (i+1)*bin_size:
                    break
                result[i/num_of_buckets,(i+1)/num_of_buckets] = None
    
    for i in range(i,num_of_buckets): ### All buckets after the last the genomic window are filled with Nones.
        result[i/num_of_buckets,(i+1)/num_of_buckets] = (j,k) if j != k else None
        j = k 
    return result


def bucketing(likelihoods,info,num_of_buckets,z,ratio):
    """ Bucketing is applied by aggregating the mean LLR per window across 
        consecutive windows. The boundaries of the bins as well as the mean LLR
        and the standard-error per bin are returned. """
    
    if ratio in info['statistics']['LLRs_per_genomic_window']:    
        LLR_stat = info['statistics']['LLRs_per_genomic_window'][ratio]
    else:
        _ = {}
        LLR_stat = {window:  mean_and_var([*starmap(LLR, ((_[ratio[0]], _[ratio[1]]) for _['monosomy'], _['disomy'], _['SPH'], _['BPH'] in L))])
                           for window,L in likelihoods.items()}
    
    K,M,V = tuple(LLR_stat.keys()), *zip(*LLR_stat.values())

    buckets = coordinates(windows=K,chr_id=info['chr_id'],num_of_buckets=num_of_buckets)
    
    X = [*buckets]
    Y = [statistics.mean(M[P[0]:P[1]]) if P else None for P in buckets.values()]
    E = [std_of_mean(V[P[0]:P[1]]) if P else None for P in buckets.values()] 
    
    return X,Y,E

def collect_data(criterias, num_of_buckets, ratio, work_dir):
    """ Iterates over all the data files in the folder and creates a dictionary
    the list all the files that fit the criterias and gives their analysis 
    (via the bucketing function). """
    
    import glob
    effective_ratio = tuple('BPH' if i=='transitions' else i for i in ratio)
    filenamesB = glob.glob(work_dir + f"simulated.{ratio[0]:s}*{criterias['depth']:.3f}*LLR.p")
    filenamesA = glob.glob(work_dir + f"simulated.{ratio[1]:s}*{criterias['depth']:.3f}*LLR.p")
    result = {r: {'chr'+str(i): {} for i in [*range(1,23)]+['X','Y'] } for r in ratio}
    for filename in filenamesA+filenamesB:
        likelihoods, info = load_llr(filename)
        if likelihoods==None: continue
        subinfo = {x: info.get(x,None) for x in criterias.keys()}
        ### print(subinfo)
        
        scenario = info.get('scenario','') 
        if criterias==subinfo and scenario in ratio:
            #show_info(filename,info)
            chr_id = info['chr_id']
            result[scenario][chr_id][filename] = {bucket: {'mean': mean, 'std': std} 
                                               for (bucket,mean,std) in zip(*bucketing(likelihoods,info,num_of_buckets=num_of_buckets[chr_id],z=1,ratio=effective_ratio))}
    return result

def prediction_rates(B, A, positive, thresholds, num_of_buckets):
    """ Creates a nested dictionary that lists chromosomes, buckets and
    thresholds. For each entry the dictionary gives the false and true positive
    rates. """ 
    
    prediction_rates = {'chr'+str(i): {} for i in [*range(1,23)]+['X','Y'] }
    
    for chr_id in prediction_rates:        
        print(len(A[chr_id]),len(B[chr_id]))
        for i in range(num_of_buckets[chr_id]):
            bucket = (i/num_of_buckets[chr_id],(i+1)/num_of_buckets[chr_id])
            prediction_rates[chr_id][bucket] = {}
            for z in thresholds:
                true_B = [file[bucket]['mean'] > z * file[bucket]['std'] for file in B[chr_id].values() if file[bucket]['mean']!=None]
                false_B = [file[bucket]['mean'] > z * file[bucket]['std'] for file in A[chr_id].values() if file[bucket]['mean']!=None]
                
                false_A = [file[bucket]['mean'] < -z * file[bucket]['std'] for file in B[chr_id].values() if file[bucket]['mean']!=None]
                true_A = [file[bucket]['mean'] < -z * file[bucket]['std'] for file in A[chr_id].values() if file[bucket]['mean']!=None]
                
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
                
                prediction_rates[chr_id][bucket][z] = (FPR,TPR)
                
    return prediction_rates

def average_prediction_rates(prediction_rates):
    f = lambda x: next(iter(x.values()))
    thresholds = f(f(prediction_rates))
    
    ###g = lambda x: {*zip(('FPR','TPR'), zip(*x))}
    g = lambda x: tuple(map(statistics.mean, zip(*x)))
    result = {z: g((prediction_rates[chr_id][bucket][z] for chr_id in prediction_rates for bucket in prediction_rates[chr_id] if z in prediction_rates[chr_id][bucket])) for z in thresholds }
        
    return result


def plot_ROC_curve(x):
    """ Plots the ROC curve. """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(*zip(*sorted(x.values())))
    #ax.legend()
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


def main(RATIO,SP,matched,mixed,C,depth,buckets_in_chr21):
    
    num_of_buckets = lambda n: {'chr'+str(i): n*chr_length('chr'+str(i))//chr_length('chr21') for i in [*range(1,23)]+['X','Y']}
    Z = [i/50 for i in [*range(-1800,-300,300)]+[*range(-300,300)]+[*range(300,1800,300)]  ]
    DATA = collect_data(criterias = configuration(C) , num_of_buckets = num_of_buckets(buckets_in_chr21), ratio = RATIO, work_dir  = f'/mybox/simulations/results{mixed:s}_{SP:s}/')
    R = prediction_rates(B = DATA[RATIO[0]], A = DATA[RATIO[1]], positive = 'both', thresholds = Z, num_of_buckets = num_of_buckets(buckets_in_chr21))
    T = average_prediction_rates(R)
    N = {'chr'+str(i): {RATIO[0]: len(DATA[RATIO[0]]['chr'+str(i)]), RATIO[1]: len(DATA[RATIO[1]]['chr'+str(i)])} for i in [*range(1,23)]+['X','Y']}
    with open(f'PREDICTED_RATES_FOR_{RATIO[0]:s}_vs_{RATIO[1]:s}_{matched:s}_{SP:s}_{depth:s}x_{buckets_in_chr21:d}BUCKETS.p' , "wb") as f:
        pickle.dump(R, f, protocol=4)
        pickle.dump(configuration(C), f, protocol=4)
        pickle.dump(N, f, protocol=4)
        pickle.dump(T, f, protocol=4)
                

if __name__ == "__main__":
    #RATIO = ('BPH','disomy')
    proc = []
    for RATIO in (('BPH','SPH'),('BPH','disomy'),('monosomy','disomy')):
        for SP in ('EUR','EAS','SAS','AMR','AFR'):
            for matched,mixed in (('MISMATCHED','_mixed'),('MATCHED','')):
                for C,depth in (('C0','0.01'),('C2','0.05'),('C1','0.1')):        
                    for buckets_in_chr21 in (5,10):
                        #if matched!='MISMATCHED' and depth!='0.01': continue
                        #main(RATIO,SP,matched,mixed,C,depth,buckets_in_chr21)
                        #buckets_in_chr21 = 5
                        try:
                            p = Process(target=main,args=(RATIO,SP,matched,mixed,C,depth,buckets_in_chr21))
                            p.start()
                            proc.append(p)
                        except:
                            print('caution: a process failed!')
    for p in proc:
        try:
            p.join()
        except:
            None
                        
        
        
                    
        #plot_ROC_curve(R['chr21'])

else:
    print("The module ROC_curve was imported.")

#from pathlib import Path
#LLR_files = [str(j) for j in Path('/home/ariad/Dropbox/postdoc_JHU/Tools/origin_V2/results_ZOUVES/').rglob('*LLR.p')]
#for filename in LLR_files:
#    LLR_dict, info = load_llr(filename); plot_single_case(LLR_dict, info, N=10, save=filename.replace('LLR.p','png'))


### END OF FILE ###
