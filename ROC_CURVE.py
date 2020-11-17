#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Daniel Ariad (daniel@ariad.org)
Aug 31, 2020
"""

import pickle, statistics, collections

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
    with open(filename, 'rb') as f:
        LLR_dict = pickle.load(f)
        info = pickle.load(f)
    
    print('\nFilename: %s' % filename)
    print('Depth: %.2f, Number of LD blocks: %d, Fraction of LD blocks with a negative LLR: %.3f' % (info['depth'], info['statistics']['num_of_LD_blocks'],info['statistics']['fraction_of_negative_LLRs']))
    print('Mean and standard error of meaningful reads per LD block: %.1f, %.1f.' % (info['statistics']['reads_mean'], info['statistics'].get('reads_std',info['statistics'].get('reads_var'))))
    print('Mean LLR: %.3f, Standard error of the mean LLR: %.3f' % ( info['statistics']['mean'], info['statistics']['std']))
    print('Calculation was done in %.3f sec.' % info['runtime'])

    return LLR_dict, info

def rate(x):
    """ Return the fraction of True elements. """
    return x.count(True)/len(x)

def confidence(LLR_dict,info,N,z):
    TEST = [num_of_reads>info['max_reads']  
            for num_of_reads in info['statistics']['reads_per_LDblock_dict'].values() if num_of_reads>=info['min_reads'] ]    
    
    WEIGHTS = [min(num_of_reads-1,info['max_reads'])   
               for num_of_reads in info['statistics']['reads_per_LDblock_dict'].values() if num_of_reads>=info['min_reads']]    
    
    LLR_stat = {block: mean_and_var(LLRs)  
                for block,LLRs in LLR_dict.items() if None not in LLRs}
    
    K,M,V = tuple(LLR_stat.keys()), *zip(*LLR_stat.values())
            
    i = lambda j: j*(len(V)//N)
    f = lambda j: (j+1)*(len(V)//N) if j!=N-1 else len(K)-1
    
    x = lambda p,q: .5*(K[p][0] + K[q][-1])
    y = lambda p,q: statistics.mean(M[p:q])
    e = lambda p,q: z * (std_of_mean(V[p:q]) if rate(TEST[p:q])>0.85 else jackknife_std(M[p:q],WEIGHTS[p:q]))
    
    X,Y,E = ([func(i(j),f(j)) for j in range(N)] for func in (x,y,e))

    return X,Y,E

def build_confidence_dict(criterias, num_of_buckets):
    import glob
    filenames = glob.glob("ROC/*.LLR.p")
    result = {'SPH': {}, 'BPH': {}}
    for filename in filenames:
        LLR_dict, info = load_llr(filename)
        subinfo = {x: info.get(x,None) for x in ('chr_id','depth','read_length','block_size','min_reads','minimal_score','max_reads','min_HF')}
        if criterias==subinfo:
            if (info.get('scenario',None)=='BPH' and info.get('recombination_spot',None)==1.00) or info.get('scenario',None)=='SPH':
                scenario = 'SPH'
            elif info.get('scenario',None)=='BPH' and info.get('recombination_spot',None)==0.00: 
                scenario = 'BPH'
            else:
                continue
            #print(scenario)
            result[scenario][filename] = tuple({'mean': mean, 'std': std} 
                                               for (pos,mean,std) in zip(*confidence(LLR_dict,info,N=num_of_buckets,z=1)))
    return result

def build_ROC_curve(criterias, positive, thresholds, num_of_buckets):
    SPH, BPH = build_confidence_dict(criterias, num_of_buckets).values()
    result = {}
    for bucket in range(num_of_buckets):
        result[bucket] = {}
        for z in thresholds:
            true_BPH = [file[bucket]['mean'] > z * file[bucket]['std'] for file in BPH.values()]
            false_BPH = [file[bucket]['mean'] > z * file[bucket]['std'] for file in SPH.values()]
            
            false_SPH = [file[bucket]['mean'] < -z * file[bucket]['std'] for file in BPH.values()]
            true_SPH = [file[bucket]['mean'] < -z * file[bucket]['std'] for file in SPH.values()]
            
            if positive == 'both':
                TPR = 0.5 * (true_BPH.count(1)/len(true_BPH) + true_SPH.count(1)/len(true_SPH))
                FPR = 0.5 * (false_BPH.count(1)/len(false_BPH) + false_SPH.count(1)/len(false_SPH))
            elif positive == 'SPH':
                TPR = true_SPH.count(1)/len(true_SPH)
                FPR = false_SPH.count(1)/len(false_SPH)
            elif positive == 'BPH':
                TPR = true_BPH.count(1)/len(true_BPH)
                FPR = false_BPH.count(1)/len(false_BPH)
            else:
                break
            
            result[bucket][z] = (FPR,TPR)
            
    return result
    
def plot(x,num_of_buckets):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    for i in range(num_of_buckets):
        ax.scatter(*zip(*x[i].values()), label=f'bin {i:d}', s=num_of_buckets-i)        
    ax.legend()
    ax.grid(True)
    plt.show()

if __name__ == "__main__":
    C0 = {'chr_id': 'chr21',
         'depth': 0.05,
         'read_length': 35,
         'block_size': 0,
         'min_reads': 4,
         'max_reads': 16,
         'minimal_score': 2,
         'min_HF': 0.15}

    C1 = {'chr_id': 'chr21',
          'depth': 0.1,
          'read_length': 35,
          'block_size': 100000.0,
          'min_reads': 6,
          'max_reads': 14,
          'minimal_score': 2,
          'min_HF': 0.15}

    C2 = {'chr_id': 'chr21',
          'depth': 0.5,
          'read_length': 35,
          'block_size': 100000.0,
          'min_reads': 6,
          'max_reads': 14,
          'minimal_score': 2,
          'min_HF': 0.15}

    Z = [i/300 for i in range(1200)]
    R = build_ROC_curve(criterias = C0, positive = 'both', thresholds = Z, num_of_buckets = 10)
    plot(R, num_of_buckets = 10)
else:
    print("The module ROC_curve was imported.")

### END OF FILE ###
        