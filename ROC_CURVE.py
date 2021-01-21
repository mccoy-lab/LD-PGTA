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

def confidence(info,N,z,ratio):
    """ Binning is applied by aggregating the mean LLR of a window across N
        consecutive windows. The boundaries of the bins as well as the mean LLR
        and the standard-error per bin are returned. """
        
    LLR_stat = info['statistics']['LLRs_per_genomic_window'][ratio]
     
    K,M,V = tuple(LLR_stat.keys()), *zip(*LLR_stat.values())
            
    i = lambda j: j*(len(V)//N)
    f = lambda j: (j+1)*(len(V)//N) if j!=N-1 else len(K)
    
    x = lambda p,q: (K[p][0],K[q-1][-1])
    y = lambda p,q: statistics.mean(M[p:q])
    e = lambda p,q: z * std_of_mean(V[p:q]) 
    
    X,Y,E = ([func(i(j),f(j)) for j in range(N)] for func in (x,y,e))

    return X,Y,E

def build_confidence_dict(criterias, num_of_buckets, ratio, work_dir):
    """ Iterates over all the data files in the folder and creates a dictionary
    the list all the files that fit the criterias and gives their analysis 
    (via the confidence function). """
    
    import glob
    filenames = glob.glob(work_dir + '*.LLR.p')
    result = {r: {} for r in ratio}
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
            result[scenario][filename] = tuple({'mean': mean, 'std': std} 
                                               for (pos,mean,std) in zip(*confidence(info,N=num_of_buckets,z=1,ratio=ratio)))
    return result

def build_ROC_curve(criterias, positive, ratio, thresholds, num_of_buckets, work_dir):
    """ Creates a nested dictionary that lists bins and thresholds and gives 
        the false and true positive rates. """ 
    
    DATA = build_confidence_dict(criterias, num_of_buckets, ratio, work_dir)
    B,A = DATA[ratio[0]], DATA[ratio[1]]
    
    print(len(A),len(B))
    result = {}
    for bucket in range(num_of_buckets):
        result[bucket] = {}
        for z in thresholds:
            true_B = [file[bucket]['mean'] > z * file[bucket]['std'] for file in B.values()]
            false_B = [file[bucket]['mean'] > z * file[bucket]['std'] for file in A.values()]
            
            false_A = [file[bucket]['mean'] < -z * file[bucket]['std'] for file in B.values()]
            true_A = [file[bucket]['mean'] < -z * file[bucket]['std'] for file in A.values()]
            
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
            
            result[bucket][z] = (FPR,TPR)
            
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

def plot_single_case(info,N,**kwargs):
    """ Plots the mean LLR vs. chromosomal position """
    import matplotlib.pyplot as plt
    save = kwargs.get('save', '')
    
    # Create lists for the plot    
    X,Y,E = confidence(info,N,z=1)
    C = [(x[1]+x[0])/2 for x in X]
    widths = [x[1]-x[0] for x in X]
    #X_boundaries = tuple(k for j in range(N+1) for k in (K[j*a][0], K[min((j+1)*a,len(K)-1)][-1]))
    X_ticks = [x[0] for x in X] + [X[-1][1]]
    X_labels = [('%.2f' % (j/chr_length(info['chr_id']))) for j in X_ticks] 
 
    # Build the plot
    fig, ax = plt.subplots(1, 1, figsize=(16, 9))  # setup the plot
    fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.07) 
    ax.bar(C, Y, yerr=E, align='center', alpha=0.5, ecolor='black', capsize=10, width=widths, color=[plt.cm.tab20(i) for i in range(N)])
    ax.set_ylabel('log-likelihood BPH/SPH ratio')
    ax.set_xlabel('Normalized chromosome position')
    ax.set_xticks(X_ticks)
    ax.set_xticklabels(X_labels)
    ax.set_title('.'.join(save.split('.')[:-2]))
    #ax.xaxis.grid(True)
    #ax.yaxis.grid(True)

    for j in range(N):
        plt.text(C[j], .5*(Y[j]-Y[j]/abs(Y[j])*E[j]), '%.2f\u00B1%.2f'% (Y[j],E[j]), ha='center', va='center',color='black',fontsize=8)
    
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
    


def plot_ROC_curve_band(R,num_of_buckets):
    """ Plots the ROC curve of all the bins as band. """
    from operator import itemgetter
    from collections import defaultdict
    from statistics import median, mean, pstdev
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    W = {i:{} for i in R}
    
    for i in R:
        for x,y in R[i].values():
            x_rounded = round(x,ndigits=2)
            y_max = W[i].get(x_rounded,0.0)
            if y>=y_max: W[i][x_rounded] = y
    
    Q = defaultdict(list)
    for bucket in W.values():
        for x,y in bucket.items():
            Q[x].append(y)
    
    X,Y,Ymin,Ymax = zip(*((x,median(y),min(y),max(y)) for x,y in sorted(Q.items(),key=itemgetter(0))))
    #X,Y,Ymin,Ymax = zip(*((x,mean(y),mean(y)-pstdev(y),mean(y)+pstdev(y)) for x,y in sorted(Q.items(),key=itemgetter(0))))
    
    
    for i in range(num_of_buckets):
        ax.scatter(*zip(*R[i].values()), label=f'bin {i:d}', s=num_of_buckets-i)        
    ax.plot([0]+[*X]+[1],[0]+[*Y]+[1],c='black')
    ax.fill_between(X, Ymin, Ymax, alpha=0.2)
    ax.legend()
    ax.grid(True)
    plt.show()

def configuration(C):
    C0 = {'chr_id': 'chr21',
         'depth': 0.01,
         'read_length': 36,
         'window_size': 0,
         'min_reads': 6,
         'max_reads': 4,
         'minimal_score': 2,
         'min_HF': 0.05}


    return locals()[f'{C:s}']

if __name__ == "__main__":    
    Z = [i/300 for i in range(-1200,1200)]
    R = build_ROC_curve(criterias = configuration('C0'), positive = 'both', ratio=('BPH','SPH'), thresholds = Z, num_of_buckets = 1, work_dir = 'results_EAS/')
    plot_ROC_curve(R)

else:
    print("The module ROC_curve was imported.")

#from pathlib import Path
#LLR_files = [str(j) for j in Path('/home/ariad/Dropbox/postdoc_JHU/Tools/origin_V2/results_ZOUVES/').rglob('*LLR.p')]
#for filename in LLR_files:
#    LLR_dict, info = load_llr(filename); plot_single_case(LLR_dict, info, N=10, save=filename.replace('LLR.p','png'))


### END OF FILE ###
