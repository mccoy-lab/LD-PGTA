#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

NO_NAME

Builds a dictionary that lists linkage disequilibrium (LD) blocks that contain
at least two reads and gives the associated log-likelihood BPH/SPH ratio (LLR).
BPH (Both Parental Homologs) correspond to the presence of three unmatched
haplotypes, while SPH (Single Parental Homolog) correspond to chromosome gains
involving identical homologs.

Daniel Ariad (daniel@ariad.org)
May 11st, 2020
"""
import collections,time,pickle,statistics,argparse,re,sys,itertools

from LLR_CALCULATOR import read_impute2
from LLR_CALCULATOR import wrapper_func_of_create_LLR_2 as get_LLR

def overlaps(obs_filename, leg_filename):
    """ Given an observation table, the fraction of overlapping reads
        is calculated; Altough the reads overlap in position, they do not
        necessary have common alleles. """
    
    leg_tab = read_impute2(leg_filename, filetype='legend')
    obs_tab = pickle.load(open(obs_filename, 'rb'))
    
    aux_dict = collections.defaultdict(list)
    for (pos, ind, read_id, base) in obs_tab:
        if base in leg_tab[ind][2:]:
            aux_dict[pos].append(read_id)   
            
    X = [set(x) for x in aux_dict.values()]

    for i in range(len(X)-1,-1,-1):
        for j in range(i-1,-1,-1):
            if not X[i].isdisjoint(X[j]):
                X[j].update(X.pop(i))
                break
    
    OVERLAPS = sum(len(x) for x in X if len(x)>1)
    TOTAL = len(set(itertools.chain.from_iterable(aux_dict.values())))
    return OVERLAPS/TOTAL
            

def mean_and_std(sample):
    """ Calculates the mean and standard deviation of normally distributed
        random variables. """
    mean = statistics.mean(sample)
    std = statistics.pstdev(sample, mu=mean)/len(sample)**.5
    return mean, std

def jackknifing(sample,weights):
    """ Given sample elements and the weight of each element, the jackknife
    estimator, the jackknife standard error and the Jackknife bias (RB)
    are calculated. More information about delete-m jackknife for unequal m
    can be found in F.M.Busing et al. (1999), [DOI:10.1023/A:1008800423698]. """
    
    N = len(sample)
    t0 = sum(sample) / N
    H = [1/w for w in weights]
    ###H = [N for w in weights]
    T = [sum(sample[:i]+sample[i+1:])/(N-1) for i in range(N)] 
    pseudo_values = [h*t0-(h-1)*t for t,h in zip(T,H)] 
    jackknife_estimator = sum((p/h for p,h in zip(pseudo_values,H)))
    jackknife_variance = sum(((p-jackknife_estimator)**2/(h-1) for p,h in zip(pseudo_values,H)))/N
    jackknife_bias = (N-1)*t0-sum((t*(h-1)/h for t,h in zip(T,H)))
    return jackknife_estimator, jackknife_variance**.5 , jackknife_bias

def build_aux_dict(obs_tab,leg_tab):
    """ Returns a dictionary that lists chromosome positions of SNPs and gives
        a dictionary that lists the observed bases at these postions. The
        nested dictionary gives all the read IDs that contain the observed base
        at the specific chromosome position. """
    
    nested = lambda: collections.defaultdict(list)
    aux_dict = collections.defaultdict(nested)

    for (pos, ind, read_id, base) in obs_tab:
        if base in leg_tab[ind][2:]:
            aux_dict[pos][base].append(read_id)  
       
    return aux_dict 


def group_alleles(aux_dict,positions,**thresholds):
    """ For each chromosome position within the tuple positions, all the
    observed alleles together with their associted read id are extracted from
    the dictionary aux_dict, which was created by build_aux_dict. 
    Then, all the alleles are grouped into haplotypes (essetially tuples),
    according to the reads they origined from. All the haplotypes together form
    a tuple, named HAPLOTYPES, which is sorted in descending order according to
    the number of alleles in each haplotype. Moreover, only the 16 longest
    haplotypes are kept in HAPLOTYPES.
    In addition, overlapping haplotypes, which have common alleles, are
    identified. For each group of overlapping haplotypes a set of indices is
    formed; The haplotypes are indexed by their position in the tuple
    HAPLOTYPES. All these sets are kept within a tuple named OVERLAPPING_HAPLOTYPES.
    Finally, the tuples HAPLOTYPES and OVERLAPPING_HAPLOTYPES are returned. """
    
    reads = collections.defaultdict(list)
    for pos in positions:
        for base in aux_dict[pos]:
            for read_id in aux_dict[pos][base]:
                reads[read_id].append((pos,base))
    
    ######################## APPLY THRESHOLDS AND FILTERS #####################    
    if thresholds.get('min_alleles_per_block',None):
        total_num_of_alleles = sum(len(alleles) for alleles in reads.values())
        if total_num_of_alleles<thresholds['min_alleles_per_block']: return None
    
    if thresholds.get('min_alleles_per_read',None):    
        for read_id,alleles in tuple(reads.items()):
            if len(alleles)<thresholds['min_alleles_per_read']: del reads[read_id]
    
    if thresholds.setdefault('min_reads_per_block',2):        
        if len(reads)<thresholds['min_reads_per_block']: return None
    ###########################################################################
    
    
    read_IDs, HAPLOTYPES = zip(*sorted(reads.items(), key=lambda x: len(x[1]), reverse=True)[:16])

    
    X = [{read_IDs.index(read_id) for read_id in aux_dict[pos][allele] if read_id in read_IDs}
                for pos in positions
                    for allele in aux_dict[pos]]
    
    for i in range(len(X)-1,-1,-1):
        for j in range(i-1,-1,-1):
            if not X[i].isdisjoint(X[j]):
                X[j].update(X.pop(i))
                break
    
    OVERLAPPING_HAPLOTYPES = tuple(x for x in X if len(x)>1)
            
    return HAPLOTYPES, OVERLAPPING_HAPLOTYPES
    
def build_blocks_dict(positions,block_size,offset):
    """ Returns a dictionary that lists blocks and gives all the SNP positions
        within them."""
    
    a = positions[0]-(block_size-1)+offset
    b = positions[-1]+(block_size-1)
    boundaries = tuple(range(int(a), int(b), int(block_size)))
    blocks = [(i,j-1) for i,j in zip(boundaries,boundaries[1:])]
    
    blocks_dict = collections.defaultdict(list)
    
    blocks_iterator = iter(blocks)
    block = next(blocks_iterator, None)
    for p in positions:
        while block:
            if block[0]<=p<=block[1]:
                blocks_dict[block].append(p)
                break
            block = next(blocks_iterator, None)
    return blocks_dict   

def main(obs_filename,leg_filename,hap_filename,block_size,offset,output_filename,**thresholds):
    """ Returns a dictionary that lists the boundaries of approximately
    independent blocks of linkage disequilibrium (LD). For each LD block it
    gives the associated log-likelihood BPH/SPH ratio (LLR)."""
    
    a = time.time()

    hap_tab = read_impute2(hap_filename, filetype='hap')
    leg_tab = read_impute2(leg_filename, filetype='legend')
    
    with open(obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        info = pickle.load(f)
      
    aux_dict = build_aux_dict(obs_tab, leg_tab)
        
    blocks_dict = build_blocks_dict(tuple(aux_dict.keys()),block_size,offset)
    blocks_dict_grouped = {block: group_alleles(aux_dict,positions,**thresholds) for block,positions in blocks_dict.items()}
    LLR = get_LLR(obs_tab, leg_tab, hap_tab)
    LLR_dict = {block: LLR(*arg) if arg!=None else None for block,arg in blocks_dict_grouped.items()}
    
    population = tuple(value for value in LLR_dict.values() if value!=None)    
    mean, std = mean_and_std(population)
    
    w0 = [len(blocks_dict[key]) for key,value in LLR_dict.items() if value!=None]
    W = sum(w0)
    weights = [i/W for i in w0]
    jk_mean, jk_std, jk_bias = jackknifing(population,weights)
    
    num_of_LD_blocks = len(population)
    fraction_of_negative_LLRs = sum([1 for i  in population if i<0])/len(population)
    
    print('Mean: %.3f, STD: %.3f' % (mean, std))
    print('Jackknife estimator: %.3f, Jackknife standard error: %.3f, Jackknife bias: %.3f' % (jk_mean, jk_std, jk_bias))
    print('Number of LD blocks: %d, Fraction of LD blocks with a negative LLR: %.3f' % (num_of_LD_blocks,fraction_of_negative_LLRs))
    
    info.update({'block_size': block_size, 'offset': offset, 'thresholds': thresholds})
    info['statistics'] = {'mean': mean, 'std': std, 'jk_mean': jk_mean,
                          'jk_std': jk_std, 'jk_bias': jk_bias,
                          'num_of_LD_blocks': num_of_LD_blocks,
                          'fraction_of_negative_LLRs': fraction_of_negative_LLRs} 
    
    if output_filename!=None:
        default_filename = re.sub('(.*)OBS','\\1LLR', obs_filename.split('/')[-1],1)    
        output_filename = default_filename if output_filename=='' else output_filename 
        with open( output_filename, "wb") as f:
            pickle.dump(obs_tab, f, protocol=4)
            pickle.dump(info, f, protocol=4)    
    
    b = time.time()
    print('Done calculating LLRs for all the LD block in %.3f sec.' % ((b-a)))
    return LLR_dict, info

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
    description='Builds a dictionary that lists linkage disequilibrium (LD) '
                'blocks that contain at least two reads and gives the '
                'associated log-likelihood BPH/SPH ratio (LLR). BPH (Both '
                'Parental Homologs) correspond to the presence of three '
                'unmatched haplotypes, while SPH (Single Parental Homolog) '
                'correspond to chromosome gains involving identical homologs.')
    parser.add_argument('obs_filename', metavar='observations_filename', type=str, 
                        help='A pickle file created by MAKE_OBS_TAB, containing base observations at known SNP positions.')
    parser.add_argument('leg_filename', metavar='legend_filename', type=str, 
                        help='IMPUTE2 legend file')
    parser.add_argument('hap_filename', metavar='haplotype_filename', type=str, 
                        help='IMPUTE2 haplotype file')
    parser.add_argument('-s', '--block-size', type=int, 
                        metavar='INT', default='100000', 
                        help='Specifies the typical size of the LD block. The default value is 10^5.')
    parser.add_argument('-f', '--offset', type=int, 
                        metavar='INT', default=0,
                        help='Shifts all the LD blocks by the requested base pairs. The default value is 0.')
    parser.add_argument('-a', '--min-alleles-per-block', type=int, metavar='INT', default=0,
                        help='Takes into consideration only LD blocks with at least INT alleles. The default value is 0.')
    parser.add_argument('-b', '--min-alleles-per-read', type=int, metavar='INT', default=0,
                        help='Uses only reads that contain at least INT alleles. The default value is 0.')
    parser.add_argument('-c', '--min-reads-per-block', type=int, metavar='INT', default=2,
                        help='Only LD blocks with at least INT reads are taken into account. The default value is 2.'
                            'Filters a,b and c are applied sequentially, one after the other.')
    parser.add_argument('-o', '--output-filename', type=str, metavar='output_filename',  default='',
                        help='The output filename. The default is the input filename with the extension \".OBS.p\" replaced by \".LLR#.p\".')
    args = parser.parse_args()
    
    
    LLR_dict, info = main(**vars(args))
            
    sys.exit(0)
else: 
    print("Executed when imported")

"""     
if __name__ == "__main__": 
    print("Executed when invoked directly")
    a = time.time()
    
    args = dict(obs_filename = 'results/mixed2haploids.X0.1.SRR10393062.SRR151495.0-2.hg38.OBS.p',
                hap_filename = '../build_reference_panel/ref_panel.HapMix.hg38.BCFtools/chr21_HapMix_panel.hap',
                leg_filename = '../build_reference_panel/ref_panel.HapMix.hg38.BCFtools/chr21_HapMix_panel.legend',
                block_size = 1e5,
                offset = 0,
                output_filename = None)
    
    #LLR_dict, info = main(**args) 
    #models_filename = 'MODELS.p'
    
    
    hap_tab = read_impute2(args['hap_filename'], filetype='hap')
    leg_tab = read_impute2(args['leg_filename'], filetype='legend')
    
    with open(args['obs_filename'], 'rb') as f:
        obs_tab = pickle.load(f)
        info = pickle.load(f)
      
    aux_dict = build_aux_dict(obs_tab, leg_tab)
    
    args = {'positions': tuple(aux_dict.keys()),
            'block_size': 1e5,
            'offset': 0 } 
    
    blocks_dict = build_blocks_dict(**args)
    
    thresholds = dict(min_alleles_per_block = 1,
                      min_alleles_per_read = 1,
                      min_reads_per_block = 4)
    
    blocks_dict_grouped = {block: group_alleles(aux_dict,positions,**thresholds) for block,positions in blocks_dict.items()}
    LLR = get_LLR(obs_tab, leg_tab, hap_tab)
    LLR_dict = {block: LLR(*arg) if arg!=None else None for block,arg in blocks_dict_grouped.items()}
    
    population = tuple(value for value in LLR_dict.values() if value!=None)    
    mean, std = mean_and_std(population)
    
    w0 = [len(blocks_dict[key]) for key,value in LLR_dict.items() if value!=None]
    W = sum(w0)
    weights = [i/W for i in w0]
    jk_mean, jk_std, jk_bias = jackknifing(population,weights)
    
    num_of_LD_blocks = len(population)
    fraction_of_negative_LLRs = sum([1 for i  in population if i<0])/len(population)
    
    print('Mean: %.3f, STD: %.3f' % (mean, std))
    print('Jackknife estimator: %.3f, Jackknife standard error: %.3f, Jackknife bias: %.3f' % (jk_mean, jk_std, jk_bias))
    print('Number of LD blocks: %d, Fraction of LD blocks with a negative LLR: %.3f' % (num_of_LD_blocks,fraction_of_negative_LLRs))
    
        
    #A = group_alleles(aux_dict,positions[:20])
    
    b = time.time()

    print('Done in %.3f sec.' % ((b-a)))
"""
