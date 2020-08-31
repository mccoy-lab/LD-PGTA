#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANEUPLOIDY_TEST
Builds a dictionary that lists linkage disequilibrium (LD) blocks that contain
at least two reads and gives the associated log-likelihood BPH/SPH ratio (LLR).
BPH (Both Parental Homologs) correspond to the presence of three unmatched
haplotypes, while SPH (Single Parental Homolog) correspond to chromosome gains
involving identical homologs.
Daniel Ariad (daniel@ariad.org)
Aug 31, 2020
"""

import collections, time, pickle, statistics, argparse, re, sys, operator, heapq, itertools, functools, random
from MAKE_OBS_TAB import read_impute2
from LLR_CALCULATOR import wrapper_func_of_create_LLR as get_LLR

def mean_and_var(sample):
    """ Calculates the mean and standard deviation of normally distributed
        random variables. """
    mean = statistics.mean(sample)
    var = statistics.pvariance(sample, mu=mean)
    return mean, var

def build_reads_dict(obs_tab,leg_tab):
    """ Returns a dictionary that lists read IDs of reads that overlap with
        SNPs and gives the alleles in each read. """

    reads = collections.defaultdict(list)

    for (pos, ind, read_id, base) in obs_tab:
        if base in leg_tab[ind][2:]:
            reads[read_id].append((pos,base))

    return reads

def build_weight_dict(reads_dict,obs_tab,leg_tab, hap_tab):
    """ Returns a dicitonary lists read_IDs and gives their weight. The weight of a
        read is inherited from the weight of the SNPs that are overlapped by the
        read. Thus, the weight of a read depends on its starting position and
        length, but not on the alleles that it contains. """

    def bools2int(x):
        """ Transforms a tuple/list of bools to a int. """
        return int(''.join(str(int(i)) for i in x),2)

    def weight(X,N):
        """ Returns the weight of a given haplotypes. """
        threshold = 0.05
        joint_freq = bin(functools.reduce(operator.and_,X)).count('1') / N
        result = threshold<=joint_freq<=(1-threshold)
        return result

    N = len(hap_tab[0])

    hap_dict = dict()
    for (pos, ind, read_id, base) in obs_tab:
        ref, alt = leg_tab[ind][2:]
        if (base==alt or base==ref) and (0.01<hap_tab[ind].count(1)/N<0.99): # if a biallelic SNP admits an allele with a frequency close to zero then this SNP would not contribute significantly to the weight of the read. However, including these SNPs in the calculation of the weight would prolong the calculation time. Thus, this function determines if the SNP should be included the calculation of the weight. 
            hap_dict[pos] = (bools2int(hap_tab[ind]), bools2int(map(operator.not_,hap_tab[ind])))
 
    weight_dict = dict()
    for read_id in reads_dict:
        hap = (hap_dict[pos] for pos,base in reads_dict[read_id] if pos in hap_dict)
        weight_dict[read_id] = sum(weight(C,N) for C in itertools.product(*hap) if len(C)!=0)

    return weight_dict

def pick_reads(reads_dict,weight_dict,read_IDs,min_reads,max_reads):
    """ Draws up to max_reads reads from a given LD block; Only reads with a
        weight larger than one are sampled. In addition, if the number of reads
        in a given LD block is less than the minimal requirment then the block
        would not be considered."""

    read_IDs_filtered = tuple(read_ID for read_ID in read_IDs if weight_dict[read_ID]>1)
    drawn_read_IDs = random.sample(read_IDs_filtered, min(max_reads,len(read_IDs_filtered)))
    haplotypes = tuple(reads_dict[read_ID] for read_ID in drawn_read_IDs)
    result = haplotypes if len(haplotypes) >= max(2,min_reads) else None
    return result
    
def iter_blocks(obs_tab,leg_tab,block_size,offset):
    """ Returns an iterator over the LD blocks together with read IDs of
        the reads that overlap with SNPs in the block. """

    aux_dict = collections.defaultdict(list) ### aux_dict is a dictionary that lists chromosome positions of SNPs and gives a list of read IDs for all the reads that overlap with the SNP.  

    for (pos, ind, read_id, base) in obs_tab:
        if base in leg_tab[ind][2:]:
            aux_dict[pos].append(read_id)
            
    first, *_, last = iter(aux_dict)
    boundaries = tuple(range(int(first+offset), int(last), int(block_size)))
    blocks_iterator = ((i,j-1) for i,j in zip(boundaries,boundaries[1:]))
   
    readIDs_in_block = set()
    block = next(blocks_iterator, None)
    
    for pos in aux_dict:
        if pos<boundaries[0]: continue
        while block:
            if block[0]<=pos<=block[1]:
                for read_IDs in aux_dict[pos]: readIDs_in_block.add(read_IDs)
                break
            yield (block, readIDs_in_block)
            block, readIDs_in_block = next(blocks_iterator, None), set()

def aneuploidy_test(obs_filename,leg_filename,hap_filename,block_size,subsamples,offset,min_reads,max_reads,output_filename):
    """ Returns a dictionary that lists the boundaries of approximately
    independent blocks of linkage disequilibrium (LD). For each LD block it
    gives the associated log-likelihood BPH/SPH ratio (LLR)."""

    a = time.time()
    random.seed(a=None, version=2) #I should set a=None after finishing to debug the code.

    with open(obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        info = pickle.load(f)
        
    hap_tab = read_impute2(hap_filename, filetype='hap')
    leg_tab = read_impute2(leg_filename, filetype='leg')

    reads_dict = build_reads_dict(obs_tab,leg_tab)
    weight_dict = build_weight_dict(reads_dict,obs_tab,leg_tab,hap_tab)
    LLR = get_LLR(obs_tab, leg_tab, hap_tab, 'MODELS/MODELS18D.p' if max_reads>16 else 'MODELS/MODELS16D.p')
    LLR_dict = collections.defaultdict(list)

    for k in range(subsamples):
        sys.stdout.write('\r')
        sys.stdout.write(f"[{'=' * int(k+1):{subsamples}s}] {int(100*(k+1)/subsamples)}% ")
        sys.stdout.flush()
        
        blocks_dict_picked = {block: pick_reads(reads_dict,weight_dict,read_IDs,min_reads,max_reads)
                           for block,read_IDs in iter_blocks(obs_tab,leg_tab,block_size,offset)}

        for block,haplotypes in blocks_dict_picked.items():
            LLR_dict[block].append(LLR(*haplotypes) if haplotypes!=None else None)
 
    LLR_stat = {block: mean_and_var(LLRs) for block,LLRs in LLR_dict.items() if None not in LLRs}
    
    M, V = zip(*LLR_stat.values())
    mean, std = sum(M)/len(M), sum(V)**.5/len(V) #The standard deviation is calculated according to the Bienaymé formula.

    num_of_LD_blocks = len(M)
    fraction_of_negative_LLRs = sum([1 for i in M if i<0])/len(M)
    
    print('\nFilename: %s' % obs_filename)
    print('Depth: %.2f, Number of LD blocks: %d, Fraction of LD blocks with a negative LLR: %.3f' % (info['depth'], num_of_LD_blocks,fraction_of_negative_LLRs))
    print('Mean: %.3f, Standard error: %.3f' % ( mean, std))
    
    info.update({'block_size': block_size,
                 'offset': offset,
                 'min_reads': min_reads,
                 'max_reads': max_reads,
                 'runtime': time.time()-a})

    info['statistics'] = {'mean': mean, 'std': std,
                          'num_of_LD_blocks': num_of_LD_blocks,
                          'fraction_of_negative_LLRs': fraction_of_negative_LLRs}

    if output_filename!=None:
        default_filename = re.sub('(.*)obs','\\1LLR', obs_filename.split('/')[-1],1)
        output_filename = default_filename if output_filename=='' else output_filename
        with open( output_filename, "wb") as f:
            pickle.dump(LLR_dict, f, protocol=4)
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
    parser.add_argument('obs_filename', metavar='OBS_FILENAME', type=str,
                        help='A pickle file created by MAKE_OBS_TAB, containing base observations at known SNP positions.')
    parser.add_argument('leg_filename', metavar='LEG_FILENAME', type=str,
                        help='IMPUTE2 legend file')
    parser.add_argument('hap_filename', metavar='HAP_FILENAME', type=str,
                        help='IMPUTE2 haplotype file')
    parser.add_argument('-b', '--block-size', type=int,
                        metavar='INT', default='100000',
                        help='Specifies the typical size of the LD block. The default value is 10^5.')
    parser.add_argument('-s', '--subsamples', type=int,
                        metavar='INT', default='32',
                        help='Sets the number of subsamples per LD block. The default value is 32.')
    parser.add_argument('-o', '--offset', type=int,
                        metavar='INT', default=0,
                        help='Shifts all the LD blocks by the requested base pairs. The default value is 0.')
    parser.add_argument('-m', '--min-reads', type=int, metavar='INT', default=2,
                        help='Takes into account only LD blocks with at least INT reads, admitting non-zero weight. The default value is 2.')
    parser.add_argument('-M', '--max-reads', type=int, metavar='INT', default=16,
                        help='Selects up to INT reads from each LD blocks. The default value is 16.')
    parser.add_argument('-O', '--output-filename', type=str, metavar='output_filename',  default='',
                        help='The output filename. The default is the input filename with the extension \".obs.p\" replaced by \".LLR.p\".')
    args = parser.parse_args()


    LLR_dict, info = aneuploidy_test(**vars(args))

    sys.exit(0)
else:
    print("The module ANEUPLOIDY_TEST was imported.")

### END OF FILE ###



"""
if __name__ == "__main__":
    print("Executed when invoked directly")
    a = time.time()
    random.seed(a=None, version=2) #I should set a=None after finishing to debug the code.

    args = dict(obs_filename = 'results_EUR/mixed2haploids.X0.5.SRR10393062.SRR151495.0-2.hg38.obs.p',
                hap_filename = '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.hap',
                leg_filename = '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.legend',
                block_size = 1e5,
                subsamples = 32,
                offset = 0,
                min_reads = 2,
                max_reads = 16,
                output_filename = None)

        with open(obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        info = pickle.load(f)
        
    hap_tab = read_impute2(hap_filename, filetype='hap')
    leg_tab = read_impute2(leg_filename, filetype='leg')

    reads_dict = build_reads_dict(obs_tab,leg_tab)
    weight_dict = build_weight_dict(reads_dict,obs_tab,leg_tab,hap_tab)
    LLR = get_LLR(obs_tab, leg_tab, hap_tab, 'MODELS/MODELS18D.p' if max_reads>16 else 'MODELS/MODELS16D.p')
    LLR_dict = collections.defaultdict(list)

    for k in range(subsamples):
        sys.stdout.write('\r')
        sys.stdout.write(f"[{'=' * int(k+1):{subsamples}s}] {int(100*(k+1)/subsamples)}% ")
        sys.stdout.flush()
        
        blocks_dict_picked = {block: pick_reads(reads_dict,weight_dict,read_IDs,min_reads,max_reads)
                           for block,read_IDs in iter_blocks(obs_tab,leg_tab,block_size,offset)}

        for block,haplotypes in blocks_dict_picked.items():
            LLR_dict[block].append(LLR(*haplotypes) if haplotypes!=None else None)
 
    LLR_stat = {block: mean_and_std(LLRs) for block,LLRs in LLR_dict.items() if None not in LLRs}
    
    M, V = zip(*LLR_stat.values())
    mean, std = sum(M)/len(M), sum(V)**.5/len(V) #The standard deviation is calculated according to the Bienaymé formula.

    LLR_jack = {block: jackknife(LLRs,[1/len(LLRs)]*len(LLRs)) for block,LLRs in LLR_dict.items() if None not in LLRs}
    M_jack, V_jack = zip(*LLR_jack.values())
    jk_mean, jk_std = sum(M_jack)/len(M_jack), sum(V_jack)**.5/len(V_jack)


    num_of_LD_blocks = len(M_jack)
    fraction_of_negative_LLRs = sum([1 for i  in M_jack if i<0])/len(M_jack)
    
    print('\nFilename: %s' % obs_filename)
    print('Depth: %.2f, Number of LD blocks: %d, Fraction of LD blocks with a negative LLR: %.3f' % (info['depth'], num_of_LD_blocks,fraction_of_negative_LLRs))
    print('Mean: %.3f, Standard error: %.3f, Jackknife standard error: %.3f' % ( mean, std, jk_std))
    
    info.update({'block_size': block_size,
                 'offset': offset,
                 'min_reads': min_reads,
                 'max_reads': max_reads,
                 'runtime': time.time()-a})

    info['statistics'] = {'mean': mean, 'std': std, 'jk_std': jk_std, 
                          'num_of_LD_blocks': num_of_LD_blocks,
                          'fraction_of_negative_LLRs': fraction_of_negative_LLRs}

    if output_filename!=None:
        default_filename = re.sub('(.*)obs','\\1LLR', obs_filename.split('/')[-1],1)
        output_filename = default_filename if output_filename=='' else output_filename
        with open( output_filename, "wb") as f:
            pickle.dump(LLR_dict, f, protocol=4)
            pickle.dump(info, f, protocol=4)

    b = time.time()
    print('Done calculating LLRs for all the LD block in %.3f sec.' % ((b-a)))
"""
