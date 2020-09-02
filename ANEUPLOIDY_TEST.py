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
    """ Calculates the mean and the sample standard deviation. """
    mean = statistics.mean(sample)
    var = statistics.variance(sample, xbar=mean)
    return mean, var

def build_reads_dict(obs_tab,leg_tab):
    """ Returns a dictionary that lists read IDs of reads that overlap with
        SNPs and gives the alleles in each read. """

    reads = collections.defaultdict(list)

    for (pos, ind, read_id, base) in obs_tab:
        if base in leg_tab[ind][2:]:
            reads[read_id].append((pos,base))

    return reads

def build_rank_dict(reads_dict,obs_tab,leg_tab, hap_tab):
    """ Returns a dicitonary lists read_IDs and gives their rank. The rank of a
        read is inherited from the rank of the SNPs that are overlapped by the
        read. Thus, the rank of a read depends on its starting position and
        length, but not on the alleles that it contains. """

    def bools2int(x):
        """ Transforms a tuple/list of bools to a int. """
        return int(''.join(str(int(i)) for i in x),2)

    def rank(X,N):
        """ Returns the rank of a given haplotypes. """
        threshold = 0.05
        joint_freq = bin(functools.reduce(operator.and_,X)).count('1') / N
        result = threshold<=joint_freq<=(1-threshold)
        return result

    N = len(hap_tab[0])

    hap_dict = dict()
    for (pos, ind, read_id, base) in obs_tab:
        ref, alt = leg_tab[ind][2:]
        if (base==alt or base==ref) and (0.01<hap_tab[ind].count(1)/N<0.99): # if a biallelic SNP admits an allele with a frequency close to zero then this SNP would not contribute significantly to the rank of the read. However, including these SNPs in the calculation of the rank would prolong the calculation time. Thus, this function determines if the SNP should be included the calculation of the rank. 
            hap_dict[pos] = (bools2int(hap_tab[ind]), bools2int(map(operator.not_,hap_tab[ind])))
 
    rank_dict = dict()
    for read_id in reads_dict:
        hap = (hap_dict[pos] for pos,base in reads_dict[read_id] if pos in hap_dict)
        rank_dict[read_id] = sum(rank(C,N) for C in itertools.product(*hap) if len(C)!=0)

    return rank_dict

def pick_reads(reads_dict,rank_dict,read_IDs,min_reads,max_reads):
    """ Draws up to max_reads reads from a given LD block; Only reads with a
        rank larger than one are sampled. In addition, if the number of reads
        in a given LD block is less than the minimal requirment then the block
        would not be considered."""

    read_IDs_filtered = tuple(read_ID for read_ID in read_IDs if rank_dict[read_ID]>1)
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
    rank_dict = build_rank_dict(reads_dict,obs_tab,leg_tab,hap_tab)
    blocks_dict = {block: read_IDs for block,read_IDs in iter_blocks(obs_tab,leg_tab,block_size,offset)}
    LLR = get_LLR(obs_tab, leg_tab, hap_tab, 'MODELS/MODELS18D.p' if max_reads>16 else 'MODELS/MODELS16D.p')
    LLR_dict = collections.defaultdict(list)

    for k in range(subsamples):
        sys.stdout.write('\r')
        sys.stdout.write(f"[{'=' * int(k+1):{subsamples}s}] {int(100*(k+1)/subsamples)}% ")
        sys.stdout.flush()
        
        blocks_dict_picked = {block: pick_reads(reads_dict,rank_dict,read_IDs,min_reads,max_reads)
                           for block,read_IDs in blocks_dict.items()}

        for block,haplotypes in blocks_dict_picked.items():
            LLR_dict[block].append(LLR(*haplotypes) if haplotypes!=None else None)
 
    LLR_stat = {block: mean_and_var(LLRs) for block,LLRs in LLR_dict.items() if None not in LLRs}
    
    M, V = zip(*LLR_stat.values())
    mean, std = sum(M)/len(M), sum(V)**.5/len(V) #The standard deviation is calculated according to the Bienaymé formula.

    num_of_LD_blocks = len(M)
    fraction_of_negative_LLRs = sum([1 for i in M if i<0])/len(M)
    
    
    reads_per_LDblock = [sum(1 for read_ID in read_IDs if rank_dict[read_ID]>1) for read_IDs in blocks_dict.values()]
    reads_per_LDblock_filtered = [i for i in reads_per_LDblock if i>1]
    reads_mean = statistics.mean(reads_per_LDblock_filtered)
    reads_var = statistics.pstdev(reads_per_LDblock_filtered, mu=reads_mean)
    
    print('\nFilename: %s' % obs_filename)
    print('Depth: %.2f, Number of LD blocks: %d, Fraction of LD blocks with a negative LLR: %.3f' % (info['depth'], num_of_LD_blocks,fraction_of_negative_LLRs))
    print('Mean and standard error of meaningful reads per LD block: %.1f, %.1f.' % (reads_mean, reads_var))
    print('Mean LLR: %.3f, Standard error of the mean LLR: %.3f' % ( mean, std))
    
    info.update({'block_size': block_size,
                 'offset': offset,
                 'min_reads': min_reads,
                 'max_reads': max_reads,
                 'runtime': time.time()-a})

    info['statistics'] = {'mean': mean, 'std': std,
                          'num_of_LD_blocks': num_of_LD_blocks,
                          'fraction_of_negative_LLRs': fraction_of_negative_LLRs,
                          'reads_mean': reads_mean, 'reads_var': reads_var}

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
                        help='Takes into account only LD blocks with at least INT reads, admitting non-zero rank. The default value is 2.')
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

    args = dict(obs_filename = 'results_EUR/mixed3haploids.X0.10.HG00096A.HG00096B.HG00097A.chr21.recomb.0.00.obs.p',
                hap_filename = '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.hap',
                leg_filename = '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.legend',
                block_size = 1e5,
                subsamples = 32,
                offset = 0,
                min_reads = 2,
                max_reads = 16,
                output_filename = None)

    with open(args['obs_filename'], 'rb') as f:
        obs_tab = pickle.load(f)
        info = pickle.load(f)
        
    hap_tab = read_impute2(args['hap_filename'], filetype='hap')
    leg_tab = read_impute2(args['leg_filename'], filetype='leg')

    reads_dict = build_reads_dict(obs_tab,leg_tab)
    rank_dict = build_rank_dict(reads_dict,obs_tab,leg_tab,hap_tab)
    blocks_dict = {block: read_IDs for block,read_IDs in iter_blocks(obs_tab,leg_tab,args['block_size'],args['offset'])}
    
    LLR = get_LLR(obs_tab, leg_tab, hap_tab, 'MODELS/MODELS18D.p' if args['max_reads']>16 else 'MODELS/MODELS16D.p')
    LLR_dict = collections.defaultdict(list)

    for k in range(args['subsamples']):
        sys.stdout.write('\r')
        sys.stdout.write(f"[{'=' * int(k+1):{args['subsamples']}s}] {int(100*(k+1)/args['subsamples'])}% ")
        sys.stdout.flush()
        
        blocks_dict_picked = {block: pick_reads(reads_dict,rank_dict,read_IDs,args['min_reads'],args['max_reads'])
                           for block,read_IDs in blocks_dict.items()}

        for block,haplotypes in blocks_dict_picked.items():
            LLR_dict[block].append(LLR(*haplotypes) if haplotypes!=None else None)
 
    LLR_stat = {block: mean_and_var(LLRs) for block,LLRs in LLR_dict.items() if None not in LLRs}
    
    M, V = zip(*LLR_stat.values())
    mean, std = sum(M)/len(M), sum(V)**.5/len(V) #The standard deviation is calculated according to the Bienaymé formula.

    num_of_LD_blocks = len(M)
    fraction_of_negative_LLRs = sum([1 for i in M if i<0])/len(M)
    
    
    reads_per_LDblock = [sum(1 for read_ID in read_IDs if rank_dict[read_ID]>1) for read_IDs in blocks_dict.values()]
    reads_per_LDblock_filtered = [i for i in reads_per_LDblock if i>1]
    reads_mean = statistics.mean(reads_per_LDblock_filtered)
    reads_var = statistics.pstdev(reads_per_LDblock_filtered, mu=reads_mean)
    
    print('\nFilename: %s' % args['obs_filename'])
    print('Depth: %.2f, Number of LD blocks: %d, Fraction of LD blocks with a negative LLR: %.3f' % (info['depth'], num_of_LD_blocks,fraction_of_negative_LLRs))
    print('Mean and standard error of meaningful reads per LD block: %.1f, %.1f.' % (reads_mean, reads_var))
    print('Mean LLR: %.3f, Standard error of the mean LLR: %.3f' % ( mean, std))
    
    info.update({'block_size': args['block_size'],
                 'offset': args['offset'],
                 'min_reads': args['min_reads'],
                 'max_reads': args['max_reads'],
                 'runtime': time.time()-a})

    info['statistics'] = {'mean': mean, 'std': std,
                          'num_of_LD_blocks': num_of_LD_blocks,
                          'fraction_of_negative_LLRs': fraction_of_negative_LLRs,
                          'reads_mean': reads_mean, 'reads_var': reads_var}

    if args['output_filename']!=None:
        default_filename = re.sub('(.*)obs','\\1LLR', args['obs_filename'].split('/')[-1],1)
        output_filename = default_filename if args['output_filename']=='' else args['output_filename']
        with open( output_filename, "wb") as f:
            pickle.dump(LLR_dict, f, protocol=4)
            pickle.dump(info, f, protocol=4)
    
    
    ##filename = 'results_EUR/mixed3haploids.X0.10.HG00096A.HG00096B.HG00097A.chr21.recomb.%.2f.LLR_100x.p' % (i * 0.1)
    ##with open(filename, 'rb') as f:
    ##    LLR_dict = pickle.load(f)
    ##    info = pickle.load(f)
    ##info['statistics']['reads_mean']=reads_mean
    ##info['statistics']['reads_var']=reads_var
    ##with open( filename, "wb") as f:
    ##    pickle.dump(LLR_dict, f, protocol=4)
    ##    pickle.dump(info, f, protocol=4)
    ##print(i)
    
    b = time.time()
    print('Done calculating LLRs for all the LD block in %.3f sec.' % ((b-a)))
"""
