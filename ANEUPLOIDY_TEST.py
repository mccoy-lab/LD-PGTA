#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANEUPLOIDY_TEST

Builds a dictionary that lists genomic windows that contain at least three
reads and gives the likelihoods to observese these reads under various
aneuploidy landscapes --- i.e., monosomy, disomy, SPH and BPH.

BPH (Both Parental Homologs) correspond to the presence of three unmatched
haplotypes, while SPH (Single Parental Homolog) correspond to chromosome gains
involving identical homologs.
Daniel Ariad (daniel@ariad.org)
Sep 1, 2020
"""

import collections, time, pickle, argparse, re, sys, random, os, bz2, gzip, platform
from HOMOGENOUES_MODELS import homogeneous
from RECENT_ADMIXTURE_MODELS import recent_admixture
from DISTANT_ADMIXTURE_MODELS import distant_admixture


from itertools import product
from functools import reduce
from operator import and_, attrgetter
from statistics import mean, variance, pstdev
from math import log

leg_tuple = collections.namedtuple('leg_tuple', ('chr_id', 'pos', 'ref', 'alt')) #Encodes the rows of the legend table
sam_tuple = collections.namedtuple('sam_tuple', ('sample_id', 'group1', 'group2', 'sex')) #Encodes the rows of the samples table
obs_tuple = collections.namedtuple('obs_tuple', ('pos', 'read_id', 'base')) #Encodes the rows of the observations table
comb_tuple = collections.namedtuple('comb_tuple', ('ref','alt','hap'))
likelihoods_tuple = collections.namedtuple('likelihoods_tuple', ('monosomy', 'disomy', 'SPH', 'BPH'))

### Getting a function to count non-zero bits in positive integer.
try:
    if platform.python_implementation()=='PyPy':
        from POPCOUNTS import get_popcount
        popcount = get_popcount(64)
    else:
        from gmpy2 import popcount
except Exception as err: 
    print(err)
    popcount = lambda x: bin(x).count('1')

try:
    from math import comb
except ImportError:
    print('caution: cound not import comb from the math module.')
    def comb(n, k):
        """ Return the number of ways to choose k items from n items without repetition and without order. """
        if not 0 <= k <= n:
            return 0
        b = 1
        for t in range(min(k, n-k)):
            b *= n
            b //= t+1
            n -= 1
        return b
    
def mean_and_var(x):
    """ Calculates the mean and variance. """
    cache = tuple(x)
    m = mean(cache)
    var = variance(cache, xbar=m)
    return m, var

def mean_and_std(x):
    """ Calculates the mean and population standard deviation. """
    cache = tuple(x)
    m = mean(cache)
    std = pstdev(cache, mu=m)
    return m, std

def mean_and_std_of_mean_of_rnd_var(A):
    """ Calculates the mean and population standard deviation of the mean of random variables.
        Each row of A represents a random variable, with observations in the columns."""
    if type(A)==dict:
        A = tuple(tuple(i) for i in A.values()) 
    
    M, N = len(A), len(A[0])
    mu = sum(sum(likelihoods_in_window)/N for likelihoods_in_window in A)
    arg = ((sum(sampled_likelihoods) - mu)**2 for sampled_likelihoods in zip(*A))
    std = (sum(arg) / (N - 1))**.5 / M
    mean = mu / M
    return mean, std

### DEPRECATED. ###
#def summarize(M,V):
#    """ Calculates chromosome-wide statistics of the LLRs """
#    result =  {'mean': mean(M),
#               'std_of_mean': sum(V)**.5/len(V),  #The standard deviation is calculated according to the Bienaymé formula.
#               'fraction_of_negative_LLRs': [i<0 for i in M].count(1)/len(M)}
#    return result

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

def invert(x,n):
    """ Inverts the bits of a positive integer. """
    return  x ^ ((1 << n) - 1)

class supporting_dictionaries:
    """ This class provides dictionaries that were created by various 
        intersections and unions of the observations table with the reference panel. """
    
    def __init__(self, obs_tab, leg_tab, hap_tab, number_of_haplotypes, min_HF):
        
        self.number_of_haplotypes = number_of_haplotypes
        self.min_HF = min_HF
    
        self.combined = self.build_combined_dict(obs_tab, leg_tab, hap_tab)
        self.reads = self.build_reads_dict(obs_tab, self.combined)
        self.score = self.build_score_dict(self.reads, self.combined ,number_of_haplotypes, min_HF)
        self.overlaps = self.build_overlaps_dict(obs_tab, self.combined)
       
        self.positions = (*self.overlaps.keys(),) ### The intersections of positions within the observation table and the reference panel.
        
    def build_combined_dict(self, obs_tab, leg_tab, hap_tab):
        """  Returns a dictionary that lists chromosomal positions of SNPs and
        gives their associated reference alleles, alternative alleles and reference
        panel. """
        cache = {leg.pos: comb_tuple(leg.ref,leg.alt,hap) 
                        for leg,hap in zip(leg_tab, hap_tab)}
        combined = {obs.pos: cache.pop(obs.pos) for obs in obs_tab 
                        if obs.pos in cache}
        return combined
    
    def build_reads_dict(self, obs_tab, combined):
        """ Returns a dictionary that lists read IDs of reads that overlap with
            SNPs and gives the alleles in each read. """
    
        reads = collections.defaultdict(list)
    
        for pos, read_id, base in obs_tab:
            if pos in combined and (base==combined[pos].ref or base==combined[pos].alt):
                reads[read_id].append((pos,base))
    
        return reads

    def build_overlaps_dict(self, obs_tab, combined):
        """ Returns a dictionary that lists chromosome positions of SNPs and gives a
        list of read IDs for all the reads that overlap with the SNP. """
        overlaps_dict = collections.defaultdict(list)
        for pos, read_id, base in obs_tab:
            if pos in combined and (base==combined[pos].ref or base==combined[pos].alt):
                overlaps_dict[pos].append(read_id)
        return overlaps_dict
    
    def build_score_dict(self, reads_dict, combined, number_of_haplotypes, min_HF):
        """ Returns a dicitonary lists read_IDs and gives their score. The scoring
        algorithm scores each read according to the number of differet haplotypes
        that the reference panel supports at the chromosomal region that overlaps
        with the read. Only bialleic SNP with a minor allele frequeancy above
        0.01 are considered for the calculation, since they are unlikely to affect
        the score. In addition, only haplotypes with a frequnecy between min_HF
        and 1-min_HF add to the score of a read. """
    
        N = number_of_haplotypes
        b = (1 << number_of_haplotypes) - 1 #### equivalent to int('1'*number_of_haplotypes,2)
    
        score_dict = dict()
        for read_id in reads_dict:
            haplotypes = ((combined[pos].hap, combined[pos].hap ^ b)
                              for pos,base in reads_dict[read_id]
                                  if 0.01 <= popcount(combined[pos].hap)/N <= 0.99)  #Include only biallelic SNPs with MAF of at least 0.01. Also, ^b flips all bits of the binary number, hap_tab[ind] using bitwise xor operator.
    
            score_dict[read_id] = sum(min_HF <= popcount(reduce(and_,hap))/N <= (1-min_HF)
                                      for hap in product(*haplotypes) if len(hap)!=0)
    
        return score_dict
    
def build_gw_dict(obs, window_size, offset, min_reads, max_reads, min_score):
    """ Returns a dictionary the lists genomic windows and gives the read IDs
        of reads that overlap with SNPs in the genomic window. Only reads with
        a score larger than one are considered. """

    max_dist = 100000 #maximal distance between consecutive observed alleles.
    max_win_size = 350000 #maximal genomic window size
    initial_win_size = 10000 #initial genomic window size

    adaptive, window_size = (False, int(window_size)) if window_size else (True, initial_win_size)
    offset = int(offset)

    first, last = obs.positions[0] + offset, obs.positions[-1] + window_size
    a, b  = first, first+window_size
    readIDs_in_window = set()
    
    gw_dict = {}
    
    for pos, overlapping_reads in obs.overlaps.items():
        if pos<first: continue
        while b<last:
            if a<=pos<b:
                readIDs_in_window.update(read_ID for read_ID in overlapping_reads if min_score<=obs.score[read_ID])
                break
            elif adaptive and 0<len(readIDs_in_window)<min_reads and b-pos<=max_dist and b-a<=max_win_size:
                b += 10000
            else:
                gw_dict[a,b-1] = tuple(readIDs_in_window) #the genomic window includes both endpoints.
                a, b = b, b+window_size
                readIDs_in_window = set()
    
    return gw_dict

def pick_reads(reads_dict,read_IDs,max_reads):
    """ Draws up to max_reads reads from a given genomic window. The reads are
        randomly sampled without replacement to meet the assumptions of the
        statistical models."""

    drawn_read_IDs = random.sample(read_IDs, min(len(read_IDs)-1,max_reads))
    haplotypes = tuple(reads_dict[read_ID] for read_ID in drawn_read_IDs)
    return haplotypes

def effective_number_of_subsamples(num_of_reads,min_reads,max_reads,subsamples):
    """ Ensures that the number of requested subsamples is not larger than the
    number of unique subsamples. In addition, it checks that a genomic window
    contains a minimal number of reads (min_reads) that overlap with known 
    SNPs. """

    if  min_reads <= num_of_reads > max_reads:
        eff_subsamples = min(comb(num_of_reads,max_reads),subsamples)
    elif min_reads <= num_of_reads <= max_reads:
        eff_subsamples = min(num_of_reads,subsamples)
    else:
        eff_subsamples = 0

    return eff_subsamples

def bootstrap(obs_tab, leg_tab, hap_tab, sam_tab, number_of_haplotypes,
              ancestral_makeup, models_dict, window_size, subsamples, offset, 
              min_reads, max_reads, min_score, min_HF):
    """ Applies a bootstrap approach in which: (i) the resample size is smaller
    than the sample size and (ii) resampling is done without replacement. """

    min_reads == max(min_reads,3) # Due to the bootstrap approach, min_reads must be at least 3.
    max_reads == max(max_reads,2) # Our statistical models require at least 2 reads.


    obs = supporting_dictionaries(obs_tab, leg_tab, hap_tab, number_of_haplotypes, min_HF)
    
    genomic_windows = build_gw_dict(obs, window_size, offset, min_reads, max_reads, min_score)
    
    if type(ancestral_makeup)==set and len(ancestral_makeup)==1:
        print('Assuming one ancestral population: %s.' % next(iter(ancestral_makeup)))
        examine = homogeneous(obs_tab, leg_tab, hap_tab, sam_tab, models_dict, number_of_haplotypes)
    elif type(ancestral_makeup)==set and len(ancestral_makeup)==2:
        print('Assuming recent-admixture between %s and %s.' % tuple(ancestral_makeup))
        examine = recent_admixture(obs_tab, leg_tab, hap_tab, sam_tab, models_dict, number_of_haplotypes)
    elif type(ancestral_makeup)==dict:
        print('Assuming the following ancestral makeup:', ", ".join("{:.1f}% {}".format(100*v, k) for k, v in ancestral_makeup.items()))       
        examine = distant_admixture(obs_tab, leg_tab, hap_tab, sam_tab, models_dict, number_of_haplotypes, ancestral_makeup)
    else:
        raise Exception('error: unsupported ancestral-makeup.')

    likelihoods = {}
                    
    for k,(window,read_IDs) in enumerate(genomic_windows.items()):

        effN = effective_number_of_subsamples(len(read_IDs),min_reads,max_reads,subsamples)
        if effN: ### Adds a genomic window only if it contains enough reads for sampling.
            cache = []
            for i in range(effN):
                sys.stdout.write(f"\r{chr(10633)+' ':s}{chr(10495) * (32*(k+1)//len(genomic_windows)) + chr(10240+(256-effN+i)%256):{33}s}{' '+chr(10634):s} {int(100*(k+1)/len(genomic_windows))}%"); sys.stdout.flush()
                sampled_reads = pick_reads(obs.reads,read_IDs,max_reads)
                cache.append(examine.get_likelihoods(*sampled_reads))
            likelihoods[window] = tuple(cache)
    return likelihoods, genomic_windows, examine.fraction_of_matches

def statistics(likelihoods,genomic_windows):
    """ Compares likelihoods of different aneuploidy scenarios and extracts
    useful information about the genmoic windows. """

    if likelihoods:
        window_size_mean, window_size_std = mean_and_std((j-i+1 for (i,j) in likelihoods))
        reads_mean, reads_std = mean_and_std((len(read_IDs) for window,read_IDs in genomic_windows.items() if window in likelihoods))
        num_of_windows = len(likelihoods)

        pairs = (('BPH','SPH'), ('disomy','monosomy'), ('BPH','disomy'), ('disomy','SPH'), ('SPH','monosomy'));
        
        LLRs = {(i,j): 
                {window: tuple(LLR(attrgetter(i)(l), attrgetter(j)(l)) for l in likelihoods_in_window)
                           for window,likelihoods_in_window in likelihoods.items()} 
                                                            for i,j in pairs}
            
        LLRs_per_genomic_window = {pair: 
                {window: mean_and_var(LLRs_in_window) for window, LLRs_in_window in LLRs[pair].items()}
                                                            for pair in pairs}
            
        cache = {pair: mean_and_std_of_mean_of_rnd_var(LLRs[pair]) for pair in pairs}
        
        LLRs_per_chromosome = {pair: {'mean_of_mean': cache[pair][0], 
                                      'std_of_mean': cache[pair][1], 
                                      'fraction_of_negative_LLRs': mean(m<0 for m,v in LLRs_per_genomic_window[pair].values())}
                                           for pair in pairs}
        
        #LLRs_per_chromosome = {pair: summarize(*zip(*LLRs_per_genomic_window[pair].values()))
        #                              for pair in pairs}
        
        result = {'num_of_windows': num_of_windows,
                  'reads_mean': reads_mean,
                  'reads_std': reads_std,
                  'window_size_mean': window_size_mean,
                  'window_size_std': window_size_std,
                  'LLRs_per_genomic_window': LLRs_per_genomic_window,
                  'LLRs_per_chromosome': LLRs_per_chromosome}
    else:
        result = {'num_of_windows': 0}
    return result

def print_summary(obs_filename,info):
    S = info['statistics']
    ancestral_makeup = ", ".join("{:.1f}% {}".format(100*v, k) for k, v in info['ancestral_makeup'].items()) if type(info['ancestral_makeup'])==dict else ', '.join(info['ancestral_makeup'])
    print('\nFilename: %s' % obs_filename)
    print('\nSummary statistics')
    print('------------------')    
    print('Chromosome ID: %s, Depth: %.2f.' % (info['chr_id'],info['depth']))
    print('Number of genomic windows: %d, Mean and standard error of genomic window size: %d, %d.' % (S.get('num_of_windows',0),S.get('window_size_mean',0),S.get('window_size_std',0)))
    print('Mean and standard error of meaningful reads per genomic window: %.1f, %.1f.' % (S.get('reads_mean',0), S.get('reads_std',0)))
    print('Ancestral makeup: %s, Fraction of alleles matched to the reference panel: %.3f.' % (ancestral_makeup, info['statistics']['matched_alleles']))

    if S.get('LLRs_per_chromosome',None):
        for (i,j), L in S['LLRs_per_chromosome'].items():
            print(f"--- Chromosome-wide LLR between {i:s} and {j:s} ----")
            print(f"Mean LLR: {L['mean_of_mean']:.3f}, Standard error of the mean LLR: {L['std_of_mean']:.3f}")
            print(f"Fraction of genomic windows with a negative LLR: {L['fraction_of_negative_LLRs']:.3f}")

def save_results(likelihoods,info,compress,obs_filename,output_filename,output_dir):
    """ Saves the likelihoods together with information about the chromosome
        number, depth of coverage, ancestry, statistics of the genomic windows
        and flags that were used. Also, data compression is supported in gzip
        and bzip2 formats. """

    Open = {'bz2': bz2.open, 'gz': gzip.open}.get(compress, open)
    ext = ('.'+compress) * (compress in ('bz2','gz'))
    obs_filename_stripped =  obs_filename.rsplit('/', 1).pop()
    default_output_filename = re.sub('(.*)obs.p(.*)',f'\\1LLR.p{ext:s}', obs_filename_stripped, 1)
    if output_filename=='': 
        output_filename = default_output_filename
    else:
        output_filename += ext
    
    output_dir = re.sub('/$','',output_dir)+'/' #undocumented option
    if output_dir!='' and not os.path.exists(output_dir): os.makedirs(output_dir)
    
    with Open(output_dir + output_filename, "wb") as f:
        pickle.dump(likelihoods, f, protocol=4)
        pickle.dump(info, f, protocol=4)
    return output_dir + output_filename

def aneuploidy_test(obs_filename,leg_filename,hap_filename,samp_filename,
                    ancestral_makeup,window_size,subsamples,offset,min_reads,
                    max_reads,min_score,min_HF,output_filename,compress,**kwargs):
    """ Returns a dictionary that lists the boundaries of approximately
    independent genomic windows. For each genomic window it gives the
    likelihood of four scenarios, namely, monosomy, disomy, SPH and BPH.
    It also returns a dictionary with various run parameters. """

    time0 = time.time()

    random.seed(a=kwargs.get('seed', None), version=2) #I should make sure that a=None after finishing to debug the code.
    path = os.path.realpath(__file__).rsplit('/', 1)[0] + '/MODELS/'
    models_filename = kwargs.get('model', path + ('MODELS18.p' if max_reads>16 else ('MODELS16.p' if max_reads>12 else 'MODELS12.p')))
        
    load = lambda filename: {'bz2': bz2.open, 'gz': gzip.open}.get(filename.rsplit('.',1)[1], open)  #Adjusts the opening method according to the file extension.

    open_hap = load(hap_filename)
    with open_hap(hap_filename,'rb') as hap_in:
        hap_tab, number_of_haplotypes = pickle.load(hap_in)

    open_leg = load(leg_filename)
    with open_leg(leg_filename,'rb') as leg_in:
        leg_tab = pickle.load(leg_in)


    open_samp = load(samp_filename)
    with open_samp(samp_filename,'rb') as samp_in:
        sam_tab = pickle.load(samp_in)

    open_obs = load(obs_filename)
    with open_obs(obs_filename, 'rb') as obs_in:
        obs_tab = pickle.load(obs_in)
        info = pickle.load(obs_in)

    open_model = load(models_filename)
    with open_model(models_filename, 'rb') as model_in:
        models_dict = pickle.load(model_in)
        
    ancestry = set(ancestral_makeup.keys()) if type(ancestral_makeup)==dict else set(ancestral_makeup)
    assert ancestry=={row.group2 for row in sam_tab}, 'error: the given ancestral makup does match the populations (group2) in the reference panel.'

    likelihoods, genomic_windows, matched_alleles = bootstrap(obs_tab, leg_tab, hap_tab, sam_tab, number_of_haplotypes, ancestral_makeup, models_dict, window_size, subsamples, offset, min_reads, max_reads, min_score, min_HF)

    some_statistics = {'matched_alleles': matched_alleles,
                       'runtime': time.time()-time0}


    info.update({'ancestral_makeup': ancestral_makeup,
                 'window_size': window_size,
                 'subsamples': subsamples,
                 'offset': offset,
                 'min_reads': min_reads,
                 'max_reads': max_reads,
                 'min_score': min_score,
                 'min_HF': min_HF,
                 'statistics': {**statistics(likelihoods,genomic_windows), **some_statistics}
                 })
    
    if output_filename!=None:
        save_results(likelihoods,info,compress,obs_filename,output_filename,kwargs.get('output_dir', 'results'))

    print_summary(obs_filename,info)

    time1 = time.time()
    print('Done calculating LLRs for all the genomic windows in %.3f sec.' % ((time1-time0)))

    return likelihoods, info

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
    description='Builds a dictionary that lists genomic windows that contain'
                'at least two reads and gives the associated log-likelihood '
                'BPH/SPH ratio (LLR). BPH (Both Parental Homologs) correspond'
                'to the presence of three unmatched haplotypes, while SPH'
                '(Single Parental Homolog) correspond to chromosome gains'
                'involving identical homologs.')
    parser.add_argument('obs_filename', metavar='OBS_FILENAME', type=str,
                        help='A observations file created by MAKE_OBS_TAB.')
    parser.add_argument('leg_filename', metavar='LEG_FILENAME', type=str,
                        help='A legend file of the reference panel.')
    parser.add_argument('hap_filename', metavar='HAP_FILENAME', type=str,
                        help='A haplotype file of the reference panel.')
    parser.add_argument('samp_filename', metavar='SAMP_FILENAME', type=str,
                        help='A samples file of the reference panel.')
    parser.add_argument('ancestral_makeup', metavar='ANCESTRAL_MAKEUP', nargs='+',
                        help='Assume an ancestral makeup:\n'
                        'a. For non-admixtures the argument consists a single superpopulation, e.g., EUR. \n'
                        'b. For recent admixtures the argument consists two superpopulations, e.g., EUR EAS. \n'
                        'c. For distant admixtures the argument consists of the superpoplations and their proportions, e.g, EUR 0.8 EAS 0.1 SAS 0.1.\n')
    parser.add_argument('-w', '--window-size', type=int,
                        metavar='INT', default='100000',
                        help='Specifies the size of the genomic window. The default value is 100 kbp. When given a zero-size genomic window, it adjusts the size of the window to include min-reads reads.')
    parser.add_argument('-s', '--subsamples', type=int,
                        metavar='INT', default='32',
                        help='Sets the number of subsamples per genomic window. The default value is 32.')
    parser.add_argument('-o', '--offset', type=int, metavar='INT', default=0,
                        help='Shifts all the genomic windows by the requested base pairs. The default value is 0.')
    parser.add_argument('-m', '--min-reads', type=int, metavar='INT', default=6,
                        help='Takes into account only genomic windows with at least INT reads, admitting non-zero score. The minimal value is 3, while the default is 6.')
    parser.add_argument('-M', '--max-reads', type=int, metavar='INT', default=4,
                        help='Selects up to INT reads from each genomic windows in each bootstrap sampling. The minimal value is 2, while the default value is 4.')
    parser.add_argument('-F', '--min-HF', type=float, metavar='FLOAT', default=0.05,
                        help='Only haplotypes with a frequnecy between FLOAT and 1-FLOAT add to the score of a read. The default value is 0.05.')
    parser.add_argument('-S', '--min-score', type=int, metavar='INT', default=2,
                        help='Consider only reads that reach the minimal score. The default value is 2.')
    parser.add_argument('-O', '--output-filename', type=str, metavar='output_filename',  default='',
                        help='The output filename. The default is the input filename with the extension \".obs.p\" replaced by \".LLR.p\".')
    parser.add_argument('-C', '--compress', metavar='gz/bz2/unc', type=str, default='unc',  choices=['gz','bz2','unc'],
                        help='Output compressed via gzip, bzip2 or uncompressed. Default is uncompressed.')
    args = vars(parser.parse_args())

    strings_even = all(i.isalpha() for i in args['ancestral_makeup'][0::2])
    strings_odd = all(i.isalpha() for i in args['ancestral_makeup'][1::2])
    floats_odd = all(i.replace('.','',1).isdigit() for i in args['ancestral_makeup'][1::2])
    if strings_even and strings_odd and len(args['ancestral_makeup']) in (1,2):
        args['ancestral_makeup'] = set(args['ancestral_makeup'])
    elif strings_even and floats_odd:
        args['ancestral_makeup'] = dict(zip(args['ancestral_makeup'][0::2],map(float,args['ancestral_makeup'][1::2])))
    else:
        raise Exception('error: Invalid argument supplied for ancestral makeup.')

    LLR_dict, info = aneuploidy_test(**args)

    sys.exit(0)
else:
    print("The module ANEUPLOIDY_TEST was imported.")

### END OF FILE ###


"""
from ANEUPLOIDY_TEST import aneuploidy_test
def aneuploidy_test_demo(obs_filename='SWI-L-10-27-May-2020_S38.chr6.obs.p.bz2',chr_id='chr6',sp='EAS_EUR',output_path='results/',ref_panel_path='../build_reference_panel'):
    args = dict(obs_filename = output_path + obs_filename,
                hap_filename = f'{ref_panel_path:s}/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.hap.gz',
                leg_filename = f'{ref_panel_path:s}/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.legend.gz',
                sam_filename = f'{ref_panel_path:s}/samples_per_panel/{sp:s}_panel.samples',
                window_size = 0,
                subsamples = 100,
                offset = 0,
                min_reads = 18,
                max_reads = 12,
                min_HF = 0.05,
                min_score = 2,
                output_dir = output_path,
                output_filename = '',
                compress = 'bz2',
                seed=0)



    LLR_dict, info = aneuploidy_test(**args)
    return LLR_dict, info"
"""
