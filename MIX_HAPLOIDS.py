#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

Simulates observed bases at known SNP positions from an aneuploid cell, based
on mixtures of haploid sequences. 

The simulation supports three scenarios: (a) trisomy with two matched haplotypes
(‘single parental homolog’; SPH), (b) three unmatched haplotypes (‘both parental
homologs’; BPH), (c) trisomy with recombination that is charecterized by a
transition between SPH to BPH along the chromosome.   

MIX_HAPLOIDS

Daniel Ariad (daniel@ariad.org)
Aug 14st, 2020

"""
import pickle, time, random, operator, collections, warnings, argparse, sys
from random import choices, randrange

warnings.formatwarning = lambda message, category, filename, lineno, file=None, line=None: 'Caution: %s\n' % message

def compare_haploids(obs_filename1, obs_filename2):
    """ Given observation tables of twp haploid sequences with a depth coverage
        of at least 1X, the fraction of exclusive alleles is calculated. """
        
    time0 = time.time()
    obs_tab1 = pickle.load(open(obs_filename1, 'rb'))
    obs_tab2 = pickle.load(open(obs_filename2, 'rb'))
    alleles = collections.defaultdict(list)
    for pos, impute2_ind, reads_id, obs_base in obs_tab1:
        alleles[pos].append(obs_base)
    for pos, impute2_ind, reads_id, obs_base in obs_tab2:
        alleles[pos].append(obs_base)
    sample = tuple(x[0]!=x[1] for x in alleles.values() if len(x)>1)
    result = sum(sample)/len(sample)
    time1 = time.time()
    print('Done in %.2f sec.' % (time1-time0))
    return result

###############################################################################

def chr_length(chr_id):
    """ Return the chromosome length for a given chromosome, based on the reference genome hg38.""" 
    #The data of chromosome length was taken from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38
    length_dict = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
                  'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422,
                  'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
                  'chr16': 90338345, 'chr17':  83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
                  'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}
    return length_dict[chr_id]

def number_of_reads(chr_id,reads_length,depth):
    """ Calculates the the number of reads for a given coverage depth and reads length."""
    number_of_fragments = depth * chr_length(chr_id) // reads_length
    return int(number_of_fragments)
    
def build_dictionary_of_reads(obs_tab,chr_id,read_length,offset):
    """ Cuts the sequence of chromosome chr_id into simulated reads of length
        read_length. Then, creates a dictionary that matches a simulated read 
        with the observed alleles at SNPs that overlap with the read."""
    
    positions = next(zip(*obs_tab))    
    a = positions[0]-(read_length-1)+offset
    b = positions[-1]+(read_length-1)
    boundaries = tuple(range(a, b, read_length))
    reads = [(i,j-1) for i,j in zip(boundaries,boundaries[1:])]
    obs_dict = collections.defaultdict(list)
    reads_iterator = iter(reads)
    read = next(reads_iterator, None)
    for i, pos in enumerate(positions):
        while read:
            if read[0]<=pos<=read[1]:
                obs_dict[read].append(obs_tab[i])
                break
            read = next(reads_iterator, None)
    return obs_dict

def sort_obs_tab(obs_dict, handle_multiple_observations):
    obs_tab = list()
    
    for rows in obs_dict.values():
        if len(rows)==1:
            obs_tab.extend(rows)
        else:
            warnings.warn('Multiple reads were found to overlap at one or more SNP positions.')
            if handle_multiple_observations=='all':
                obs_tab.extend(rows)
            elif handle_multiple_observations=='first':
                if len(rows)>0: obs_tab.append(rows[0])
            elif handle_multiple_observations=='random':
                if len(rows)>0: obs_tab.append(random.choice(rows))
            elif handle_multiple_observations=='skip':
                pass
            else:
                raise Exception('error: handle_multiple_observations only supports the options \"skip\", \"all\", \"first\" and \"random\".')
  
    obs_tab_sorted =  sorted(obs_tab, key=operator.itemgetter(0))
    
    return obs_tab_sorted

def build_cache(obs_tabs,chr_id, read_length):
    """ Caches the dictionary of reads. """
    
    cache = [{} for _ in range(len(obs_tabs))]
    for i,obs_tab in enumerate(obs_tabs):
        for offset in range(read_length):
            cache[i].update(build_dictionary_of_reads(obs_tab,chr_id,read_length,offset)) 

    return cache

def build_obs_dict(cache, chr_id, read_length, depth, recombination_spot):
    """ Mixes simulated reads to mimic Illumina dye sequencing of a trisomic cell. """ 
    
    num_of_reads = number_of_reads(chr_id,read_length,depth)
         
    L = len(cache)
        
    obs_dict = collections.defaultdict(list)
    
    for i in range(num_of_reads):
        p = randrange(chr_length(chr_id))+1
        read_boundaries = (p,p+read_length-1)
        if L==1:
            W = [1,]
        elif L==2:
            W = [2,1]
        elif L==3:
            W = [1,1,1] if p > recombination_spot * chr_length(chr_id) else [2,1,0]
        else:
            raise Exception('Error: given more than three OBS files.')   
           
        rnd = choices(range(len(W)), weights=W, k=1)[0]
        reads_id = '%d.%d.%s.%d' % (read_boundaries[0],read_boundaries[1],chr(65+rnd),i) 
        for pos, impute2_ind, _, obs_base in cache[rnd].get(read_boundaries,[]):
            obs_dict[pos].append((pos, impute2_ind, reads_id, obs_base))            

    return obs_dict

def MixHaploids(obs_filenames, read_length, depth, **kwargs):
    """ Given N observation tables of haploid sequences, an observation
        table that depicts a trisomic cell is created. """    
        
    time0 = time.time()
    random.seed(a=None, version=2) #I should set a=None after finishing to debug the code.        
    
    handle_multiple_observations = kwargs.get('handle_multiple_observations','all')
    rs = kwargs.get('recombination_spots', 0)
    recombination_spots = [rs] if type(rs) in (float,int) else rs
    work_dir = kwargs.get('work_dir', '')
    given_output_filename = kwargs.get('output_filename','')
    
    work_dir += '/' if len(work_dir)!=0 and work_dir[-1]!='/' else ''
    
    obs_tabs, info_dicts = [], []
    for filename in obs_filenames:
        with open(work_dir + filename, 'rb') as f:
            obs_tabs.append(pickle.load(f))
            info_dicts.append(pickle.load(f))
        
    chr_id = info_dicts[0]['chr_id'] 

    if not all(info['chr_id']==chr_id for info in info_dicts):
        raise Exception('Error: the chr_id differs from one OBS file to another.')   
    
    cache = build_cache(obs_tabs,chr_id, read_length)
    
    for u, recombination_spot in enumerate(recombination_spots,start=1):
        
        obs_dict = build_obs_dict(cache, chr_id, read_length, depth, recombination_spot)
        
        info = info_dicts[0]
        info.update({'depth': depth, 'read_length': read_length})
        
        obs_tab_sorted = sort_obs_tab(obs_dict, handle_multiple_observations)
        
        info['handle_multiple_observations_when_mixing'] = handle_multiple_observations
        
        if given_output_filename!=None:
            MIX = '.'.join(filename.strip().split('/')[-1].strip().split('.')[0] for filename in obs_filenames)
            RECOMB = f'.recomb.{recombination_spot:.2f}' if len(obs_filenames)==3 else ''
            default_output_filename = f'mixed{len(obs_filenames):d}haploids.X{depth:.2f}.' + MIX + '.' + chr_id + RECOMB + '.obs.p'    
            output_filename = default_output_filename if given_output_filename=='' else f'{u:d}.' + given_output_filename 
            with open(  work_dir + output_filename , 'wb' ) as f:
                    pickle.dump( obs_tab_sorted, f, protocol=4)
                    pickle.dump( info, f, protocol=4)
        
        sys.stdout.write('\r')
        sys.stdout.write(f"[{'=' * int(u):{len(recombination_spots)}s}] {int(100*u/len(recombination_spots))}% ")
        sys.stdout.flush()
        
    time1 = time.time()
    print('Done simulating the observations table of a trisomic cell in %.2f sec.' % (time1-time0))
    return tuple(obs_tab_sorted), info

def MixHaploids2(*obs_filenames, read_length, depth, **kwargs):
    return MixHaploids(obs_filenames, read_length, depth, **kwargs)

if __name__ == "__main__":     
    parser = argparse.ArgumentParser(
        description='Simulates observed bases at known SNP positions from a trisomic cell.')
    parser.add_argument('obs_filenames', metavar='OBS_FILENAME', type=str, nargs='+',
                        help='Pickle files created by MAKE_OBS_TAB, containing base observations at known SNP positions.')
    parser.add_argument('-d', '--depth', type=float, 
                        metavar='FLOAT', default=0.1, 
                        help='The average coverage for the whole chromosome.  Default value 0.1')
    parser.add_argument('-l', '--read-length', type=int, 
                        metavar='INT', default=150,
                        help='The number of base pairs (bp) sequenced from a DNA fragment. Default value 150.')
    parser.add_argument('-r', '--recombination-spots', type=list,
                        metavar='FLOAT/LIST', default=0,
                        help='Introduces a transion between SPH and BPH along the chromosome. '
                             'The location of the recombination spot is determined by a fraction, ranging between 0 to 1. '
                             'Giving 0 and 1 as arguemnts means having the SPH and BPH senarios along the entire chromosome, respectively. '
                             'In addition, giving a list of fractions, e.g. 0.2,0.4,0.6, would create a batch of simulations.')
    parser.add_argument('-u', '--handle-multiple-observations', type=str, 
                        metavar='all/first/random/skip', default='all', 
                        help='We expect to observe at most a single base per SNP. When encountering '
                             'an exception the default behavior is to keep all the alleles. However, a '
                             'few alternative options to handle multiple observations are available: ' 
                             '(a) take the first observed base, (b) pick randomly an observed base'
                             'and (c) skip the SNP.')
    parser.add_argument('-o', '--output-filename', metavar='OUTPUT_FILENAME', type=str, default='',
                        help='Output filename. The default filename is a combination of both obs filenames.')    
    kwargs = vars(parser.parse_args())
    kwargs['recombination_spots'] = [float(i.strip()) for i in ''.join(kwargs['recombination_spots']).split(',')]
    MixHaploids(**kwargs)    
    sys.exit(0)
