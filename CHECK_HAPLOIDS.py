#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

CHECK_HAPLOIDS

Daniel Ariad (daniel@ariad.org)
May 15st, 2020

"""
import pickle, time, random, operator, collections, warnings, argparse, sys

warnings.formatwarning = lambda message, category, filename, lineno, file=None, line=None: 'Caution: %s\n' % message

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
        with the observed bases at SNPs within the read."""
    
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

def build_obs_dict(obs_tab, info, read_length, depth):
    """ Mixes simulated reads to mimic the outcome of an Illumina dye sequencing
        for an aneuploid cell with two matched haplotypes (‘single parental
        homolog’; SPH). """ 
    
    random.seed(a=None, version=2) #I should set a=None after finishing to debug the code.
    
    chr_id = info['chr_id']
    
    num_of_reads = number_of_reads(chr_id,read_length,depth)
    
    cache = dict()
    
    for offset in range(read_length):
        cache.update(build_dictionary_of_reads(obs_tab,chr_id,read_length,offset)) 
    
    obs_dict = collections.defaultdict(list)
    
    for i in range(num_of_reads):
        p = random.randrange(chr_length(chr_id))+1
        read_boundaries = (p,p+read_length-1)
        reads_id = '%d.%d.%d' % (read_boundaries[0],read_boundaries[1],i) 
        for pos, impute2_ind, _, obs_base in cache.get(read_boundaries,[]):
            obs_dict[pos].append((pos, impute2_ind, reads_id, obs_base))            

    info.update({'depth': depth, 'read_length': read_length})

    return obs_dict, info


def check_haploid(obs_filename, read_length, depth, handle_multiple_observations, **kwargs):
    """ Wraps build_obs_dict. In addition, copies the values of obs_dict into
        the list obs_tab, while handling multiple observations at SNP positions.
        Then stores obs_tab as well as a dictionary with all the arguments that
        were used. """
    
    time0 = time.time()
    
    with open('results/'+obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        info = pickle.load(f)
 
    obs_dict, info = build_obs_dict(obs_tab, info, read_length, depth)
            
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
    
    info['handle_multiple_observations_when_mixing'] = handle_multiple_observations
    
    if kwargs.get('save',True):
        default_output_filename = ('check_haploid.X%f' % depth).rstrip('0')+'.'+obs_filename  
        output_filename = default_output_filename if kwargs.get('output_filename','')=='' else kwargs.get('output_filename','') 
        with open(  'results/' + output_filename , "wb" ) as f:
                pickle.dump( obs_tab_sorted, f )
                pickle.dump( info, f )    
        
    time1 = time.time()
    print('Done in %.2f sec.' % (time1-time0))
    return tuple(obs_tab_sorted), info

