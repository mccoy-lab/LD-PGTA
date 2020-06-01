#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

Simulates observed bases at known SNP positions from an aneuploid cell
with two matched haplotypes (‘single parental homolog’; SPH), by mixing
two haploids. 

MIX2HAPLOIDS

Daniel Ariad (daniel@ariad.org)
May 15st, 2020

"""
import pickle, time, random, operator, collections, warnings, argparse, sys

warnings.formatwarning = lambda message, category, filename, lineno, file=None, line=None: 'Caution: %s\n' % message

###############################################################################

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
    
def check_haploid(obs_filename, read_length, depth, handle_multiple_observations, **kwargs):
    """ Given an observation table of a haploid sequence with a depth coverage
        of at least 1X, a new observation table is created to mimic the outcome
        of an Illumina dye sequencing."""
    
    time0 = time.time()
    
    with open('results/'+obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        info = pickle.load(f)
 
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

def build_obs_dict(obs_tab1, obs_tab2, info, read_length, depth):
    """ Mixes simulated reads to mimic the outcome of an Illumina dye sequencing
        for an aneuploid cell with two matched haplotypes (‘single parental
        homolog’; SPH). """ 
    
    random.seed(a=None, version=2) #I should set a=None after finishing to debug the code.
    
    chr_id = info['chr_id']
    
    num_of_reads = number_of_reads(chr_id,read_length,depth)
    
    cache1, cache2 = dict(), dict()
    
    for offset in range(read_length):
        cache1.update(build_dictionary_of_reads(obs_tab1,chr_id,read_length,offset)) 
        cache2.update(build_dictionary_of_reads(obs_tab2,chr_id,read_length,offset)) 
    
    obs_dict = collections.defaultdict(list)
    
    for i in range(num_of_reads):
        p = random.randrange(chr_length(chr_id))+1
        read_boundaries = (p,p+read_length-1)
        rand_dict, ID = (cache1,'A') if random.randrange(3)==0 else (cache2,'B')  
        reads_id = '%d.%d.%s.%d' % (read_boundaries[0],read_boundaries[1],ID,i) 
        for pos, impute2_ind, _, obs_base in rand_dict.get(read_boundaries,[]):
            obs_dict[pos].append((pos, impute2_ind, reads_id, obs_base))            

    info.update({'depth': depth, 'read_length': read_length})

    return obs_dict, info


def mix2haploids(obs_filename1, obs_filename2, read_length, depth, handle_multiple_observations, **kwargs):
    """ Wraps build_obs_dict. In addition, copies the values of obs_dict into
        the list obs_tab, while handling multiple observations at SNP positions.
        Then stores obs_tab as well as a dictionary with all the arguments that
        were used. """
    
    time0 = time.time()
    
    with open('results/'+obs_filename1, 'rb') as f:
        obs_tab1 = pickle.load(f)
        info1 = pickle.load(f)
    
    with open('results/'+obs_filename2, 'rb') as f:
        obs_tab2 = pickle.load(f)
        info2 = pickle.load(f)
   
    if info1!=info2:
        raise Exception("Error: mismatch between the details of %s and %s." % (obs_filename1, obs_filename2))    
    
    obs_dict, info = build_obs_dict(obs_tab1, obs_tab2, info1, read_length, depth)
            
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
        default_output_filename = ('mixed2haploids.X%f' % depth).rstrip('0')+'.'+obs_filename1.lstrip().split('.')[0]+'.'+obs_filename2    
        output_filename = default_output_filename if kwargs.get('output_filename','')=='' else kwargs.get('output_filename','') 
        with open(  'results/' + output_filename , "wb" ) as f:
                pickle.dump( obs_tab_sorted, f, protocol=4)
                pickle.dump( info, f, protocol=4)    
        
    time1 = time.time()
    print('Done simulating the observations table of an aneuploid cell with SPH in %.2f sec.' % (time1-time0))
    return tuple(obs_tab_sorted), info
    
if __name__ == "__main__":     
    parser = argparse.ArgumentParser(
        description='Simulates observed bases at known SNP positions from an aneuploid cell '
                    'with two matched haplotypes (‘single parental homolog’; SPH), by mixing '
                    'two haploids.')
    parser.add_argument('obs_filename1', metavar='OBS_FILENAME1', type=str, 
                        help='OBS file')
    parser.add_argument('obs_filename2', metavar='OBS_FILENAME2', type=str, 
                        help='OBS file')
    parser.add_argument('-d', '--depth', type=float, 
                        metavar='FLOAT', default=0.01, 
                        help='The average coverage for a whole genome.  Default value 0.01')
    parser.add_argument('-l', '--read-length', type=int, 
                        metavar='INT', default=150,
                        help='The number of base pairs (bp) sequenced from a DNA fragment. Default value 150.')
    parser.add_argument('-u', '--handle-multiple-observations', type=str, 
                        metavar='all/first/random/skip', default='all', 
                        help='We expect to observe at most a single base per SNP. When encountering '
                             'an exception the default behavior is to keep all the alleles. However, a '
                             'few alternative options to handle multiple observations are avaible: ' 
                             '(a) take the first observed base, (b) pick randomly an observed base'
                             'and (c) skip the SNP.')
    parser.add_argument('-o', '--output-filename', metavar='OUTPUT_FILENAME', type=str, default='',
                        help='Output filename. The default filename is a combination of both OBS filenames.')    
    mix2haploids(**vars(parser.parse_args()))    
    sys.exit(0)
