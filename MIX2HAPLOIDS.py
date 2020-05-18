#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

MIX2HAPLOIDS

Daniel Ariad (daniel@ariad.org)
May 15st, 2020

"""
import pickle, time, random, operator, collections, math, warnings

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

def mix2haploids(obs_filename1, obs_filename2, read_length, depth):
    """ Simulates observed bases at known SNP positions from an aneuploid cell
        with two matched haplotypes (‘single parental homolog’; SPH), by mixing
        two haploids. """ 
    
    time0 = time.time()
    random.seed(a=None, version=2) #I should set a=None after finishing to debug the code.
    
    with open('results/'+obs_filename1, 'rb') as f:
        obs_tabA = pickle.load(f)
        infoA = pickle.load(f)
    
    with open('results/'+obs_filename2, 'rb') as f:
        obs_tabB = pickle.load(f)
        infoB = pickle.load(f)
   
    if infoA!=infoB:
        raise Exception("Error: mismatch between the details of %s and %s." % (obs_filename1, obs_filename2))
    else:
        chr_id = infoA['chr_id']
    
    num_of_reads = number_of_reads(chr_id,read_length,depth)
    
    obs_dictA, obs_dictB = dict(), dict()
    
    for offset in range(read_length):
        obs_dictA.update(build_dictionary_of_reads(obs_tabA,chr_id,read_length,offset)) 
        obs_dictB.update(build_dictionary_of_reads(obs_tabB,chr_id,read_length,offset)) 
    
    obs_dict = collections.defaultdict(list)
    
    for i in range(num_of_reads):
        p = random.randrange(chr_length(chr_id))+1
        read_boundaries = (p,p+read_length-1)
        rand_dict = obs_dictA if random.randrange(3)==0 else obs_dictB  
        
        for pos, impute2_ind, _, obs_base in rand_dict.get(read_boundaries,[]):
            reads_id = '%d%d' % read_boundaries            
            if pos in obs_dict:
                print(pos)
                warnings.warn('Multiple reads were found to overlap at one or more SNP positions.')
            obs_dict[pos].append((pos, impute2_ind, reads_id, obs_base))            
            
    obs_tab =  sorted((row for row in obs_dict.values()), key=operator.itemgetter(0))
    info = {**infoA, **{'depth': depth, 'read_length': read_length}}
    
    output_filename  = ('mixed2haploids.X%f' % depth).rstrip('0')+'.'+obs_filename1.lstrip().split('.')[0]+'.'+obs_filename2    
    with open(  'results/' + output_filename , "wb" ) as f:
            pickle.dump( obs_tab, f )
            pickle.dump( info, f )    
    time1 = time.time()
    print('Done in %.2f sec.' % (time1-time0))
    return obs_tab, info

