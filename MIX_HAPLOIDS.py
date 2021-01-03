#!/usr/bin/env python3
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
Dec 29th, 2020

"""
import pickle, time, random, operator, collections, warnings, argparse, sys, os
from random import choices, randrange
from itertools import islice

warnings.formatwarning = lambda message, category, filename, lineno, file=None, line=None: 'Caution: %s\n' % message

class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter): 
    pass

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
    
def build_dictionary_of_reads(obs_tab,chr_id,read_length): 
    """ Returns a dictionary that lists all the possible intervals of length
        read_length in the chromosome chr_id. For each iterval it gives   
        a range indices of obs_tab, which correspond to alleles that overlap
        with the interval. """
    
    positions = next(zip(*obs_tab))        
    obs_dict = {}
    for offset in range(read_length):
        reads_iterator = ((i,i+read_length) for i in range(positions[0]-offset,positions[-1]+read_length,read_length))
        read = next(reads_iterator, None)
        for i, pos in enumerate(positions):
            while read:
                if read[0]<=pos<read[1]:
                    obs_dict.setdefault(read,[i,0])[1] = i+1
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

def build_obs_dict(cache, obs_tabs, chr_id, read_length, depth, scenario, recombination_spot):
    """ Mixes simulated reads to DNA sequencing of various aneuploidy landscapes. """ 
    
    num_of_reads = number_of_reads(chr_id,read_length,depth)
    L = len(cache)
        
    obs_dict = collections.defaultdict(list)
    
    for i in range(num_of_reads):
        p = randrange(chr_length(chr_id))+1
        read_boundaries = (p,p+read_length)
        if scenario=='monosomy':
            W = [1,] + [0] * (L - 1)
        elif scenario=='disomy':
            W = [1,1] + [0] * (L - 2)
        elif scenario=='SPH':
            W = [2,1] + [0] * (L - 2)
        elif scenario=='BPH':
            inequality = p > recombination_spot * chr_length(chr_id)
            W = [1,1,1] + [0] * (L - 3) if inequality else [2,1,] + [0] * (L - 2)
        else:
            raise Exception('error: undefined scenario.')   
           
        rnd = choices(range(len(W)), weights=W, k=1)[0]
        reads_id = '%d.%d.%s.%d' % (read_boundaries[0],read_boundaries[1]-1,chr(65+rnd),i) 
        A,B = cache[rnd].get(read_boundaries, (None,None))
        if A!=None:
            for pos, impute2_ind, _, obs_base in islice(obs_tabs[rnd],A,B):
                obs_dict[pos].append((pos, impute2_ind, reads_id, obs_base))            

    return obs_dict

def senarios_iter(sc, rs):
    """ Iterates over the different scenarios, taking into account a possible
    list of recombination spots. """
    
    recombination_spots = {rs} if type(rs) in (float, int) else {*rs}
    scenarios = {sc} if type(sc) is str else {*sc}

    for s in {*scenarios}: 
        if s=='BPH':
            for b in recombination_spots:
                yield s,b
        else:
            yield s, None

def save_results(obs_tab_sorted,info,ind,recombination_spot,given_output_filename,**kwargs):
    """ Saves the simulated observation table togther with the 
        supplementary information. """
    
    suffix = f'.{ind:d}' if ind else ''
    rs = f'.rs{recombination_spot:.2f}' if info['scenario']=='BPH' else ''
    
    default_output_filename = f"simulated.{info['scenario']:s}.{info['chr_id']:s}.x{info['depth']:.3f}.{'.'.join(info['sample_ids']):s}{rs:s}.obs.p"
    output_filename = default_output_filename if given_output_filename=='' else given_output_filename.rsplit('/', 1).pop()+suffix
    
    output_dir = kwargs.get('output_dir', 'results')
    output_dir += '/' if output_dir[-1:]!='/' else ''
    if output_dir!='' and not os.path.exists(output_dir): os.makedirs(output_dir)
    
    with open(  output_dir + output_filename , 'wb' ) as f:
            pickle.dump( obs_tab_sorted, f, protocol=4)
            pickle.dump( info, f, protocol=4)    
    
    return output_dir + output_filename
    
def MixHaploids(obs_filenames, read_length, depth, scenarios, **kwargs):
    """ Given N observation tables of haploid sequences, an observation
        table that depicts a chromosomal aneuploidy is created. """    
        
    time0 = time.time()
    random.seed(a=None, version=2) #I should set a=None after finishing to debug the code.        
    
    handle_multiple_observations = kwargs.get('handle_multiple_observations','all')
    rs = kwargs.get('recombination_spots', 0)
    given_output_filename = kwargs.get('output_filename','')
        
    obs_tabs, info_dicts = [], []
    for filename in obs_filenames:
        with open(filename, 'rb') as f:
            obs_tabs.append(pickle.load(f))
            info_dicts.append(pickle.load(f))
        
    chr_id = info_dicts[0]['chr_id'] 

    if not all(info['chr_id']==chr_id for info in info_dicts):
        raise Exception('Error: the chr_id differs from one OBS file to another.')   
    
    cache = [build_dictionary_of_reads(obs_tab,chr_id,read_length) for i, obs_tab in enumerate(obs_tabs)]
    
    number_of_required_obs_files = {'monosomy': 1, 'disomy': 2, 'SPH': 2, 'BPH': 3}
    
    output_filenames = []
    
    cases = tuple(senarios_iter(scenarios, rs))

    for ind, (scenario, recombination_spot) in enumerate(cases, start=1): 
        
        if len(obs_filenames) < number_of_required_obs_files[scenario]:
            raise Exception(f'error: The {scenario:s} scenario requires at least {number_of_required_obs_files[scenario]:d} observation files.') 
            
        obs_dict = build_obs_dict(cache, obs_tabs, chr_id, read_length, depth, scenario, recombination_spot)
        obs_tab_sorted = sort_obs_tab(obs_dict, handle_multiple_observations)
        
        sample_ids = [info_dicts[i].get('sample_id',obs_filenames[i].strip().rsplit('/',1).pop()[:-6])
                          + info_dicts[i].get('haplotype','')
                                   for i in range(number_of_required_obs_files[scenario])] 
        
        info = {'chr_id': chr_id,
                'depth': depth,
                'read_length': read_length,
                'scenario': scenario,
                'recombination_spot': recombination_spot,
                'sample_ids': sample_ids,
                'handle_multiple_observations_when_mixing': handle_multiple_observations}
        
        if given_output_filename!=None:
            fn = save_results(obs_tab_sorted,info,ind,recombination_spot,given_output_filename,**kwargs)
            output_filenames.append(fn)
        
        sys.stdout.write(f"\r[{'=' * int(ind):{len(cases)}s}] {int(100*ind/len(cases))}% "); sys.stdout.flush()
        
    time1 = time.time()
    print(f'\nDone simulating the observations table of a trisomic cell in {time1-time0:.2f} sec.')
    return output_filenames

def MixHaploids_wrapper(*obs_filenames, read_length, depth, scenarios, **kwargs):
    return MixHaploids(obs_filenames, read_length, depth, scenarios, **kwargs)

if __name__ == "__main__":       
    parser = argparse.ArgumentParser(
        description="Simulates an observation table of various aneuploidy landscapes:\n"
                    "\t(1) Monosomy - a single copy of a chromosome pair.\n"
                    "\t(2) Disomy - two unmatched haplotypes.\n"
                    "\t(3) SPH (`single parental homolog') - trisomy with two matched haplotypes.\n"
                    "\t(4) BPH (`both parental homologs') - trisomy with three unmatched haplotypes.\n"
                    "\t(5) Meiotic trisomy - trisomy with recombination, charecterized by a transition between SPH to BPH.\n", formatter_class=Formatter)
    parser.add_argument('obs_filenames', metavar='OBS_FILENAME', type=str, nargs='+',
                        help='Pickle files created by SIMULATE_HAPLOIDS, containing base observations at known SNP positions.')
    parser.add_argument('-d', '--depth', type=float, 
                        metavar='FLOAT', default=0.1, 
                        help='The average coverage for the whole chromosome.  Default value 0.1')
    parser.add_argument('-l', '--read-length', type=int, 
                        metavar='INT', default=150,
                        help='The number of base pairs (bp) sequenced from a DNA fragment. Default value 150.')
    parser.add_argument('-s', '--scenarios', type=list,
                        metavar='monosomy/disomy/SPH/BPH', default='BPH',
                        help="The simulation supports four scenarios: monosomy/disomy/SPH/BPH."
                             "Giving a list of scenarios, e.g. SPH,BPH would create a batch of simulations.")
    parser.add_argument('-r', '--recombination-spots', type=list,
                        metavar='FLOAT', default=0,
                        help='Relevent only for the BPH scenario. Introduces a transion between SPH and BPH along the chromosome. '
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
