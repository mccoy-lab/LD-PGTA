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
Jan 14th, 2021

"""
import argparse, sys, os, collections, bz2, gzip
from random import choices, randrange, seed
from collections import defaultdict
from operator import itemgetter
from time import time
from pickle import load, dump

obs_tuple = collections.namedtuple('obs_tuple', ('pos', 'read_id', 'base')) #Encodes the rows of the observations table

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

def inBPHregion(x,transitions):
    """ Checks if x is in a BPH region. """
    l = sorted(transitions[1:]+(x,)).index(x)
    return (transitions[0]=='BPH') ^ (l%2)

def build_obs_tab(obs_dicts, chr_id, read_length, depth, scenario, transitions):
    """ Mixes reads of DNA sequencing to simulate various aneuploidy landscapes
        in non-admixtures as well as recent-admixtures. """

    num_of_reads = number_of_reads(chr_id,read_length,depth)
    L = len(obs_dicts)

    obs_tab = list()
    dx, odd = divmod(read_length, 2)
    for i in range(num_of_reads):
        p = randrange(chr_length(chr_id))+1
        read_boundaries = (p-dx,p+dx+odd)
        if scenario=='monosomy':
            W = [1,] + [0] * (L - 1)
        elif scenario=='disomy':
            W = [1,1] + [0] * (L - 2)
        elif scenario=='SPH':
            W = [2,1] + [0] * (L - 2)
        elif scenario=='BPH':
            W = [1,1,1] + [0] * (L - 3)
        elif scenario=='transitions':
            ####(BPH,0.3,0.7)
            x = inBPHregion(p/chr_length(chr_id),transitions)
            ####print(p/chr_length(chr_id),transitions,x)
            W = [1,1,1] + [0] * (L - 3) if x else [2,1,] + [0] * (L - 2)
        else:
            raise Exception('error: undefined scenario.')

        rnd = choices(range(len(W)), weights=W, k=1)[0]
        reads_id = '%d.%d.%s.%d' % (read_boundaries[0],read_boundaries[1]-1,chr(65+rnd),i)

        obs_tab.extend(obs_tuple(pos, reads_id, obs_dicts[rnd][pos]) for pos in range(*read_boundaries) if pos in obs_dicts[rnd])
        #### obs_dicts[rnd][pos] is the observed base
    obs_tab.sort(key=itemgetter(0))

    return obs_tab

def build_obs_tab_distant(obs_dicts, chr_id, read_length, depth, scenario, transitions, proportions):
    """ Mixes reads of DNA sequencing to simulate various aneuploidy landscapes in distant admixtures. """

    num_of_reads = number_of_reads(chr_id,read_length,depth)
    regions = 40
    obs_tab = list()
    dx, odd = divmod(read_length, 2)

    number_of_haplotypes = {'monosomy':1, 'disomy':2, 'BPH': 3, 'SPH': 3, 'transitions': 3}
    rsize = chr_length(chr_id)//regions
    
    
    rnd_SPH = lambda: sum(choices([[2,0],[0,2]], weights=proportions) + choices([[1,0],[0,1]], weights=proportions),start=[])
    rnd_BPH = lambda: sum(choices([[1,0],[0,1]], weights=proportions, k=number_of_haplotypes[scenario]),start=[])
    if scenario=='SPH':
        configurations = {(i*rsize,(i+1)*rsize): rnd_SPH() for i in range(regions)}
    elif scenario in {'monosomy', 'disomy', 'BPH'}:
        configurations = {(i*rsize,(i+1)*rsize): rnd_BPH() for i in range(regions)}
        
    elif scenario=='transitions':
        configurations = {(i*rsize,(i+1)*rsize): 
                              rnd_BPH() if inBPHregion((i+0.5)/regions,transitions) else rnd_SPH() 
                                  for i in range(regions)}
        
    else:
        raise Exception('error: undefined scenario.')

    for (start,stop),W in configurations.items():
        for k in range(num_of_reads//regions):
            p = randrange(start,stop) + 1
            read_boundaries = (p-dx,p+dx+odd)
            rnd = choices(range(len(W)), weights=W, k=1)[0]
            reads_id = '%d.%d.%s.%d' % (read_boundaries[0],read_boundaries[1]-1,chr(65+rnd),k)

            obs_tab.extend(obs_tuple(pos, reads_id, obs_dicts[rnd][pos]) for pos in range(*read_boundaries) if pos in obs_dicts[rnd])
            #### obs_dicts[rnd][pos] is the observed base
    obs_tab.sort(key=itemgetter(0))

    return obs_tab

def senarios_iter(scenarios, list_of_transitions):
    """ Iterates over the different scenarios, taking into account a batch of
        transitions. """

    for s in scenarios:
        if s=='transitions':
            for b in list_of_transitions:
                yield s,b
        else:
            yield s, None

def save_results(obs_tab,info,ind,transitions,given_output_filename,output_dir,compress):
    """ Saves the simulated observation table togther with the
        supplementary information.  Also, data compression is supported in gzip
        and bzip2 formats. """

    Open = {'bz2': bz2.open, 'gz': gzip.open}.get(compress, open)
    ext = ('.'+compress) * (compress in ('bz2','gz'))

    suffix = f'.{ind:d}' if ind else ''
    T = '.trans_'+transitions[0]+'_'+'_'.join(f'{i:.1f}' for i in transitions[1:]) if info['scenario']=='transitions' else ''

    default_output_filename = f"simulated.{info['scenario']:s}.{info['chr_id']:s}.x{info['depth']:.3f}.{'.'.join(info['sample_ids']):s}{T:s}.obs.p{ext:s}"
    output_filename = default_output_filename if given_output_filename=='' else given_output_filename.rsplit('/', 1).pop()+suffix

    output_dir += '/' if output_dir[-1:]!='/' else ''
    if output_dir!='' and not os.path.exists(output_dir): os.makedirs(output_dir)

    with Open(  output_dir + output_filename , 'wb' ) as f:
            dump(obs_tab, f, protocol=4)
            dump(info, f, protocol=4)

    return output_dir + output_filename

def MixHaploids(obs_filenames, read_length, depth, scenarios, **kwargs):
    """ Given N observation tables of haploid sequences, an observation
        table that depicts a chromosomal aneuploidy is created. """

    time0 = time()
    seed(a=None, version=2) #I should set a=None after finishing to debug the code.

    list_of_transitions = kwargs.get('transitions', [])
    given_output_filename = kwargs.get('output_filename','')
    output_dir = kwargs.get('output_dir', 'results')
    distant_admixture = kwargs.get('distant_admixture', [])
    compress = kwargs.get('compress', None)



    obs_dicts, info_dicts = [], []
    for filename in obs_filenames:
        with open(filename, 'rb') as f:
            obs_dict = {pos: obs_base for (pos, reads_id, obs_base) in load(f)}
            obs_dicts.append(obs_dict)
            info_dicts.append(load(f))

    chr_id = info_dicts[0]['chr_id']

    if not all(info['chr_id']==chr_id for info in info_dicts):
        raise Exception('Error: the chr_id differs from one OBS file to another.')


    number_of_required_obs_files = {'monosomy': 1, 'disomy': 2, 'SPH': 2, 'BPH': 3, 'transitions': 3}
    output_filenames = []

    cases = tuple(senarios_iter(scenarios, list_of_transitions))


    if distant_admixture:
        print('mode: distant admixture (each homolog is a mosaic of haplotypes, where each haplotype is associated with one out of two possible ancestries).')
    else:
        print('mode: non-admixed/recent-admixture (each homolog is associated with a specific ancestry).')

    for ind, (scenario, transitions) in enumerate(cases, start=1):

        if len(obs_filenames) < number_of_required_obs_files[scenario] * (1+(distant_admixture!=[])):
            raise Exception(f'error: The {scenario:s} scenario requires at least {number_of_required_obs_files[scenario]*(1+(distant_admixture==[]))} observation files.')

        if distant_admixture==[]:
            obs_tab = build_obs_tab(obs_dicts, chr_id, read_length, depth, scenario, transitions)
        else:
            obs_tab = build_obs_tab_distant(obs_dicts, chr_id, read_length, depth, scenario, transitions, distant_admixture)

        sample_ids = [info_dicts[i].get('sample_id',obs_filenames[i].strip().rsplit('/',1).pop()[:-6])
                          + info_dicts[i].get('haplotype','')
                                   for i in range(number_of_required_obs_files[scenario])]

        info = {'chr_id': chr_id,
                'depth': depth,
                'read_length': read_length,
                'scenario': scenario,
                'transitions': transitions,
                'sample_ids': sample_ids,
                'handle-multiple-observations': 'all',
                'distant': distant_admixture}

        if given_output_filename!=None:
            fn = save_results(obs_tab,info,ind,transitions,given_output_filename,output_dir,compress)
            output_filenames.append(fn)

        sys.stdout.write(f"\r[{'=' * int(ind):{len(cases)}s}] {int(100*ind/len(cases))}% "); sys.stdout.flush()

    time1 = time()
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
                        help='Pickle files created by EXTRACT_GENOTYPES, containing base observations at known SNP positions. '
                             'When simulating SPH, the first given haploid would serve as the duplicated homolog.')
    parser.add_argument('-d', '--depth', type=float,
                        metavar='FLOAT', default=0.1,
                        help='The average coverage for the whole chromosome.  Default value 0.1')
    parser.add_argument('-l', '--read-length', type=int,
                        metavar='INT', default=36,
                        help='The number of base pairs (bp) sequenced from a DNA fragment. Default value 36.')
    parser.add_argument('-s', '--scenarios', type=str, nargs='+',
                        metavar='monosomy/disomy/SPH/BPH/transitions', default='disomy', choices=['monosomy','disomy','SPH','BPH','transitions'],
                        help="The simulation supports five scenarios: monosomy/disomy/SPH/BPH/transitions. Default scenario is disomy."
                             "Giving a list of scenarios, e.g. \"SPH BPH\" would create a batch of simulations."
                             "In batch mode, the first observation table would be used to simulate monosomy (and the first two tables for disomy).")
    parser.add_argument('-o', '--output-filename', metavar='OUTPUT_FILENAME', type=str,
                        help='Output filename. The default filename is a combination of both obs filenames.')
    parser.add_argument('-c', '--compress', metavar='gz/bz2/unc', type=str, default='unc',  choices=['gz','bz2','unc'],
                        help='Output compressed via gzip, bzip2 or uncompressed. Default is uncompressed.')
    parser.add_argument('-t', '--transitions', type=str, nargs='+', metavar='STR,FLOAT,...,FLOAT',
                        help='Relevant only for the transitions scenario. Introduces transitions between SPH and BPH along the chromosome. '
                             'The locations of the transition is determined by a fraction of chromosome length, ranging between 0 to 1. '
                             'For example a BPH-SPH-BPH transition that equally divides the chromosomes is exressed as BPH,0.333,0.666 and,'
                             'similarly, a SPH-BPH transition at the middle of the chromosome is expressed as SPH,0.5. '
                             'In addition, giving a list of cases, e.g. \"SPH,0.2 SPH,0.4 SPH,0.6\" would create a batch of three simulations. ')
    parser.add_argument('-a', '--distant-admixture', metavar='FLOAT FLOAT', nargs='+', default=[],
                        help='Assume a distant admixture with a certain ancestry proportion, e.g, AFR 0.8 EUR 0.2. '
                             'In addition, the order of observation tables that are given as arguments is important; '
                             'Odd positions are associated with population 1, while even positions with population 2. '
                             'For example, in order to simulate a SPH case the observation tables should be given '
                             'as follows: \"python MIX_HAPLOIDS -s SPH HAPLOID1_AFR.obs.p HAPLOID2_EUR.obs.p HAPLOID3_AFR.obs.p HAPLOID4_EUR.obs.p\". '
                             'When simulating SPH, the first two observation tables would be associated with the duplicated homolog.')

    kwargs = vars(parser.parse_args())
    kwargs['transitions'] = [[float(i) if i.isdigit() else i for i in j.split(',')] for j in kwargs.get('transitions','')]
    MixHaploids(**kwargs)
    sys.exit(0)
