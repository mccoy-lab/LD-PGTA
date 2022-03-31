#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXTRACT_GENOTYPES

Simulates an observation table, obs_tab of a haploid, using phased genotypes from a LD-PGTA reference panel.

Daniel Ariad (daniel@ariad.org)
Jan 13, 2021
"""
import pickle, os, sys, time, argparse, random, gzip, collections

leg_tuple = collections.namedtuple('leg_tuple', ('chr_id', 'pos', 'ref', 'alt')) #Encodes the rows of the legend table
sam_tuple = collections.namedtuple('sam_tuple', ('sample_id', 'group1', 'group2', 'sex')) #Encodes the rows of the samples table
obs_tuple = collections.namedtuple('obs_tuple', ('pos', 'read_id', 'base')) #Encodes the rows of the observations table
    
def get_haplotypes(sample_filename, hap_filename, sample_id, group2):
    """ Extracts haplotypes that correspond to a specific sample ID. """

    with gzip.open(sample_filename, 'rb') as sam_in:
        SAM = pickle.load(sam_in)[group2]
    
    samples = [s.sample_id for s in SAM]
    
    if sample_id in samples:
        ind = len(samples) - samples.index(sample_id)
    else:
        raise Exception('Error: sample_id not found.')
    
    bits2tuple = {0: (0,0), 1: (0,1), 2: (1,0), 3: (1,1)}
    with gzip.open(hap_filename,'rb') as hap_in:
        hap_tab = pickle.load(hap_in)[group2]
    shift_count =  2*(ind-1)
    result = [bits2tuple[(h >> shift_count) & 0b11] for h in hap_tab]

    return result


def extract(leg_filename,hap_filename,samp_filename,chr_id,sample_id,group2,**kwargs):
    """ Builds an observation tables of effective haploids by extracting 
        phased genotypes from a LD-PGTA reference panel. """

    a = time.time()
    random.seed(None,version=2)

    genotypes = kwargs.get('genotypes', 'AB')

    output_dir = kwargs.get('output_dir', '')
    if output_dir!='' and not os.path.exists(output_dir): os.makedirs(output_dir)
    output_dir += '/' if output_dir[-1:]!='/' else ''

    haplotypes = get_haplotypes(samp_filename, hap_filename, sample_id, group2)
    
    with gzip.open(leg_filename,'rb') as leg_in:
        legend = pickle.load(leg_in)

    info = {'chr_id': chr_id,
            'depth': 1,
            'read_length': 1,
            'sample_id': sample_id}

    if genotypes in ('A','AB'):
        obs_tab1 = tuple(obs_tuple(pos, 'XXX', alt if allele1 else ref)
                             for (chrID,pos,ref,alt),(allele1,allele2) in zip(legend,haplotypes)
                                 if chr_id==chrID)
    
        with open(output_dir+sample_id+'A.%s.hg38.obs.p' % chr_id, 'wb') as binfile:
            info1 = {**info, 'haplotype': 'A'}
            pickle.dump(obs_tab1, binfile, protocol=4)
            pickle.dump(info1 , binfile, protocol=4)

    if genotypes in ('B','AB'):
        obs_tab2 = tuple(obs_tuple(pos, 'XXX', alt if allele2 else ref)
                            for (chrID,pos,ref,alt),(allele1,allele2) in zip(legend,haplotypes)
                                if chr_id==chrID)
    
        with open(output_dir+sample_id+'B.%s.hg38.obs.p' % chr_id, 'wb') as binfile:
            info2 = {**info, 'haplotype': 'B'}
            pickle.dump(obs_tab2, binfile, protocol=4)
            pickle.dump(info2, binfile, protocol=4)

    b = time.time()
    print('Done in %.3f sec.' % ((b-a)))

    return 0

if __name__ == "__main__":

    parser = argparse.ArgumentParser( description='Simulates two observation tables of haploids, using phased genotypes from a LD-PGTA reference panel. ')
    
    parser.add_argument('leg_filename', metavar='legend_filename', type=str,
                        help='IMPUTE2 legend file.')
    parser.add_argument('hap_filename', metavar='haplotypes_filename', type=str,
                        help='IMPUTE2 haplotypes file.')
    parser.add_argument('samp_filename', metavar='samples_filename', type=str,
                        help='IMPUTE2 samples file.')
    parser.add_argument('chr_id', metavar='chromosomeID', type=str,
                        help='Chromosome ID.')
    parser.add_argument('sample_id', metavar='sampleID', type=str,
                        help='Sample ID.')
    parser.add_argument('group2', metavar='group2', type=str,
                        help='The group2/superpopulation associated with Sample ID.')
    parser.add_argument('-g', '--genotypes', metavar='A/B/AB', type=str, default='AB',
                        help='Which of the individual\'s haplotypes should be used. For each specified haplotype, one haploid would be genereated. Default is both (AB).')


    args = parser.parse_args()
    sys.exit(extract(**vars(args)))


def test():
    if sys.platform == "linux" or sys.platform == "linux2":
        print('Detected OS: linux.')
        HOME='home'
    elif sys.platform == "darwin":
        print('Detected OS: macOS.')
        HOME='Users'
    elif sys.platform == "win32":
        print('Detected OS: windows.')
        HOME='????????'
        
    sample_id = 'NA06985'
    chr_id = 'chr21'
    sp = 'EUR'
    ref_path = f'/{HOME:s}/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/reference_panels/'
    leg_filename = ref_path + f'{sp:s}_panel/{chr_id:s}_{sp:s}_panel.legend.gz'
    hap_filename = ref_path + f'{sp:s}_panel/{chr_id:s}_{sp:s}_panel.hap.gz'
    samp_filename = ref_path + f'{sp:s}_panel/{sp:s}_panel.sample.gz'

    work_dir='results_TEMP'
    haplotypes = get_haplotypes(samp_filename, hap_filename, sample_id, sp)
    #return extract(leg_filename,hap_filename,samp_filename,chr_id,sample_id,output_dir=work_dir)
