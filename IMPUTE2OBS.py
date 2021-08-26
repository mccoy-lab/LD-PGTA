#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IMPUTE2OBS

Simulates an observation table, obs_tab of a haploid, using phased genotypes from IMPUTE2 files.

Daniel Ariad (daniel@ariad.org)
Jan 13, 2021
"""
import pickle, os, sys, time, argparse, random, gzip, collections

obs_tuple = collections.namedtuple('obs_tuple', ('pos', 'read_id', 'base')) #Encodes the rows of the observations table

def read_legend(filename):
    """ Reads an IMPUTE2 legend fileand builds a list
        of lists, containing the dataset. """

    def leg_format(line):
        rs_id, pos, ref, alt = line.strip().split()
        return ('chr'+rs_id[:2].rstrip(':'), int(pos), ref, alt)

    ###with open(filename, 'r') as impute2_in:
    with (gzip.open(filename,'rt') if filename[-3:]=='.gz' else open(filename, 'r')) as impute2_in:
        impute2_in.readline()   # Bite off the header
        result = tuple(map(leg_format,impute2_in))

    return result

def read_haplotypes(sample_filename, hap_filename, sample_id):
    """ Reads an IMPUTE2 samples file as well as an IMPUTE2 haplotypes file
        and returns a list of haplotypes associated with the sample ID. """

    with open(sample_filename, 'r') as impute2_in:
         impute2_in.readline()   # Bite off the header
         samples = [line.split(None, 1)[0] for line in impute2_in]

    if sample_id in samples:
        ind = samples.index(sample_id)
    else:
        raise Exception('Error: sample_id not found.')

    string2tuple = {'0 0':(0,0),'0 1':(0,1),'1 0':(1,0),'1 1':(1,1)}
    a,b = 4*ind, 4*ind+3
    ###with open(hap_filename, 'r') as impute2_in:
    with (gzip.open(hap_filename,'rt') if hap_filename[-3:]=='.gz' else open(hap_filename, 'r')) as impute2_in:
        result = [string2tuple[x[a:b]] for x in impute2_in]

    return result

def main(leg_filename,hap_filename,samp_filename,chr_id,sample_id,**kwargs):

    a = time.time()
    random.seed(None,version=2)

    genotypes = kwargs.get('genotypes', 'AB')

    output_dir = kwargs.get('output_dir', '')
    if output_dir!='' and not os.path.exists(output_dir): os.makedirs(output_dir)
    output_dir += '/' if output_dir[-1:]!='/' else ''

    haplotypes = read_haplotypes(samp_filename, hap_filename, sample_id)
    legend = read_legend(leg_filename)

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

    parser = argparse.ArgumentParser(
        description='Simulates two observation tables of haploids, using phased genotypes from IMPUTE2 files. '
                    'Each observation table includes SNP positions, from the same haplotype of a certain ' 
                    'individual and their corresponding line number within the IMPUTE2 legend file.')
    parser.add_argument('leg_filename', metavar='legend_filename', type=str,
                        help='IMPUTE2 legend file')
    parser.add_argument('hap_filename', metavar='haplotypes_filename', type=str,
                        help='IMPUTE2 haplotypes file')
    parser.add_argument('samp_filename', metavar='samples_filename', type=str,
                        help='IMPUTE2 samples file')
    parser.add_argument('chr_id', metavar='chromosomeID', type=str,
                        help='Chromosome ID')
    parser.add_argument('sample_id', metavar='sampleID', type=str,
                        help='Sample ID')
    parser.add_argument('-g', '--genotypes', metavar='A/B/AB', type=str, default='AB',
                        help='Which of the individual\'s haplotypes should be used. For each specified haplotype, one haploid would be genereated. Default is both (AB).')


    args = parser.parse_args()
    sys.exit(main(**vars(args)))


def test():
    sample_id = 'HG00097'
    chr_id = 'chr1'
    leg_filename = f'../build_reference_panel/AFR_EUR_panel.hg38.BCFtools/{chr_id:s}_AFR_EUR_panel.legend.gz'
    hap_filename = f'../build_reference_panel/AFR_EUR_panel.hg38.BCFtools/{chr_id:s}_AFR_EUR_panel.hap.gz'
    samp_filename = '../build_reference_panel/samples_per_panel/AFR_EUR_panel.samples'

    work_dir='results_TEMP'
    return main(leg_filename,hap_filename,samp_filename,chr_id,sample_id,output_dir=work_dir)
