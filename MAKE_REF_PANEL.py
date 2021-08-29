#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAKE_REF_PANEL

This script creates reference panels for LD-PGTA, using genotype calls in VCF 
files. The reference panels of LD-PGTA have a similar structure to the IMPUTE2
format.

Daniel Ariad (daniel@ariad.org)
AUG 27, 2022
"""

import sys, os, time, random, argparse, re, pickle, gzip, bz2, collections, itertools

leg_tuple = collections.namedtuple('leg_tuple', ('chr_id', 'pos', 'ref', 'alt')) #Encodes the rows of the legend table
sam_tuple = collections.namedtuple('sam_tuple', ('sample_id', 'group1', 'group2', 'sex')) #Encodes the rows of the samples table
obs_tuple = collections.namedtuple('obs_tuple', ('pos', 'read_id', 'base')) #Encodes the rows of the observations table
    
try:
    import pysam
except ModuleNotFoundError:
    print('Caution: The module pysam is missing.')

def read_impute2(filename,**kwargs):
    """ Reads an IMPUTE2 file format (LEGEND/HAPLOTYPE/SAMPLE) and builds a list
        of lists, containing the dataset. """

    filetype = kwargs.get('filetype', None)

    def leg_format(line):
        rs_id, pos, ref, alt = line.strip().split()
        return leg_tuple('chr'+rs_id[:2].rstrip(':'), int(pos), ref, alt)

    def sam_format(line):
          sample_id, group1, group2, sex = line.strip().split(' ')
          return sam_tuple(sample_id, group1, group2, int(sex))

    with (gzip.open(filename,'rt') if filename[-3:]=='.gz' else open(filename, 'r')) as impute2_in:
        if filetype == 'leg':
            impute2_in.readline()   # Bite off the header
            result = tuple(map(leg_format,impute2_in))

        elif filetype == 'hap':
            firstline = impute2_in.readline()   # Get first line
            a0 = int(firstline.replace(' ', ''), 2)
            a1 = (int(line.replace(' ', ''), 2) for line in impute2_in)
            hap_tab = (a0, *a1)
            number_of_haplotypes = len(firstline.strip().split())
            result = hap_tab, number_of_haplotypes

        elif filetype == 'sam':
            impute2_in.readline()   # Bite off the header
            result = tuple(map(sam_format,impute2_in))

        else:
            result = tuple(line.strip().split() for line in impute2_in)

    return result

def test_module(impute2_leg_filename, impute2_hap_filename, legend, haplotypes):
    """ Compares the IMPUTE2 reference panels to LD-PGTA reference panels. """
    impute2_leg = read_impute2(impute2_leg_filename,filetype='leg')
    impute2_hap = read_impute2(impute2_hap_filename,filetype='hap')
    print('Legend:', all(a==b for a,b in zip(impute2_leg,legend)))
    print('Haplotypes:', all(a==b for a,b in zip(impute2_hap[0],haplotypes[0])))
    return 0

def build_ref_panel(samp_filename,vcf_filename):
    """ Builds a reference panel with similar structure to the IMPUTE2 format.
        The reference panel is encoded for efficient storage and retrieval. """
    
    time0 = time.time()
       
    def sam_format(line):
          sample_id, group1, group2, sex = line.strip().split(' ')
          return sam_tuple(sample_id, group1, group2, int(sex))
    
    with (gzip.open(samp_filename,'rt') if samp_filename[-3:]=='.gz' else open(samp_filename, 'r')) as impute2_in:
        impute2_in.readline()   # Bite off the header
        SAMPLES = tuple(map(sam_format,impute2_in))
    
    vcf_in = pysam.VariantFile(vcf_filename,'r')  # auto-detect input format
    print(vcf_in.description) ### Based on the VCF header, prints a description of the VCF file. 

    SAM = [s.sample_id for s in SAMPLES if s.sample_id in vcf_in.header.samples]

    lenSAM = len(SAM) ### The number of samples that are also included in the VCF.
    
    vcf_in.subset_samples(SAM) ### Read only a subset of samples to reduce processing time and memory. Must be called prior to retrieving records.

    HAPLOTYPES = []
    LEGEND = []
    
    for record in vcf_in.fetch():
        if record.info["VT"]==('SNP',): ### ### Only encode SNPs
        
            phased = all((record.samples[sample].phased for sample in SAM))
            if not phased: continue ### Only encode phased SNPs
        
            ALLELES = tuple(itertools.chain.from_iterable((record.samples[sample].allele_indices for sample in SAM)))
            an = ALLELES.count(1)
            if an==2*lenSAM or an==0: continue ### Only encode SNPs with a non-zero minor allele count.

            LEGEND.append(leg_tuple('chr'+record.contig, record.pos, *record.alleles)) ### Add the record to the legend list        
            binary = sum(v<<i for i, v in enumerate(ALLELES[::-1])) ### Encode the alleles as bits
            HAPLOTYPES.append(binary) ### Add the record to the haplotypes list
                      
    time1 = time.time()
    print('Done building the reference panel in %.3f sec.' % (time1-time0))
    
    result = tuple(LEGEND), (tuple(HAPLOTYPES), 2*lenSAM), SAMPLES

    return result

def save_ref_panel(samp_filename, legend, haplotypes, samples):
    """ Saves the reference panel as a compressed pickle file. """
    time0 = time.time()
    base = samp_filename.rsplit('/', 1)[1].rsplit('.', 1)[0]
    with gzip.open(''.join([legend[0].chr_id,'_',base,'.legend.gz']),'wb') as f:
        pickle.dump(legend,f)
    with gzip.open(''.join([legend[0].chr_id,'_',base,'.hap.gz']),'wb') as f:
        pickle.dump(haplotypes,f)
    with gzip.open(base+'.samples.gz','wb') as f:
        pickle.dump(samples,f)
    time1 = time.time()
    print('Done saving the reference panel in %.3f sec.' % (time1-time0))
    return 0


def main(samp_filename,vcf_filename):
    """ Builds and saves the reference panel. """ 
    legend, haplotypes, samples = build_ref_panel(samp_filename,vcf_filename)
    save_ref_panel(samp_filename, legend, haplotypes, samples)
    return 0
    
        

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='creates reference panels for LD-PGTA, using phased genotypes in VCF files. '
        'The reference panels of LD-PGTA have a similar structure to the IMPUTE2 format. ')
    parser.add_argument('samp_filename', metavar='samples_filename', type=str,
                    help='IMPUTE2 samples file')
    parser.add_argument('vcf_filename', metavar='vcf_filename', type=str,
                        help='IMPUTE2 legend file')

    args = parser.parse_args()
    sys.exit(main(**vars(args)))
else:
    print('The module MAKE_REF_PANEL was imported.')

"""
if __name__ == "__main__": 
    SP ='ALL'
    samp_filename = f'/home/ariad/Dropbox/postdoc_JHU/Project1_LD-PGTA/LD-PGTA_ecosystem/build_reference_panel/samples_per_panel/{SP:s}_panel.samples'
    for i in ['X',*range(22,0,-1)]:    
        #vcf_filename = f'/home/ariad/Dropbox/postdoc_JHU/Project1_LD-PGTA/LD-PGTA_ecosystem/vcf_phase3_hg38_v2/ALL.chr{i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz'
        #legend, haplotypes, samples = build_ref_panel(samp_filename,vcf_filename)
        #save_ref_panel(legend, haplotypes, samples)
        
        impute2_leg_filename = f'/home/ariad/Dropbox/postdoc_JHU/Project1_LD-PGTA/LD-PGTA_ecosystem/build_reference_panel/{SP:s}_panel.hg38.BCFtools/chr{i}_{SP:s}_panel.legend.gz'
        impute2_hap_filename = f'/home/ariad/Dropbox/postdoc_JHU/Project1_LD-PGTA/LD-PGTA_ecosystem/build_reference_panel/{SP:s}_panel.hg38.BCFtools/chr{i}_{SP:s}_panel.hap.gz'
        
        impute2_leg = read_impute2(impute2_leg_filename,filetype='leg')
        impute2_hap = read_impute2(impute2_hap_filename,filetype='hap')
        impute2_sam = read_impute2(samp_filename,filetype='sam')
        save_ref_panel(impute2_leg, impute2_hap, impute2_sam)
        
        
        #print(test_module(impute2_leg_filename, impute2_hap_filename, legend, haplotypes))
"""        
    
    
