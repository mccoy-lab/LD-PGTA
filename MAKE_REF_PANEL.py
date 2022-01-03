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

import sys, os, time, random, argparse, re, pickle, gzip, bz2, collections, itertools, operator

leg_tuple = collections.namedtuple('leg_tuple', ('chr_id', 'pos', 'ref', 'alt')) #Encodes the rows of the legend table
sam_tuple = collections.namedtuple('sam_tuple', ('sample_id', 'group1', 'group2', 'sex')) #Encodes the rows of the samples table
obs_tuple = collections.namedtuple('obs_tuple', ('pos', 'read_id', 'base')) #Encodes the rows of the observations table


global handle_vcf

try:
    import cyvcf2
    handle_vcf = 'cyvcf2'
except ModuleNotFoundError:
    print('Caution: The module cyvcf2 is missing. Trying to use pysam instead.')
    try:
        import pysam
        handle_vcf = 'pysam'
    except ModuleNotFoundError:
        print('Caution: The module pysam is missing.')
        print('Error: Either the module cyvcf2 or pysam is required. Use cyvcf2 for faster performance.')
        exit(1)

def add_prefix(x):
    """ Adds the suffix 'chr' when missing. """
    return 'chr' + x.removeprefix('chr')


def read_impute2(filename,**kwargs):
    """ Reads an IMPUTE2 file format (LEGEND/HAPLOTYPE/SAMPLE) and builds a list
        of lists, containing the dataset. """

    filetype = kwargs.get('filetype', None)

    def leg_format(line):
        rs_id, pos, ref, alt = line.strip().split()
        return leg_tuple(add_prefix(rs_id.split(':',1)[0]), int(pos), ref, alt)

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

def load_mask_tab_fasta_gz(mask_filename):
    """ Loads accessibility masks from fasta.gz files. """
    print(f'--- Mask filename: {mask_filename:s}')
    try:
        with gzip.open(mask_filename,'rt') as f:
            f.readline() #bite header off
            result = ''.join((i.rstrip('\n') for i in f))
    except Exception as err:
        ### Checks if it is a nested gzip file ###
        if err.__str__()=="\'utf-8\' codec can\'t decode byte 0x8b in position 1: invalid start byte":
            with gzip.open(mask_filename,'rb') as f0:
                with gzip.open(f0,'rt') as f1:
                    f1.readline() #bite header off
                    result = ''.join((i.rstrip('\n') for i in f1))
        else:
            raise err

    return result

def parse_samples(samp_filename):
    """ Parses the samples file. """

    def sam_format(line):
        """ Auxaliry function for parsing a single line in the samples file. """
        sample_id, group1, group2, sex = line.strip().split(' ')
        return sam_tuple(sample_id, group1, group2, int(sex))

    with (gzip.open(samp_filename,'rt') if samp_filename[-3:]=='.gz' else open(samp_filename, 'r')) as impute2_in:
        impute2_in.readline()   # Bite off the header
        SAMPLES = tuple(map(sam_format,impute2_in))

    return SAMPLES


def build_ref_panel_via_bcftools(samp_filename,vcf_filename,mask_filename):
    """ Builds a reference panel via bcftools with similar structure to the IMPUTE2 format.
        The reference panel is encoded for efficient storage and retrieval. """
    time0 = time.time()
    import os, tempfile
    print(f'--- Samples Filename: {samp_filename:s}')
    
    tmp_dir = tempfile._get_default_tempdir()
    indv_filename = tmp_dir +'/' + next(tempfile._get_candidate_names()) 
    bcf_filename = tmp_dir +'/' + next(tempfile._get_candidate_names()) 
    impute2_leg_filename = tmp_dir +'/' + next(tempfile._get_candidate_names()) 
    impute2_hap_filename = tmp_dir +'/' + next(tempfile._get_candidate_names()) 
    impute2_samp_filename = tmp_dir +'/' + next(tempfile._get_candidate_names()) 
    
    
    pipeline = [f"sed '1d' {samp_filename:s} | cut -f 1 -d ' ' > {indv_filename:s}",
                f"bcftools view \
                {vcf_filename:s} \
                --samples-file {indv_filename:s} \
                --exclude-types indels,mnps,ref,bnd,other \
                --min-alleles 2 \
                --max-alleles 2 \
                --min-ac 1:minor \
                --phased \
                --exclude 'AN!=2*N_SAMPLES' \
                --output-file {bcf_filename:s} \
                --output-type u \
                --force-samples",
                f"bcftools convert {bcf_filename:s} --haplegendsample {impute2_hap_filename:s},{impute2_leg_filename:s},{impute2_samp_filename:s}"]
    
    try:
        for cmd in pipeline:
            print('--- Executing:',cmd)
            stream = os.popen(cmd)
            print('--- Output:',stream.read())
                      
        leg_tab = read_impute2(impute2_leg_filename,filetype='leg')
        hap_tab, number_of_haplotypes = read_impute2(impute2_hap_filename,filetype='hap')
        sam_tab = read_impute2(samp_filename,filetype='sam')
        
        if mask_filename!='':
            mask = load_mask_tab_fasta_gz(mask_filename)  
            skipped_SNPs=len(leg_tab)
            hap_tab = tuple(hap for leg, hap in zip(leg_tab,hap_tab) if mask[leg.pos-1]=='P')
            leg_tab = tuple(leg for leg in leg_tab if mask[leg.pos-1]=='P')
            skipped_SNPs-=len(leg_tab)
            print(f'--- Based on the genome accessibility mask, {skipped_SNPs:d} SNP records were skipped.')

            
    finally:
        os.remove(indv_filename)
        os.remove(bcf_filename)
        os.remove(impute2_leg_filename)
        os.remove(impute2_hap_filename)
        os.remove(impute2_samp_filename)
    
    time1 = time.time()
    print('Done building the reference panel in %.3f sec.' % (time1-time0))
    return leg_tab, (hap_tab, number_of_haplotypes), sam_tab

def build_ref_panel_via_pysam(samp_filename,vcf_filename,mask_filename):
    """ Builds a reference panel via pysam with similar structure to the IMPUTE2 format.
        The reference panel is encoded for efficient storage and retrieval. """
    time0 = time.time()
    print(f'--- Samples Filename: {samp_filename:s}')

    ACTG = set('ACTG')
    mask = load_mask_tab_fasta_gz(mask_filename) if mask_filename!='' else None


    vcf_in = pysam.VariantFile(vcf_filename,'r')  # auto-detect input format
    print(f'--- VCF Filename: {vcf_filename:s}')
    print(f'--- VCF description: {vcf_in.description:s}') ### Based on the VCF header, prints a description of the VCF file.

    SAMPLES = parse_samples(samp_filename)
    SAM = [s.sample_id for s in SAMPLES if s.sample_id in vcf_in.header.samples]

    lenSAM = len(SAM) ### The number of samples that are also included in the VCF.
    get_samples = operator.itemgetter(*SAM)


    vcf_in.subset_samples(SAM) ### Read only a subset of samples to reduce processing time and memory. Must be called prior to retrieving records.

    HAPLOTYPES = []
    LEGEND = []
    skipped_SNPs = 0

    for record in vcf_in.fetch():
        if record.info.get('VT')!=('SNP',) and not (record.alleles[0] in ACTG and record.alleles[1] in ACTG): continue
        phased = all((record.samples[sample].phased for sample in SAM))
        if not phased: continue ### Only encode phased SNPs

        ALLELES = tuple(itertools.chain.from_iterable((s.allele_indices for s in get_samples(record.samples))))
        an = ALLELES.count(1)
        if an==2*lenSAM or an==0: continue ### Only encode SNPs with a non-zero minor allele count.
        if mask!=None and mask[record.pos-1]!='P': ### Record start position on chrom/contig is 1-based inclusive.
            skipped_SNPs +=1
            continue ### Include only SNPs in regions accessible to NGS, according to accessibility masks.

        LEGEND.append(leg_tuple(add_prefix(record.contig), record.pos, *record.alleles)) ### Add the record to the legend list. pos is 1-based inclusive!
        binary = sum(v<<i for i, v in enumerate(reversed(ALLELES)) if v) ### Encode the alleles as bits
        HAPLOTYPES.append(binary) ### Add the record to the haplotypes list
    time1 = time.time()
    if mask!=None:
        print(f'--- Based on the genome accessibility mask, {skipped_SNPs:d} SNP records were skipped.')
    print(f'--- The reference panel contains {len(LEGEND):d} SNPs.')
    print('Done building the reference panel in %.3f sec.' % (time1-time0))

    result = tuple(LEGEND), (tuple(HAPLOTYPES), 2*lenSAM), SAMPLES

    return result

def build_ref_panel_via_cyvcf2(samp_filename,vcf_filename,mask_filename):
    """ Builds a reference panel via cyvcf2 with similar structure to the IMPUTE2 format.
        The reference panel is encoded for efficient storage and retrieval. """

    time0 = time.time()
    print(f'--- Samples Filename: {samp_filename:s}')
    reverse_haplotypes = operator.itemgetter(1,0)
    ACTG = set('ACTG')
    mask = load_mask_tab_fasta_gz(mask_filename) if mask_filename!='' else None

    vcf_in = cyvcf2.VCF(vcf_filename,'r', strict_gt=True)  # auto-detect input format

    print(f'--- VCF Filename: {vcf_filename:s}')

    assert len(vcf_in.seqnames)==1, 'All records in the VCF must correspond to a single chromosome.'

    SAMPLES = parse_samples(samp_filename)

    SAM = [s.sample_id for s in SAMPLES if s.sample_id in vcf_in.samples]
    lenSAM = len(SAM) ### The number of samples that are also included in the VCF.
    vcf_in.set_samples(SAM) ### Read only a subset of samples to reduce processing time and memory. Must be called prior to retrieving records.


    X = [s for s in vcf_in.samples if s in SAM]
    reversed_order = operator.itemgetter(*(X.index(s) for s in reversed(SAM))) #Returns the reversed order of samples with respect to the list SAMPLES.

    HAPLOTYPES = []
    LEGEND = []
    skipped_SNPs = 0

    for record in vcf_in():
        if record.INFO.get('VT')!='SNP' and not (record.REF in ACTG and record.ALT[0] in ACTG): continue
        if not record.gt_phases.all(): continue ### Only encode phased SNPs
        if record.num_unknown>0 or record.num_hom_ref==lenSAM or record.num_hom_alt==lenSAM: continue ### Only encode SNPs with a non-zero minor allele count.
        if mask!=None and mask[record.POS-1]!='P': # According to description of the VCF format, positions are 1-based.
            skipped_SNPs +=1
            continue ### Include only SNPs in regions accessible to NGS, according to accessibility masks.

        LEGEND.append(leg_tuple(add_prefix(record.CHROM), record.POS, record.REF, *record.ALT)) ### Add the record to the legend list
        alleles = itertools.chain.from_iterable(map(reverse_haplotypes,reversed_order(record.genotypes)))
        binary = sum(v<<i for i, v in enumerate(alleles) if v) ### Encode the alleles as bits
        HAPLOTYPES.append(binary) ### Add the record to the haplotypes list
    time1 = time.time()
    if mask!=None:
        print(f'--- Based on the genome accessibility mask, {skipped_SNPs:d} SNP records were skipped.')
    print(f'--- The reference panel contains {len(LEGEND):d} SNPs.')
    print('Done building the reference panel in %.3f sec.' % (time1-time0))

    result = tuple(LEGEND), (tuple(HAPLOTYPES), 2*lenSAM), SAMPLES

    return result

def save_ref_panel(samp_filename, legend, haplotypes, samples, output_dir):
    """ Saves the reference panel as a compressed pickle file. """
    time0 = time.time()
    if output_dir.strip()!='' and not os.path.exists(output_dir): os.makedirs(output_dir)
    path = output_dir.strip().rstrip('/')+'/' if output_dir.strip()!='' else ''
    strip_samp_filename = samp_filename.rsplit('/', 1)[1].rsplit('.', 1)[0]
    base_filename = ''.join([legend[0].chr_id,'_',strip_samp_filename])
    with gzip.open(path+base_filename+'.legend.gz','wb') as f:
        pickle.dump(legend,f)
    with gzip.open(path+base_filename+'.hap.gz','wb') as f:
        pickle.dump(haplotypes,f)
    with gzip.open(path+strip_samp_filename+'.samples.gz','wb') as f:
        pickle.dump(samples,f)
    time1 = time.time()
    print('Done saving the reference panel in %.3f sec.' % (time1-time0))
    return 0


def main(samp_filename,vcf_filename,mask,output_directory,force_module):
    """ Builds and saves the reference panel. """

    if force_module=='pysam' or handle_vcf+force_module == 'pysam':
        if 'pysam' not in sys.modules:
            global pysam; import pysam
        print('--- Creating the reference panel via the module pysam.')
        legend, haplotypes, samples = build_ref_panel_via_pysam(samp_filename,vcf_filename,mask)
    elif force_module=='cyvcf2' or handle_vcf+force_module == 'cyvcf2':
        print('--- Creating the reference panel via the module cyvcf2.')
        legend, haplotypes, samples = build_ref_panel_via_cyvcf2(samp_filename,vcf_filename,mask)
    else:
        print('--- Creating the reference panel via bcftools.')
        legend, haplotypes, samples = build_ref_panel_via_bcftools(samp_filename,vcf_filename,mask)

    save_ref_panel(samp_filename, legend, haplotypes, samples, output_directory)
    return 0



if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='creates reference panels for LD-PGTA, using phased genotypes in VCF files. '
        'The reference panels of LD-PGTA have a similar structure to the IMPUTE2 format. ')
    parser.add_argument('samp_filename', metavar='samples_filename', type=str,
                    help='IMPUTE2 samples file')
    parser.add_argument('vcf_filename', metavar='vcf_filename', type=str,
                        help='IMPUTE2 legend file')
    parser.add_argument('-m','--mask', type=str,metavar='GZIPPED_FASTA_FILENAME', default='',
                        help='An accessibility mask file in a gzipped FASTA format.'
                             'Supplying an accessibility mask file will reduce false SNPs in regions of the genome that are less accessible to NGS methods.')
    parser.add_argument('-o','--output-directory', type=str,metavar='PATH', default='',
                        help='The directory in which the reference panel would be created.')
    parser.add_argument('-f','--force-module', type=str,metavar='cyvcf2/pysam/bcftools', default='',
                        help='By deafult cyvcf2 module would be used. This allows to use pysam or bcftools instead. In order to use bcftools, it needs to be in search path.')
    args = parser.parse_args()
    sys.exit(main(**vars(args)))
else:
    print('The module MAKE_REF_PANEL was imported.')

"""
if __name__ == "__main__":
    for SP in 'EUR','EAS','SAS','AFR','AMR','AFR_EUR','EAS_EUR','SAS_EUR','EAS_SAS':
        #SP = 'ALL'
        print(SP)
        samp_filename = f'/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/reference_panels/samples_per_panel/{SP:s}_panel.samples'
        output_directory = f'/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/reference_panels/{SP:s}_panel'
        for i in ['X',*range(22,0,-1)]:
            print(i)
            if i=='X':
                vcf_filename = '/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/1000_genomes_30x_on_GRCh38_3202_samples/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz'
            else:
                vcf_filename = f'/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/1000_genomes_30x_on_GRCh38_3202_samples/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{i}.filtered.shapeit2-duohmm-phased.vcf.gz'
            mask_filename = ''# f'/home/ariad/Dropbox/postdoc_JHU/Project1_LD-PGTA/LD-PGTA_ecosystem/mask/20160622.chr{i}.mask.fasta.gz'
            legend, haplotypes, samples = build_ref_panel_via_cyvcf2(samp_filename,vcf_filename,mask_filename)
            save_ref_panel(samp_filename, legend, haplotypes, samples, output_directory)
    
            #impute2_leg_filename = f'/home/ariad/Dropbox/postdoc_JHU/Project1_LD-PGTA/LD-PGTA_ecosystem/build_reference_panel/{SP:s}_panel.hg38.BCFtools/chr{i}_{SP:s}_panel.legend.gz'
            #impute2_hap_filename = f'/home/ariad/Dropbox/postdoc_JHU/Project1_LD-PGTA/LD-PGTA_ecosystem/build_reference_panel/{SP:s}_panel.hg38.BCFtools/chr{i}_{SP:s}_panel.hap.gz'
    
            #impute2_leg = read_impute2(impute2_leg_filename,filetype='leg')
            #impute2_hap = read_impute2(impute2_hap_filename,filetype='hap')
            #impute2_sam = read_impute2(samp_filename,filetype='sam')
            #save_ref_panel(samp_filename, impute2_leg, impute2_hap, impute2_sam)
    
    
            #print(test_module(impute2_leg_filename, impute2_hap_filename, legend, haplotypes))

"""
