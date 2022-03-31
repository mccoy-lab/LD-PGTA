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

import sys, os, time, argparse, pickle, gzip, bz2, collections, itertools, operator

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

def read_impute2(impute2_leg_filename,impute2_hap_filename,impute2_samp_filename,samp_filename):
    """ Reads an IMPUTE2 file format (LEGEND/HAPLOTYPE/SAMPLE) and builds a list
        of lists, containing the dataset. """

    ### LEGEND ###
    
    LEGEND = []    
    Open = {'bz2': bz2.open, 'gz': gzip.open}.get(impute2_leg_filename.rsplit('.',1).pop(), open)
    with Open(impute2_leg_filename,'r') as impute2_in:
        impute2_in.readline()   # Bite off the header
        for line in impute2_in:
            rs_id, pos, ref, alt = line.strip().split()
            LEGEND.append(leg_tuple(rs_id.split(':',1)[0], int(pos), ref, alt))
            
    ### SAMPLE ###
    Open = {'bz2': bz2.open, 'gz': gzip.open}.get(impute2_samp_filename.rsplit('.',1).pop(), open)
    with Open(impute2_samp_filename,'r') as impute2_in:
        impute2_in.readline()   # Bite off the header
        SAMPLES_IN_VCF = {line.strip().split(' ')[0] for line in impute2_in}
    
    
    IMPUTE2_SAMPLE = []
    Open = {'bz2': bz2.open, 'gz': gzip.open}.get(samp_filename.rsplit('.',1).pop(), open)
    with Open(samp_filename,'r') as impute2_in:
        impute2_in.readline()   # Bite off the header
        for line in impute2_in:
            sample_id, group1, group2, sex = line.strip().split(' ')
            if sample_id in SAMPLES_IN_VCF:
                IMPUTE2_SAMPLE.append(sam_tuple(sample_id, group1, group2, int(sex)))
    
    SAMPLE = {group2: tuple(samples) 
                  for group2, samples in itertools.groupby(IMPUTE2_SAMPLE, key=operator.attrgetter('group2'))} #For each group2 gives all the associated sample IDs.    
        
    ### HAPLOTYPES ###
    
    GROUP2 = tuple(s.group2 for s in IMPUTE2_SAMPLE for i in range(2)) 
    HAPLOTYPES = {group2: [] for group2 in SAMPLE}
    
    Open = {'bz2': bz2.open, 'gz': gzip.open}.get(impute2_hap_filename.rsplit('.',1).pop(), open)
    with Open(impute2_hap_filename,'r') as impute2_in:
        for line in impute2_in:
            
            cache = {group2: [] for group2 in SAMPLE}
            for group2, allele in zip(GROUP2,line.replace(' ','')):
                cache[group2].append(allele=='1')
            
            for group2, alleles in cache.items():
                binary = sum(v<<i for i, v in enumerate(reversed(alleles)) if v)
                HAPLOTYPES[group2].append(binary)
    
    HAPLOTYPES = {group2: tuple(hap) for group2, hap in HAPLOTYPES.items()}
       
    HAP_length = {group2: 2*len(samples) for group2, samples in SAMPLE.items()} ### The number of samples that are associated with a particular group2/superpopulation.         
   
    result = tuple(LEGEND), HAPLOTYPES, HAP_length, SAMPLE
        
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
            print('--- Executing:'," ".join(cmd.split()))
            stream = os.popen(cmd)
            print('--- Output:',stream.read())
                      
        LEGEND, HAPLOTYPES, HAP_length, SAMPLE = read_impute2(impute2_leg_filename,impute2_hap_filename,impute2_samp_filename,samp_filename)
        
        if mask_filename!='':
            mask = load_mask_tab_fasta_gz(mask_filename)  
            skipped_SNPs=len(LEGEND)
            for group2 in HAPLOTYPES:
                HAPLOTYPES[group2] = [hap for leg, hap in zip(LEGEND,HAPLOTYPES[group2]) 
                                           if mask[leg.pos-1]=='P']
            LEGEND = tuple(leg for leg in LEGEND if mask[leg.pos-1]=='P')
            skipped_SNPs -= len(LEGEND)
            print(f'--- Based on the genome accessibility mask, {skipped_SNPs:d} SNP records were skipped.')

            
    finally:
        os.remove(indv_filename)
        os.remove(bcf_filename)
        os.remove(impute2_leg_filename)
        os.remove(impute2_hap_filename)
        os.remove(impute2_samp_filename)
    
    time1 = time.time()
    print('Done building the reference panel in %.3f sec.' % (time1-time0))
    return LEGEND, HAPLOTYPES, HAP_length, SAMPLE

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

    IMPUTE2_SAMPLE = parse_samples(samp_filename)
    sampleIDs = [s.sample_id for s in IMPUTE2_SAMPLE if s.sample_id in vcf_in.header.samples] #All sample IDs that exist in both the SAMPLE file and the VCF file.
    
    SAMPLES = {group2: tuple(samples) 
               for group2, samples in itertools.groupby(IMPUTE2_SAMPLE, key=operator.attrgetter('group2'))} #For each group2 gives all the associated sample IDs.    
    
    
    HAP_length = {group2: 2*len(samples) for group2, samples in SAMPLES.items()} ### The number of samples that are associated with a particular group2/superpopulation.
            
    get_samples = {group2: operator.itemgetter(*(s.sample_id for s in SAMPLES[group2])) 
                       for group2 in SAMPLES}
    
    vcf_in.subset_samples(sampleIDs) ### Read only a subset of samples to reduce processing time and memory. Must be called prior to retrieving records.

    HAPLOTYPES = collections.defaultdict(list)
    LEGEND = []
    skipped_SNPs = 0

    for record in vcf_in.fetch():
        if len(record.alleles)==2 and not (record.alleles[0] in ACTG and record.alleles[1] in ACTG): 
            continue ### Keep only records of biallelic SNPs.
        
        if not all((info.phased for info in record.samples.values())): 
            continue ### Only encode phased SNPs

        an = sum(sum(info.allele_indices) for info in record.samples.values())
        if an==2*len(record.samples) or an==0: continue ### Only encode SNPs with a non-zero minor allele count.
        
        if mask!=None and mask[record.pos-1]!='P': ### Record start position on chrom/contig is 1-based inclusive.
            skipped_SNPs +=1
            continue ### Include only SNPs in regions accessible to NGS, according to accessibility masks.

        LEGEND.append(leg_tuple(record.contig, record.pos, *record.alleles)) ### Add the record to the legend list. pos is 1-based inclusive!
        
        for group2, get in get_samples.items():    
            ALLELES = [allele for info in get(record.samples) for allele in info.allele_indices]
            binary = sum(v<<i for i, v in enumerate(reversed(ALLELES)) if v) ### Encode the alleles as bits
            HAPLOTYPES[group2].append(binary) ### Add the record to the haplotypes list
    
    HAPLOTYPES = {k:tuple(v) for k,v in HAPLOTYPES.items()} #Transfer deafultdict into dict and lists into tuples.        
    
    time1 = time.time()
    if mask!=None:
        print(f'--- Based on the genome accessibility mask, {skipped_SNPs:d} SNP records were skipped.')
    print(f'--- The reference panel contains {len(LEGEND):d} SNPs.')
    print('Done building the reference panel in %.3f sec.' % (time1-time0))

    result = tuple(LEGEND), HAPLOTYPES, HAP_length, SAMPLES

    return result

def build_ref_panel_via_cyvcf2(samp_filename,vcf_filename,mask_filename):
    """ Builds a reference panel via cyvcf2 with similar structure to the IMPUTE2 format.
        The reference panel is encoded for efficient storage and retrieval. """

    time0 = time.time()
    print(f'--- Samples Filename: {samp_filename:s}')
    print(f'--- VCF Filename: {vcf_filename:s}')

    get_alleles = operator.itemgetter(0,1)
    ACTG = set('ACTG')
    mask = load_mask_tab_fasta_gz(mask_filename) if mask_filename!='' else None

    vcf_in = cyvcf2.VCF(vcf_filename,'r', strict_gt=True)  # auto-detect input format
    
    assert len(vcf_in.seqnames)==1, 'All records in the VCF must correspond to a single chromosome.'

    REQUESTED_SAMPLES = parse_samples(samp_filename)   
    IMPUTE2_SAMPLE = [s for s in REQUESTED_SAMPLES if s.sample_id in vcf_in.samples] #The requested samples that also exist in the VCF file.
    
    MISSING_SAMPLES = [s.sample_id for s in set(REQUESTED_SAMPLES)-set(IMPUTE2_SAMPLE)] #The requested samples that were missing from the VCF file.
    if MISSING_SAMPLES != []:
        print(f"The following samples were missing from the VCF file: {', '.join(MISSING_SAMPLES):s} ")
    
    sampleIDs = [s.sample_id for s in IMPUTE2_SAMPLE] #All sample IDs that exist in both the SAMPLE file and the VCF file.
    
    len_sampleIDs = len(sampleIDs)
    
    SAMPLES = {group2: tuple(samples) 
               for group2, samples in itertools.groupby(IMPUTE2_SAMPLE, key=operator.attrgetter('group2'))} #For each group2 gives all the associated sample IDs.    
    
    
    HAP_length = {group2: 2*len(samples) for group2, samples in SAMPLES.items()} ### The number of samples that are associated with a particular group2/superpopulation.
            


    vcf_in.set_samples(sampleIDs) ### Read only a subset of samples to reduce processing time and memory. Must be called prior to retrieving records.

    get_samples = {group2: operator.itemgetter(*(vcf_in.samples.index(s.sample_id) for s in SAMPLES[group2])) 
                       for group2 in SAMPLES}
    
    HAPLOTYPES = collections.defaultdict(list)
    LEGEND = []
    skipped_SNPs = 0

    for record in vcf_in():
        if  not (record.is_snp and len(record.ALT)==1 and record.REF in ACTG and record.ALT[0] in ACTG): 
            continue  ### Keep only records of biallelic SNPs
        
        if not record.gt_phases.all(): 
            continue ### Only encode phased SNPs
        
        if record.num_unknown>0 or record.num_hom_ref==len_sampleIDs or record.num_hom_alt==len_sampleIDs: 
            continue ### Only encode SNPs with a non-zero minor allele count.
        
        if mask!=None and mask[record.POS-1]!='P': # According to description of the VCF format, positions are 1-based.
            skipped_SNPs +=1
            continue ### Include only SNPs in regions accessible to NGS, according to accessibility masks.

        LEGEND.append(leg_tuple(record.CHROM, record.POS, record.REF, *record.ALT)) ### Add the record to the legend list
        
        for group2, get_group2 in get_samples.items(): 
            ALLELES = [*itertools.chain.from_iterable(map(get_alleles,get_group2(record.genotypes)))]
            binary = sum(v<<i for i, v in enumerate(reversed(ALLELES)) if v) ### Encode the alleles as bits
            HAPLOTYPES[group2].append(binary) ### Add the record to the haplotypes list
    
    HAPLOTYPES = {k:tuple(v) for k,v in HAPLOTYPES.items()} #Transfer deafultdict into dict and lists into tuples.        
    
    time1 = time.time()
    if mask!=None:
        print(f'--- Based on the genome accessibility mask, {skipped_SNPs:d} SNP records were skipped.')
    print(f'--- The reference panel contains {len(LEGEND):d} SNPs.')
    print('Done building the reference panel in %.3f sec.' % (time1-time0))

    result = tuple(LEGEND), HAPLOTYPES, HAP_length, SAMPLES

    return result

def remove_duplicates(legend, haplotypes):
    """ In some VCF files, multi-allelic SNPs are splited into separate rows
        that share the same chromosomal position. After the spliting the rows 
        of the multiallelic SNPs have the same representation as biallelic SNP,
        so these multiallelic SNPs are not removed when filtering multi-allelic
        SNPs."""
        
    
    indices = []
    for index, previous_leg, next_leg, leg in zip(range(len(legend)), legend[-1:]+legend[:-1], legend[1:]+legend[:1], legend):
        if leg.pos != previous_leg.pos and leg.pos != next_leg.pos:
            indices.append(index)
            
    get_indices = operator.itemgetter(*indices)
    
    legend_filtered = get_indices(legend)
    haplotypes_filtered = {group2: get_indices(hap) for group2, hap in haplotypes.items()}
            
    
    return legend_filtered, haplotypes_filtered
            

def save_ref_panel(samp_filename, legend, haplotypes, number_of_haplotypes, samples, output_dir):
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
        pickle.dump(number_of_haplotypes,f)
    with gzip.open(path+strip_samp_filename+'.sample.gz','wb') as f:
        pickle.dump(samples,f)
    time1 = time.time()
    print('Done saving the reference panel in %.3f sec.' % (time1-time0))
    return 0


def main(samp_filename,vcf_filename,ignore_duplicates,mask,output_directory,force_module):
    """ Builds and saves the reference panel. """

    if force_module=='pysam' or handle_vcf+force_module == 'pysam':
        if 'pysam' not in sys.modules:
            global pysam; import pysam
        print('--- Creating the reference panel via the module pysam.')
        legend, haplotypes, number_of_haplotypes, samples = build_ref_panel_via_pysam(samp_filename,vcf_filename,mask)
    elif force_module=='cyvcf2' or handle_vcf+force_module == 'cyvcf2':
        print('--- Creating the reference panel via the module cyvcf2.')
        legend, haplotypes, number_of_haplotypes, samples = build_ref_panel_via_cyvcf2(samp_filename,vcf_filename,mask)
    else:
        print('--- Creating the reference panel via bcftools.')
        legend, haplotypes, number_of_haplotypes, samples = build_ref_panel_via_bcftools(samp_filename,vcf_filename,mask)

    if ignore_duplicates:    
        print('--- Ignoring multiple records with the same chromosomal position.')
        legend, haplotypes = remove_duplicates(legend, haplotypes)

    save_ref_panel(samp_filename, legend, haplotypes, number_of_haplotypes, samples, output_directory)
    return 0

def test(samp_filename,vcf_filename,mask):
    """ This function creates a reference panel via the modules
        pysam, cyvcf2 and bcftools and checkes that they agree
        with each other."""

    global pysam; import pysam    
    global cyvcf2; import cyvcf2
    
    legend = {}
    haplotypes = {}
    samples = {}
    
    print('--- Creating the reference panel via the module pysam.')
    legend['pysam'], haplotypes['pysam'], samples['pysam'] = build_ref_panel_via_pysam(samp_filename,vcf_filename,mask)
    print('--- Creating the reference panel via the module cyvcf2.')
    legend['cyvcf2'], haplotypes['cyvcf2'], samples['cyvcf2'] = build_ref_panel_via_cyvcf2(samp_filename,vcf_filename,mask)
    print('--- Creating the reference panel via bcftools.')
    legend['bcftools'], haplotypes['bcftools'], samples['bcftools'] = build_ref_panel_via_bcftools(samp_filename,vcf_filename,mask)
    print("legend['pysam']==legend['cyvcf2']:", legend['pysam']==legend['cyvcf2'])
    print("legend['pysam']==legend['bcftools']:", legend['pysam']==legend['bcftools'])
    print("legend['cyvcf2']==legend['bcftools']:", legend['cyvcf2']==legend['bcftools'])
    print("haplotypes['pysam']==haplotypes['cyvcf2']:", haplotypes['pysam']==haplotypes['cyvcf2'])
    print("haplotypes['pysam']==haplotypes['bcftools']:", haplotypes['pysam']==haplotypes['bcftools'])
    print("haplotypes['cyvcf2']==haplotypes['bcftools']:", haplotypes['cyvcf2']==haplotypes['bcftools'])
    print("samples['pysam']==samples['cyvcf2']:", samples['pysam']==samples['cyvcf2'])
    print("samples['pysam']==samples['bcftools']:", samples['pysam']==samples['bcftools'])
    print("samples['cyvcf2']==samples['bcftools']:", samples['cyvcf2']==samples['bcftools'])
    
    return 0

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='creates reference panels for LD-PGTA, using phased genotypes in VCF files. '
        'The reference panels of LD-PGTA have a similar structure to the IMPUTE2 format. ')
    parser.add_argument('samp_filename', metavar='samples_filename', type=str,
                    help='IMPUTE2 sample file')
    parser.add_argument('vcf_filename', metavar='vcf_filename', type=str,
                        help='IMPUTE2 legend file')
    parser.add_argument('-i','--ignore-duplicates', action='store_true', default=True,
                        help='Ignore multiple records with the same chromosomal position.')
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


#samp_filename = '/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/reference_panels/samples_per_panel/EAS_EUR_panel.sample'
#vcf_filename = '/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/1000_genomes_30x_on_GRCh38_3202_samples/CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.filtered.shapeit2-duohmm-phased.vcf.gz'
#test(samp_filename,vcf_filename,mask='')

"""
if __name__ == "__main__":
    from multiprocessing import Process
    process = 32
    proc = []
    force_module = ''
    for SP in {'EUR','EAS','SAS','AFR','AMR','EAS_EUR','EUR_SAS','EAS_SAS',
               'AFR_EUR','AFR_AMR','AMR_EUR','AMR_EAS','AMR_SAS','AFR_SAS',
               'AFR_EAS','AMR_EUR_SAS','EAS_EUR_SAS','AMR_EAS_SAS','AFR_AMR_EUR',
               'AFR_AMR_SAS','AFR_AMR_EAS'}:
        #SP = 'ALL'
        print(SP)
        samp_filename = f'/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/reference_panels/samples_per_panel/{SP:s}_panel.sample'
        output_directory = f'/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/reference_panels/{SP:s}_panel'
        for i in ['X',*range(22,0,-1)]:
            print(i)
            if i=='X':
                #vcf_filename = '/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/1000_genomes_30x_on_GRCh38_3202_samples/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz'
                vcf_filename = f'/home/ariad/Dropbox/postdoc_JHU/Project1_LD-PGTA/LD-PGTA_ecosystem/vcf_phase3_hg38_v2/ALL.chrX.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz'
            else:
                #vcf_filename = f'/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/1000_genomes_30x_on_GRCh38_3202_samples/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{i}.filtered.shapeit2-duohmm-phased.vcf.gz'
                vcf_filename = f'/home/ariad/Dropbox/postdoc_JHU/Project1_LD-PGTA/LD-PGTA_ecosystem/vcf_phase3_hg38_v2/ALL.chr{i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz'
            mask_filename = ''# f'/home/ariad/Dropbox/postdoc_JHU/Project1_LD-PGTA/LD-PGTA_ecosystem/mask/20160622.chr{i}.mask.fasta.gz'
            
            #main(samp_filename,vcf_filename,mask_filename,output_directory,force_module)
            try:
                p = Process(target=main,args=(samp_filename,vcf_filename,mask_filename,output_directory,force_module))
                p.start()
                proc.append(p)
                
            except Exception as error:
                print('caution: a process failed!')
                print(error)
            
            while( len(proc)==process ):
                for p in proc:
                    if not p.is_alive():
                        p.close()
                        proc.remove(p)
    for p in proc:
        p.join()            
"""     
"""
for SP in {'EUR','EAS','SAS','AFR','AMR','EAS_EUR','EUR_SAS','EAS_SAS',
           'AFR_EUR','AFR_AMR','AMR_EUR','AMR_EAS','AMR_SAS','AFR_SAS',
           'AFR_EAS','AMR_EUR_SAS','EAS_EUR_SAS','AMR_EAS_SAS','AFR_AMR_EUR',
           'AFR_AMR_SAS','AFR_AMR_EAS'}:    

    for i in ['X',*range(22,0,-1)]:
        print(SP,i)


        hap_filename = f'/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/reference_panels_NYGC/{SP}_panel/chr{i}_{SP}_panel.hap.gz'
        leg_filename = f'/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/reference_panels_NYGC/{SP}_panel/chr{i}_{SP}_panel.legend.gz'
        
        load = lambda filename: {'bz2': bz2.open, 'gz': gzip.open}.get(filename.rsplit('.',1)[1], open)  #Adjusts the opening method according to the file extension.
        
        open_hap = load(hap_filename)
        open_leg = load(leg_filename)

        with open_hap(hap_filename,'rb') as hap_in:
            hap_tab_per_group = pickle.load(hap_in)
            number_of_haplotypes_per_group = pickle.load(hap_in)
        
        with open_leg(leg_filename,'rb') as leg_in:
            leg_tab = pickle.load(leg_in)
                    
        leg_tab, hap_tab_per_group = remove_duplicates(leg_tab, hap_tab_per_group)
        
        with open_hap(hap_filename,'wb') as hap_in:
            pickle.dump(hap_tab_per_group,hap_in)
            pickle.dump(number_of_haplotypes_per_group,hap_in)
        
        with open_leg(leg_filename,'wb') as leg_in:
            pickle.dump(leg_tab,leg_in)
"""