#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SIMULATE_OBS_TAB

Simulates an observation table, obs_tab of a haploid, using phased genotypes from VCF files.

Daniel Ariad (daniel@ariad.org)
Aug 31, 2020
"""
import pickle, sys, os, time, argparse, random, string

def get_random_string(length):
    letters = string.ascii_lowercase
    result = ''.join(random.choice(letters) for i in range(length))
    return result

def get_alleles_from_bcftools(vcf_file,chr_id,sample_id,bcftools_dir):
    """ Runs bcftools as a subprocess and returns the output """
    
    if not os.path.exists('tmp'): os.makedirs('tmp')
    tmp_filename = 'tmp/'+get_random_string(8)+'.tmp'
    ins_dir = bcftools_dir if bcftools_dir=='' or bcftools_dir[-1]=='/' else bcftools_dir+'/'
    cmd1 = '%sbcftools view %s --samples %s --targets %s --phased --exclude-types indels,mnps,ref,bnd,other --output-type v' % (ins_dir,vcf_file,sample_id,chr_id[3:])
    cmd2 = '%sbcftools query -f \'%%CHROM\\t%%POS\\t[%%TGT]\\n\' - > %s' % (ins_dir,tmp_filename)
    run = cmd1 + ' | ' + cmd2
    print('Calling bcftools:\n%s' % run)
    os.system(run)
    return tmp_filename 


def read_impute2(impute2_filename,**kwargs):
    """ Reads an IMPUTE2 file format (SAMPLE/LEGEND/HAPLOTYPE) and builds a list
        of lists, containing the dataset. """

    filetype = kwargs.get('filetype', None)
    with open(impute2_filename, 'r') as impute2_in:
        if filetype == 'legend':
            impute2_in.readline()   # Bite off the header
            def parse(x):
                y=x.strip().split()
                y[0] = 'chr'+y[0].split(':')[0]
                y[1]=int(y[1])
                return y
        elif filetype == 'hap':
            def parse(x):
                return [i=='1' for i in x.strip().split()]
        else:
            def parse(x):
                return x.strip().split() # Trailing whitespaces stripped from the ends of the string. Then it splits the string into a list of words.

        impute2_tab = [parse(line) for line in impute2_in]
    return impute2_tab

def parse(x):
    """ Parses the output of the bcftools query. """ 
    
    line = x.strip().split('\t')
    result = ['chr'+line[0],int(line[1])]+line[-1].strip().split('|')
    return result


def main(vcf_filename,leg_filename,chr_id,sample_id,bcftools_dir,**kwargs):
    
    a = time.time()
    random.seed(None,version=2)
    output_dir = kwargs.get('output_dir', '')
    output_dir += '/' if len(output_dir)!=0 and output_dir[-1]!='/' else ''
    
    
    
    tmp_filename = get_alleles_from_bcftools(vcf_filename,chr_id,sample_id,bcftools_dir)

    with open(tmp_filename, 'r') as txtfile:
        tab = [parse(line) for line in txtfile]
    os.remove(tmp_filename)

    REF = {pos:ref for chr_id,pos,ref,alt in tab}
    ALT = {pos:alt for chr_id,pos,ref,alt in tab}

    impute2_tab = read_impute2(leg_filename,filetype='legend')

    obs_tab1 = tuple((pos, impute2_index, 'XXX', REF[pos])
                     for impute2_index,(chr_id,pos,ref,alt) in enumerate(impute2_tab)
                     if pos in REF)

    obs_tab2 = tuple((pos, impute2_index, 'XXX', ALT[pos])
                     for impute2_index,(chr_id,pos,ref,alt) in enumerate(impute2_tab)
                     if pos in ALT)

    info = {'chr_id': chr_id}
    
    with open(output_dir+sample_id+'A.%s.hg38.obs.p' % chr_id, 'wb') as binfile:
        pickle.dump(obs_tab1, binfile, protocol=4)
        pickle.dump(info, binfile, protocol=4)

    with open(output_dir+sample_id+'B.%s.hg38.obs.p' % chr_id, 'wb') as binfile:
        pickle.dump(obs_tab2, binfile, protocol=4)
        pickle.dump(info, binfile, protocol=4)

    b = time.time()
    print('Done in %.3f sec.' % ((b-a)))
    
    return 0

if __name__ == "__main__": 
    
    parser = argparse.ArgumentParser(
        description='Simulates two observation tables of a haploid, using phased genotypes from VCF files.'
                    'Each observation table includes SNP positions, alleles in a genotype of an individual'
                    'and their corresponding line number within the IMPUTE2 legend file.')
    parser.add_argument('vcf_filename', metavar='vcf_filename', type=str, 
                        help='VCF file')
    parser.add_argument('leg_filename', metavar='legend_filename', type=str, 
                        help='IMPUTE2 legend file')
    parser.add_argument('chr_id', metavar='chromosomeID', type=str,
                        help='Chromosome ID')
    parser.add_argument('sample_id', metavar='sampleID', type=str,
                        help='Sample ID')
    parser.add_argument('-s', '--bcftools-dir', type=str,
                        default='', metavar='BCFTOOLS_DIR',
                        help='The directory where bcftools are installed.')
    
    args = parser.parse_args()
    sys.exit(main(**vars(args)))

def test():
    bcftools_dir = ''
    sample_id = 'HG00097'
    chr_id = 'chr21'
    #leg_filename = '../build_reference_panel/ref_panel.ALL.hg38.BCFtools/chr21_ALL_panel.legend'
    #leg_filename = '../build_reference_panel/ref_panel.Ashkenazi.hg38.BCFtools/chr21_Ashkenazi_panel.legend'
    #leg_filename = '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.legend'
    leg_filename = '../build_reference_panel/ref_panel.COMMON.hg38.BCFtools/%s_COMMON_panel.legend' % chr_id
    #leg_filename = '../build_reference_panel/ref_panel.TEST.hg38.BCFtools/chr21_TEST_panel.legend'
    vcf_filename = '../vcf_phase3_hg38_v2/ALL.%s.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz' % chr_id
    work_dir='results_COMMON'
    return main(vcf_filename,leg_filename,chr_id,sample_id,bcftools_dir,output_dir=work_dir)