#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

PIPELINE

Daniel Ariad (daniel@ariad.org)
Feb 27, 2020

"""
import time, pickle, re, sys
from multiprocessing import Process
from MIX_HAPLOIDS import MixHaploids2
from SIMULATE_OBS_TAB import main as simulate

def runInParallel(*fns):
    proc = []
    for fn in fns:
        try:  
            p = Process(target=fn)
            p.start()
            proc.append(p)
        except:
            print('caution: a process failed!')
    for p in proc:
        try:
          p.join()
        except:
            None

def make_obs_tab_demo2(bam_filename,legend_filename,handle):
    from MAKE_OBS_TAB import retrive_bases
    args = {'bam_filename': '../For_Rajiv/'+bam_filename,
            'legend_filename': '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/'+legend_filename,
            'max_depth': 0,
            'min_bq': 30,
            'min_mq': 30,
            'handle_multiple_observations': handle,
            'fasta_filename': '',#'../genome_ref_hg38/hg38.fa', 
            'output_filename': ''}
    
    args['output_filename'] = 'results_EUR/'+re.sub('.bam$','',args['bam_filename'].strip().split('/')[-1])+'.obs.p'

    result = retrive_bases(**args)
    
    return result
 
def aneuploidy_test_demo(obs_filename,chr_id):
    from ANEUPLOIDY_TEST import aneuploidy_test
    args = dict(obs_filename = 'results_COMMON/ABC.obs.p',
                hap_filename = '../build_reference_panel/ref_panel.COMMON.hg38.BCFtools/%s_COMMON_panel.hap' % chr_id,
                leg_filename = '../build_reference_panel/ref_panel.COMMON.hg38.BCFtools/%s_COMMON_panel.legend' % chr_id,
                block_size = 15e4,
                adaptive = True,
                subsamples = 300,
                offset = 0,
                min_reads = 2,
                max_reads = 16,
                output_filename = None)
    args['obs_filename'] = 'results_COMMON/' + obs_filename
    args['output_filename'] = 'results_COMMON/'+re.sub('(.*)obs','\\1LLR_x300_adaptive_LDblocks', obs_filename.split('/')[-1],1)    
    LLR_dict, info = aneuploidy_test(**args)
    return LLR_dict, info
   
def make_obs_tab_demo(bam_filename,legend_filename,handle):
    from MAKE_OBS_TAB import retrive_bases
    args = {'bam_filename': '../BAMs_hg38/'+bam_filename,
            'legend_filename': '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/'+legend_filename,
            'max_depth': 0,
            'min_bq': 30,
            'min_mq': 30,
            'handle_multiple_observations': handle,
            'fasta_filename': '',#'../genome_ref_hg38/hg38.fa', 
            'output_filename': ''}
    
    args['output_filename'] = 'results_EUR/'+re.sub('.bam$','',args['bam_filename'].strip().split('/')[-1])+'.obs.p'

    result = retrive_bases(**args)
    
    return result

def make_simulated_obs_tab():
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
    return simulate(vcf_filename,leg_filename,chr_id,sample_id,bcftools_dir,output_dir=work_dir)
        
if __name__ == "__main__":
    #A = MixHaploids2('HG00096A.chr13.hg38.obs.p', 'HG00096B.chr13.hg38.obs.p', 'HG00097A.chr13.hg38.obs.p', read_length=150, depth=0.20, work_dir='results_EUR', recombination_spots=[i/10 for i in range(11)])
    #A = MixHaploids2('HG00096A.chr6.hg38.obs.p', 'HG00096B.chr6.hg38.obs.p', 'HG00097A.chr6.hg38.obs.p', read_length=36, depth=0.01, work_dir='results_EUR', recombination_spots=[i/10 for i in range(11)])
    #A = MixHaploids2('HG00096A.chr21.hg38.obs.p', 'HG00096B.chr21.hg38.obs.p', 'HG00097A.chr21.hg38.obs.p', read_length=150, depth=0.05, work_dir='results_EUR', recombination_spots=[i/10 for i in range(11)])
    
    #A = MixHaploids2('HG00096A.chr21.hg38.obs.p', 'HG00096B.chr21.hg38.obs.p', 'HG00097A.chr21.hg38.obs.p', read_length=150, depth=0.50, work_dir='results_COMMON', recombination_spots=[i/10 for i in range(11)])
    
    
    #for i in range(0,11):
    #   print('Recombination spot: %.2f' % (i * 0.1))
    #    LLR_dict, info = aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.recomb.%.2f.obs.p' % (i * 0.1))       

    
    filenames = ('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.%.2f.obs.p' % (i * 0.1) for i in range(11))
    functions = (eval('lambda: aneuploidy_test_demo(\'%s\',\'%s\')' % (f,'chr21')) for f in filenames)
    runInParallel(*functions)
    
    #filenames = ('mixed3haploids.X0.10.HG00096A.HG00096B.HG00097A.chr21.recomb.%.2f.obs.p' % (i * 0.1) for i in range(11))
    #functions = (eval('lambda: aneuploidy_test_demo(\'%s\',\'%s\')' % (f,'chr21')) for f in filenames)
    #runInParallel(*functions)
    
    #filenames = ('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.%.2f.obs.p' % (i * 0.1) for i in range(11))
    #functions = (eval('lambda: aneuploidy_test_demo(\'%s\',\'%s\')' % (f,'chr21')) for f in filenames)
    #runInParallel(*functions)
        
    #filenames = ('mixed3haploids.X0.01.HG00096A.HG00096B.HG00097A.chr6.recomb.%.2f.obs.p' % (i * 0.1) for i in range(5,11))
    #functions = (eval('lambda: aneuploidy_test_demo2(\'%s\',\'%s\')' % (f,'chr6')) for f in filenames)
    #runInParallel(*functions)
    
    #LLR_dict, info = aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.0.00.obs.p','chr21')
    #LLR_dict, info = aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.1.00.obs.p','chr21')
    
    
    #A = aneuploidy_test_demo2('11909FA_2.merged.hg38.obs.p','chr6')
    #A = aneuploidy_test_demo2('11694FA_3.merged.hg38.obs.p','chr13')
    
    
    #A = aneuploidy_test_demo2('11909FA_2_T1.chr6.hg38.obs.p','chr6')
    #A = aneuploidy_test_demo2('11909FA_2_T2.chr6.hg38.obs.p','chr6')
    #A = aneuploidy_test_demo2('11909FA_2_I1.chr6.hg38.obs.p','chr6')
    #A = aneuploidy_test_demo2('11909FA_2_I2.chr6.hg38.obs.p','chr6')
    #A = aneuploidy_test_demo2('11694FA_3_T1.chr13.hg38.obs.p','chr13')
    #A = aneuploidy_test_demo2('11694FA_3_T2.chr13.hg38.obs.p','chr13')
    #A = aneuploidy_test_demo2('11694FA_3_I1.chr13.hg38.obs.p','chr13')

    #sys.exit(0)
    #make_obs_tab_demo('SRR10965088.hg38.bam','chr21_EUR_panel.legend')
    
    #make_obs_tab_demo2('chr6/11909FA_2.merged.hg38.bam','chr6_EUR_panel.legend','all')
    #make_obs_tab_demo2('chr13/11694FA_3.merged.hg38.bam','chr13_EUR_panel.legend','all')
    
    #make_obs_tab_demo('SRR6676163.hg38.bam','chr21_MixHap_panel.legend','all')
    #make_obs_tab_demo('SRR151495.0-2.hg38.bam','chr21_HapMix_panel.legend','random')
    #make_obs_tab_demo('SRR10393062.hg38.bam','chr21_HapMix_panel.legend','random')
    #make_obs_tab_demo('SRR8257209.hg38.bam','chr21_HapMix_panel.legend','random')
    #A = make_obs_tab_demo('HG00096.0.5.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG00096.HG00097.0.5.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG00096.0.4.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG00096.HG00097.0.4.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG00096.0.3.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG00096.HG00097.0.3.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG00096.0.2.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG00096.HG00097.0.2.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG00096.0.1.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG00096.HG00097.0.1.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG01880.0.5.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG01880.HG01885.0.5.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG01880.0.4.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG01880.HG01885.0.4.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG01880.0.3.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG01880.HG01885.0.3.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG01880.0.2.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG01880.HG01885.0.2.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG01880.0.1.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG01880.HG01885.0.1.BPH.bam','chr21_EUR_panel.legend','all')

    
    #A = MixHaploids('HG00096.A.hg38.obs.p', 'HG00096.B.hg38.obs.p', 'HG00097.A.hg38.obs.p', read_length=150, depth=0.5, work_dir='results_EUR', recombination_spot=0)
    #B = MixHaploids('HG00096.A.hg38.obs.p', 'HG00097.A.hg38.obs.p', read_length=150, depth=0.5, work_dir='results_EUR', recombination_spot=0)
   
    #for i in range(11):
    #    A = MixHaploids2('HG00096A.hg38.obs.p', 'HG00096B.hg38.obs.p', 'HG00097A.hg38.obs.p', read_length=36, depth=0.01, work_dir='results_EUR', recombination_spot=i/10)
        #B = MixHaploids2('HG00096A.hg38.obs.p', 'HG00097A.hg38.obs.p', read_length=36, depth=0.01, work_dir='results_EUR', recombination_spot=i/10)
    
    #LLR_dict, info = aneuploidy_test_demo('mixed2haploids.X0.50.HG00096A.HG00097A.recomb.0.90.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('mixed2haploids.X0.50.HG00096A.HG00097A.recomb.0.00.obs.p')

    
    
    pass
    
    