#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

PIPELINE

Daniel Ariad (daniel@ariad.org)
Feb 27, 2020

"""
import time, pickle, re, sys
from multiprocessing import Process

from MIX2HAPLOIDS import mix2haploids
from MIX2HAPLOIDS import check_haploid
from MIX2HAPLOIDS import mix3haploids
from ANEUPLOIDY_TEST import aneuploidy_test

def LDblockHIST(x):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.hist(x,bins=int(len(x)**.5),histtype='step', linewidth=2.2, label='LLR distribution accross LD blocks')
    ax.set_xlabel('Aggregated log-likelihood ratio')
    ax.set_ylabel('Counts')
    ax.set_title('LLR distribution accross LD blocks' )
    #ax.legend()
    plt.show()
    
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

def aneuploidy_test_demo(obs_filename):
    from ANEUPLOIDY_TEST import aneuploidy_test
    args = dict(obs_filename = 'results_EUR/HG00096.HG00097.0.5.BPH.obs.p',
                hap_filename = '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.hap',
                leg_filename = '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.legend',
                block_size = 1e5,
                offset = 0,
                min_reads = 2,
                max_reads = 16,
                output_filename = None)
    args['obs_filename'] = 'results_EUR/' + obs_filename
    args['output_filename'] = 'results_EUR/'+re.sub('(.*)obs','\\1LLR', obs_filename.split('/')[-1],1)    
    LLR_dict, info = aneuploidy_test(**args)
    
    return LLR_dict, info
        
if __name__ == "__main__":
    LLR_dict, info = aneuploidy_test_demo('mixed2haploids.X0.5.HG00096.HG00096.B.hg38.obs.p')
    LLR_dict, info = aneuploidy_test_demo('mixed3haploids.X0.5.HG00096.HG00096.HG00097.A.hg38.obs.p')


    
    
    
    #sys.exit(0)
    #make_obs_tab_demo('SRR10965088.hg38.bam','chr21_EUR_panel.legend')
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

    
    #A = mix3haploids('HG00096.A.hg38.obs.p', 'HG00096.B.hg38.obs.p', 'HG00097.A.hg38.obs.p', 150, 0.5, 'all', 'results_EUR', '')
    #B = mix2haploids('HG00096.A.hg38.obs.p', 'HG00096.B.hg38.obs.p', 150, 0.5, 'all', 'results_EUR', '')


    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.01, 'all', 'results_HapMix', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.02, 'all', 'results_HapMix', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.03, 'all', 'results_HapMix', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.04, 'all', 'results_HapMix', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.05, 'all', 'results_HapMix', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.06, 'all', 'results_HapMix', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.07, 'all', 'results_HapMix', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.08, 'all', 'results_HapMix', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.09, 'all', 'results_HapMix', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.1, 'all', 'results_HapMix', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.2, 'all', 'results_HapMix', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.3, 'all', 'results_HapMix', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.4, 'all', 'results_HapMix', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.5, 'all', 'results_HapMix', '')

    #LLR_dict, info = aneuploidy_test_demo('mixed2haploids.X0.3.SRR10393062.SRR151495.0-2.hg38.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('mixed3haploids.X0.3.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('mixed2haploids.X0.4.SRR10393062.SRR151495.0-2.hg38.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('mixed3haploids.X0.4.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('mixed2haploids.X0.5.SRR10393062.SRR151495.0-2.hg38.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('mixed3haploids.X0.5.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('mixed2haploids.X0.1.SRR10393062.SRR151495.0-2.hg38.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('mixed3haploids.X0.1.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('mixed2haploids.X0.2.SRR10393062.SRR151495.0-2.hg38.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('mixed3haploids.X0.2.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('HG00096.0.5.SPH.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('HG00096.HG00097.0.5.BPH.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('HG00096.0.4.SPH.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('HG00096.HG00097.0.4.BPH.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('HG00096.0.3.SPH.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('HG00096.HG00097.0.3.BPH.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('HG00096.0.2.SPH.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('HG00096.HG00097.0.2.BPH.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('HG00096.0.1.SPH.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('HG00096.HG00097.0.1.BPH.obs.p')    
        
    #aneuploidy_test_demo('mixed2haploids.X0.02.SRR10393062.SRR151495.0-2.hg38.obs.p')

    
    #A = lambda: aneuploidy_test_demo('mixed2haploids.X0.1.SRR10393062.SRR151495.0-2.hg38.obs.p')
    #B = lambda: aneuploidy_test_demo('mixed2haploids.X0.2.SRR10393062.SRR151495.0-2.hg38.obs.p')
    #C = lambda: aneuploidy_test_demo('mixed2haploids.X0.3.SRR10393062.SRR151495.0-2.hg38.obs.p')
    #D = lambda: aneuploidy_test_demo('mixed2haploids.X0.4.SRR10393062.SRR151495.0-2.hg38.obs.p')
    #E = lambda: aneuploidy_test_demo('mixed2haploids.X0.5.SRR10393062.SRR151495.0-2.hg38.obs.p')
    #runInParallel(A,B,C,D,E)
    
    #aneuploidy_test_demo('mixed2haploids.X0.02.SRR10393062.SRR151495.0-2.hg38.obs.p')
    #aneuploidy_test_demo('mixed2haploids.X0.2.SRR10393062.SRR151495.0-2.hg38.obs.p')
    #aneuploidy_test_demo('mixed2haploids.X0.3.SRR10393062.SRR151495.0-2.hg38.obs.p')
    #aneuploidy_test_demo('mixed2haploids.X0.4.SRR10393062.SRR151495.0-2.hg38.obs.p')
    #aneuploidy_test_demo('mixed2haploids.X0.5.SRR10393062.SRR151495.0-2.hg38.obs.p')
 
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.01, 'all', 'results_HapMix', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.02, 'all', 'results_HapMix', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.03, 'all', 'results_HapMix', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.04, 'all', 'results_HapMix', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.05, 'all', 'results_HapMix', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.06, 'all', 'results_HapMix', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.07, 'all', 'results_HapMix', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.08, 'all', 'results_HapMix', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.09, 'all', 'results_HapMix', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.1, 'all', 'results_HapMix', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.2, 'all', 'results_HapMix', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.3, 'all', 'results_HapMix', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.4, 'all', 'results_HapMix', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.5, 'all', 'results_HapMix', '')

    #aneuploidy_test_demo('mixed3haploids.X0.1.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    #aneuploidy_test_demo('mixed3haploids.X0.2.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    #aneuploidy_test_demo('mixed3haploids.X0.3.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    #aneuploidy_test_demo('mixed3haploids.X0.4.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    #aneuploidy_test_demo('mixed3haploids.X0.5.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')

    #A = aneuploidy_test_demo('mixed3haploids.X0.1.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    #B = aneuploidy_test_demo('mixed3haploids.X0.2.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    #C = aneuploidy_test_demo('mixed3haploids.X0.3.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    #D = aneuploidy_test_demo('mixed3haploids.X0.4.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    #E = aneuploidy_test_demo('mixed3haploids.X0.5.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    #runInParallel(A,B,C,D,E)
  
    #A = check_haploid('SRR10393062.hg38.obs.p', 150, 0.1, 'all')
    #A = check_haploid('SRR151495.0-2.hg38.obs.p', 150, 0.1, 'all')
    #A = check_haploid('SRR151495.0-2.hg38.obs.p', 150, 0.05, 'all')
    #A = check_haploid('SRR151495.0-2.hg38.obs.p', 150, 0.15, 'all')
    #A = check_haploid('SRR151495.0-2.hg38.obs.p', 150, 0.2, 'all')



    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.1)
    #make_llr_dict_demo('mixed2haploids.X0.025.SRR10393062.SRR151495.0-2.hg38.obs.p','chr21','pair',100)
    #make_llr_dict_demo('mixed2haploids.X0.025.SRR10393062.SRR151495.0-2.hg38.obs.p','chr21','triplet',100)
    #make_llr_dict_demo('mixed2haploids.X0.1.SRR10393062.SRR151495.0-2.hg38.obs.p','chr21','quartet',100)

    #make_llr_dict_demo('SRR10965088.hg38.obs.p','chr21','quartet',50)
    #make_llr_dict_demo('SRR10393062.hg38.obs.p','chr21','quartet',25)


    pass
    
    