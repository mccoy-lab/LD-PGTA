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
            'legend_filename': '../build_reference_panel/ref_panel.HapMix_EXT.hg38.BCFtools/'+legend_filename,
            'max_depth': 0,
            'min_bq': 30,
            'min_mq': 30,
            'handle_multiple_observations': handle,
            'fasta_filename': '',#'../genome_ref_hg38/hg38.fa', 
            'output_filename': ''}
    
    args['output_filename'] = 'results_HapMix_EXT/'+re.sub('.bam$','',args['bam_filename'].strip().split('/')[-1])+'.obs.p'

    retrive_bases(**args)
    
    return 0

def aneuploidy_test_demo(obs_filename):
    args = dict(obs_filename = 'results_HapMix_EXT/mixed2haploids.X0.1.SRR10393062.SRR151495.0-2.hg38.obs.p',
                hap_filename = '../build_reference_panel/ref_panel.HapMix_EXT.hg38.BCFtools/chr21_HapMix_EXT_panel.hap',
                leg_filename = '../build_reference_panel/ref_panel.HapMix_EXT.hg38.BCFtools/chr21_HapMix_EXT_panel.legend',
                block_size = 5e4,
                offset = 0,
                min_reads = 16,
                max_reads = 16,
                output_filename = None)
    args['obs_filename'] = 'results_HapMix_EXT/' + obs_filename
    args['output_filename'] = 'results_HapMix_EXT/'+re.sub('(.*)obs','\\1LLR', obs_filename.split('/')[-1],1)    
    LLR_dict, info = aneuploidy_test(**args)
    
    return LLR_dict, info
        
if __name__ == "__main__":
    #sys.exit(0)
    #make_obs_tab_demo('SRR10965088.hg38.bam','chr21_EUR_panel.legend')
    #make_obs_tab_demo('SRR6676163.hg38.bam','chr21_MixHap_EXT_panel.legend','all')
    #make_obs_tab_demo('SRR151495.0-2.hg38.bam','chr21_HapMix_EXT_panel.legend','random')
    #make_obs_tab_demo('SRR10393062.hg38.bam','chr21_HapMix_EXT_panel.legend','random')
    #make_obs_tab_demo('SRR8257209.hg38.bam','chr21_EUR_panel.legend')
    #make_reads_bed3_demo('SRR6676163.hg38.bam','chr21_EUR_panel.legend')
    
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.01, 'all', 'results_HapMix_EXT', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.02, 'all', 'results_HapMix_EXT', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.03, 'all', 'results_HapMix_EXT', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.04, 'all', 'results_HapMix_EXT', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.05, 'all', 'results_HapMix_EXT', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.06, 'all', 'results_HapMix_EXT', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.07, 'all', 'results_HapMix_EXT', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.08, 'all', 'results_HapMix_EXT', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.09, 'all', 'results_HapMix_EXT', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.1, 'all', 'results_HapMix_EXT', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.2, 'all', 'results_HapMix_EXT', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.3, 'all', 'results_HapMix_EXT', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.4, 'all', 'results_HapMix_EXT', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.5, 'all', 'results_HapMix_EXT', '')

    aneuploidy_test_demo('mixed2haploids.X0.1.SRR10393062.SRR151495.0-2.hg38.obs.p')
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
 
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.01, 'all', 'results_HapMix_EXT', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.02, 'all', 'results_HapMix_EXT', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.03, 'all', 'results_HapMix_EXT', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.04, 'all', 'results_HapMix_EXT', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.05, 'all', 'results_HapMix_EXT', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.06, 'all', 'results_HapMix_EXT', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.07, 'all', 'results_HapMix_EXT', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.08, 'all', 'results_HapMix_EXT', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.09, 'all', 'results_HapMix_EXT', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.1, 'all', 'results_HapMix_EXT', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.2, 'all', 'results_HapMix_EXT', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.3, 'all', 'results_HapMix_EXT', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.4, 'all', 'results_HapMix_EXT', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.5, 'all', 'results_HapMix_EXT', '')

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
    
    