#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

PIPELINE

Daniel Ariad (daniel@ariad.org)
Feb 27, 2020

"""
import time, pickle, re, sys
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
    
def make_obs_tab_demo(bam_filename,legend_filename,handle):
    from MAKE_OBS_TAB import retrive_bases
    args = {'bam_filename': '../BAMs_hg38/'+bam_filename,
            'legend_filename': '../build_reference_panel/ref_panel.HapMix.hg38.BCFtools/'+legend_filename,
            'max_depth': 0,
            'min_bq': 30,
            'min_mq': 30,
            'handle_multiple_observations': handle,
            'fasta_filename': '',#'../genome_ref_hg38/hg38.fa', 
            'output_filename': ''}
    
    args['output_filename'] = 'results/'+re.sub('.bam$','',args['bam_filename'].strip().split('/')[-1])+'.obs.p'

    retrive_bases(**args)
    
    return 0

def aneuploidy_test_demo(obs_filename):
    args = dict(obs_filename = 'results/mixed2haploids.X0.1.SRR10393062.SRR151495.0-2.hg38.obs.p',
                hap_filename = '../build_reference_panel/ref_panel.HapMix.hg38.BCFtools/chr21_HapMix_panel.hap',
                leg_filename = '../build_reference_panel/ref_panel.HapMix.hg38.BCFtools/chr21_HapMix_panel.legend',
                block_size = 1e5,
                offset = 0,
                output_filename = None,
                min_frequency = None,
                min_alleles_per_block = None,
                min_alleles_per_read = 1,
                min_reads_per_block = 2)
    args['obs_filename'] = 'results/' + obs_filename
    args['output_filename'] = 'results/'+re.sub('(.*)obs','\\1LLR', obs_filename.split('/')[-1],1)    
    LLR_dict, info = aneuploidy_test(**args)
    
    return LLR_dict, info
        
if __name__ == "__main__":
    #sys.exit(0)
    #make_obs_tab_demo('SRR10965088.hg38.bam','chr21_EUR_panel.legend')
    #make_obs_tab_demo('SRR6676163.hg38.bam','chr21_MixHap_panel.legend','all')
    #make_obs_tab_demo('SRR151495.0-2.hg38.bam','chr21_HapMix_panel.legend','random')
    #make_obs_tab_demo('SRR10393062.hg38.bam','chr21_HapMix_panel.legend','random')
    #make_obs_tab_demo('SRR8257209.hg38.bam','chr21_EUR_panel.legend')
    #make_reads_bed3_demo('SRR6676163.hg38.bam','chr21_EUR_panel.legend')
    
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.01, 'all', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.02, 'all', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.03, 'all', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.04, 'all', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.05, 'all', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.06, 'all', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.07, 'all', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.08, 'all', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.09, 'all', '')
    #A = mix3haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 'HG002_NA24385_A.hg38.obs.p', 150, 0.1, 'all', '')

    aneuploidy_test_demo('mixed3haploids.X0.01.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    aneuploidy_test_demo('mixed3haploids.X0.02.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    aneuploidy_test_demo('mixed3haploids.X0.03.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    aneuploidy_test_demo('mixed3haploids.X0.04.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    aneuploidy_test_demo('mixed3haploids.X0.05.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    aneuploidy_test_demo('mixed3haploids.X0.06.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    aneuploidy_test_demo('mixed3haploids.X0.07.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    aneuploidy_test_demo('mixed3haploids.X0.08.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    aneuploidy_test_demo('mixed3haploids.X0.09.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')
    aneuploidy_test_demo('mixed3haploids.X0.1.SRR10393062.SRR151495.HG002_NA24385_A.hg38.obs.p')

    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.01, 'all', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.02, 'all', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.03, 'all', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.04, 'all', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.05, 'all', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.06, 'all', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.07, 'all', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.08, 'all', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.09, 'all', '')
    #A = mix2haploids('SRR10393062.hg38.obs.p', 'SRR151495.0-2.hg38.obs.p', 150, 0.1, 'all', '')

    aneuploidy_test_demo('mixed2haploids.X0.01.SRR10393062.SRR151495.0-2.hg38.obs.p')
    aneuploidy_test_demo('mixed2haploids.X0.02.SRR10393062.SRR151495.0-2.hg38.obs.p')
    aneuploidy_test_demo('mixed2haploids.X0.03.SRR10393062.SRR151495.0-2.hg38.obs.p')
    aneuploidy_test_demo('mixed2haploids.X0.04.SRR10393062.SRR151495.0-2.hg38.obs.p')
    aneuploidy_test_demo('mixed2haploids.X0.05.SRR10393062.SRR151495.0-2.hg38.obs.p')
    aneuploidy_test_demo('mixed2haploids.X0.06.SRR10393062.SRR151495.0-2.hg38.obs.p')
    aneuploidy_test_demo('mixed2haploids.X0.07.SRR10393062.SRR151495.0-2.hg38.obs.p')
    aneuploidy_test_demo('mixed2haploids.X0.08.SRR10393062.SRR151495.0-2.hg38.obs.p')
    aneuploidy_test_demo('mixed2haploids.X0.09.SRR10393062.SRR151495.0-2.hg38.obs.p')
    aneuploidy_test_demo('mixed2haploids.X0.1.SRR10393062.SRR151495.0-2.hg38.obs.p')
    
    
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
    
    