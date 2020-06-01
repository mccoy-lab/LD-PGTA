#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

PIPELINE

Daniel Ariad (daniel@ariad.org)
Feb 27, 2020

"""
import time, pickle, re
from MAKE_OBS_TAB import retrive_bases
from MIX2HAPLOIDS import mix2haploids
from MIX2HAPLOIDS import check_haploid



def make_obs_tab_demo(bam_filename,legend_filename,handle):   
    args = {'bam_filename': '../BAMs_hg38/'+bam_filename,
            'legend_filename': '../build_reference_panel/ref_panel.HapMix.hg38.BCFtools/'+legend_filename,
            'max_depth': 0,
            'min_bq': 30,
            'min_mq': 30,
            'handle_multiple_observations': handle,
            'fasta_filename': '',#'../genome_ref_hg38/hg38.fa', 
            'output_filename': '',
            'save': False}
    
    obs_tab, info  = retrive_bases(**args)
    
    output_filename = re.sub('.bam$','',args['bam_filename'].strip().split('/')[-1])+'.obs.p'
    
    with open( 'results/'+output_filename, "wb") as f:
        pickle.dump(obs_tab, f)
        pickle.dump(info, f)

    return 0

if __name__ == "__main__":
    #make_obs_tab_demo('SRR10965088.hg38.bam','chr21_EUR_panel.legend')
    #make_obs_tab_demo('SRR6676163.hg38.bam','chr21_EUR_panel.legend','all')
    #make_obs_tab_demo('SRR151495.0-2.hg38.bam','chr21_HapMix_panel.legend','first')
    #make_obs_tab_demo('SRR10393062.hg38.bam','chr21_HapMix_panel.legend','first')
    #make_obs_tab_demo('SRR8257209.hg38.bam','chr21_EUR_panel.legend')
    #make_reads_bed3_demo('SRR6676163.hg38.bam','chr21_EUR_panel.legend')
    
    #A = mix2haploids('SRR10393062.hg38.OBS.p', 'SRR151495.0-2.hg38.OBS.p', 150, 0.001 'all')
    A = mix2haploids('SRR10393062.hg38.OBS.p', 'SRR151495.0-2.hg38.OBS.p', 150, 0.01, 'all')
    A = mix2haploids('SRR10393062.hg38.OBS.p', 'SRR151495.0-2.hg38.OBS.p', 150, 0.1, 'all')
    A = mix2haploids('SRR10393062.hg38.OBS.p', 'SRR151495.0-2.hg38.OBS.p', 150, 0.02, 'all')
    A = mix2haploids('SRR10393062.hg38.OBS.p', 'SRR151495.0-2.hg38.OBS.p', 150, 0.03, 'all')
    A = mix2haploids('SRR10393062.hg38.OBS.p', 'SRR151495.0-2.hg38.OBS.p', 150, 0.04, 'all')
    A = mix2haploids('SRR10393062.hg38.OBS.p', 'SRR151495.0-2.hg38.OBS.p', 150, 0.05, 'all')
    A = mix2haploids('SRR10393062.hg38.OBS.p', 'SRR151495.0-2.hg38.OBS.p', 150, 0.06, 'all')
    A = mix2haploids('SRR10393062.hg38.OBS.p', 'SRR151495.0-2.hg38.OBS.p', 150, 0.07, 'all')
    A = mix2haploids('SRR10393062.hg38.OBS.p', 'SRR151495.0-2.hg38.OBS.p', 150, 0.08, 'all')
    A = mix2haploids('SRR10393062.hg38.OBS.p', 'SRR151495.0-2.hg38.OBS.p', 150, 0.09, 'all')

    
    #A = check_haploid('SRR10393062.hg38.OBS.p', 150, 0.1, 'all')
    #A = check_haploid('SRR151495.0-2.hg38.OBS.p', 150, 0.1, 'all')
    #A = check_haploid('SRR151495.0-2.hg38.OBS.p', 150, 0.05, 'all')
    #A = check_haploid('SRR151495.0-2.hg38.OBS.p', 150, 0.15, 'all')
    #A = check_haploid('SRR151495.0-2.hg38.OBS.p', 150, 0.2, 'all')




    #make_llr_dict_demo('SRR10965088.hg38.OBS.p','chr21','pair',25)
    #make_llr_dict_demo('SRR6676163.hg38.OBS.p','chr21','pair',25)   
    #make_llr_dict_demo('SRR151495.0-2.hg38.OBS.p','chr21','pair',10)
    #make_llr_dict_demo('SRR10393062.hg38.OBS.p','chr21','pair',10)
    #make_llr_dict_demo('SRR8257209.hg38.OBS.p','chr21','pair',25)
    #make_llr_dict_demo('SRR151495.0-2.hg38.OBS.p','chr21','triplet', 25)
    #make_llr_dict_demo('SRR10393062.hg38.OBS.p','chr21','triplet', 25)
    #make_llr_dict_demo('SRR8257209.hg38.OBS.p','chr21','triplet',25)
    #make_llr_dict_demo('SRR6676163.hg38.OBS.p','chr21','triplet',100)
    #make_llr_dict_demo('SRR10965088.hg38.OBS.p','chr21','quartet',15)
    #make_llr_dict_demo('SRR151495.0-2.hg38.OBS.p','chr21','quartet', 20)
    #make_llr_dict_demo('SRR10393062.hg38.OBS.p','chr21','quartet', 20)
    #make_llr_dict_demo('SRR8257209.hg38.OBS.p','chr21','quartet',20)
    #make_llr_dict_demo('SRR6676163.hg38.OBS.p','chr21','quartet',20)

    
    #make_llr_dict_demo('mixed2haploids.X0.001.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','pair',100)
    #make_llr_dict_demo('mixed2haploids.X0.005.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','pair',100)
    #make_llr_dict_demo('mixed2haploids.X0.01.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','pair',100)
    #make_llr_dict_demo('mixed2haploids.X0.05.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','pair',100)
    #make_llr_dict_demo('mixed2haploids.X0.1.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','pair',100)
    #make_llr_dict_demo('mixed2haploids.X0.5.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','pair',100)
    #make_llr_dict_demo('mixed2haploids.X1.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','pair',100)
    #make_llr_dict_demo('mixed2haploids.X10.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','pair',100)
    #make_llr_dict_demo('mixed2haploids.X0.001.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','triplet',100)
    #make_llr_dict_demo('mixed2haploids.X0.005.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','triplet',100)
    #make_llr_dict_demo('mixed2haploids.X0.01.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','triplet',100)
    #make_llr_dict_demo('mixed2haploids.X0.05.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','triplet',100)
    #make_llr_dict_demo('mixed2haploids.X0.1.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','triplet',100)
    #make_llr_dict_demo('mixed2haploids.X0.5.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','triplet',100)
    #make_llr_dict_demo('mixed2haploids.X1.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','triplet',100)
    #make_llr_dict_demo('mixed2haploids.X10.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','triplet',100)
    #make_llr_dict_demo('mixed2haploids.X0.05.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','quartet',100)

    #A = mix2haploids('SRR10393062.hg38.OBS.p', 'SRR151495.0-2.hg38.OBS.p', 150, 0.1)
    #make_llr_dict_demo('mixed2haploids.X0.025.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','pair',100)
    #make_llr_dict_demo('mixed2haploids.X0.025.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','triplet',100)
    #make_llr_dict_demo('mixed2haploids.X0.1.SRR10393062.SRR151495.0-2.hg38.OBS.p','chr21','quartet',100)

    #make_llr_dict_demo('SRR10965088.hg38.OBS.p','chr21','quartet',50)
    #make_llr_dict_demo('SRR10393062.hg38.OBS.p','chr21','quartet',25)


    pass
    
    