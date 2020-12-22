#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

PIPELINE

Daniel Ariad (daniel@ariad.org)
Feb 27, 2020

"""
import time, pickle, re, sys
from multiprocessing import Process

def runInParallel(*fns,**kwargs):
    proc = []
    for fn in fns:
        try:
            p = Process(target=fn,args=kwargs.get('args',tuple()))
            p.start()
            proc.append(p)
        except:
            print('caution: a process failed!')
    for p in proc:
        try:
          p.join()
        except:
            None

def aneuploidy_test_demo(obs_filename,chr_id,sp):
    from ANEUPLOIDY_TEST import aneuploidy_test
    args = dict(obs_filename = '/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/' + obs_filename,
                hap_filename = f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.hap',
                leg_filename = f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.legend',
                window_size = 0,
                subsamples = 100,
                offset = 0,
                min_reads = 3,
                max_reads = 5,
                min_HF = 0.05,
                minimal_score = 2,
                output_dir = '/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/', 
                output_filename = '')
    LLR_dict, info = aneuploidy_test(**args)
    return LLR_dict, info
   
def make_obs_tab_demo(bam_filename,chr_id,sp):
    from MAKE_OBS_TAB import retrive_bases
    args = {'bam_filename': '/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/'+bam_filename,
            'output_dir': '/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/',
            'legend_filename': f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.legend',
            'max_depth': 0,
            'min_bq': 30,
            'min_mq': 30,
            'handle_multiple_observations': 'all',
            'fasta_filename': '',#'../genome_ref_hg38/hg38.fa', 
            'output_filename': ''}
    
    result = retrive_bases(**args)
    
    return result

if __name__ == "__main__":
    make_obs_tab_demo('11522FA-AP925_3.bam', 'chr11', 'SAS')
    aneuploidy_test_demo('11522FA-AP925_3.chr11.obs.p', 'chr11', 'SAS')
    pass
    
    