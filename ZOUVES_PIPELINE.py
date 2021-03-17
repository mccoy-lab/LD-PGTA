#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

ZOUVES_PIPELINE

Daniel Ariad (daniel@ariad.org)
Nov 18, 2020

"""
import time, re, pickle, os
from multiprocessing import Process

def make_obs_tab_demo(bam_filename,chr_id,sp):
    from MAKE_OBS_TAB import retrive_bases
    args = {'bam_filename': '/home/ariad/Dropbox/postdoc_JHU/LD-PGTA_ecosystem/LD-PGTA_V2/results_ZOUVES/'+bam_filename,
            'output_dir': '/home/ariad/Dropbox/postdoc_JHU/LD-PGTA_ecosystem/LD-PGTA_V2/results_ZOUVES/',
            'legend_filename': f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.legend',
            'max_depth': 0,
            'min_bq': 30,
            'min_mq': 30 if 'chrX'!=chr_id!='chrY' else 0,
            'handle_multiple_observations': 'all',
            'fasta_filename': '',#'../genome_ref_hg38/hg38.fa', 
            'output_filename': ''}
    
    result = retrive_bases(**args)
    
    return result
    
def aneuploidy_test_demo(obs_filename,chr_id,sp):
    from ANEUPLOIDY_TEST import aneuploidy_test
    args = dict(obs_filename = '/home/ariad/Dropbox/postdoc_JHU/LD-PGTA_ecosystem/LD-PGTA_V2/results_ZOUVES/' + obs_filename,
                hap_filename = f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.hap',
                leg_filename = f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.legend',
                window_size = 0,
                subsamples = 100,
                offset = 0,
                min_reads = 6,
                max_reads = 4,
                min_HF = 0.05,
                minimal_score = 2,
                output_dir = '/home/ariad/Dropbox/postdoc_JHU/LD-PGTA_ecosystem/LD-PGTA_V2/results_ZOUVES/', 
                output_filename = '')
    LLR_dict, info = aneuploidy_test(**args)
    return LLR_dict, info

if __name__ == "__main__":
    CHECK = False
    with open('/home/ariad/Dropbox/postdoc_JHU/Zouves-BlueFuse/Play/all.p', 'rb') as f:
        db_TEST = pickle.load(f)
   
    DONE = []
    ERRORS = []
    output_dir = '/home/ariad/Dropbox/postdoc_JHU/LD-PGTA_ecosystem/LD-PGTA_V2/results_ZOUVES/'
    for case in db_TEST:
        if case not in DONE:
            bam_filename = case['filename']
            print(case['filename'])
            sp = case['sp']
            proc = []
            try:
                for chr_num in case['chr_num']:
                    chr_id = 'chr'+str(chr_num)
                    obs_filename = re.sub('.bam$','',bam_filename.split('/')[-1]) + f'.{chr_id:s}.obs.p'
                    LLR_filename = re.sub('.bam$','',bam_filename.split('/')[-1]) + f'.{chr_id:s}.LLR.p'
                    
                    if not os.path.isfile(output_dir+obs_filename) or not CHECK:   
                        make_obs_tab_demo(case['filename'],chr_id, sp)
                    else:
                        print(f'{obs_filename:s} already exists.')
                    
                    if not os.path.isfile(output_dir+LLR_filename) or not CHECK: 
                        #aneuploidy_test_demo(obs_filename, chr_id, sp)
                        p = Process(target=aneuploidy_test_demo,args=(obs_filename, chr_id, sp))
                        p.start()
                        proc.append(p)
                    else:
                        print(f'{LLR_filename:s} already exists.')
            except Exception as error: 
                print('ERROR: ', error)
                if os.path.isfile(output_dir+obs_filename): os.remove(output_dir+obs_filename)
                if os.path.isfile(output_dir+LLR_filename): os.remove(output_dir+LLR_filename)
                ERRORS.append((bam_filename.strip().split('/')[-1],error))
                
            for p in proc:
                try: p.join()
                except: None
        DONE.append(bam_filename.strip().split('/')[-1])

    print(ERRORS)                
else:
    print("The module BUILD_SIMULATED_SEQUENCES was imported.")
    
