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
    args = {'bam_filename': '/home/ariad/Dropbox/postdoc_JHU/RPL_merged/'+bam_filename,
            'output_dir': '/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_RPL/',
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
    args = dict(obs_filename = '/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_RPL/' + obs_filename,
                hap_filename = f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.hap',
                leg_filename = f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.legend',
                window_size = 0,
                subsamples = 100,
                offset = 0,
                min_reads = 12,
                max_reads = 8,
                min_HF = 0.05,
                minimal_score = 2,
                output_dir = '/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_RPL/', 
                output_filename = '')
    LLR_dict, info = aneuploidy_test(**args)
    return LLR_dict, info

if __name__ == "__main__":
    
    
    filenames = ['2019-12-23_RIG003-2_S11.bam', '2020-05-22_RPL005-4_S8.bam',
                 '2019-12-23_RIG003-3_S4.bam',   '2020-05-22_RPL006-1_S9.bam',
                 '2019-12-23_RIG003-4_S8.bam',   '2020-05-22_RPL006-2_S10.bam',
                 '2020-05-22_RPL004-1_S1.bam',   '2020-05-22_RPL006-3_S11.bam',
                 '2020-05-22_RPL004-2_S2.bam',   '2020-05-22_RPL006-4_S12.bam',
                 '2020-05-22_RPL004-3_S3.bam',   '2020-05-22_RPL007-1_S13.bam',
                 '2020-05-22_RPL004-4_S4.bam',   '2020-05-22_RPL007-2_S14.bam',
                 '2020-05-22_RPL005-1_S5.bam',   '2020-05-22_RPL007-3_S15.bam',
                 '2020-05-22_RPL005-2_S6.bam',   '2020-05-22_RPL007-4_S16.bam',
                 '2020-05-22_RPL005-3_S7.bam']
 
    DONE = []
    ERRORS = []
    output_dir = '/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_RPL/'
    for bam_filename in filenames:
        if bam_filename not in DONE:
            print(bam_filename)
            sp = 'EUR'
            proc = []
            try:
                for chr_num in range(1,23):
                    chr_id = 'chr'+str(chr_num)
                    obs_filename = re.sub('.bam$','',bam_filename.split('/')[-1]) + f'.{chr_id:s}.obs.p'
                    LLR_filename = re.sub('.bam$','',bam_filename.split('/')[-1]) + f'.{chr_id:s}.LLR.p'
                    
                    if not os.path.isfile(output_dir+obs_filename):   
                        make_obs_tab_demo(bam_filename,chr_id, sp)
                    else:
                        print(f'{obs_filename:s} already exists.')
                    
                    if not os.path.isfile(output_dir+LLR_filename): 
                        #aneuploidy_test_demo(obs_filename, chr_id, sp)
                        p = Process(target=aneuploidy_test_demo,args=(obs_filename, chr_id, sp))
                        p.start()
                        proc.append(p)
                    else:
                        print(f'{LLR_filename:s} already exists.')
                        
                    ###make_obs_tab_demo(bam_filename,chr_id, sp)
                    ####aneuploidy_test_demo(obs_filename, chr_id, sp)
                    ###p = Process(target=aneuploidy_test_demo,args=(obs_filename, chr_id, sp))
                    ###p.start()
                    ###proc.append(p)
            except Exception as error: 
                print('ERROR: ', error)
                if os.path.isfile(output_dir+obs_filename): os.remove(output_dir+obs_filename)
                if os.path.isfile(output_dir+LLR_filename): os.remove(output_dir+LLR_filename)
                ERRORS.append((bam_filename.strip().split('/')[-1],error))
                
            for p in proc:
                try: p.join()
                except: None
        DONE.append(bam_filename)

    print(ERRORS)                
else:
    print("The module BUILD_SIMULATED_SEQUENCES was imported.")
    
