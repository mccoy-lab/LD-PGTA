#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

ZOUVES_PIPELINE

Daniel Ariad (daniel@ariad.org)
Nov 18, 2020

"""
import time, re

def make_obs_tab_demo(bam_filename,chr_id,sp):
    from MAKE_OBS_TAB import retrive_bases
    args = {'bam_filename': '/home/ariad/Dropbox/postdoc_JHU/Tools/origin_V2/results_ZOUVES/'+bam_filename,
            'output_dir': '/home/ariad/Dropbox/postdoc_JHU/Tools/origin_V2/results_ZOUVES/',
            'legend_filename': f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.legend',
            'max_depth': 0,
            'min_bq': 30,
            'min_mq': 30,
            'handle_multiple_observations': 'all',
            'fasta_filename': '',#'../genome_ref_hg38/hg38.fa', 
            'output_filename': ''}
    
    result = retrive_bases(**args)
    
    return result
    
def aneuploidy_test_demo(obs_filename,chr_id,sp):
    from ANEUPLOIDY_TEST import aneuploidy_test
    args = dict(obs_filename = '/home/ariad/Dropbox/postdoc_JHU/Tools/origin_V2/results_ZOUVES/' + obs_filename,
                hap_filename = f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.hap',
                leg_filename = f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.legend',
                window_size = 0,
                subsamples = 100,
                offset = 0,
                min_reads = 3,
                max_reads = 8,
                min_HF = 0.05,
                minimal_score = 2,
                output_dir = '/home/ariad/Dropbox/postdoc_JHU/Tools/origin_V2/results_ZOUVES/', 
                output_filename = '')
    LLR_dict, info = aneuploidy_test(**args)
    return LLR_dict, info


if __name__ == "__main__":
    db = [{'filename': '10469FA-A9DBT_4.bam', 'sp': 'EUR', 'chr_num': (19,) },
          {'filename': '10534FA-AFFRU_23.bam', 'sp': 'AMR', 'chr_num': (19,) },
          {'filename': '10568FA-AFFAN_4.bam', 'sp': 'AMR', 'chr_num': (7,20) },
          {'filename': '10846FA-AFPAB_2.bam', 'sp': 'EAS', 'chr_num': (22,) },
          {'filename': '10969FA-AJ470_2.bam', 'sp': 'EAS', 'chr_num': (16,) },
          {'filename': '12699FA-B8F4K_4.bam', 'sp': 'AMR', 'chr_num': (16,) },
          {'filename': '12853FA-AWKB8_1.bam', 'sp': 'EAS', 'chr_num': (16,) },
          {'filename': '13782FA-C6KKV_20.bam', 'sp': 'SAS', 'chr_num': (13,) },
          {'filename': '14286FA-CFHN7_14.bam', 'sp': 'EAS', 'chr_num': (2,3,17,20) },
          {'filename': '14799FA-CRGVJ_8.bam', 'sp': 'EUR', 'chr_num': (8,9,18) },
          {'filename': '14994FA-CWLYV_2.bam', 'sp': 'EUR', 'chr_num': (16,22) },
          {'filename': '14214FA-CFP2Y_4.bam', 'sp': 'EUR', 'chr_num': (16,22) },
          {'filename': '14220FA-CFP2Y_1.bam', 'sp': 'EAS', 'chr_num': (19,) },
          {'filename': 'RES-11878-AU3TL_2I.bam', 'sp': 'SAS', 'chr_num': (22,) },
          {'filename': '20150608-22-AFDL2_NegC.bam', 'sp': 'EUR', 'chr_num': (15,) },
          {'filename': '10967FA-AJ470_F6.bam', 'sp': 'EAS', 'chr_num': (9,11,22) },
          {'filename': '13948FA-C7R7Y_30.bam', 'sp': 'EUR', 'chr_num': (5,) },
          {'filename': '13935FA-C7R7Y_35.bam', 'sp': 'EUR', 'chr_num': (22,) },
          {'filename': '12880FA-AWK94_3.bam', 'sp': 'EUR', 'chr_num': (13,) },
          {'filename': '11433FA-ALNVD_3.bam', 'sp': 'SAS', 'chr_num': (5,10,16) },
          {'filename': '11952FA-AP1KV_1.bam', 'sp': 'EUR', 'chr_num': (18,) },
          {'filename': '12940FA-BK3D2_13.bam', 'sp': 'EUR', 'chr_num': (22,) },
          {'filename': '13935FA-C7R7Y_28.bam', 'sp': 'EUR', 'chr_num': (16,) },
          {'filename': '14214FA-CFP2Y_3.bam', 'sp': 'EUR', 'chr_num': (16,) }]
    
    
    for case in db:
        bam_filename = case['filename']
        sp = case['sp']
        for chr_num in case['chr_num']:
            chr_id = f'chr{chr_num:d}'
            #make_obs_tab_demo(case['filename'],chr_id,sp)
            obs_filename = re.sub('.bam$','',bam_filename.strip().split('/')[-1]) + f'.{chr_id:s}.obs.p'
            aneuploidy_test_demo(obs_filename, chr_id,sp)  
else:
    print("The module BUILD_SIMULATED_SEQUENCES was imported.")