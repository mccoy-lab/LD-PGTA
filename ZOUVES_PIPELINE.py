#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

ZOUVES_PIPELINE

Daniel Ariad (daniel@ariad.org)
Nov 18, 2020

"""
import time, re, pickle

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


if __name__ == "__main__":
    db_BPH = [{'filename': '10469FA-A9DBT_4.bam', 'sp': 'EUR', 'chr_num': (19,) },
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
    
    db_DIPLOID = [  {'filename': '10523FA-AFFRU_3.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10523FA-AFFRU_4.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10560FA-AFFPH_3.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10675FA-BJNTV_2.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10675FA-BJNTV_3.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10675FA-BJNTV_4.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10675FA-BJNTV_5.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10686FA-AFFPE_9.bam', 'sp': 'SAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10846FA-AFPAB_3.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10871FA-AJ3U4_12.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10951FA-AJ3WW_3.bam', 'sp': 'SAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10969FA-AJ470_2.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '11522FA-AP925_3.bam', 'sp': 'SAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '11578FA-AR3WC_3.bam', 'sp': 'AMR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '11598FA-AP923_9.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '12662FA-B5Y5R_1.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '12662FA-B5Y5R_3.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '12699FA-B8F4K_4.bam', 'sp': 'AMR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '12789FA-AWL1L_12.bam', 'sp': 'AFR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '12962FA-BK2G8_6.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '14529FA-CM2GK_2.bam', 'sp': 'AMR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': 'GP-CWFRM_8.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': 'MZ-AFFC4_1.bam', 'sp': 'AMR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '13068FA-BK2G5_23.bam', 'sp': 'AMR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10675FA-BJNTV_3c.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': 'MZ-AFFC4_2.bam', 'sp': 'SAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10964FA-AJ470_1.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '13121FA-BK23M_23.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '13086FA-BK2G5_6.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10668FA-AFDL2_2.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '11550FA-AP91V_5.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10722FA-AFFCT_3a.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '12055FA-ANTJ1_15.bam', 'sp': 'SAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '12454FA-AW7BB_2.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10967FA-AJ470_F10.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '11946FA-AR452_9.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '11550FA-AP91V_4.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '13744FA-C4RPY_2.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '13086FA-BK2G5_8.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '10658FA-AFDL2_F10.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '14220FA-CFP2Y_1.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '12446FA-AU3UC_29.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '14212FA-CFP2Y_5.bam', 'sp': 'AMR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '11946FA-AR452_1.bam', 'sp': 'EAS', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '11944FA-AR452_13.bam', 'sp': 'EUR', 'chr_num': [11,17,19,20,21,22]},
                    {'filename': '11511FA-AP91V_8.bam', 'sp': 'AMR', 'chr_num': [11,17,19,20,21,22]}]
    
    db_AAAAA = [{'filename': '12751FA-AWL31_14.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]}]
    
    with open('/home/ariad/Dropbox/postdoc_JHU/BlueFuse/Play/diploid_females.p', 'rb') as f:
        db_TEST = pickle.load(f)
   
    ###DONE = []
    for case in db_TEST:
        if case not in DONE:
            bam_filename = case['filename']
            print(case['filename'])
            sp = case['sp']
            try:
                for chr_num in case['chr_num']:
                    chr_id = f'chr{chr_num:d}'
                    make_obs_tab_demo(case['filename'],chr_id,sp)
                    obs_filename = re.sub('.bam$','',bam_filename.strip().split('/')[-1]) + f'.{chr_id:s}.obs.p'
                    aneuploidy_test_demo(obs_filename, chr_id,sp)
            except:
                continue
        DONE.append(case)
                
else:
    print("The module BUILD_SIMULATED_SEQUENCES was imported.")
    
