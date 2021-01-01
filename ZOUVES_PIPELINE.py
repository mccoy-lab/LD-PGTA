#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

ZOUVES_PIPELINE

Daniel Ariad (daniel@ariad.org)
Nov 18, 2020

"""
import time, re, pickle, os

def make_obs_tab_demo(bam_filename,chr_id,sp):
    from MAKE_OBS_TAB import retrive_bases
    args = {'bam_filename': '/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/'+bam_filename,
            'output_dir': '/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/',
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
    args = dict(obs_filename = '/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/' + obs_filename,
                hap_filename = f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.hap',
                leg_filename = f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.legend',
                window_size = 0,
                subsamples = 100,
                offset = 0,
                min_reads = 10,
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
    
    db_DIPLOID = [  {'filename': '10523FA-AFFRU_3.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10523FA-AFFRU_4.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10560FA-AFFPH_3.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10675FA-BJNTV_2.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10675FA-BJNTV_3.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10675FA-BJNTV_4.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10675FA-BJNTV_5.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10686FA-AFFPE_9.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10846FA-AFPAB_3.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10871FA-AJ3U4_12.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10951FA-AJ3WW_3.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10969FA-AJ470_2.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '11522FA-AP925_3.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '11578FA-AR3WC_3.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '11598FA-AP923_9.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '12662FA-B5Y5R_1.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '12662FA-B5Y5R_3.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '12699FA-B8F4K_4.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '12789FA-AWL1L_12.bam', 'sp': 'AFR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '12962FA-BK2G8_6.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '14529FA-CM2GK_2.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': 'GP-CWFRM_8.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': 'MZ-AFFC4_1.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '13068FA-BK2G5_23.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10675FA-BJNTV_3c.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': 'MZ-AFFC4_2.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10964FA-AJ470_1.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '13121FA-BK23M_23.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '13086FA-BK2G5_6.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10668FA-AFDL2_2.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '11550FA-AP91V_5.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10722FA-AFFCT_3a.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '12055FA-ANTJ1_15.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '12454FA-AW7BB_2.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10967FA-AJ470_F10.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '11946FA-AR452_9.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '11550FA-AP91V_4.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '13744FA-C4RPY_2.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '13086FA-BK2G5_8.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '10658FA-AFDL2_F10.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '14220FA-CFP2Y_1.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '12446FA-AU3UC_29.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '14212FA-CFP2Y_5.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '11946FA-AR452_1.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '11944FA-AR452_13.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                    {'filename': '11511FA-AP91V_8.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}]
    
    db_DIPLOID_noisy =  [{'filename': '10469FA-A9DBT_5.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '10560FA-AFFPH_5.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '11095FA-AK9FG_5.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '11422FA-APA42_4.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '11517FA-AP925_2.bam', 'sp': 'SAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '11594FA-AP923_15.bam', 'sp': 'AMR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '12853FA-AWKB8_1.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '12862FA-AWKB8_1.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '14529FA-CM2GK_3.bam', 'sp': 'AMR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '14529FA-CM2GK_4.bam', 'sp': 'AMR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '14953FA-CWLHB_6.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '14974FA-CWFV7_13.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '14974FA-CWFV7_14.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '14992FA-CWLYV_1.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '15022FA-J2MJB_11.bam', 'sp': 'AMR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '15042FA-J2JJH_4.bam', 'sp': 'AMR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': 'GP-CWFG9_1.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': 'GP-CWGJL_3.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': 'GP-CWGJL_6.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '11946FA-AR452_7.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': 'RES-11892-AU3TL_5T.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '12446FA-AU3UC_32.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '10722FA-AFFCT_5.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '13068FA-BK2G5_27.bam', 'sp': 'AMR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '11944FA-AR452_15.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '11557FA-AP2W6_4.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '12446FA-AU3UC_19.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '11557FA-AP2W6_1.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)}, 
                         {'filename': '11143FA-AL3GC_27.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '11494FA-AP9FW_6.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '14217FA-CFP2Y_3.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '11180FA-AK9G0_8.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '13126FA-B6RLR_12944FA-1.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '13948FA-C7R7Y_23.bam', 'sp': 'EUR', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '12881FA-AWK94_5.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)},
                         {'filename': '12446FA-AU3UC_28.bam', 'sp': 'EAS', 'chr_num': (1, 2, 3, 11, 17, 21)}]
    
    db_HAPLOIDS = [{'filename': '12751FA-AWL31_14.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12459FA-AU3UR_8.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '14338FA-CFHJ8_14.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '14715FA-CTGFC_21.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '14080FA-CBRH3_24.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '14260FA-CFPGD_19.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '14451FA-CFN45_14.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '12751FA-B5Y5C_14.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '11968FA-AP1KV_6.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '13885FA-C6JPJ_10.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12459FA-ARHR1_1.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13402FA-BK37C_22.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '12459FA-AU3UM_12.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '14417FA-CFMY3_11.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12456FA-ARHR1_16.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13474FA-BRCNL_11.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '13335FA-BK7V7_16.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '13198FA-BK2FN_12.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '11143FA-AK9EU_26.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '13327FA-BK7V7_21.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': '12150FA-AW823_LB-8-AMP-17.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13213FA-BK2FN_20.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13213FA-BK2FN_2.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13297FA-B6TD2_12.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '11946FA-AR48V_9.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   {'filename': 'RES-11733-AW727_LB-8N-AMP16.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13578FA-BRD32_12.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '11946FA-AR48V_3.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13597FA-C476J_35.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13659FA-BJKF9_10.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12073FA-ANTJ1_4.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12394FA-AU3TP_2.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '13297FA-B6TD2_17.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12459FA-AU3UM_8.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '11257FA-ANDVJ_11.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '11826FA-AR21R_28.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12155FA-AW829_LB-7-AMP17.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12149FA-AW823_LB-16-AMP17.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']},
                   #{'filename': '12150FA-AW823_LB-8-AMP-17.bam', 'sp': 'EUR','chr_num': [*range(1,23)]+['X']}
                   ]
    
    
    db_Ydisomy = [{'filename': '12663FA-B8DPD_11.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']},
                  {'filename': '13068FA-BK2G5_23.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']},
                  {'filename': '12405FA-AU3TR_3.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}]
    
    db_TRIPLOIDS_maybe = [{'filename': '13331FA-BK8NW_2.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14096FA-CBRH3_6.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11749FA-AR3WC_3.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14088FA-C7T6M_15.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12431FA-B22P2_14.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12279FA-ATVDW_6.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13809FA-C6JBM_9.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14119FA-CB9JP_2.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12149FA-AW829_11-AMP16.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12149FA-AW738_4.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14149FA-CB9JP_14.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11494FA-AP9FW_2.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11624FA-ANTHL_6.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13718FA-C4T43_8.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14277FA-CFKCF_20.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14767FA-CTFFN_17.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11270-RES-C7RGH_1T1.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13528FA-BRCJ2_7.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10672FA-AFFPR_7.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13642FA-BJKDV_4.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14473FA-CHTK6_25.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10405FA-AAAM9_2.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12610FA-AV2K5_4.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10667FA-AFFCM_2.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12682FA-B8FT9_2.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11288FA-ALNLD_9.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13700FA-C4T36_9.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11627FA-ANTHL_6.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13809FA-C6JBM_10.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13515FA-BRCJ2_7.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13534FA-BRCJ2_7.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11232FA-AK90E_5.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12279FA-AY8J4_6.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12716FA-AYW5E_12.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12103FA-AW60G_9.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11171FA-AK9FM_19.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14236FA-CFJKV_1.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14328FA-CFHJ8_14.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13564FA-BRCN9_11.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13098FA-B6MTL_3.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12742FA-B8FHD_5.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11510FA-APA12_5.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12566FA-AV2NV_1.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13175FA-BJNV7_3.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13175FA-B6RJM_3.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13762FA-C4T43_3.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10864FA-AFPAN_1.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14444FA-CFN45_17.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10837FA-AFPAJ_1.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14731FA-CTGFC_1.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12103FA-AW60G_5.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12279FA-ATVDW_18.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10771FA-AFPNC_3.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11270FA-ALNKA_21.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11947FA-AR45Y_8.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11070FA-AJY20_3.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11573FA-AP9A4_4.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13617FA-C22MH_9.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14570FA-CLPKK_6-10925.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11826FA-AR21R_25.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13288FA-BJNVW_6.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14473FA-CHTK6_26.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12568FA-AV2NV_2.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11021FA-AJK7A_22.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10635FA-AFDK3_13.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12289FA-AW60A_26.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12279FA-ATVDW_7.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10671FA-AFFPE_24.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12258FA-AU3DR_9.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11270FA-ALNKA_11.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13692FA-C4VK7_15.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11261FA-AL2TY_10.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10417FA-AAALW_7.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11913FA-APHAJ_14.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10749FA-AFP9W_14.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10940FA-AJ47P_5.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12831FA-AWK85_10.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11279FA-ALNLD_23.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13845FA-C6JP9_10.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '14149FA-CB9JP_4.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11541FA-AP99W_1.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11657FA-AP1KJ_4.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13699FA-C4T36_1.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12073FA-AW12E_7.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12946FA-AV39W_12.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11059FA-AJY20_27.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13441FA-BRMJT_20.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13885FA-C6JPJ_8.bam', 'sp': 'SAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13143FA-BJNV7_32.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13025FA-BJNV6_4.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13062FA-BJNTW_4.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12258FA-AU3DR_7.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12478FA-ARHL5_3.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12896FA-AWK94_1.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10882FA-AJ3U4_3.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '12478FA-ARHL5_5.bam', 'sp': 'EAS', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13510FA-C4BW2_2.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '13334FA-BK8HK_3.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '10851FA-AFP9Y_28.bam', 'sp': 'EUR', 'chr_num': [*range(1,23)]+['X']}, {'filename': '11947FA-AP1GW_10.bam', 'sp': 'AMR', 'chr_num': [*range(1,23)]+['X']}]
    
    with open('/home/ariad/Dropbox/postdoc_JHU/BlueFuse/Play/diploid_females.p', 'rb') as f:
        db_TEST = pickle.load(f)
   
    DONE = []
    ERRORS = []
    output_dir = '/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_ZOUVES/'
    for case in db_TEST[::-1]:
        if case not in DONE:
            bam_filename = case['filename']
            print(case['filename'])
            sp = case['sp']
            try:
                for chr_num in case['chr_num']:
                    chr_id = 'chr'+str(chr_num)
                    obs_filename = re.sub('.bam$','',bam_filename.split('/')[-1]) + f'.{chr_id:s}.obs.p'
                    LLR_filename = re.sub('.bam$','',bam_filename.split('/')[-1]) + f'.{chr_id:s}.LLR.p'
                    
                    if not os.path.isfile(output_dir+obs_filename): 
                        make_obs_tab_demo(case['filename'],chr_id, sp)
                    else:
                        print(f'{obs_filename:s} already exists.')
                    
                    if not os.path.isfile(output_dir+LLR_filename): 
                        aneuploidy_test_demo(obs_filename, chr_id, sp)
                    else:
                        print(f'{LLR_filename:s} already exists.')
            except Exception as error: 
                print('ERROR: ', error)
                if os.path.isfile(output_dir+obs_filename): os.remove(output_dir+obs_filename)
                if os.path.isfile(output_dir+LLR_filename): os.remove(output_dir+LLR_filename)
                ERRORS.append((bam_filename.strip().split('/')[-1],error))
                continue
        DONE.append(bam_filename.strip().split('/')[-1])

    print(ERRORS)                
else:
    print("The module BUILD_SIMULATED_SEQUENCES was imported.")
    
