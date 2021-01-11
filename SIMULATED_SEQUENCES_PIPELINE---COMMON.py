#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

BUILD_SIMULATED_SEQUENCES

Daniel Ariad (daniel@ariad.org)
Nov 18, 2020

"""
import time, pickle, re, sys
from random import sample, choices, seed
from multiprocessing import Process
from MIX_HAPLOIDS import MixHaploids2
from SIMULATE_OBS_TAB import main as simulate

def read_ref(filename):
    with open(filename, 'r') as data_in:
        tab = tuple(str(line.replace('\n','')) for line in data_in)
    return tab

def runInParallel(*fns,**kwargs):
    proc = []
    for fn in fns:
        try:
            p = Process(target=fn,args=kwargs.get('args',tuple()))
            p.start()
            proc.append(p)
            time.sleep(5)
        except:
            print('caution: a process failed!')
    for p in proc:
        try:
          p.join()
        except:
            None
   
def aneuploidy_test_demo(obs_filename,chr_id,sp,model,min_reads,max_reads,output_dir):
    from ANEUPLOIDY_TEST import aneuploidy_test
    args = dict(obs_filename = f'results_{sp:s}/ABC.obs.p',
                hap_filename = f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.hap',
                leg_filename = f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.legend',
                window_size = 0,
                subsamples = 100,
                offset = 0,
                min_reads = min_reads, #3,
                max_reads = max_reads, #8,
                min_HF = 0.05,
                minimal_score = 2,
                output_dir = output_dir, #f'results_{sp:s}/',
                output_filename = '')
                #model = model)
    args['obs_filename'] = obs_filename #f'results_{sp:s}/' + obs_filename
    LLR_dict, info = aneuploidy_test(**args)
    return LLR_dict, info

def make_simulated_obs_tab(sample_id,sp,chr_id,output_dir):
    bcftools_dir = '' #'../bcftools-1.10.2/bin'
    #sample_id = 'HG00096'
    #chr_id = 'chr21'
    #leg_filename = '../build_reference_panel/ref_panel.ALL.hg38.BCFtools/chr21_ALL_panel.legend'
    #leg_filename = '../build_reference_panel/ref_panel.Ashkenazi.hg38.BCFtools/chr21_Ashkenazi_panel.legend'
    #leg_filename = '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.legend'
    leg_filename = f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.legend'
    #leg_filename = '../build_reference_panel/ref_panel.TEST.hg38.BCFtools/chr21_TEST_panel.legend'
    vcf_filename = f'../vcf_phase3_hg38_v2/ALL.{chr_id:s}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz'
    #output_dir  = f'results_{sp:s}/'
    return simulate(vcf_filename,leg_filename,chr_id,sample_id,bcftools_dir,output_dir=output_dir)

def main(depth,sp,ind,chr_id,read_length,min_reads,max_reads):
    ###depth = 0.5
    ###sp = 'EUR'
    ###chr_id = 'chr21'
    work_dir = f'results_{sp:s}/'
    #####################
    seed(None, version=2)
    work_dir = work_dir.rstrip('/') + '/' if len(work_dir)!=0 else ''
    INDIVIDUALS = read_ref(f'../build_reference_panel/{ind:s}_panel.txt') #{sp:s}_panel.txt')
    A = sample(INDIVIDUALS,k=3)
    B = choices(['A','B'],k=3)
    C = [i+j for i,j in zip(A,B)]
    #for a in A: make_simulated_obs_tab(a,sp,chr_id,work_dir)
    func = (eval(f'lambda: make_simulated_obs_tab(\'{a:s}\', \'{sp:s}\', \'{chr_id:s}\', \'{work_dir:s}\')') for a in A)
    runInParallel(*func)
    MixHaploids2(f'{work_dir:s}{C[0]:s}.chr21.hg38.obs.p', f'{work_dir:s}{C[1]:s}.chr21.hg38.obs.p', f'{work_dir:s}{C[2]:s}.chr21.hg38.obs.p', read_length=read_length, depth=depth, output_dir=work_dir, recombination_spots=[0.00,1.00])
    filenames = (f'mixed3haploids.X{depth:.2f}.{C[0]:s}.{C[1]:s}.{C[2]:s}.chr21.recomb.{i:.2f}.obs.p' for i in (0,1))
    func2 = (eval(f'lambda: aneuploidy_test_demo(\'{work_dir:s}{f:s}\',\'{chr_id:s}\',\'{sp:s}\',\'MODELS/MODELS16D.p\',{min_reads:d},{max_reads:d},\'{work_dir:s}\')') for f in filenames)
    runInParallel(*func2)
    return 0


if __name__ == "__main__":
    sp='COMMON'
    chr_id='chr21'
    read_length = 35
    
    #depth=0.01
    #min_reads,max_reads = 3,4
    #for ind in ('EUR','AFR','AMR','SAS','EAS','ALL'):
    #    print(ind)
    #    runInParallel(*(main for _ in range(14)),args=(depth,sp,ind,chr_id,read_length,min_reads,max_reads))
    #depth=0.02
    #min_reads,max_reads = 3,3
    #for ind in ('EUR','AFR','AMR','SAS','EAS','ALL'):
    #    print(ind)
    #    runInParallel(*(main for _ in range(14)),args=(depth,sp,ind,chr_id,read_length,min_reads,max_reads))
    #depth=0.02
    #min_reads,max_reads = 3,4
    #for ind in ('EUR','AFR','AMR','SAS','EAS','ALL'):
    #    print(ind)
    #    runInParallel(*(main for _ in range(12)),args=(depth,sp,ind,chr_id,read_length,min_reads,max_reads))
    #depth=0.02
    #min_reads,max_reads = 3,5
    #for ind in ('EUR','AFR','AMR','SAS','EAS','ALL'):
    #    print(ind)    
    #    runInParallel(*(main for _ in range(12)),args=(depth,sp,ind,chr_id,read_length,min_reads,max_reads))
    min_reads,max_reads = 3,6
    depth=0.1
    for ind in ('EUR','AFR','AMR','SAS','EAS','ALL'):
        print(ind)
        runInParallel(*(main for _ in range(6)),args=(depth,sp,ind,chr_id,read_length,min_reads,max_reads))
    #min_reads,max_reads = 3,14
    #depth=0.25
    #for ind in ('EUR','AFR','AMR','SAS','EAS','ALL'):
    #    runInParallel(*(main for _ in range(12)),args=(depth,sp,ind,chr_id,read_length,min_reads,max_reads))
    pass
else:
    print("The module BUILD_SIMULATED_SEQUENCES was imported.")
