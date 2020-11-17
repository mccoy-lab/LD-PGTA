#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

BUILD_SIMULATED_SEQUENCES

Daniel Ariad (daniel@ariad.org)
Feb 27, 2020

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

def aneuploidy_test_demo(obs_filename,chr_id,sp,model):
    from ANEUPLOIDY_TEST import aneuploidy_test
    args = dict(obs_filename = f'results_{sp:s}/ABC.obs.p',
                hap_filename = f'../build_reference_panel/ref_panel.{sp:s}.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.hap',
                leg_filename = f'../build_reference_panel/ref_panel.{sp:s}.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.legend',
                window_size = 0,
                subsamples = 100,
                offset = 0,
                min_reads = 4,
                max_reads = 18,
                min_HF = 0.15,
                minimal_score = 2,
                output_dir = f'results_{sp:s}/',
                output_filename = '',
                model = model)
    args['obs_filename'] = f'results_{sp:s}/' + obs_filename
    LLR_dict, info = aneuploidy_test(**args)
    return LLR_dict, info

def make_simulated_obs_tab(sample_id,sp):
    bcftools_dir = ''
    #sample_id = 'HG00096'
    chr_id = 'chr21'
    #leg_filename = '../build_reference_panel/ref_panel.ALL.hg38.BCFtools/chr21_ALL_panel.legend'
    #leg_filename = '../build_reference_panel/ref_panel.Ashkenazi.hg38.BCFtools/chr21_Ashkenazi_panel.legend'
    #leg_filename = '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.legend'
    leg_filename = f'../build_reference_panel/ref_panel.{sp:s}.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.legend'
    #leg_filename = '../build_reference_panel/ref_panel.TEST.hg38.BCFtools/chr21_TEST_panel.legend'
    vcf_filename = f'../vcf_phase3_hg38_v2/ALL.{chr_id:s}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz'
    work_dir = f'results_{sp:s}'
    return simulate(vcf_filename,leg_filename,chr_id,sample_id,bcftools_dir,output_dir=work_dir)

def main():
    depth = 0.01
    seed(None, version=2)
    INDIVIDUALS = read_ref('/home/ariad/Dropbox/postdoc_JHU/Tools/build_reference_panel/EUR_panel.txt')
    for r in range(10):
        A = sample(INDIVIDUALS,k=3)
        B = choices(['A','B'],k=3)
        C = [i+j for i,j in zip(A,B)]
        func = [eval(f'lambda: make_simulated_obs_tab(\'{a:s}\',\'EUR\')') for a in A]
        runInParallel(*func)
        D = MixHaploids2(f'{C[0]:s}.chr21.hg38.obs.p', f'{C[1]:s}.chr21.hg38.obs.p', f'{C[2]:s}.chr21.hg38.obs.p', read_length=35, depth=depth, work_dir='results_EUR', recombination_spots=[0.00,1.00])
        filenames = (f'mixed3haploids.X{depth:.2f}.{C[0]:s}.{C[1]:s}.{C[2]:s}.chr21.recomb.{i:.2f}.obs.p' for i in (0,1))
        func2 = (eval('lambda: aneuploidy_test_demo(\'%s\',\'%s\',\'EUR\',\'MODELS/MODELS16D.p\')' % (f,'chr21')) for f in filenames)
        runInParallel(*func2)
    return 0


if __name__ == "__main__":
    runInParallel(*[main for _ in range(4)])