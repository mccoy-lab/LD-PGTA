#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

BUILD_SIMULATED_SEQUENCES

Daniel Ariad (daniel@ariad.org)
Dec 30, 2020

"""
import time, pickle, re, sys, random
from random import sample, choices, seed
from multiprocessing import Process
from MIX_HAPLOIDS import MixHaploids_wrapper
from IMPUTE2OBS import main as simulate
from os import remove
import os

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

def transitions(chr_id):
    x = int(chr_id[3:]) if chr_id[3:].isnumeric() else chr_id[3:]
    if type(x) is int and 1<=x<=6:
        #BPH-SPH-BPH-SPH
        result = ('BPH',random.uniform(0,.25),random.uniform(.5,.75),random.uniform(.75,1))
    elif type(x) is int and 7<=x<=12:
        #BPH-SPH-BPH
        result = ('BPH',random.uniform(0,.333),random.uniform(.666,1))
    elif x=='X' or (type(x) is int and 13<=x<=22):
        #SPH-BPH
        result = ('SPH',random.uniform(.5,1))
    else:
        result = ('SPH',1)
    return (result,)
    
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

#def make_simulated_obs_tab(sample_id,sp,chr_id,genotypes,output_dir):
#    """ Based on VCF2OBS """    
#    bcftools_dir = '' #'../bcftools-1.10.2/bin'
#    #sample_id = 'HG00096'
#    #chr_id = 'chr21'
#    leg_filename = f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/{chr_id:s}_{sp:s}_panel.legend'
#    vcf_filename = f'../vcf_phase3_hg38_v2/ALL.{chr_id:s}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz'
#    #output_dir  = f'results_{sp:s}/'
#    return simulate(vcf_filename,leg_filename,chr_id,sample_id,bcftools_dir,genotypes=genotypes,output_dir=output_dir)

def make_simulated_obs_tab(sample_id,sp,chr_id,genotypes,output_dir):
    """ Based on IMPUTE2OBS """
    path = f'../build_reference_panel/{sp:s}_panel.hg38.BCFtools/'
    #path = f'../build_reference_panel/ref_panel.{sp:s}.hg38.BCFtools/'
    leg_filename = path + f'{chr_id:s}_{sp:s}_panel.legend'
    hap_filename = path + f'{chr_id:s}_{sp:s}_panel.hap'
    samp_filename = path + f'{chr_id:s}_{sp:s}_panel.samples'
    return simulate(leg_filename,hap_filename,samp_filename,chr_id,sample_id,genotypes=genotypes,output_dir=output_dir)

def main(depth,sp,chr_id,read_length,min_reads,max_reads):
    ###depth = 0.5
    ###sp = 'EUR'
    ###chr_id = 'chr21'
    work_dir = f'results_{sp:s}/' #'results_EAS' #
    #####################
    seed(None, version=2)
    work_dir = work_dir.rstrip('/') + '/' if len(work_dir)!=0 else ''
    INDIVIDUALS = read_ref(f'../build_reference_panel/{sp:s}_panel.txt') #EAS_panel.txt') 
    A = sample(INDIVIDUALS,k=3)
    B = choices(['A','B'],k=3)
    C = [i+j for i,j in zip(A,B)]
    #for a,b in zip(A,B): make_simulated_obs_tab(a,sp,chr_id,b,work_dir)
    func = (eval(f"lambda: make_simulated_obs_tab('{a:s}', '{sp:s}', '{chr_id:s}', '{b:s}', '{work_dir:s}')") for a,b in zip(A,B))
    runInParallel(*func)
    filenames = MixHaploids_wrapper(f'{work_dir:s}{C[0]:s}.{chr_id:s}.hg38.obs.p', f'{work_dir:s}{C[1]:s}.{chr_id:s}.hg38.obs.p', f'{work_dir:s}{C[2]:s}.{chr_id:s}.hg38.obs.p', read_length=read_length, depth=depth, scenarios=('monosomy','disomy','SPH','BPH','transitions'), transitions=transitions(chr_id),output_dir=work_dir)
    
    for c in C: os.remove(f'./{work_dir:s}{c:s}.{chr_id:s}.hg38.obs.p')
    
    #filenames = ["/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_EUR/simulated.SPH.chr21.x0.010.NA20536B.NA20536A.obs.p",
    #             "/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_EUR/simulated.BPH.chr21.x0.010.NA20536B.NA20536A.NA20802B.rs0.00.obs.p",
    #             "/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_EUR/simulated.monosomy.chr21.x0.010.NA20536B.obs.p",
    #             "/home/ariad/Dropbox/postdoc_JHU/origin_ecosystem/origin_V2/results_EUR/simulated.disomy.chr21.x0.010.NA20536B.NA20536A.obs.p"]
    print(filenames)
    func2 = (eval(f"lambda: aneuploidy_test_demo('{f:s}','{chr_id:s}','{sp:s}','MODELS/MODELS16.p',{min_reads:d},{max_reads:d},'{work_dir:s}')") for f in filenames)
    runInParallel(*func2)

    return 0


if __name__ == "__main__":
    depth=0.01
    sp='AFR'
    chr_id='chr21'
    read_length = 36
    min_reads,max_reads = 6,4
    for n in ([*range(1,23)]+['X'])*20:
        chr_id = 'chr' + str(n)
        runInParallel(*([main]*6),args=(depth,sp,chr_id,read_length,min_reads,max_reads) )
    #main(depth,sp,chr_id,read_length,min_reads,max_reads)
    pass
else:
    print("The module BUILD_SIMULATED_SEQUENCES was imported.")
