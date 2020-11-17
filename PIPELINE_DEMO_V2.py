#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

PIPELINE

Daniel Ariad (daniel@ariad.org)
Feb 27, 2020

"""
import time, pickle, re, sys
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
                window_size = 10e4,
                subsamples = 100,
                offset = 0,
                min_reads = 6,
                max_reads = 14,
                min_HF = 0.15,
                minimal_score = 2,
                output_filename = None,
                model = model)
    print(model)
    #M = model.replace('.', '/').split('/')[-2]
    args['obs_filename'] = f'results_{sp:s}/' + obs_filename
    args['output_filename'] = f'results_{sp:s}/'+re.sub('(.*)obs','\\1LLR123', obs_filename.split('/')[-1],1)    
    LLR_dict, info = aneuploidy_test(**args)
    return LLR_dict, info
   
def make_obs_tab_demo(bam_filename,legend_filename,handle):
    from MAKE_OBS_TAB import retrive_bases
    args = {'bam_filename': '../BAMs_hg38/'+bam_filename,
            'legend_filename': '../build_reference_panel/ref_panel.EUR.hg38.BCFtools/'+legend_filename,
            'max_depth': 0,
            'min_bq': 30,
            'min_mq': 30,
            'handle_multiple_observations': handle,
            'fasta_filename': '',#'../genome_ref_hg38/hg38.fa', 
            'output_filename': ''}
    
    args['output_filename'] = 'results_EUR/'+re.sub('.bam$','',args['bam_filename'].strip().split('/')[-1])+'.obs.p'

    result = retrive_bases(**args)
    
    return result

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
        
if __name__ == "__main__":
    from random import sample, choices, seed
    seed(None, version=2)
    INDIVIDUALS = read_ref('/home/ariad/Dropbox/postdoc_JHU/Tools/build_reference_panel/EUR_panel.txt')
    for r in range(10):
        A = sample(INDIVIDUALS,k=3)
        B = choices(['A','B'],k=3) 
        C = [i+j for i,j in zip(A,B)]    
        func = [eval(f'lambda: make_simulated_obs_tab(\'{a:s}\',\'EUR\')') for a in A] 
        runInParallel(*func)    
        D = MixHaploids2(f'{C[0]:s}.chr21.hg38.obs.p', f'{C[1]:s}.chr21.hg38.obs.p', f'{C[2]:s}.chr21.hg38.obs.p', read_length=35, depth=0.50, work_dir='results_EUR', recombination_spots=[0.00,1.00])        
        filenames = (f'mixed3haploids.X0.50.{C[0]:s}.{C[1]:s}.{C[2]:s}.chr21.recomb.{i:.2f}.obs.p' for i in (0,1))
        func2 = (eval('lambda: aneuploidy_test_demo(\'%s\',\'%s\',\'EUR\',\'MODELS/MODELS16D.p\')' % (f,'chr21')) for f in filenames)
        runInParallel(*func2)
    
        
    
    
    #make_simulated_obs_tab('HG00096','COMMON')
    #make_simulated_obs_tab('HG00097','COMMON')
    #make_simulated_obs_tab('HG00357','COMMON')
    #make_simulated_obs_tab('NA20524','COMMON')
    
    #A = MixHaploids2('HG00096A.chr13.hg38.obs.p', 'HG00096B.chr13.hg38.obs.p', 'HG00097A.chr13.hg38.obs.p', read_length=150, depth=0.20, work_dir='results_EUR', recombination_spots=[i/10 for i in range(11)])
    #A = MixHaploids2('HG00096A.chr6.hg38.obs.p', 'HG00096B.chr6.hg38.obs.p', 'HG00097A.chr6.hg38.obs.p', read_length=36, depth=0.01, work_dir='results_EUR', recombination_spots=[i/10 for i in range(11)])
    #A = MixHaploids2('HG00096A.chr21.hg38.obs.p', 'HG00096B.chr21.hg38.obs.p', 'HG00097A.chr21.hg38.obs.p', read_length=150, depth=0.05, work_dir='results_EUR', recombination_spots=[i/10 for i in range(11)])
    
    #A = MixHaploids2('HG00096A.chr21.hg38.obs.p', 'HG00357B.chr21.hg38.obs.p', 'NA20524A.chr21.hg38.obs.p', read_length=150, depth=0.01, work_dir='results_EUR', recombination_spots=[i/10 for i in range(11)])
    #A = MixHaploids2('HG00096A.chr21.hg38.obs.p', 'HG00357B.chr21.hg38.obs.p', 'NA20524A.chr21.hg38.obs.p', read_length=150, depth=0.05, work_dir='results_EUR', recombination_spots=[i/10 for i in range(11)])

    #A = MixHaploids2('HG00096A.chr21.hg38.obs.p', 'HG00357B.chr21.hg38.obs.p', 'NA20524A.chr21.hg38.obs.p', read_length=35, depth=0.01, work_dir='results_EUR', recombination_spots=[i/10 for i in range(11)])        
    #A = MixHaploids2('HG00096A.chr21.hg38.obs.p', 'HG00357B.chr21.hg38.obs.p', 'NA20524A.chr21.hg38.obs.p', read_length=35, depth=0.05, work_dir='results_EUR', recombination_spots=[i/10 for i in range(11)])        
    #A = MixHaploids2('HG00096A.chr21.hg38.obs.p', 'HG00357B.chr21.hg38.obs.p', 'NA20524A.chr21.hg38.obs.p', read_length=35, depth=0.10, work_dir='results_EUR', recombination_spots=[i/10 for i in range(11)])
    #A = MixHaploids2('HG00096A.chr21.hg38.obs.p', 'HG00357B.chr21.hg38.obs.p', 'NA20524A.chr21.hg38.obs.p', read_length=35, depth=0.50, work_dir='results_EUR', recombination_spots=[i/10 for i in range(11)])
    
    

    #A = lambda: aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.0.00.obs.p','chr21','EUR','MODELS/MODELS16D.p') 
    #B = lambda: aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.1.00.obs.p','chr21','EUR','MODELS/MODELS16D.p')             
    #C = lambda: aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00357B.NA20524A.chr21.recomb.0.00.obs.p','chr21','EUR','MODELS/MODELS16D.p') 
    #D = lambda: aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00357B.NA20524A.chr21.recomb.1.00.obs.p','chr21','EUR','MODELS/MODELS16D.p')             
  
    #runInParallel(A,B,C,D)
    
    #A = lambda: aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00357B.NA20524A.chr21.recomb.0.00.obs.p','chr21','MODELS_CLASSIC.p') 
    #B = lambda: aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00357B.NA20524A.chr21.recomb.1.00.obs.p','chr21','MODELS_CLASSIC.p')      
    #A = lambda: aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.0.00.obs.p','chr21','MODELS_6.p') 
    #B = lambda: aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.1.00.obs.p','chr21','MODELS_6.p')     
    #A = lambda: aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.0.00.obs.p','chr21','MODELS_5.p') 
    #B = lambda: aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.1.00.obs.p','chr21','MODELS_5.p') 
    #aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.0.00.obs.p','chr21','MODELS_4.p') 
    #aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.1.00.obs.p','chr21','MODELS_4.p') 
    #aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.0.00.obs.p','chr21','MODELS_3.p') 
    #aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.1.00.obs.p','chr21','MODELS_3.p') 
    #aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.0.00.obs.p','chr21','MODELS_2.p') 
    #aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.1.00.obs.p','chr21','MODELS_2.p') 
    #aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.0.00.obs.p','chr21','MODELS_1.p') 
    #aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.1.00.obs.p','chr21','MODELS_1.p') 
    #C = lambda: aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.0.00.obs.p','chr21','MODELS_CLASSIC.p')
    #D = lambda: aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.1.00.obs.p','chr21','MODELS_CLASSIC.p')
    #runInParallel(A,B)
    
    #aneuploidy_test_demo('HG00097A.chr21.hg38.obs.p','chr21','MODELS_CLASSIC.p')
    

    #for i in range(0,11):
    #   print('Recombination spot: %.2f' % (i * 0.1))
    #    LLR_dict, info = aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.recomb.%.2f.obs.p' % (i * 0.1))       

    
    #filenames = ('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.%.2f.obs.p' % (i * 0.1) for i in range(11))
    #functions = (eval('lambda: aneuploidy_test_demo(\'%s\',\'%s\')' % (f,'chr21')) for f in filenames)
    #runInParallel(*functions)
    
    #filenames = ('mixed3haploids.X0.10.HG00096A.HG00096B.HG00097A.chr21.recomb.%.2f.obs.p' % (i * 0.1) for i in range(11))
    #functions = (eval('lambda: aneuploidy_test_demo(\'%s\',\'%s\')' % (f,'chr21')) for f in filenames)
    #runInParallel(*functions)
    
    #filenames = ('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.%.2f.obs.p' % (i * 0.1) for i in range(11))
    #functions = (eval('lambda: aneuploidy_test_demo(\'%s\',\'%s\')' % (f,'chr21')) for f in filenames)
    #runInParallel(*functions)
        
    #filenames = ('mixed3haploids.X0.01.HG00096A.HG00096B.HG00097A.chr6.recomb.%.2f.obs.p' % (i * 0.1) for i in range(5,11))
    #functions = (eval('lambda: aneuploidy_test_demo2(\'%s\',\'%s\')' % (f,'chr6')) for f in filenames)
    #runInParallel(*functions)
    
    #LLR_dict, info = aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.0.00.obs.p','chr21')
    #LLR_dict, info = aneuploidy_test_demo('mixed3haploids.X0.50.HG00096A.HG00096B.HG00097A.chr21.recomb.1.00.obs.p','chr21')
    
    
    #A = aneuploidy_test_demo2('11909FA_2.merged.hg38.obs.p','chr6')
    #A = aneuploidy_test_demo2('11694FA_3.merged.hg38.obs.p','chr13')
    
    
    #A = aneuploidy_test_demo2('11909FA_2_T1.chr6.hg38.obs.p','chr6')
    #A = aneuploidy_test_demo2('11909FA_2_T2.chr6.hg38.obs.p','chr6')
    #A = aneuploidy_test_demo2('11909FA_2_I1.chr6.hg38.obs.p','chr6')
    #A = aneuploidy_test_demo2('11909FA_2_I2.chr6.hg38.obs.p','chr6')
    #A = aneuploidy_test_demo2('11694FA_3_T1.chr13.hg38.obs.p','chr13')
    #A = aneuploidy_test_demo2('11694FA_3_T2.chr13.hg38.obs.p','chr13')
    #A = aneuploidy_test_demo2('11694FA_3_I1.chr13.hg38.obs.p','chr13')

    #sys.exit(0)
    #make_obs_tab_demo('SRR10965088.hg38.bam','chr21_EUR_panel.legend')
    
    #make_obs_tab_demo2('chr6/11909FA_2.merged.hg38.bam','chr6_EUR_panel.legend','all')
    #make_obs_tab_demo2('chr13/11694FA_3.merged.hg38.bam','chr13_EUR_panel.legend','all')
    
    #make_obs_tab_demo('SRR6676163.hg38.bam','chr21_MixHap_panel.legend','all')
    #make_obs_tab_demo('SRR151495.0-2.hg38.bam','chr21_HapMix_panel.legend','random')
    #make_obs_tab_demo('SRR10393062.hg38.bam','chr21_HapMix_panel.legend','random')
    #make_obs_tab_demo('SRR8257209.hg38.bam','chr21_HapMix_panel.legend','random')
    #A = make_obs_tab_demo('HG00096.0.5.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG00096.HG00097.0.5.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG00096.0.4.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG00096.HG00097.0.4.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG00096.0.3.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG00096.HG00097.0.3.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG00096.0.2.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG00096.HG00097.0.2.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG00096.0.1.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG00096.HG00097.0.1.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG01880.0.5.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG01880.HG01885.0.5.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG01880.0.4.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG01880.HG01885.0.4.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG01880.0.3.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG01880.HG01885.0.3.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG01880.0.2.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG01880.HG01885.0.2.BPH.bam','chr21_EUR_panel.legend','all')
    #A = make_obs_tab_demo('HG01880.0.1.SPH.bam','chr21_EUR_panel.legend','all')
    #B = make_obs_tab_demo('HG01880.HG01885.0.1.BPH.bam','chr21_EUR_panel.legend','all')

    
    #A = MixHaploids('HG00096.A.hg38.obs.p', 'HG00096.B.hg38.obs.p', 'HG00097.A.hg38.obs.p', read_length=150, depth=0.5, work_dir='results_EUR', recombination_spot=0)
    #B = MixHaploids('HG00096.A.hg38.obs.p', 'HG00097.A.hg38.obs.p', read_length=150, depth=0.5, work_dir='results_EUR', recombination_spot=0)
   
    #for i in range(11):
    #    A = MixHaploids2('HG00096A.hg38.obs.p', 'HG00096B.hg38.obs.p', 'HG00097A.hg38.obs.p', read_length=36, depth=0.01, work_dir='results_EUR', recombination_spot=i/10)
        #B = MixHaploids2('HG00096A.hg38.obs.p', 'HG00097A.hg38.obs.p', read_length=36, depth=0.01, work_dir='results_EUR', recombination_spot=i/10)
    
    #LLR_dict, info = aneuploidy_test_demo('mixed2haploids.X0.50.HG00096A.HG00097A.recomb.0.90.obs.p')
    #LLR_dict, info = aneuploidy_test_demo('mixed2haploids.X0.50.HG00096A.HG00097A.recomb.0.00.obs.p')

    
    
    pass
    
    