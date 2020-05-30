#!/usr/local/python/intelpython3/bin/python
# -*- coding: utf-8 -*-
"""
MAKE_OBS_TAB

This script extracts single base observations at SNP positions from a given
sequence. It requires an aligned and sorted BAM file with the sequence, as well
as an IMPUTE2 legend format, which contains the SNPs positions.
The observed alleles, together with their chromosome position and line number
in the legend file, are organized in a table.

Daniel Ariad (daniel@ariad.org)
May 11st, 2020
"""
import sys, os, time, random, warnings, argparse, re, pickle, pysam

warnings.formatwarning = lambda message, category, filename, lineno, file=None, line=None: 'Caution: %s\n' % message

def read_impute2(impute2_filename,**kwargs):
    """ Reads an IMPUTE2 file format (SAMPLE/LEGEND/HAPLOTYPE) and builds a list
        of lists, containing the dataset. """
    
    filetype = kwargs.get('filetype', None)
    with open(impute2_filename, 'r') as impute2_in:
        if filetype == 'legend':
            impute2_in.readline()   # Bite off the header
            def parse(x): 
                y=x.strip().split()
                y[0] = 'chr'+y[0].split(':')[0]
                y[1]=int(y[1])
                return y
        elif filetype == 'hap':
            def parse(x):
                return [i=='1' for i in x.strip().split()]
        else:
            def parse(x):
                return x.strip().split() # Trailing whitespaces stripped from the ends of the string. Then it splits the string into a list of words.
       
        impute2_tab = [parse(line) for line in impute2_in]
    return impute2_tab

time0 = time.time()

def retrive_bases(bam_filename,legend_filename,fasta_filename,handle_multiple_observations,min_bq,min_mq,max_depth,**kwargs):
    """ Retrives observed bases from known SNPs position. """
    time0 = time.time()

    if not os.path.isfile(bam_filename): raise Exception('Error: BAM file does not exist.')
    if not os.path.isfile(legend_filename): raise Exception('Error: LEGEND file does not exist.')
    if fasta_filename!='' and not os.path.isfile(fasta_filename): raise Exception('Error: FASTA file does not exist.')

    obs_tab = list()
    
    try:
        genome_reference = pysam.FastaFile(fasta_filename) if fasta_filename!='' else None
        samfile = pysam.AlignmentFile(bam_filename, 'rb' )
        leg_tab = read_impute2(legend_filename, filetype='legend')

        if next(zip(*leg_tab)).count(leg_tab[0][0])!=len(leg_tab):
            raise Exception('Error: Unsuitable legend file. All SNP positions should refer to the same chr_id.')
        else:
            chr_id = leg_tab[0][0]
        
        arg = {'contig': chr_id,              # The chr_id of the considered chromosome.
               'start': leg_tab[0][1]-1,      # The first snp in chr_id.
               'end': leg_tab[-1][1],         # The last snp in chr_id.
               'truncate': True,              # By default, the samtools pileup engine outputs all reads overlapping a region. If truncate is True and a region is given, only columns in the exact region specificied are returned.
               'max_depth': max_depth,        # Maximum read depth permitted. The default limit is ‘8000’.
               'stepper': 'samtools',         # The following arguments all pertain to the samtools stepper:
               'min_base_quality': min_bq,    # Minimum base quality. Bases below the minimum quality will not be output.
               'min_mapping_quality': min_mq, # Only use reads above a minimum mapping quality. The default is 0.
               'ignore_overlaps': True,       # If set to True, detect if read pairs overlap and only take the higher quality base.
               'ignore_orphans': True,        # Ignore orphans (paired reads that are not in a proper pair).
               'fastafile': genome_reference, # FastaFile object of a reference sequence.    
               'compute_baq': True}           # By default, performs re-alignment computing per-Base Alignment Qualities (BAQ), if a reference sequence is given.'
                                 
        leg_tab_iterator = enumerate(leg_tab)        
        pos = 0         
        
        for pileupcolumn in samfile.pileup(**arg):
            
            while pileupcolumn.pos > pos-1: 
                impute2_index, (chr_id,pos,ref,alt) = next(leg_tab_iterator)  
            
            if pileupcolumn.pos == pos-1:
                
                rows = [(pos, impute2_index, pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]) for pileupread in pileupcolumn.pileups if pileupread.query_position!=None] # query_position is None if the base on the padded read is a deletion or a skip (e.g. spliced alignment). 
            
                if pileupcolumn.get_num_aligned()==1:
                    obs_tab.extend(rows)
                else:
                    warnings.warn('Multiple reads were found to overlap at one or more SNP positions.')
                    if handle_multiple_observations=='all':
                        obs_tab.extend(rows)
                    elif handle_multiple_observations=='first':
                        if len(rows)>0: obs_tab.append(rows[0])
                    elif handle_multiple_observations=='random':
                        if len(rows)>0: obs_tab.append(random.choice(rows))
                    elif handle_multiple_observations=='skip':
                        pass
                    else:
                        raise Exception('error: handle_multiple_observations only supports the options \"skip\", \"all\", \"first\" and \"random\".')
    
        info = {'redo-BAQ': fasta_filename=='',
                'handle-multiple-observations' : handle_multiple_observations,
                'min-bq': min_bq,
                'min-mq' :  min_mq,
                'max-depth' :  max_depth,
                'chr_id': chr_id} 

        if kwargs.get('save',True):
            default_output_filename = re.sub('.bam$','',bam_filename.strip().split('/')[-1])+'.obs.p'
            output_filename = default_output_filename if kwargs.get('output_filename','')=='' else kwargs.get('output_filename','') 
            with open( output_filename, "wb") as f:
                pickle.dump(obs_tab, f)
                pickle.dump(info, f)    
    
    finally:
        if genome_reference!=None: genome_reference.close()
        samfile.close()
    
    time1 = time.time()
    print('Done building the observations table in %.2f sec.' % (time1-time0))
    return tuple(obs_tab), info

if __name__ == "__main__": 
    
    parser = argparse.ArgumentParser(
        description='Builds a table of single base observations at known SNP positions.')
    parser.add_argument('bam_filename', metavar='BAM_FILENAME', type=str, 
                        help='BAM file')
    parser.add_argument('legend_filename', metavar='LEGEND_FILENAME', type=str, 
                        help='IMPUTE2 legend file')
    parser.add_argument('-f','--fasta_filename', type=str,metavar='FASTA_FILENAME', default='',
                        help='The faidx-indexed reference file in the FASTA format. ' 
                             'Supplying a reference file will reduce false SNPs caused by misalignments using the Base Alignment Quality (BAQ) method described in the paper “Improving SNP discovery by base alignment quality”, Heng Li, Bioinformatics, Volume 27, Issue 8.')
    parser.add_argument('-u', '--handle-multiple-observations', type=str, 
                        metavar='all/first/random/skip', default='all', 
                        help='We expect to observe at most a single base per SNP. When encountering an exception the default behavior is to keep all the observations.'
                             'However, a few alternative options to handle multiple observations are available: (a) take the first observed base, (b) pick randomly an observed base and (c) skip the observed bases.')
    parser.add_argument('-b', '--min-bq', type=int, 
                        metavar='INT', default=30, 
                        help='Minimum base quaility for observations. Default value is 30.')
    parser.add_argument('-m', '--min-mq', type=int, 
                        metavar='INT', default=30,
                        help='Minimum mapping quaility for observations. Default value 30.')
    parser.add_argument('-d', '--max-depth', type=int, 
                        metavar='INT', default=0,
                        help='Maximum depth coverage to be considered (inclusive). Default value is 0, effectively removing the depth limit.')
    parser.add_argument('-o', '--output-filename', metavar='OUTPUT_FILENAME', type=str, default='',
                        help='Output filename. The default filename is the same as the BAM filename, but with an extension of .obs.p')
    
    retrive_bases(**vars(parser.parse_args()))
    sys.exit(0)