#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAKE_OBS_TAB

This script extracts single base observations at SNP positions from a given
sequence. It requires an aligned and sorted BAM file with the sequence, as well
as an IMPUTE2 legend format, which contains the SNPs positions.
The observed alleles, together with their associated read ID, chromosome
position and line number in the legend file, are organized in a table.

Daniel Ariad (daniel@ariad.org)
Jan 3rd, 2021
"""
import sys, os, time, random, warnings, argparse, re, pickle, gzip, bz2

try:
    import pysam
except ModuleNotFoundError:
    print('caution: the module pysam is missing.')
  
warnings.formatwarning = lambda message, category, filename, lineno, file=None, line=None: 'Caution: %s\n' % message

def read_impute2(filename,**kwargs):
    """ Reads an IMPUTE2 file format (LEGEND/HAPLOTYPE/SAMPLE) and builds a list
        of lists, containing the dataset. """
    
    filetype = kwargs.get('filetype', None)
    
    def leg_format(line):
        rs_id, pos, ref, alt = line.strip().split()
        return ('chr'+rs_id[:2].rstrip(':'), int(pos), ref, alt)  
    
    def sam_format(line):
          sample_id, group1, group2, sex = line.strip().split(' ')
          return (sample_id, group1, group2, int(sex))  
      
    with (gzip.open(filename,'rt') if filename[-3:]=='.gz' else open(filename, 'r')) as impute2_in:
        if filetype == 'leg': 
            impute2_in.readline()   # Bite off the header
            result = tuple(map(leg_format,impute2_in))
            
        elif filetype == 'hap':
            firstline = impute2_in.readline()   # Get first line
            a0 = int(firstline.replace(' ', ''), 2)
            a1 = (int(line.replace(' ', ''), 2) for line in impute2_in)
            hap_tab = (a0, *a1)
            number_of_haplotypes = len(firstline.strip().split())
            result = hap_tab, number_of_haplotypes
            
        elif filetype == 'sam': 
            impute2_in.readline()   # Bite off the header
            result = tuple(map(sam_format,impute2_in))
        
        else:
            result = tuple(line.strip().split() for line in impute2_in)
    
    return result 

def save_obs(obs_tab,info,compress,bam_filename,output_filename,output_dir):
    """ Saves the observations table together with information about
        the chromosome number, depth of coverage, and flags that were used.
        Also, data compression is supported in gzip and bzip2 formats. """
        
    Open = {'bz2': bz2.open, 'gz': gzip.open}.get(compress, open)
    ext = ('.'+compress) * (compress in ('bz2','gz'))
    default_output_filename = re.sub('.bam$',f".{info['chr_id']:s}.obs.p{ext:s}",bam_filename.strip().split('/')[-1])
    output_filename = default_output_filename * (output_filename=='')
    output_dir = re.sub('/$','',output_dir)+'/' #undocumented option
    if output_dir!='' and not os.path.exists(output_dir): os.makedirs(output_dir)
    with Open(output_dir + output_filename, "wb") as f:
        pickle.dump(obs_tab, f, protocol=4)
        pickle.dump(info, f, protocol=4)       
    return output_dir + output_filename

def retrive_bases(bam_filename,legend_filename,fasta_filename,handle_multiple_observations,min_bq,min_mq,max_depth,output_filename,compress,**kwargs):
    """ Retrives observed bases from known SNPs position. """
    time0 = time.time()
    random.seed(a=None, version=2) #I should set a=None after finishing to debug the code.

    if not os.path.isfile(bam_filename): raise Exception('Error: BAM file does not exist.')
    if not os.path.isfile(legend_filename): raise Exception('Error: LEGEND file does not exist.')
    if fasta_filename!='' and not os.path.isfile(fasta_filename): raise Exception('Error: FASTA file does not exist.')

    obs_tab = list()
    
    try:
        genome_reference = pysam.FastaFile(fasta_filename) if fasta_filename!='' else None
        samfile = pysam.AlignmentFile(bam_filename, 'rb' )
        leg_tab = read_impute2(legend_filename, filetype='leg')

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
            
            while pileupcolumn.pos > pos-1:  ### Chromosomal position starts from 1 in bam files, while it starts from 0 in obs files.
                impute2_index, (chr_id,pos,ref,alt) = next(leg_tab_iterator)  
            
            if pileupcolumn.pos == pos-1:
                
                rows = [(pos, impute2_index, pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]) 
                        for pileupread in pileupcolumn.pileups if pileupread.query_position!=None] # query_position is None if the base on the padded read is a deletion or a skip (e.g. spliced alignment). 
                
                if pileupcolumn.get_num_aligned()==1:
                    obs_tab.extend(rows) #Each tuple in the list has the form (chromosomal posistion, associated line number in the legend file, reads id, observed allele)
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
    
        info = {'redo-BAQ': fasta_filename!='',
                'handle-multiple-observations' : handle_multiple_observations,
                'min-bq': min_bq,
                'min-mq' :  min_mq,
                'max-depth' :  max_depth,
                'chr_id': chr_id,
                'depth': len(obs_tab)/len(leg_tab)} 
        
        if output_filename!=None:
            save_obs(obs_tab,info,compress,bam_filename,output_filename,kwargs.get('output_dir', 'results'))
    
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
    parser.add_argument('legend_filename', metavar='LEG_FILENAME', type=str, 
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
                        help='Output filename. The default filename is the same as the BAM filename, but with an extension of .chr_id.obs.p')
    parser.add_argument('-c', '--compress', metavar='gz/bz2/unc', type=str, default='unc',
                        help='Output compressed via gzip, bzip2 or uncompressed. Default is uncompressed.')
    
    retrive_bases(**vars(parser.parse_args()))
    sys.exit(0)
else:
    print('The module MAKE_OBS_TAB was imported.')