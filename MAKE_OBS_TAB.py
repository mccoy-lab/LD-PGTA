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

def retrive_bases(bam_filename,legend_filename,fasta_filename,handle_multiple_observations,min_bq,min_mq,max_depth,**more_kwargs):
    """ Retrives observed bases from known SNPs position. """
    
    if not os.path.isfile(bam_filename): raise Exception('Error: BAM file does not exist.')
    if not os.path.isfile(legend_filename): raise Exception('Error: LEGEND file does not exist.')
    if fasta_filename!='' and not os.path.isfile(fasta_filename): raise Exception('Error: FASTA file does not exist.')

    obs_tab = list()
    
    try:
        genome_reference = pysam.FastaFile(fasta_filename) if fasta_filename!='' else None
        samfile = pysam.AlignmentFile(bam_filename, 'rb' )
        leg_tab = read_impute2(legend_filename, filetype='legend')

        
        kwarg = {'contig': leg_tab[0][0],       # The chr_id of the considered chromosome.
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
        
        for pileupcolumn in samfile.pileup(**kwarg):
            
            while pileupcolumn.pos > pos-1: 
                i, (chr_id,pos,ref,alt) = next(leg_tab_iterator)  
            
            if pileupcolumn.pos == pos-1:
                
                rows = [(chr_id, pos, (ref,alt), i, pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]) for pileupread in pileupcolumn.pileups if pileupread.query_position!=None] # query_position is None if the base on the padded read is a deletion or a skip (e.g. spliced alignment). 
            
                if pileupcolumn.get_num_aligned()==1:
                    obs_tab.extend(rows)
                else:
                    warnings.warn('Multiple reads were found to overlap at one or more SNP positions.')
                    if handle_multiple_observations=='all':
                        obs_tab.extend(rows)
                    elif handle_multiple_observations=='first' and len(rows)>0:
                        obs_tab.append(rows[0])
                    elif handle_multiple_observations=='random' and len(rows)>0:
                        obs_tab.append(random.choice(rows))
                    elif handle_multiple_observations=='skip':
                        pass
                    else:
                        raise Exception('error: handle_multiple_observations only supports the options \"skip\", \"all\", \"first\" and \"random\".')
        
    finally:
        if genome_reference!=None: genome_reference.close()
        samfile.close()
    
    return tuple(obs_tab)

if __name__ == "__main__": 
    time0 = time.time()
    
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
                        metavar='all/first/random/skip', default='skip', 
                        help='We expect to observe at most a single base per SNP. When encountering an exception the default behavior is to skip the SNP. '
                             'However, a few alternative options to handle multiple observations are avaible: (a) take the first observed base, (b) pick randomly an observed base and (c) keep all the observed bases.')
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
    
    args = parser.parse_args()
    
    obs_tab  = retrive_bases(**vars(args))
    
    info = {'redo-BAQ': args.fasta_filename=='', 'handle-multiple-observations' : args.handle_multiple_observations, 'min-bq': args.min_bq, 'min-mq' :  args.min_mq, 'max-depth' :  args.max_depth} 

    default_output_filename = re.sub('.bam$','',args.bam_filename.strip().split('/')[-1])+'.obs.p'
    output_filename = default_output_filename if args.output_filename=='' else args.output_filename 
    
    with open( output_filename, "wb") as f:
        pickle.dump(obs_tab, f)
        pickle.dump(info, f)
     
    time1 = time.time()
    print('Done in %.2f sec.' % (time1-time0))
    sys.exit(0)
    

######################
#if __name__ == "__main__": 
#    print("Executed when invoked directly")
#    bam_filename = '../BAMs_hg38/SRR6676163.hg38.bam'
#    fasta_filename = None#'../genome_ref_hg38/hg38.fa'
#    legend_filename = '../make_reference_panel/ref_panel.EUR.hg38.BCFtools/chr21_EUR_panel.legend'
#    handle_multiple_observations = 'skip'
#    min_bq = 30
#    min_mq = 30
#    max_depth = 0
#    result = retrive_bases(bam_filename,legend_filename,fasta_filename,handle_multiple_observations,min_bq,min_mq,max_depth)
#    #test()
#    time1 = time.time()
#    print('Done in %.3f sec.' % (time1-time0))
#else: 
#    print("Executed when imported")
#
#
#def pileup(genome_reference,samfile,chr_id,start_position,end_position):
#    additional_arguments = {'truncate': True,               # By default, the samtools pileup engine outputs all reads overlapping a region. If truncate is True and a region is given, only columns in the exact region specificied are returned.
#                            'max_depth': 0,                 # Maximum read depth permitted. The default limit is ‘8000’.
#                            'stepper': 'samtools',          # The following arguments all pertain to the samtools stepper:
#                            'min_base_quality': 30,         # Minimum base quality. Bases below the minimum quality will not be output.
#                            'min_mapping_quality': 30,      # Only use reads above a minimum mapping quality. The default is 0.
#                            'ignore_overlaps': True,        # If set to True, detect if read pairs overlap and only take the higher quality base.
#                            'ignore_orphans': True,         # Ignore orphans (paired reads that are not in a proper pair).
#                            'compute_baq': False,            # Re-alignment computing per-Base Alignment Qualities (BAQ). The default is to do re-alignment. 
#                            'fastafile': genome_reference}  # Requires a FastaFile object of a reference sequence. If none is present, no realignment will be performed.                 
#
#    for pileupcolumn in samfile.pileup(chr_id,start_position,end_position, **additional_arguments):
#        print ("\ncoverage at base %s = %s" %
#               (pileupcolumn.pos, pileupcolumn.get_num_aligned()))
#        for pileupread in pileupcolumn.pileups:
#            position = pileupread.query_position
#            #print(position)
#            if position!=None: # Query position is None if the base on the padded read is a deletion or a skip (e.g. spliced alignment). 
#                print ('\tbase in read %s = %s' %
#                      (pileupread.alignment.query_name,
#                       pileupread.alignment.query_sequence[position]))
#                
#def test():
#    try:
#        genome_reference = None #pysam.FastaFile('/Users/ariad/Dropbox/postdoc_JHU/Tools/genome_ref_hg38/hg38.fa')
#        samfile = pysam.AlignmentFile('../BAMs_hg38/SRR6676163.hg38.bam', 'rb' )
#        #pileup(genome_reference, samfile,'chr21', 5033830, 5033850)
#        pileup(genome_reference, samfile,'chr21', 38119107, 38119108)
#    finally:
#        if genome_reference!=None: genome_reference.close()
#        samfile.close()
#
########################    

