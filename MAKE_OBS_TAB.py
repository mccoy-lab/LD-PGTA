#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

MAKE_OBS_TAB

This script extracts single base observations at SNPs positions from a given
sequence. It requires an aligned and sorted BAM file with the sequence, as well
as an IMPUTE2 legend format, which contains the SNPs positions.
The observed alleles, together with their chromosome position and line number in
the legend file, are organized in a table.

THIS SCRIPT IS BASED ON MAKE_OBS_TABLE FROM THE TILDE PROJECT (https://github.com/svohr/tilde).

Daniel Ariad (daniel@ariad.org)
April 1st, 2020

"""

import sys, os, subprocess, time, tempfile, pickle, argparse, re, random

def read_impute2_legend(legend_filename):
    """ Reads in the SNPs from the impute legend file and builds a list of
        tuples, containing snp positions, reference and alternate alleles. """
    
    leg_in = open(legend_filename, 'r')
    snp_positions = list()
    leg_in.readline()   # Bite off the header
    for line in leg_in:
        pos = line.rstrip().split()[1] # Trailing whitespaces stripped from the right of the string. Then it splits the string into a list of words.
        snp_positions.append(int(pos))
    leg_in.close()

    return snp_positions


def generate_snps_position_list(snp_positions, chrm_id):
    """ Creates a position list file of our SNPs, which would be used by
        samtools mpileup to find observed bases. 
        The position list file contains two columns for chromosome and position.
        The columns are TAB-separated and countings starts from 1. """ 
    
    pos_file = tempfile.NamedTemporaryFile(delete=False, mode='w') 

    for pos in snp_positions:
        pos_file.write('%s\t%d\n' % (chrm_id, pos)) #Enter data into file.
    
    pos_file.close()
    
    return pos_file.name

def parse_bases_field(ref, bases_str):
    """ Parses the bases field of mpileup output and returns the observed
        bases at this position. Positions where indels occur are not
        considered. """
    
    dictionary = {'.*[*+-].*': '',  #(1)  Skip deletions\insertions and remove sequence end mark.
                  '\\^.': '',       #(2)  Remove sequence start mark followed by mapping quality.
                  '\$': '',         #(3)  Remove sequence end mark
                  '[,.]': ref}      #(4)  Replace with matching reference.
    
    for pattern, substitute in dictionary.items():
        bases_str = re.sub(pattern, substitute, bases_str)
    
    bases_str = bases_str.upper()
    
    return bases_str
     
def extract_bases_from_pileup(mpu_output, max_cov):
    """ Extracts, from the pileup format data, the observed bases at SNPs
        positions. If more than one base is present, one is chosen at random."""
        
    extracted = []
    for line in mpu_output.splitlines():
        fields = line.rstrip().split('\t')
        pos, ref, cov, bases = int(fields[1]), fields[2], int(fields[3]), fields[4]
        if cov > 0  and (max_cov == 0 or cov <= max_cov):
            reads = parse_bases_field(ref, bases)
            if len(set(reads))==1:  #Only non overlaping reads are considered.
                extracted.append((pos,reads[0]))
                
    return extracted    

def get_snps_samtools(bam_filename, pos_filename, 
                      min_mapq, min_baseq, max_cov, fasta_ref, verbose, samtools_dir):
    """ Runs samtools mpileup as a subprocess and returns the output """
    
    samtools_inst_dir = samtools_dir if samtools_dir=='' or samtools_dir[-1]=='/' else samtools_dir+'/'

    samtools_args = [samtools_inst_dir+'samtools', 'mpileup',
                     '-q', str(min_mapq), #Minimum mapping quality for an alignment to be used.
                     '-Q', str(min_baseq), #Minimum base quality for a base to be considered.
                     '-d', str(max_cov), #At a position, read maximally INT reads per input file. 
                     '-l', pos_filename] #Position list file containing a list of regions or sites where pileup or BCF should be generated.
    
    samtools_args.append(bam_filename) if fasta_ref == '' else samtools_args.extend(['-f', fasta_ref, bam_filename])
        
    if verbose: 
        sys.stderr.write('Calling samtools:\n  %s\n' % ' '.join(samtools_args))
    
    mpileup = subprocess.Popen(samtools_args, 
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    out, err = mpileup.communicate()
    return out, err 



def main(bam, impute_leg, chrm, min_mq, min_bq, max_cov, fasta_ref, verbose, samtools_dir, write_to_stdout): 
    """ Builds a table of SNP positions and base observations from mapped reads. """
    
    if not os.path.isfile(bam): raise Exception('error: BAM file does not exist.')
    if not os.path.isfile(impute_leg): raise Exception('error: LEGEND file does not exist.')
    
    try:
        snp_positions = read_impute2_legend(impute_leg) 
        pos_filename = generate_snps_position_list(snp_positions, chrm)
        mpu_out, mpu_err = get_snps_samtools(bam, pos_filename, 
                                             min_mq, min_bq, 
                                             max_cov, fasta_ref, verbose, samtools_dir) 
    finally:
        os.remove(pos_filename)
    
    extracted = extract_bases_from_pileup(mpu_out.decode('UTF-8'), max_cov)   
    
    snp_positions_iter, obs_tab, i0 = iter(snp_positions), list(), 0
    
    for (pos1, read) in extracted:
        for i, pos2 in enumerate(snp_positions_iter,start=i0):
            if pos1 == pos2:
                obs_tab.append([i, pos1, read])
                i0 = i + 1
                break
                        
    return obs_tab
    
if __name__ == "__main__": 
    time0 = time.time()
    
    parser = argparse.ArgumentParser(
        description='Builds a table of SNP positions, single base observations from mapped reads and '
                    'their corresponding line number within the IMPUTE2 legend file.')
    parser.add_argument('bam', metavar='BAM_filename', type=str, 
                        help='BAM file')
    parser.add_argument('impute_leg', metavar='legend_filename', type=str, 
                        help='IMPUTE2 legend file')
    parser.add_argument('chrm', metavar='chromosomeID', type=str,
                        help='Chromosome ID')
    parser.add_argument('-q', '--min-bq', type=int, 
                        metavar='BASEQ', default=30, 
                        help='Minimum base quaility for observations. Default value is 30.')
    parser.add_argument('-Q', '--min-mq', type=int, 
                        metavar='MAPQ', default=30,
                        help='Minimum mapping quaility for observations. Default value 30.')
    parser.add_argument('-c', '--max-cov', type=int, 
                        metavar='COV', default=0,
                        help='Maximum coverage to allow (inclusive). Default value is 0.')
    parser.add_argument('-f', '--fasta-ref', type=str,
                        default='', metavar='FASTA_REF',
                        help='The faidx-indexed reference file in the FASTA format. ' 
                             'The file can be optionally compressed by bgzip. '
                             'Supplying a reference file will reduce false SNPs caused by misalignments using the Base Alignment Quality (BAQ) method described in the paper “Improving SNP discovery by base alignment quality”, Heng Li, Bioinformatics, Volume 27, Issue 8.')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print debug information')
    parser.add_argument('-s', '--samtools-dir', type=str,
                        default='', metavar='SAMTOOLS_DIR',
                        help='The directory where samtools are installed.')
    parser.add_argument('-w', '--write-to-stdout', action='store_true',
                        help='Write output to stdout (standard output).')
    
    args = parser.parse_args()
    obs_tab = main(**vars(args))
    
    info = {'chr' : args.chrm, 'min-bq': args.min_bq, 'min-mq' :  args.min_mq, 'max-cov' :  args.max_cov} 

    if args.write_to_stdout:
        tab_out = sys.stdout
        tab_out.write('SNPs indices, SNPs Positions, Observed Bases\n')
        for (ind, pos, read) in obs_tab: 
            tab_out.write('%d\t%d\t%s\n' % (ind, pos, read))

    else:
        EXT = '.OBS[BAQ].p' if args.fasta_ref != '' else '.OBS.p'
        output_filename = re.sub('.bam$','',args.bam.split('/')[-1])+EXT
        with open( output_filename, "wb") as f:
            pickle.dump(obs_tab, f)
            pickle.dump(info, f)
     
    time1 = time.time()
    print('Done in %.2f sec.' % (time1-time0))
    sys.exit(0)