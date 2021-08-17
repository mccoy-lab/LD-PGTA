#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 21:45:36 2021

@author: ariad
"""
import collections, time, pickle, argparse, re, sys, random, os, bz2, gzip

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
    description='Builds a dictionary that lists genomic windows that contain'
                'at least two reads and gives the associated log-likelihood '
                'BPH/SPH ratio (LLR). BPH (Both Parental Homologs) correspond'
                'to the presence of three unmatched haplotypes, while SPH'
                '(Single Parental Homolog) correspond to chromosome gains'
                'involving identical homologs.')
    parser.add_argument('obs_filename', metavar='OBS_FILENAME', type=str,
                        help='A pickle file created by MAKE_OBS_TAB, containing base observations at known SNP positions.')
    parser.add_argument('leg_filename', metavar='LEG_FILENAME', type=str,
                        help='IMPUTE2 legend file')
    parser.add_argument('hap_filename', metavar='HAP_FILENAME', type=str,
                        help='IMPUTE2 haplotype file')
    parser.add_argument('sam_filename', metavar='SAM_FILENAME', type=str,
                        help='IMPUTE2 samples file')
    parser.add_argument('-a', '--admixture', metavar='STR FLOAT', type=str, nargs=2, default=['None','-1'],
                        help='Assume an admixture with a certain ancestry proportion, e.g, EUR 0.8.')
    parser.add_argument('-w', '--window-size', type=int,
                        metavar='INT', default='100000',
                        help='Specifies the size of the genomic window. The default value is 100 kbp. When given a zero-size genomic window, it adjusts the size of the window to include min-reads reads.')
    parser.add_argument('-s', '--subsamples', type=int,
                        metavar='INT', default='32',
                        help='Sets the number of subsamples per genomic window. The default value is 32.')
    parser.add_argument('-o', '--offset', type=int,
                        metavar='INT', default=0,
                        help='Shifts all the genomic windows by the requested base pairs. The default value is 0.')
    parser.add_argument('-m', '--min-reads', type=int, metavar='INT', default=6,
                        help='Takes into account only genomic windows with at least INT reads, admitting non-zero score. The minimal value is 3, while the default is 6.')
    parser.add_argument('-M', '--max-reads', type=int, metavar='INT', default=4,
                        help='Selects up to INT reads from each genomic windows in each bootstrap sampling. The minimal value is 2, while the default value is 4.')
    parser.add_argument('-F', '--min-HF', type=int, metavar='FLOAT', default=0.05,
                        help='Only haplotypes with a frequnecy between FLOAT and 1-FLOAT add to the score of a read. The default value is 0.05.')
    parser.add_argument('-S', '--min-score', type=int, metavar='INT', default=16,
                        help='Consider only reads that reach the minimal score. The default value is 2.')
    parser.add_argument('-O', '--output-filename', type=str, metavar='output_filename',  default='',
                        help='The output filename. The default is the input filename with the extension \".obs.p\" replaced by \".LLR.p\".')
    parser.add_argument('-C', '--compress', metavar='gz/bz2/unc', type=str, default='unc',  choices=['gz','bz2','unc'],
                        help='Output compressed via gzip, bzip2 or uncompressed. Default is uncompressed.')
    args = parser.parse_args()

    print(args.__dict__)
    print(vars(args))

    sys.exit(0)