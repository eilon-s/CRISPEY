#!/usr/bin/env python2

"""
Eilon Sharon April 2016

Extract guides for SNPs

"""


from __future__ import division

import sys
import os
import getopt
import argparse
import itertools
  
  
import math
import re
import numpy as np
import pandas as pd
import random


sys.path.append(os.path.expanduser('~/bin/python/crispr'))
from extract_guides_functions import extract_guides_for_snps, cal_guide_features, design_donor_for_SNP_guides

################################################################################################################

def main():
  
    parser = argparse.ArgumentParser("Extracts guides for SNPs")

    # input files

    parser.add_argument("input_snps_vcf_filename", help= "SNPs VCF file")


    parser.add_argument("input_genome_fasta_filename", \
                      help="genome fasta file. If bowtie2 is used for offtargets mapping (default), bowtie2 index should exist with similar prefix")
    
    parser.add_argument("input_SNP_donor_design_table_filename", help='input SNP donor design table filename')

    #output files
    
    parser.add_argument("output_SNP_table_filename", help="out SNPs table")

    parser.add_argument("output_guides_table_filename", help="out guides table")
    
    
    parser.add_argument("output_guides_with_features_table_filename", help='output file for the guides table')
    
    parser.add_argument("output_donor_table_filename", help='output file for the donor sequences table')
    
    
    # options

    parser.add_argument("-p", "--pam", dest='PAM_seq', type=str, default='GA',
                      help='2 letters pam sequence')
    
    parser.add_argument("-g", "--guide_length", dest='guide_length', default=20, type=int,
                      help=('guide length (20)'))
    
    parser.add_argument("-d", "--donor_length", dest='donor_length', default=100, type=int,
                      help=('sequences that will make the construct non functional'))
    
    parser.add_argument("-m", "--edit_max_distance_from_PAM5prime", dest='edit_max_distance_from_PAM5prime', default=17,  type=int,
                      help=('max position in guide of the edit'))
    
    parser.add_argument("-e", "--excluded_seqs", dest='excluded_seqs', default=['TTTTTTTTTT', 'AAAAAAAAAA'], nargs='*',
                      help=('sequences that will make the construct non functional'))


    parser.add_argument("-c", "--min_dist_cut_to_donor_edge", dest='min_dist_cut_to_donor_edge', type=int, default=30,
                      help='path to bowtie2 executable')
    

    args = parser.parse_args()

    # running the function
    
    print('extracting guides')
    extract_guides_for_snps(input_snps_vcf_filename = args.input_snps_vcf_filename, input_genome_fasta_filename = args.input_genome_fasta_filename, 
                            output_SNP_table_filename = args.output_SNP_table_filename, output_guides_table_filename = args.output_guides_table_filename,
                            pam_seqs = [args.PAM_seq], guide_length = args.guide_length, edit_max_distance_from_PAM5prime = args.edit_max_distance_from_PAM5prime)
    

    print('adding features to guides')
    BOWTIE_exe = "bowtie2"
    out_guides_with_features_df = cal_guide_features(args.output_guides_table_filename, args.input_genome_fasta_filename, 
                                                 args.output_guides_with_features_table_filename, 
                                                 args.PAM_seq, args.excluded_seqs, args.donor_length,BOWTIE_exe = BOWTIE_exe)
    
    print('designing donors')
    out_guide_donor_df = design_donor_for_SNP_guides(input_design_table_filename = args.input_SNP_donor_design_table_filename, 
                            input_SNP_table_filename = args.output_SNP_table_filename, input_guides_with_features_table_filename = args.output_guides_with_features_table_filename,
                            input_genome_fasta_filename = args.input_genome_fasta_filename,
                            donor_length = args.donor_length, excluded_seqs = args.excluded_seqs, min_dist_cut_to_donor_edge = args.min_dist_cut_to_donor_edge,
                            output_donor_table_filename = args.output_donor_table_filename)

    
    
if __name__ == '__main__':
  main()

