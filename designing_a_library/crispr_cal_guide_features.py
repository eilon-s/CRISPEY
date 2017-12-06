#!/usr/bin/env python2


"""
Eilon Sharon April 2016

Calculate features for guides

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
from extract_guides_functions import cal_guide_features,which




def main():
  
    parser = argparse.ArgumentParser("""Calculates features for guides. For example: 
                                    1. Azimuth score 
                                    2. Off targets 
                                    3. excluded sequences yes/no 
                                    4. is in nuclear chromosomes. 
                                    The following columns will be added to the table: 
                                    1. Azimuth
                                    2. guide_map_mismatch_0 
                                    3. guide_map_mismatch_1 
                                    4. guide_map_mismatch_2 
                                    5. guide_map_mismatch_3
                                    6. contain_excluded_sequences
                                    7. is_nuclear_chromosome
                                    8. dna_around_guide""")


    parser.add_argument("input_guides_table_filename", \
                      help= \
                        """
                        Guide table without features column that contains the columns (order is not important): 
                        1. SNP_id 
                        2. SNP_chr_pos 
                        3. SNP_pos_in_guide 
                        4. guide_id 
                        5. chrom 
                        6. guide0_chr_pos
                        7. guide_cut_chr_pos
                        8. guide_strand 
                        9. upstream_or_downstream (downstream = + strand, upstream = - strand. This is for backward compatability)
                        10. guide_noPAM
                        11. guide_PAM_m4p3 
                        (the input can be the output of crispr_extract_guide_for_snps)""")

    parser.add_argument("input_genome_fasta_filename", \
                      help="genome fasta file. If bowtie2 is used for offtargets mapping (default), bowtie2 index should exist with similar prefix")


    parser.add_argument("-p", "--pam", dest='PAM_seq', type=str, default='GA',
                      help='2 letters pam sequence')

    parser.add_argument("-e", "--excluded_seqs", dest='excluded_seqs', default=['TTTT'], nargs='*',
                      help=('sequences that will make the construct non functional'))
    
    parser.add_argument("-l", "--donor_length", dest='donor_length', default=100, type=int,
                      help=('sequences that will make the construct non functional'))


    parser.add_argument("-b", "--BOWTIE_exe", dest='BOWTIE_exe', type=str, default=which('bowtie2'),
                      help='path to bowtie2 executable')


    parser.add_argument("-o", "--output_guides_with_features_table_filename", dest='output_guides_with_features_table_filename', type=str, default='',
                      help='output file for the guides table')

    args = parser.parse_args()

    # running the function
    out_guides_with_features_df = cal_guide_features(args.input_guides_table_filename, 
                                                   args.input_genome_fasta_filename, 
                                                   args.output_guides_with_features_table_filename, 
                                                   PAM_seq = args.PAM_seq, 
                                                   excluded_seqs = args.excluded_seqs, donor_length = args.donor_length,
                                                   BOWTIE_exe = args.BOWTIE_exe)


if __name__ == '__main__':
  main()



