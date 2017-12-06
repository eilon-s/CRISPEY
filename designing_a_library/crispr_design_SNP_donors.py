#!/usr/bin/env python2

"""
Eilon Sharon April 2016

Designs donor sequences for guides that target SNPs according to a design table

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
from extract_guides_functions import cal_guide_features




def main():
  
    parser = argparse.ArgumentParser("Designs donor sequences for guides that target SNPs according to a design table")


    parser.add_argument("input_SNP_donor_design_table_filename", help= "SNP donor design table")
    
    parser.add_argument("input_SNP_table_filename", help= "SNPs table (processed from VCF)")
    
    parser.add_argument("input_guides_table_filename", help= "SNP guides table (optionaly - with features)")

    parser.add_argument("input_genome_fasta_filename", \
                      help="genome fasta file. If bowtie2 is used for offtargets mapping (default), bowtie2 index should exist with similar prefix")



    parser.add_argument("-e", "--excluded_seqs", dest='excluded_seqs', default=['TTTT'], nargs='*',
                      help=('sequences that will make the construct non functional'))
    
    parser.add_argument("-l", "--donor_length", dest='donor_length', default=100, type=int,
                      help=('sequences that will make the construct non functional'))


    parser.add_argument("-m", "--min_dist_cut_to_donor_edge", dest='min_dist_cut_to_donor_edge', type=int, default=30,
                      help='path to bowtie2 executable')


    parser.add_argument("-o", "--output_donor_table_filename", dest='output_donor_table_filename', type=str, default='',
                      help='output file for the donor sequences table')

    args = parser.parse_args()

    # running the function
    out_guide_donor_df = design_donor_for_SNP_guides(input_design_table_filename = args.input_SNP_donor_design_table_filename, 
                            input_SNP_table_filename = args.input_SNP_table_filename, input_guides_with_features_table_filename = args.input_guides_table_filename,
                            input_genome_fasta_filename = args.input_genome_fasta_filename,
                            donor_length = args.donor_length, excluded_seqs = args.excluded_seqs, min_dist_cut_to_donor_edge = args.min_dist_cut_to_donor_edge,
                            output_donor_table_filename = args.output_donor_table_filename)



if __name__ == '__main__':
  main()


