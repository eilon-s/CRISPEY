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
from extract_guides_functions import extract_guides_for_snps

################################################################################################################

def main():
  
    parser = argparse.ArgumentParser("Extracts guides for SNPs")


    parser.add_argument("input_snps_vcf_filename", help= "SNPs VCF file")


    parser.add_argument("input_genome_fasta_filename", \
                      help="genome fasta file. If bowtie2 is used for offtargets mapping (default), bowtie2 index should exist with similar prefix")

    parser.add_argument("output_SNP_table_filename", help="out SNPs table")

    parser.add_argument("output_guides_table_filename", help="out guides table")

    parser.add_argument("-p", "--pam", dest='PAM_seq', type=str, default='GA',
                      help='2 letters pam sequence')
    
    parser.add_argument("-l", "--guide_length", dest='guide_length', default=20, type=int,
                      help=('guide length (20)'))
    
    parser.add_argument("-m", "--edit_max_distance_from_PAM5prime", dest='edit_max_distance_from_PAM5prime', default=17,  type=int,
                      help=('max position in guide of the edit'))


    parser.add_argument("-o", "--output_guides_with_features_table_filename", dest='output_guides_with_features_table_filename', type=str, default='',
                      help='output file for the guides table')

    args = parser.parse_args()

    # running the function
    extract_guides_for_snps(input_snps_vcf_filename = args.input_snps_vcf_filename, input_genome_fasta_filename = args.input_genome_fasta_filename, 
                            output_SNP_table_filename = args.output_SNP_table_filename, output_guides_table_filename = args.output_guides_table_filename,
                            pam_seqs = [args.PAM_seq], guide_length = args.guide_length, edit_max_distance_from_PAM5prime = args.edit_max_distance_from_PAM5prime)
    

if __name__ == '__main__':
  main()

