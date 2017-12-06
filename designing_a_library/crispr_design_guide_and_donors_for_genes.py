#!/usr/bin/env python2

"""
Eilon Sharon April 2016

Design guides and donor sequence as specified in a design table
"""


from __future__ import division


#import platform
#print(platform.python_version())


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


#sys.path.append(os.path.expanduser('~/anaconda/lib/python2.7/site-packages'))


sys.path.append(os.path.expanduser('~/bin/python/crispr'))
from crispr_design_gene_guides import design_guide_and_donor_for_genes
from extract_guides_functions import which

# ~/bin/python/crispr/crispr_design_guide_and_donors_for_genes.py /Users/eilon/Dropbox/postdoc/crispr/Output/debug/design_guides_in_genes_edited.txt /Users/eilon/Dropbox/postdoc/crispr/Input/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208_noFasta.gff /Users/eilon/Dropbox/postdoc/crispr/01-GATK_Bloom_RM_Data/S288C_R64_REAN_140405.fa -g /Users/eilon/Dropbox/postdoc/crispr/Output/debug/debug_cds_GA18bp_guides_d2.tab -d /Users/eilon/Dropbox/postdoc/crispr/Output/debug/debug_cds_GA18bp_donors_d2.tab

#########################################################################################

def main():
  
  
  parser = argparse.ArgumentParser("Design guides and donor sequence as specified in a design table")
  
  parser.add_argument("input_design_table_filename", type=str,  help="Input table of the duide and donor sequences design. \
                      Columns:gene_id,guide_num,donor_mut_type,num_donor_variants,mut_pos_in_guide,donor_length,donor_seq_offsets")
  parser.add_argument("input_annot_gff_filename", help="SGD genes gff file")
  parser.add_argument("input_genome_fasta_filename", help="genome fasta file. If bowtie2 is used for offtargets mapping (default), bowtie2 index should exist with similar prefix")
  
  
  parser.add_argument("-p", "--pam", dest='PAM_seq', type=str, default='GA',
                      help='2 letters pam sequence')
  
  parser.add_argument("-e", "--excluded_seqs", dest='excluded_seqs', default=['TTTT'], nargs='*',
                      help=('sequences that will make the construct non functional'))
  
  parser.add_argument("-b", "--BOWTIE_exe", dest='BOWTIE_exe', type=str, default=which('bowtie2'),
                      help='path to bowtie2 executable')
  
  parser.add_argument("-f", "--off_targets_min_mismatch", dest='off_targets_min_mismatch', type=int, default=10,
                      help='min number of mismatches in allowed off targets (>3 means no off targets are allowed)')
  
  parser.add_argument("-a", "--min_azimuth_score", dest='min_azimuth_score', type=float, default=0.5,
                      help='min Azimuth score allowed, guides below the score will be filtered out')
  
  
  parser.add_argument("-s", "--sort_by", dest='sort_by', type=str, default='5to3',
                      help='by what method to sort the guides (currently only 5to3 is implemented)')
  
  parser.add_argument("-g", "--output_guides_df_filename", dest='output_guides_df_filename', type=str, default='',
                      help='output file for the guides table')
  
  parser.add_argument("-d", "--output_donors_df_filename", dest='output_donors_df_filename', type=str, default='',
                      help='output file for the donors table')
  

  args = parser.parse_args()
    
  
  out_guide_df, out_guide_donor_df = design_guide_and_donor_for_genes(args.input_design_table_filename, args.input_annot_gff_filename, args.input_genome_fasta_filename, \
                                      PAM_seq = args.PAM_seq, excluded_seqs = args.excluded_seqs, BOWTIE_exe = args.BOWTIE_exe, \
                                      off_targets_min_mismatch = args.off_targets_min_mismatch, min_azimuth_score = args.min_azimuth_score, sort_by = args.sort_by, \
                                      output_guides_df_filename = args.output_guides_df_filename, output_donors_df_filename = args.output_donors_df_filename)


if __name__ == '__main__':
  main()



