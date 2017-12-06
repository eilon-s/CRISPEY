#!/usr/bin/env python2

"""
Eilon Sharon April 2016

Designs oligo sequences from guide and donor sequences

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
from extract_guides_functions import generate_oligo_from_guide_and_donors



def main():
  
    parser = argparse.ArgumentParser("Designs oligo sequences form guides and donors using a design table")


    parser.add_argument("input_oligo_design_table_filename", help= "Oligo design table. 'Name' column must contain 'guide' and 'donor',	pos_start & pos_end are zero based. 'DNA_seq' contains the sequence")
    
    parser.add_argument("input_guide_table_filename", help= "input guides table file name")
    
    parser.add_argument("input_donor_table_filename", help= "input donor table file name")

    parser.add_argument("input_genome_fasta_filename", \
                      help="genome fasta file. If bowtie2 is used for offtargets mapping (default), bowtie2 index should exist with similar prefix")


    parser.add_argument("-o", "--output_donor_table_filename", dest='output_donor_table_filename', type=str, default='',
                      help='output file for oligo sequences table')

    args = parser.parse_args()

    # running the function
    generate_oligo_from_guide_and_donors(input_oligo_design_table_filename,
                                         input_guide_table_filename, input_donor_table_filename,
                                         input_guide_iloc = None, input_donor_iloc = None,
                                         output_oligos_table_filename = output_donor_table_filename)


if __name__ == '__main__':
  main()


