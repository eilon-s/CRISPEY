#!/usr/bin/env python2

"""
Eilon Sharon April 2016

Draw SNP-Guide general statistics figures
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
from crispr_draw_func import draw_snp_guide_stats




def main():
  
    parser = argparse.ArgumentParser("Draw SNP-guide statistics figures")


    parser.add_argument("input_SNP_table_filename", help= "SNPs table")
    parser.add_argument("input_guides_with_features_table_filename", help= "Guides with features table")
    
    parser.add_argument("-d", "--input_donors_table_filename", dest='input_donors_table_filename', type=str, default='',
                      help='input donors table')
    
    parser.add_argument("output_filenames_prefix", help= "Output prefix for the figures")

    args = parser.parse_args()

    # running the function
    draw_design_stats(input_SNP_table_filename = args.input_SNP_table_filename,
                      input_guides_with_features_table_filename = args.input_guides_with_features_table_filename,
                      input_donors_table_filename = args.input_donors_table_filename,
                      output_filenames_prefix = args.output_filenames_prefix)


if __name__ == '__main__':
  main()

