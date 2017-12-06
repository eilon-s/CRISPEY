#!/usr/bin/env python2

"""
Eilon Sharon April 2016

Merging SNPs with QTLs (taken for example from Bloom et al.)

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
  
    parser = argparse.ArgumentParser("Merging SNPs with QTLs")


    parser.add_argument("input_SNP_table_filename", help= "input SNP table file name")
    
    parser.add_argument("input_QTLs_filename", help= "input QTL table file name")
    
    parser.add_argument("input_selected_QTLs", help= "input selected QTL file name")
    
    parser.add_argument("input_known_genes_with_large_effect", help= "input a table of genes with known large effect")


    parser.add_argument("-o1", "--output_SNP_withQTL_table_filename", dest='output_SNP_withQTL_table_filename', type=str, default='',
                      help='output the SNPs table with additional columns for QTL info table file name (the main output)')

    parser.add_argument("-o2", "--output_selected_qtls_table_filename", dest='output_selected_qtls_table_filename', type=str, default='',
                      help='output selected QTL table file name')
    
    
    parser.add_argument("-o3", "--output_snp_in_selected_QTLs_filename", dest='output_snp_in_selected_QTLs_filename', type=str, default='',
                      help='output SNPs in selected QTLs table file name')
    
    parser.add_argument("-o4", "--output_SNP_merged_with_QTL_table_filename", dest='output_SNP_merged_with_QTL_table_filename', type=str, default='',
                      help='SNPs merged with QTLs cartesian product table file name')
    
    

    args = parser.parse_args()

    # running the function
    merge_SNPs_with_QTLs(input_SNP_table_filename, input_QTLs_filename, input_QTLs_H_filename,
                         input_selected_QTLs, input_known_genes_with_large_effect,
                         output_selected_qtls_table_filename, output_snp_in_selected_QTLs_filename,
                         output_SNP_merged_with_QTL_table_filename, output_SNP_withQTL_table_filename)


if __name__ == '__main__':
  main()




