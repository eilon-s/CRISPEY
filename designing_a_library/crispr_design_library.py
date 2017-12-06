#!/usr/bin/env python2

# coding: utf-8


# Eilon Sharon April 2016

#get_ipython().magic(u'matplotlib inline')

# pip install roman -t /home/eilon/software/python_modules

from __future__ import division

import sys
import os

import getopt
import argparse

import random
import pandas as pd
import numpy as np


import matplotlib
#matplotlib.style.use('ggplot')
import matplotlib.pyplot as plt


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import pairwise2


sys.path.append(os.path.expanduser('~/software/python_modules'))

sys.path.append(os.path.expanduser('~/bin/python/crispr'))
from extract_guides_functions import which, extract_guides_for_snps, cal_guide_features, \
                                        generate_qtl_gene_guide_design_table, \
                                        design_donor_for_SNP_guides, \
                                        generate_oligo_from_guide_and_donors, \
                                        parse_essential_genes_to_design_table, \
                                        rank_and_filter_SNP_guides, \
                                        merge_SNPs_with_QTLs, \
                                        filter_1pbMut_donor_table, \
                                        write_output_oligos, \
                                        write_small_control_oligo_set
                        
from crispr_draw_func import draw_design_stats
from crispr_design_gene_guides import design_guide_and_donor_for_genes


# setting the random seed
random.seed(1)

###############################################################################################################
# input parameters
###############################################################################################################

parser = argparse.ArgumentParser("Design crispr library")

parser.add_argument("lib_name", help= "lib_name can be pol2,pol3,noPol")
args = parser.parse_args()
lib_name = args.lib_name


###############################################################################################################
# input parameters
###############################################################################################################

guide_length = 20
edit_max_distance_from_PAM5prime = 17

PAM_seq = 'GG' # pam_seqs[0] 
# the code for desinging guides for SNPs supports multiple PAM sequences
# pam_seqs = [PAM_seq] # ["GG"]



#lib_name = 'pol2'
# lib_name = 'pol3'
#lib_name = 'noPol'

agilent_homopolymer_max_len = 10


excluded_seqs_pol2 = ['TATATA','TATGTA', 'A' * agilent_homopolymer_max_len, 'C' * agilent_homopolymer_max_len, 'G' * agilent_homopolymer_max_len, 'T' * agilent_homopolymer_max_len]
excluded_seqs_pol3 = ['TTTTT', 'A' * agilent_homopolymer_max_len, 'C' * agilent_homopolymer_max_len, 'G' * agilent_homopolymer_max_len, 'T' * agilent_homopolymer_max_len]
excluded_seqs_noPol = ['A' * agilent_homopolymer_max_len, 'C' * agilent_homopolymer_max_len, 'G' * agilent_homopolymer_max_len, 'T' * agilent_homopolymer_max_len]
excluded_seqs_bothPol = ['TTTTT','TATATA','TATGTA', 'A' * agilent_homopolymer_max_len, 'C' * agilent_homopolymer_max_len, 'G' * agilent_homopolymer_max_len, 'T' * agilent_homopolymer_max_len]

excluded_seqs_SNPs = []
excluded_seqs_Genes = []
if lib_name.startswith('pol2'):    # (lib_name in  ['pol2']):
    excluded_seqs_SNPs = excluded_seqs_pol2
    excluded_seqs_Genes = excluded_seqs_pol2
elif lib_name.startswith('pol3'):  #(lib_name in  ['pol3']):
    excluded_seqs_SNPs = excluded_seqs_pol3
    excluded_seqs_Genes = excluded_seqs_pol3
elif lib_name.startswith('noPol'):  #(lib_name in  ['noPol']):
    excluded_seqs_SNPs = excluded_seqs_noPol
    excluded_seqs_Genes = excluded_seqs_noPol
elif lib_name.startswith('bothPolForGeneNonForSNP'):  #(lib_name in  ['noPol']):
    excluded_seqs_SNPs = excluded_seqs_noPol
    excluded_seqs_Genes = excluded_seqs_bothPol
else:
    raise ValueError('Unknown excluded sequences for the library name: %s' % (lib_name))
   


donor_length = 100
min_dist_cut_to_donor_edge = 30

#BOWTIE_exe = "/Users/eilon/software/bowtie2-2.2.8/bowtie2"
BOWTIE_exe = "bowtie2"

design_guide_for_dubious_gene_in_qtls_flag = False
design_guide_for_essential_dubious_genes_flag = False
design_essential_genes_rand_seed = 14


min_ok_Azimuth_score_guides_in_genes = 0.5
min_ok_Azimuth_score_SNP_guides = 0

off_targets_min_mismatch_gene_guides = 10
off_targets_min_mismatch_SNP_guides = 1

max_num_of_essential_genes_to_select = 100
qtl_trait_for_essential_genes = ['YPD']


# for guides in essential genes
max_guide_pos_frac_in_gene_essential = 0.66

# for 1bpMut
#K_donors_for_each_guide_pos_1bpMut = 100
#max_guide_pos_frac_in_gene_1bpMut = 0.53
K_donors_for_each_guide_pos_1bpMut = 150
max_guide_pos_frac_in_gene_1bpMut = 0.6

oligo_production_batch_size = 2000

small_control_set_guide_num = 50 
small_control_set_editing_position = -3

###############################################################################################################

# directories
working_dir = "/science/eilon/crispr"
output_directory = working_dir + "/Output"
output_lib_dir = output_directory + "/lib_design"
output_fig_dir = output_lib_dir + "/Figures"

input_directory = working_dir + "/Input"
lib_design_input_dir = input_directory + "/lib_design"
input_Bloom_dir = input_directory + "/Bloom"

os.chdir(working_dir)
print "Current directory: %s" % os.getcwd()

###############################################################################################################

# input files

# from Carlo's analysis
input_by_rm_snps_vcf = working_dir + "/01-GATK_Bloom_RM_Data/GATK_RM_Vars_SCC30SEC10_LowQualFilt_141022.vcf"

# http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/
input_genome_fasta_filename = working_dir + "/01-GATK_Bloom_RM_Data/S288C_R64_REAN_140405.fa"


# annotations files are from # http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/ 
# version 61_1_1
# cat Input/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208.gff |  head -n 16425 > Input/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208_noFasta.gff
input_annot_gff_filename = input_directory + "/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208_noFasta.gff"


#####################################################
# input - essential genes
#####################################################
# esential genes used to be from:
# http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt
# essential_genes_filename = os.path.expanduser("~/Dropbox/postdoc/crispr/Input/essential/crispr_analyze_guides.txt")
# replaced by: http://www.yeastgenome.org/phenotype/inviable/overview

essential_genes_filename = input_directory +  "/essential/inviable_annotations.txt"


#####################################################
# input - Bloom QTL
#####################################################

input_known_genes_with_large_effect = input_directory + "/gene_with_known_large_effect.txt"

# selected 8 bloom QTL
input_bloom_selected_QTLs =  lib_design_input_dir + "/selected_bloom_traits.txt"

input_QTLs_filename   = input_Bloom_dir + "/nature11867-s4.csv"
input_QTLs_H_filename = input_Bloom_dir + "/nature11867-s2.txt"


#####################################################
# donor design for SNPs
#####################################################

# SNP guides donor design table 
input_SNP_donor_design_table_filename = lib_design_input_dir + "/design_donor_for_snps.txt"

#####################################################
# guide and donor design table for genes- selected QTL and essential genes in selected QTL 
#####################################################

# design table for QTL genes
input_qtl_gene_design_table_filename = lib_design_input_dir + "/design_guides_in_qtl_genes.txt"

# desing table for essential genes
input_guide_design_cols_for_essential_genes_filename = lib_design_input_dir + "/design_guides_in_essential_genes.txt"

# desing table for essential genes 1bpMut
input_guide_design_cols_for_essential_genes_1bpMut_filename = lib_design_input_dir + "/design_guides_in_essential_genes_1bpMut.txt"

# desing table for essential genes stop2aa
input_guide_design_cols_for_essential_genes_stop2aa_filename = lib_design_input_dir + "/design_guides_in_essential_genes_stop2aa.txt"


#####################################################
# oligo design table - contains information about the other segments in the oligo 
#####################################################

input_oligo_design_table_filename = lib_design_input_dir + "/oligo_design_table.txt"

#####################################################
# input output sets by priority.
# NOTICE - only sets within this list will be outputed
#####################################################
input_sets_by_priority_filename = lib_design_input_dir + "/sets_by_priority.txt"


###############################################################################################################
# output files (no need to pass as argument -  depends on the output directory, PAM sequence and max edit distance)
###############################################################################################################


output_files_uniq_str = lib_name + "_" + PAM_seq + "_" + str(edit_max_distance_from_PAM5prime) + "bp"


# log
output_log_stdout_filename = output_lib_dir + "/log_" + output_files_uniq_str + "_stdout.txt"
output_log_stderr_filename = output_lib_dir + "/log_" + output_files_uniq_str + "_stderr.txt"

################ genes in SNPs

# step: extract all guides for SNPs
output_SNP_table_filename =    output_lib_dir + "/all_SNPs_" + output_files_uniq_str + "_SNP.tab"
output_guides_table_filename = output_lib_dir + "/all_SNPs_" + output_files_uniq_str + "_GUIDE.tab"

# step: calculate features for SNPs
output_guides_with_features_table_filename = output_lib_dir + "/all_SNPs_" + output_files_uniq_str + "_GUIDE_withFeatures.tab"

# step: design SNP donor sequences
output_SNP_donor_table_filename = output_lib_dir + "/all_SNPs_" + output_files_uniq_str + "_DONOR.tab"

# step: rank guides for each SNP and add pass filter column
output_guides_with_features_and_rank_table_filename = output_lib_dir + "/all_SNPs_" + output_files_uniq_str + "_GUIDE_withFeatures_withRank.tab"


# step: merge SNPs with Bloom QTLs
output_selected_qtls_table_filename = output_lib_dir + "/selected_QTLs_" + output_files_uniq_str + ".tab"
output_SNP_merged_with_QTL_table_filename = output_lib_dir + "/all_SNPs_" + output_files_uniq_str + "_SNP_merged_with_QTL.tab"
output_SNP_withQTL_table_filename = output_lib_dir + "/all_SNPs_" + output_files_uniq_str + "_SNP_withQTL.tab"
output_snp_in_selected_QTLs_filename = output_lib_dir + "/SNPs_in_selected_QTLs_" + output_files_uniq_str + ".tab"

# step: create oligos from guide and donor tables for SNPs
output_oligos_for_SNPs_in_selected_qtls_df_filename = output_lib_dir + "/all_SNPs_in_selected_QTLs_" + output_files_uniq_str + "_OLIGO.tab"
output_oligos_for_SNPs_not_in_selected_qtls_df_filename = output_lib_dir + "/all_SNPs_not_in_selected_QTLs_" + output_files_uniq_str + "_OLIGO.tab"

################ genes in QTLs 


# step: generate design table for genes in QTLs
output_qtl_uniq_gene_guides_design_table_filename = output_lib_dir + "/design_guides_in_selected_qtl_genes_uniq_" + output_files_uniq_str + ".tab"
output_qtl_nonuniq_gene_guides_design_table_filename = output_lib_dir + "/design_guides_in_selected_qtl_genes_nonuniq_" + output_files_uniq_str + ".tab"
output_qtl_genes_that_were_filtered_out_filename = output_lib_dir + "/genes_in_qtls_that_were_filtered_out_" + output_files_uniq_str + ".tab"

# step: designing guides and donors for genes in QTLs
output_guides_for_genes_in_QTLs_df_filename = output_lib_dir + "/selected_QTL_genes_" + output_files_uniq_str + "_GUIDE.tab"
output_donors_for_genes_in_QTLs_df_filename = output_lib_dir + "/selected_QTL_genes_" + output_files_uniq_str + "_DONOR.tab"

# step: create oligos from guide and donor tables for Genes in QTL
output_oligos_for_genes_in_selected_qtls_df_filename= output_lib_dir + "/selected_QTL_genes_" + output_files_uniq_str + "_OLIGO.tab"

################ essential overlap with eql - for part of the edits

# step: create design table for guides and donors in essential genes by loading the genes and merging with the design table 
output_essential_gene_design_table_filename = output_lib_dir + "/essential_genes_design_table_" + output_files_uniq_str + ".tab"

# step designing guides and donors for genes in QTLs
output_guides_for_essential_genes_df_filename = output_lib_dir + "/essential_genes_" + output_files_uniq_str + "_GUIDE.tab"
output_donors_for_essential_genes_df_filename = output_lib_dir + "/essential_genes_" + output_files_uniq_str + "_DONOR.tab"
output_oligos_for_essential_genes_df_filename= output_lib_dir + "/essential_genes_" + output_files_uniq_str + "_OLIGO.tab"

# step: create design table for guides and donors in essential genes by loading the genes and merging with the design table 
output_essential_gene_1bpMut_design_table_filename = output_lib_dir + "/essential_genes_1bpMut_design_table.tab"

# step designing guides and donors for genes in QTLs
output_guides_for_essential_genes_1bpMut_df_filename = output_lib_dir + "/essential_genes_1bpMut_" + output_files_uniq_str + "_GUIDE.tab"
output_donors_for_essential_genes_1bpMut_unfiltered_df_filename = output_lib_dir + "/essential_genes_1bpMut_" + output_files_uniq_str + "_unfilt_DONOR.tab"
output_donors_for_essential_genes_1bpMut_df_filename = output_lib_dir + "/essential_genes_1bpMut_" + output_files_uniq_str + "_DONOR.tab"
output_oligos_for_essential_genes_1bpMut_df_filename= output_lib_dir + "/essential_genes_1bpMut_" + output_files_uniq_str + "_OLIGO.tab"

# step: create design table for guides and donors in essential genes by loading the genes and merging with the design table 
output_essential_gene_stop2aa_design_table_filename = output_lib_dir + "/essential_genes_stop2aa_design_table_" + output_files_uniq_str + ".tab"

# step designing 1bpMut guides and donors for genes 
output_guides_for_essential_genes_stop2aa_df_filename = output_lib_dir + "/essential_genes_stop2aa_" + output_files_uniq_str + "_GUIDE.tab"
output_donors_for_essential_genes_stop2aa_df_filename = output_lib_dir + "/essential_genes_stop2aa_" + output_files_uniq_str + "_DONOR.tab"
output_oligos_for_essential_genes_stop2aa_df_filename= output_lib_dir + "/essential_genes_stop2aa_" + output_files_uniq_str + "_OLIGO.tab"

# joining the oligos to one table

output_oligo_for_production_nonuniq_filename = output_lib_dir + "/by2rm_lib_nonuniq_" + output_files_uniq_str + "_OLIGO.txt"
output_oligo_for_production_nonuniq_with_align_filename = output_lib_dir + "/by2rm_lib_nonuniq_" + output_files_uniq_str + "_OLIGO.txt"

output_oligo_for_production_uniq_filename = output_lib_dir + "/by2rm_lib_uniq_" + output_files_uniq_str + "_OLIGO.txt"
output_oligo_for_production_uniq_batch_prefix_filename = output_lib_dir + "/by2rm_lib_uniq_" + output_files_uniq_str + "_"



################ draw stat

# step: draw SNP-guides statistics
output_snp_guides_fig_pref = output_fig_dir + "/by2rm_" + output_files_uniq_str + "_"



###############################################################################################################
# opening log files
###############################################################################################################

org_stdout = sys.stdout
org_stderr = sys.stderr

sys.stdout = open(output_log_stdout_filename, 'w')
sys.stderr = open(output_log_stderr_filename, 'w')

print "========================================================================="
print "=========================" + " Starting... " + "========================="
print "========================================================================="


# In[4]:

###############################################################################################################
# creating the output directories
###############################################################################################################

if not os.path.exists(output_lib_dir):
    os.makedirs(output_lib_dir)
    print "create directory:" + output_lib_dir
    
if not os.path.exists(output_fig_dir):
    os.makedirs(output_fig_dir)
    print "create directory:" + output_fig_dir



# In[5]:

###############################################################################################################
# extract guides for all SNPs
###############################################################################################################

extract_guides_for_snps(input_by_rm_snps_vcf, input_genome_fasta_filename, 
                        output_SNP_table_filename, output_guides_table_filename,
                        [PAM_seq], guide_length, edit_max_distance_from_PAM5prime)




# In[6]:

###############################################################################################################
# calculate guide features for all SNPs
###############################################################################################################

out_guides_with_features_df = cal_guide_features(output_guides_table_filename, input_genome_fasta_filename, 
                                                 output_guides_with_features_table_filename, 
                                                 PAM_seq, excluded_seqs_SNPs, donor_length,BOWTIE_exe = BOWTIE_exe)



# In[7]:

###############################################################################################################
# merge SNPs with Bloom QTLs
###############################################################################################################

snps_withQTLs_df = merge_SNPs_with_QTLs(output_SNP_table_filename, input_QTLs_filename, input_QTLs_H_filename,
                                        input_bloom_selected_QTLs, input_known_genes_with_large_effect,
                                        output_selected_qtls_table_filename, output_snp_in_selected_QTLs_filename,
                                        output_SNP_merged_with_QTL_table_filename, output_SNP_withQTL_table_filename)


# In[8]:

###############################################################################################################
# design donors for SNPs guides (for all SNPs)
###############################################################################################################

out_SNP_donor_df = design_donor_for_SNP_guides(input_SNP_donor_design_table_filename, 
                            output_SNP_withQTL_table_filename, output_guides_with_features_table_filename,
                            input_genome_fasta_filename,
                            donor_length, excluded_seqs_SNPs, min_dist_cut_to_donor_edge,
                            output_SNP_donor_table_filename)


# In[9]:

###############################################################################################################
# add filter and ranking to the SNP guides (depends on having a SNP and a donor tables) 
###############################################################################################################
out_SNP_guides_withFandR_df = rank_and_filter_SNP_guides(input_guides_with_features_table_filename = output_guides_with_features_table_filename,
                              input_SNP_table_filename = output_SNP_withQTL_table_filename,
                              input_donor_table_filename = output_SNP_donor_table_filename,
                              output_guides_with_features_and_rank_table_filename = output_guides_with_features_and_rank_table_filename,
                              off_targets_min_mismatch_SNP_guides = off_targets_min_mismatch_SNP_guides, 
                              min_ok_Azimuth_score_SNP_guides = min_ok_Azimuth_score_SNP_guides, 
                              edit_max_distance_from_PAM5prime = edit_max_distance_from_PAM5prime)


# In[10]:


###############################################################################################################
# design oligos for SNPs guides
###############################################################################################################


# TODO the filter parameter is not used - using the filter below 

snp_df = pd.read_table(output_SNP_withQTL_table_filename, sep='\t', na_values = "")
snp_guides_df = pd.read_table(output_guides_with_features_and_rank_table_filename, sep='\t', na_values = "")
snp_donors_df = pd.read_table(output_SNP_donor_table_filename, sep='\t', na_values = "")

############### indices for SNPs in selected QTLs
guides_in_sel_qtl_iloc_vec = np.array(range(snp_guides_df.shape[0]))
#guides_in_sel_qtl_iloc_vec = guides_in_sel_qtl_iloc_vec[ (~snp_guides_df['contain_excluded_sequences'].iloc[guides_in_sel_qtl_iloc_vec]).values ]

guides_in_sel_qtl_iloc_vec = guides_in_sel_qtl_iloc_vec[ (snp_guides_df['is_nuclear_chromosome'].iloc[guides_in_sel_qtl_iloc_vec]).values ]
guides_in_sel_qtl_iloc_vec = guides_in_sel_qtl_iloc_vec[ (snp_guides_df['guide_map_mismatch_0'].iloc[guides_in_sel_qtl_iloc_vec] < 1).values ]

guides_in_sel_qtl_iloc_vec = guides_in_sel_qtl_iloc_vec[ (snp_guides_df['SNP_id'].iloc[guides_in_sel_qtl_iloc_vec].isin(
        snp_df['SNP_id'][ snp_df['set_name'] ==  'selected_qtls']) ).values ]

############## indices for SNPs not in selected QTLs
guides_not_in_sel_qtl_iloc_vec = np.array(range(snp_guides_df.shape[0]))
#guides_not_in_sel_qtl_iloc_vec = guides_not_in_sel_qtl_iloc_vec[ (~snp_guides_df['contain_excluded_sequences'].iloc[guides_not_in_sel_qtl_iloc_vec]).values ]

guides_not_in_sel_qtl_iloc_vec = guides_not_in_sel_qtl_iloc_vec[ (snp_guides_df['is_nuclear_chromosome'].iloc[guides_not_in_sel_qtl_iloc_vec]).values ]
guides_not_in_sel_qtl_iloc_vec = guides_not_in_sel_qtl_iloc_vec[ (snp_guides_df['guide_map_mismatch_0'].iloc[guides_not_in_sel_qtl_iloc_vec] < 1).values ]
guides_not_in_sel_qtl_iloc_vec = guides_not_in_sel_qtl_iloc_vec[ (snp_guides_df['guide_map_mismatch_1'].iloc[guides_not_in_sel_qtl_iloc_vec] < 1).values ]
guides_not_in_sel_qtl_iloc_vec = guides_not_in_sel_qtl_iloc_vec[ (snp_guides_df['guide_map_mismatch_2'].iloc[guides_not_in_sel_qtl_iloc_vec] < 1).values ]
guides_not_in_sel_qtl_iloc_vec = guides_not_in_sel_qtl_iloc_vec[ (snp_guides_df['guide_map_mismatch_3'].iloc[guides_not_in_sel_qtl_iloc_vec] < 100).values ]


# one guide per SNP
guides_not_in_sel_qtl_iloc_vec = guides_not_in_sel_qtl_iloc_vec[ (snp_guides_df['quality_rank'].iloc[guides_not_in_sel_qtl_iloc_vec] == 0).values ]

guides_not_in_sel_qtl_iloc_vec = guides_not_in_sel_qtl_iloc_vec[ (snp_guides_df['SNP_id'].iloc[guides_not_in_sel_qtl_iloc_vec].isin(
        snp_df['SNP_id'][ ( snp_df['set_name'] ==  'not_in_sel_qtls').values]) ).values ]

# donors optimaized to exclude excluded sequences
donors_iloc_vec = np.array(range(snp_donors_df.shape[0]))
donors_iloc_vec = donors_iloc_vec[ (snp_donors_df['set_name'].iloc[donors_iloc_vec] == 'SNP_avoidExcludedSeqs').values ]



# In[11]:

# desgining oligos for SNPs in selected QTLs
SNPs_in_sel_qtls_oligo_df = generate_oligo_from_guide_and_donors(input_oligo_design_table_filename,
                                     output_guides_with_features_and_rank_table_filename, 
                                     output_SNP_donor_table_filename,
                                     input_guide_iloc = guides_in_sel_qtl_iloc_vec, input_donor_iloc = donors_iloc_vec,
                                     input_SNP_table_filename = output_SNP_withQTL_table_filename,                           
                                     output_oligos_table_filename = output_oligos_for_SNPs_in_selected_qtls_df_filename)

# designing oligos for SNPs not in selected QTLs
SNPs_not_in_sel_qtls_oligo_df = generate_oligo_from_guide_and_donors(input_oligo_design_table_filename,
                                     output_guides_with_features_and_rank_table_filename, 
                                     output_SNP_donor_table_filename,
                                     input_guide_iloc = guides_not_in_sel_qtl_iloc_vec, input_donor_iloc = donors_iloc_vec,
                                     input_SNP_table_filename = output_SNP_withQTL_table_filename,
                                     output_oligos_table_filename = output_oligos_for_SNPs_not_in_selected_qtls_df_filename)

# In[14]:

###############################################################################################################
# designing guides that target genes
###############################################################################################################


# In[15]:

###############################################################################################################
# generate design table for genes in selected QTLs
###############################################################################################################

uniq_qtl_genes_design_df, nonuniq_qtl_genes_design_df = generate_qtl_gene_guide_design_table(output_selected_qtls_table_filename, input_qtl_gene_design_table_filename,
                                         input_annot_gff_filename,
                                         design_guide_for_dubious_gene_in_qtls_flag,
                                         output_qtl_uniq_gene_guides_design_table_filename, 
                                         output_qtl_nonuniq_gene_guides_design_table_filename,
                                         output_qtl_genes_that_were_filtered_out_filename)



# In[5]:

###############################################################################################################
# design guides and donor for genes in selected QTLs
###############################################################################################################

# takes about 2 hours to run on my laptop (25 genes per minutes)

# designing guides and donor sequences for genes in QTLs
out_selected_qtl_genes_guide_df, out_selected_qtl_genes_guide_donor_df = design_guide_and_donor_for_genes(output_qtl_uniq_gene_guides_design_table_filename,
                                     input_annot_gff_filename, input_genome_fasta_filename,
                                     PAM_seq = PAM_seq, excluded_seqs = excluded_seqs_Genes, BOWTIE_exe = BOWTIE_exe,
                                     off_targets_min_mismatch = off_targets_min_mismatch_gene_guides,
                                     min_azimuth_score = min_ok_Azimuth_score_guides_in_genes, sort_by = '5to3',
                                     min_dist_cut_to_donor_edge = 30,
                                     filter_out_donors_containing_excluded_seqs = True,
                                     max_guide_pos_frac_in_gene = 1.1,
                                     output_guides_df_filename = output_guides_for_genes_in_QTLs_df_filename,
                                     output_donors_df_filename = output_donors_for_genes_in_QTLs_df_filename)



# In[6]:

###############################################################################################################
# designing oligos for genes in QTLs
###############################################################################################################
genes_in_sel_QTLs_oligo_df = generate_oligo_from_guide_and_donors(input_oligo_design_table_filename,
                                     output_guides_for_genes_in_QTLs_df_filename, 
                                     output_donors_for_genes_in_QTLs_df_filename,
                                     input_guide_iloc = None, input_donor_iloc = None,
                                     input_SNP_table_filename = None,
                                     output_oligos_table_filename = output_oligos_for_genes_in_selected_qtls_df_filename)



# In[ ]:

###############################################################################################################
# Designing guides in essential genes to examine differnt aspects of the technology
###############################################################################################################


# In[5]:

###############################################################################################################
# Creating the design table for the essential genes 
###############################################################################################################

# setting the random seed
random.seed(1)
#np.random.seed(10)
np.random.seed(11)

essential_genes_design_df = parse_essential_genes_to_design_table(essential_genes_filename, 
                                                    input_annot_gff_filename, 
                                                    output_selected_qtls_table_filename,
                                                    qtl_trait_for_essential_genes,
                                                    design_guide_for_essential_dubious_genes_flag,
                                                    max_num_of_essential_genes_to_select,
                                                    True,
                                                    input_guide_design_cols_for_essential_genes_filename,
                                                    output_essential_gene_design_table_filename,
                                                    design_essential_genes_rand_seed)


# In[8]:

###############################################################################################################
# design guides and donor for genes in essential genes
###############################################################################################################

# takes about 2 hours to run on my laptop (25 genes per minutes)

# designing guides and donor sequences for genes in QTLs
out_essential_genes_guide_df, out_essential_genes_guide_donor_df = design_guide_and_donor_for_genes(output_essential_gene_design_table_filename,
                                     input_annot_gff_filename, input_genome_fasta_filename,
                                     PAM_seq = PAM_seq, excluded_seqs = excluded_seqs_Genes, BOWTIE_exe = BOWTIE_exe,
                                     off_targets_min_mismatch = off_targets_min_mismatch_gene_guides, 
                                     min_azimuth_score = min_ok_Azimuth_score_guides_in_genes, sort_by = '5to3',
                                     min_dist_cut_to_donor_edge = min_dist_cut_to_donor_edge,
                                     filter_out_donors_containing_excluded_seqs = True,
                                     max_guide_pos_frac_in_gene = max_guide_pos_frac_in_gene_essential,
                                     output_guides_df_filename = output_guides_for_essential_genes_df_filename,
                                     output_donors_df_filename = output_donors_for_essential_genes_df_filename)


# In[7]:

###############################################################################################################
# generating oligos from guides and donor that target genes
# can get iloc vector for which rows to include in the design
###############################################################################################################

# desgining oligos for essential genes 
essential_oligo_df = generate_oligo_from_guide_and_donors(input_oligo_design_table_filename,
                                     output_guides_for_essential_genes_df_filename, 
                                     output_donors_for_essential_genes_df_filename,
                                     input_guide_iloc = None, input_donor_iloc = None,
                                     input_SNP_table_filename = None,
                                     output_oligos_table_filename = output_oligos_for_essential_genes_df_filename)




# In[6]:

###############################################################################################################
# Creating the design table for the essential genes 1bp Mutations 
# In this table the genes are not filtered to be in QTLs
###############################################################################################################

essential_genes_1bpMut_design_df =  parse_essential_genes_to_design_table(essential_genes_filename, 
                                    input_annot_gff_filename, 
                                    output_selected_qtls_table_filename,
                                    None,
                                    design_guide_for_essential_dubious_genes_flag,
                                    1000000,
                                    False,
                                    input_guide_design_cols_for_essential_genes_1bpMut_filename,
                                    output_essential_gene_1bpMut_design_table_filename,
                                    design_essential_genes_rand_seed)


# In[5]:

###############################################################################################################
# design 1bpMut guides and donor for genes in essential genes
###############################################################################################################

# takes about 2 hours to run on my laptop (25 genes per minutes)

# designing guides and donor sequences for genes in QTLs
out_essential_genes_1bpMut_guide_df, out_essential_genes_1bpMut_guide_donor_df = design_guide_and_donor_for_genes(output_essential_gene_1bpMut_design_table_filename,
                                     input_annot_gff_filename, input_genome_fasta_filename,
                                     PAM_seq = PAM_seq, excluded_seqs = excluded_seqs_Genes, BOWTIE_exe = BOWTIE_exe,
                                     off_targets_min_mismatch = off_targets_min_mismatch_gene_guides, 
                                     min_azimuth_score = min_ok_Azimuth_score_guides_in_genes, sort_by = '5to3',
                                     min_dist_cut_to_donor_edge = min_dist_cut_to_donor_edge,
                                     filter_out_donors_containing_excluded_seqs = True,
                                     max_guide_pos_frac_in_gene = max_guide_pos_frac_in_gene_1bpMut,
                                     output_guides_df_filename = output_guides_for_essential_genes_1bpMut_df_filename,
                                     output_donors_df_filename = output_donors_for_essential_genes_1bpMut_unfiltered_df_filename)


###############################################################################################################
# filter 1bpMut donors such that the will be 
###############################################################################################################

out_essential_genes_1bpMut_donor_filt_df = \
        filter_1pbMut_donor_table(output_guides_for_essential_genes_1bpMut_df_filename, 
                                  output_donors_for_essential_genes_1bpMut_unfiltered_df_filename,
                                  K_donors_for_each_guide_pos = K_donors_for_each_guide_pos_1bpMut, 
                                  min_Azimuth_score = min_ok_Azimuth_score_guides_in_genes,
                                  max_gene_pos_frac = max_guide_pos_frac_in_gene_1bpMut,
                                  output_donors_for_essential_genes_1bpMut_filt_df_filename = output_donors_for_essential_genes_1bpMut_df_filename)




###############################################################################################################
# generating oligos from guides and donor that target genes
# can get iloc vector for which rows to include in the design
###############################################################################################################

# desgining oligos for essential genes 
essential_oligo_df = generate_oligo_from_guide_and_donors(input_oligo_design_table_filename,
                                     output_guides_for_essential_genes_1bpMut_df_filename, 
                                     output_donors_for_essential_genes_1bpMut_df_filename,
                                     input_guide_iloc = None, input_donor_iloc = None,
                                     input_SNP_table_filename = None,
                                     output_oligos_table_filename = output_oligos_for_essential_genes_1bpMut_df_filename)



###############################################################################################################
# Creating the design table for the essential genes stop2aa Mutations 
# In this table the genes are not filtered to be in QTLs
###############################################################################################################

essential_genes_stop2aa_design_df = \
        parse_essential_genes_to_design_table(essential_genes_filename, 
                                              input_annot_gff_filename, 
                                              output_selected_qtls_table_filename,
                                              None,
                                              design_guide_for_essential_dubious_genes_flag,
                                              1000000,
                                              False,
                                              input_guide_design_cols_for_essential_genes_stop2aa_filename,
                                              output_essential_gene_stop2aa_design_table_filename,
                                              design_essential_genes_rand_seed)

###############################################################################################################
# design stop2aa guides and donor for genes in essential genes
###############################################################################################################

# takes about 2 hours to run on my laptop (25 genes per minutes)

# designing guides and donor sequences for genes in QTLs
out_essential_genes_stop2aa_guide_df, out_essential_genes_stop2aa_donor_df = \
    design_guide_and_donor_for_genes(output_essential_gene_stop2aa_design_table_filename,
                                     input_annot_gff_filename, input_genome_fasta_filename,
                                     PAM_seq = PAM_seq, excluded_seqs = excluded_seqs_Genes, BOWTIE_exe = BOWTIE_exe,
                                     off_targets_min_mismatch = off_targets_min_mismatch_gene_guides,
                                     min_azimuth_score = min_ok_Azimuth_score_guides_in_genes, sort_by = 'onlyStopCodon',
                                     min_dist_cut_to_donor_edge = min_dist_cut_to_donor_edge,
                                     filter_out_donors_containing_excluded_seqs = True,
                                     max_guide_pos_frac_in_gene = 1.1,
                                     output_guides_df_filename = output_guides_for_essential_genes_stop2aa_df_filename,
                                     output_donors_df_filename = output_donors_for_essential_genes_stop2aa_df_filename)


###############################################################################################################
# generating oligos from guides and donor that target genes
# can get iloc vector for which rows to include in the design
###############################################################################################################

# desgining oligos for essential genes 
essential_oligo_df = generate_oligo_from_guide_and_donors(input_oligo_design_table_filename,
                                     output_guides_for_essential_genes_stop2aa_df_filename, 
                                     output_donors_for_essential_genes_stop2aa_df_filename,
                                     input_guide_iloc = None, input_donor_iloc = None,
                                     input_SNP_table_filename = None,
                                     output_oligos_table_filename = output_oligos_for_essential_genes_stop2aa_df_filename)



###############################################################################################################
# concating oligo tables
###############################################################################################################

# loading the tables

oligo_snps_sel_df = pd.read_table(output_oligos_for_SNPs_in_selected_qtls_df_filename, sep='\t', na_values = "None")
oligo_snps_notSel_df = pd.read_table(output_oligos_for_SNPs_not_in_selected_qtls_df_filename, sep='\t', na_values = "None")
oligo_qtl_genes_df = pd.read_table(output_oligos_for_genes_in_selected_qtls_df_filename, sep='\t', na_values = "None")
oligo_ess_genes_df = pd.read_table(output_oligos_for_essential_genes_df_filename, sep='\t', na_values = "None")
oligo_1bpMut_genes_df = pd.read_table(output_oligos_for_essential_genes_1bpMut_df_filename, sep='\t', na_values = "None")
oligo_stop2aa_genes_df = pd.read_table(output_oligos_for_essential_genes_stop2aa_df_filename, sep='\t', na_values = "None")


#print oligo_snps_sel_df.columns
#print oligo_snps_notSel_df.columns
#print oligo_qtl_genes_df.columns
#print oligo_ess_genes_df.columns
#print oligo_1bpMut_genes_df.columns

oligo_colnames = ['set_name','guide_id', 'donor_id','oligo_id','donor_info_str','guide_noPAM','donor_seq','oligo_seq']

oligo_dfs = [oligo_snps_sel_df[oligo_colnames],
             oligo_snps_notSel_df[oligo_colnames],
             oligo_qtl_genes_df[oligo_colnames],
             oligo_ess_genes_df[oligo_colnames],
             oligo_1bpMut_genes_df[oligo_colnames],
             oligo_stop2aa_genes_df[oligo_colnames] ]


oligo_all_df = pd.concat(oligo_dfs, ignore_index = True)



#plt.figure();
set_sizes_df =  oligo_all_df['set_name'].value_counts().sort_index()
#set_sizes_df.plot.bar()

print "set sizes:"
print set_sizes_df
print "total # oligos:"
print oligo_all_df.shape[0]



oligo_for_order_dfs = write_output_oligos(oligo_all_df, 
                    input_sets_by_priority_filename,
                    batch_size = oligo_production_batch_size,
                    output_oligo_for_production_nonuniq_filename = output_oligo_for_production_nonuniq_filename,
                    output_oligo_for_production_uniq_filename = output_oligo_for_production_uniq_filename,
                    output_oligo_for_production_uniq_batch_prefix_filename = output_oligo_for_production_uniq_batch_prefix_filename)


print "# unique oligos:"
print oligo_for_order_dfs.shape[0]

#print oligo_all_df.head(3)

write_small_control_oligo_set(output_oligo_for_production_nonuniq_filename, 
                                output_oligo_for_production_uniq_filename, 
                                output_oligos_for_essential_genes_1bpMut_df_filename, 
                                small_control_set_guide_num,
                                small_control_set_editing_position)

print "--------------- Finished designing the library! -------------------"


# In[7]:
###############################################################################################################
# align guide to the donor to make sure the strand matching is correct
###############################################################################################################

print "# aligning guide and donor:"

def align_guide_and_donor(row):
    cur_align = pairwise2.align.localxs(row['guide_noPAM'],row['donor_seq'],-1,-1)
    if len(cur_align)>0:
        return(max(cur_align)[2])
    else:
        return(-1)
    
def align_revcomp_guide_and_donor(row):
    cur_align = pairwise2.align.localxs( str(Seq(row['guide_noPAM'],alphabet=generic_dna).reverse_complement())  ,row['donor_seq'],-1,-1)
    if len(cur_align)>0:
        return(max(cur_align)[2])
    else:
        return(-1)  


oligo_all_df['guide_donor_align_score'] = oligo_all_df.apply(align_guide_and_donor, axis=1)
oligo_all_df['revcomp_guide_donor_align_score'] = oligo_all_df.apply(align_revcomp_guide_and_donor, axis=1)
oligo_all_df['guide_and_donor_align_better_than_revcomp'] = oligo_all_df['guide_donor_align_score'] >= oligo_all_df['revcomp_guide_donor_align_score']


if output_oligo_for_production_nonuniq_with_align_filename:
    print "Saving : " + output_oligo_for_production_nonuniq_with_align_filename
    oligo_all_df.to_csv(output_oligo_for_production_nonuniq_with_align_filename, sep='\t', index = False)
    







###############################################################################################################
# draw SNP guide general statistics
###############################################################################################################

print "--------------- Drawing -------------------"

draw_design_stats(output_SNP_withQTL_table_filename, output_guides_with_features_table_filename, 
                  output_SNP_donor_table_filename, output_snp_guides_fig_pref)



# In[6]:

###############################################################################################################
# closing log files 
###############################################################################################################


print "========================================================="
print "=========================" + " Done! " + "========================="
print "========================================================="

sys.stdout = org_stdout
sys.stderr = org_stderr

print "========================================================="
print "=========================" + " Done! " + "========================="
print "========================================================="


# In[5]:



