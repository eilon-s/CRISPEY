#!/usr/bin/env python2

from __future__ import division

import sys
import os
import itertools
  
  
import math
import re
import random

import numpy as np
import pandas as pd



import matplotlib.pyplot as plt
from scipy.stats import kendalltau
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from scipy.stats import gaussian_kde

import seaborn as sns
sns.set(rc={"figure.figsize": (6, 6)})

sys.path.append(os.path.expanduser('~/bin/python/crispr'))
from extract_guides_functions import *


###############################################################################################################
# SNP - guide stats 
###############################################################################################################

###############################################################################################################
def draw_design_stats(input_SNP_table_filename,
                      input_guides_with_features_table_filename,
                      input_donors_table_filename,
                      output_filenames_prefix):
    
    print "in draw_design_stats loading tables"
    
    snps_df = pd.DataFrame(data=None)
    guides_df = pd.DataFrame(data=None)
    donors_df = pd.DataFrame(data=None)
    
    # load SNPs
    if input_SNP_table_filename:
      snps_df = pd.read_table(input_SNP_table_filename, sep='\t', na_values = "None")

    # load guides
    if input_guides_with_features_table_filename:
      guides_df = pd.read_table(input_guides_with_features_table_filename, sep='\t', na_values = "None")
    
    # load donors
    if input_donors_table_filename:
      donors_df = pd.read_table(input_donors_table_filename, sep='\t', na_values = "None")
    
    print "finished loading tables"
    
    
    ################################################################################
    # draw figures
    ################################################################################
    
    print "Start drawing plots..."
    
    draw_num_guides_per_SNP_fromGuideTable(guides_df, output_filenames_prefix)
    
    draw_num_guides_per_SNP_fromSNPtable(snps_df, output_filenames_prefix)
    
    draw_num_strands_with_guides_per_SNP(snps_df, output_filenames_prefix)
    
    draw_SNP_pos_in_best_guide_fromSNPtable(snps_df, output_filenames_prefix)
    
    draw_SNP_best_pos_guide_Azimuth_score(snps_df, guides_df, output_filenames_prefix)
    
    draw_guides_offtragets(guides_df, output_filenames_prefix)
    
    print "Done drawing plots."

###############################################################################################################
def draw_num_guides_per_SNP_fromGuideTable(guides_df, output_filenames_prefix):
    
    ################################################################################
    # Ploting statistics of the number of guides per SNP (not including Zero)
    ################################################################################
    
    guides_per_snp_df =  guides_df["SNP_id"].value_counts().to_frame("GuideCounts")


    # cdf of guides number
    ser = guides_per_snp_df['GuideCounts']
    ser = ser.sort_values(ascending=False)
    cum_dist = np.linspace(0.,1.,len(ser))
    ser_cdf = pd.Series(cum_dist, index=ser)
    sns.set(style="whitegrid", palette="deep", color_codes=True)

    fig = plt.figure(figsize=(6,6))
    ax = ser_cdf.plot(drawstyle='steps', figsize=(6,6), grid=True,xlim=(15,0))
    ax.set(xlabel='# guides', ylabel='%SNPs')
    fig.savefig( output_filenames_prefix + 'num_guides_per_snp_cdf_v1.pdf', bbox_inches='tight')

    # histograms
    with sns.axes_style("whitegrid"):
        fig = plt.figure(figsize=(6,6))
        ax = sns.distplot(guides_per_snp_df['GuideCounts'], kde=False, color='black')
        ax.set(xlabel='# guides', ylabel='#SNPs')
        fig.savefig( output_filenames_prefix + 'num_guides_per_snp.pdf', bbox_inches='tight')

        fig = plt.figure(figsize=(6,6))
        ax = sns.distplot(guides_per_snp_df['GuideCounts'], kde=False, norm_hist=True, bins=np.arange(0.5,15.5,1))
        ax.set(xlabel='# guides', ylabel='%SNPs')
        fig.savefig( output_filenames_prefix + 'num_guides_per_snp_dist.pdf', bbox_inches='tight')

###############################################################################################################
def draw_num_guides_per_SNP_fromSNPtable(snps_df, output_filenames_prefix):
    ################################################################################
    # number of guides per SNP (including zero)
    ################################################################################

    ser = snps_df['num_possible_guides']
    ser = ser.sort_values(ascending=False) #ascending=False
    cum_dist = np.linspace(0.,1.,len(ser))
    ser_cdf = pd.Series(cum_dist, index=ser)
    sns.set(style="whitegrid", palette="deep", color_codes=True)
    fig = plt.figure(figsize=(6,6))
    ax = ser_cdf.plot(drawstyle='steps', figsize=(6,6), grid=True,xlim=(15,0))
    ax.set(xlabel='# guides', ylabel='%SNPs')
    fig.savefig( output_filenames_prefix + 'num_guides_per_snp_cdf_includingZeros.pdf', bbox_inches='tight')

###############################################################################################################
def draw_num_strands_with_guides_per_SNP(snps_df, output_filenames_prefix):
    ################################################################################
    # number of SNPs with guides on both strands (histogram)
    ################################################################################

    snps_df['num_strands_with_guides'] = snps_df['pos_in_guide_downstream'].notnull().astype(int) + snps_df['pos_in_guide_upstream'].notnull().astype(int)

    print "Fraction of SNPs that have no guide" + str(sum(snps_df['num_strands_with_guides']==0)/len(snps_df['num_strands_with_guides']))
    print "Fraction of SNPs that have guides in one strand" + str(sum(snps_df['num_strands_with_guides']==1)/len(snps_df['num_strands_with_guides']))
    print "Fraction of SNPs that have guides in two strands" + str(sum(snps_df['num_strands_with_guides']==2)/len(snps_df['num_strands_with_guides']))

    with sns.axes_style("whitegrid"):
        fig = plt.figure(figsize=(6,6))
        ax = sns.distplot(snps_df['num_strands_with_guides'], kde=False, norm_hist=True, bins=np.arange(-0.5,3.5,1))
        ax.set(xlabel='# strands ', ylabel='%SNPs')
        fig.savefig( output_filenames_prefix + 'num_strands_with_guides.pdf', bbox_inches='tight')

###############################################################################################################
def draw_guides_offtragets(guides_df, output_filenames_prefix):
  
  guides_df['does_have_off_tragets'] = (guides_df['guide_map_mismatch_0'] != 0) | \
                                     (guides_df['guide_map_mismatch_1'] != 0) | \
                                     (guides_df['guide_map_mismatch_2'] != 0) | \
                                     (guides_df['guide_map_mismatch_3'] != 0) 
  

  with sns.axes_style("whitegrid"):
      fig = plt.figure(figsize=(6,6))
      ax = sns.distplot(guides_df['does_have_off_tragets'], kde=False)
      ax.set(xlabel='# strands ', ylabel='%guides')
      fig.savefig( output_filenames_prefix + 'num_off_tragets_cnt.pdf', bbox_inches='tight')
      
  summary_df = pd.DataFrame({'DoesHaveOffTargets_cnt' : [100*sum(~guides_df['does_have_off_tragets'])/guides_df.shape[0], 
                                                       100*sum(guides_df['does_have_off_tragets'])/guides_df.shape[0]],
                           'DoesHaveOffTargets' : ['No', 'Yes'] } , index = ['No','Yes'] )


  #ax = fig.add_axes([0.07,0.12,0.4,0.8])
  fig, ax = plt.subplots(1, 1, figsize=(6, 6))
  summary_df.plot.bar(ax=ax)
  ax.set(xlabel='Does have off targets', ylabel='% guides')
  ax.legend_.remove()
  fig.savefig( output_filenames_prefix + 'num_off_tragets.pdf', bbox_inches='tight')

###############################################################################################################
def draw_SNP_pos_in_best_guide_fromSNPtable(snps_df, output_filenames_prefix):

    ################################################################################
    # minimal (maximal of negative) position of SNP in guide over two sides
    ################################################################################

    snps_df['best_pos_in_guide_both'] =  snps_df[ ['pos_in_guide_downstream','pos_in_guide_upstream'] ].max(axis=1)

    with sns.axes_style("whitegrid"):
        fig = plt.figure(figsize=(6,6))
        ax = sns.distplot(snps_df['best_pos_in_guide_both'].dropna(), kde=False, norm_hist=True, bins=np.arange(-18.5,3.5,1), 
                          hist_kws={"alpha": 0.8, "color": sns.xkcd_rgb["denim blue"]}) 
        ax.set(xlabel='SNP position in guide', ylabel='%SNPs')
        fig.savefig( output_filenames_prefix + 'guide_best_pos.pdf', bbox_inches='tight')

        # number of guides per SNP (including zero)
        ser = snps_df['best_pos_in_guide_both']
        ser = ser.sort_values(ascending=False) #
        cum_dist = np.linspace(0.,1.,len(ser))
        ser_cdf = pd.Series(cum_dist, index=ser)

        fig = plt.figure(figsize=(6,6))
        ax = ser_cdf.plot(drawstyle='steps', figsize=(6,6), grid=True,xlim=(5,-20),color=sns.xkcd_rgb["denim blue"])
        ax.set(xlabel='SNP position in guide', ylabel='%SNPs')
        fig.savefig( output_filenames_prefix + 'guide_best_pos_cdf.pdf', bbox_inches='tight')

###############################################################################################################
def draw_SNP_best_pos_guide_Azimuth_score(snps_df, guides_df, output_filenames_prefix):

    ################################################################################
    # distribution of scores of best position guide
    ################################################################################

    # adding columns to the SNP tables to allow easy merging with the guide table
    snps_df['is_downstream_pos_better'] = (snps_df['pos_in_guide_downstream'] >= snps_df['pos_in_guide_upstream']) | (snps_df['pos_in_guide_upstream'].isnull())
    snps_df['best_guide_strand'] = np.where((snps_df['pos_in_guide_downstream'] >= snps_df['pos_in_guide_upstream']) | (snps_df['pos_in_guide_upstream'].isnull()),
                                            'downstream',
                                            'upstream')
    best_pos_guides = pd.merge(snps_df, guides_df, left_on = ['SNP_id','best_guide_strand', 'best_pos_in_guide_both'], right_on = ['SNP_id', 'upstream_or_downstream','SNP_pos_in_guide'])

    # distributions of scores of the guides with the best positions
    with sns.axes_style("whitegrid"):
        fig = plt.figure(figsize=(6,6))
        ax = sns.kdeplot(best_pos_guides['Azimuth'], shade=True, color="r", cumulative=True, legend=False)
        ax.set(xlabel='Azimuth', ylabel='SNPs (fraction)')
        fig.savefig( output_filenames_prefix + 'bestPos_AzimuthScore.pdf', bbox_inches='tight')
        
    
    ################################################################################
    # best score per SNP in OK [-18,2] position
    ################################################################################

    # max azimuth
    idx = guides_df.groupby(['SNP_id'])['Azimuth'].transform(max) == guides_df['Azimuth']
    guides_best_azimuth = guides_df[idx]
    # max position
    idx = guides_best_azimuth.groupby(['SNP_id'])['SNP_pos_in_guide'].transform(max) == guides_best_azimuth['SNP_pos_in_guide']
    guides_best_azimuth = guides_best_azimuth[idx]


    with sns.axes_style("whitegrid"):
        fig = plt.figure(figsize=(6,6))
        ax = sns.kdeplot(guides_best_azimuth['Azimuth'], shade=True, color="r", cumulative=True, legend=False)
        ax.set(xlabel='Azimuth', ylabel='SNPs (fraction)')
        fig.savefig( output_filenames_prefix + 'bestAzimuth_AzimuthScore.pdf', bbox_inches='tight')


    with sns.axes_style("whitegrid"):
        fig = plt.figure(figsize=(6,6))
        ax = sns.kdeplot(best_pos_guides['Azimuth'], shade=True, color="r", cumulative=True, legend=False)
        ax = sns.kdeplot(guides_best_azimuth['Azimuth'], shade=True, color="b", cumulative=True, legend=False)
        ax.set(xlabel='Azimuth', ylabel='SNPs (fraction)')
        fig.savefig( output_filenames_prefix + 'bestPosOrAzimuth_AzimuthScore.pdf', bbox_inches='tight')

###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################


###############################################################################################################

