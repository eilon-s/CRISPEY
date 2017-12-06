#!/usr/bin/env python2

from __future__ import division

import sys
import os
import warnings
import itertools

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

import numpy as np
import math
import re

import pandas as pd
import random

import ast

from roman import fromRoman

# Azimuth
sys.path.append(os.path.expanduser("~/software/Azimuth"))
import azimuth.model_comparison

sys.path.append(os.path.expanduser('~/bin/python/crispr'))
from gff2dataframe import gff2dataframe, sgd_gff2dataframe
from crispr_codon import *
from extract_guides_functions import *

####################################################################################################################
class CRISPR_Gene(object):
    """Gene hold data and do operation for a single coding gene (may have several coding exons)"""
    
    
    ###############################################################################################################
    def __init__(self, gene_name, gene_gff_df, genome_seq):
        """
        Initialize CRISPR_Gene object.
        
        Parameters:
        -----------
        gene_name : SGD name (column name in genes gff df - 'Name')
        gene_gff_df:  df of loaded gene annotation file (gff; should be filtered and chromosome fixed)
        genome_seq: genome sequences dictionary
        
        
        Attributes:
        -----------
        name : str
            gene name
        chrom : str
            gene chromosome
        strand : str ['+'|'-'] 
            gene strand
        
        coordinates are python style zero based and not including the last position
        gene_left : int
        gene_right : int
        
        cds_lefts : np.array
        cds_rights : np.array
        
        genome_seq : dict
        gene_gff_df :  pandas data frame 
        
        """
        
        # saving input parameters (pointers)
        self.name = gene_name
        self.genome_seq = genome_seq
        self.gene_gff_df = gene_gff_df
        
        this_gene_gff_df = self.get_gene_gff_df()
        
        self.strand = this_gene_gff_df['strand'].iloc[0]
        self.chrom = this_gene_gff_df['chrom'].iloc[0]
        
        self.cds_lefts = np.sort(this_gene_gff_df['start'][this_gene_gff_df['feature'].isin(['CDS'])].values.astype(int))-1
        #if (self.is_neg_strand()):
        #    self.cds_lefts = self.cds_lefts[::-1]
        
        self.cds_rights = np.sort(this_gene_gff_df['end'][this_gene_gff_df['feature'].isin(['CDS'])].values.astype(int))#-1
        #if (self.is_neg_strand()):
        #    self.cds_rights = self.cds_rights[::-1]
        
        self.cds_length_nts = sum(self.cds_rights - self.cds_lefts)
        
        self.gene_left = np.min(np.concatenate((self.cds_lefts, self.cds_rights)))
        self.gene_right = np.max(np.concatenate((self.cds_lefts, self.cds_rights)))
        
        self.CRISPR_CUT_INDEX = 3
    
    ###############################################################################################################
    def __repr__(self):
        return "CRISPR_Gene()"
    
    ###############################################################################################################
    def __str__(self):
        return "Name: " + self.name + "\n" + \
               "Strand: " + self.strand + "\n" + \
               "Chrom: " + self.chrom + "\n" + \
               "Left: " + str(self.gene_left) + "\n" + \
               "Right: " + str(self.gene_right) + "\n" + \
               "CDS lefts: " + str(self.cds_lefts) + "\n" + \
               "CDS rights: " + str(self.cds_rights) + "\n"
        
    
    ###############################################################################################################
    def get_gene_gff_df(self):
        
        this_gene_df = self.gene_gff_df[self.gene_gff_df['Name'].isin([self.name])]
        
        # gene not in gff file
        if (this_gene_df.shape[0] == 0):
            raise ValueError(self.name + " is missing from the gff file data frame") 
            
        return(this_gene_df)
    
    ###############################################################################################################
    def is_neg_strand(self):
        return (self.strand == '-')
    
    ###############################################################################################################
    def get_cds_dna(self):
        ret_gene_cds_dna = ''
        
        for left,right in zip(self.cds_lefts, self.cds_rights):
            ret_gene_cds_dna = ret_gene_cds_dna + self.genome_seq[self.chrom].seq[left:right]
        
        if (self.is_neg_strand()):
            ret_gene_cds_dna = ret_gene_cds_dna.reverse_complement()
        
        return(ret_gene_cds_dna)
    
    ###############################################################################################################
    def get_protein_seq(self):
        return (self.get_cds_dna().translate())
    
    ###############################################################################################################
    def __guide_positions_to_df(self,pam, is_neg_strand, seq_len_around_cut, excluded_seqs, \
                                guides0_chr_pos, guides_cut_chr_pos, guides_cut_gene_dna_pos, guides_exon_num):
     
        
        out_df = pd.DataFrame({"Gene" : self.name , \
                               "CDS_len_nts" : self.cds_length_nts, \
                               "chrom" : self.chrom, \
                               "gene_strand" : self.strand,  \
                               "guide0_chr_pos" : guides0_chr_pos, \
                               "guide_cut_chr_pos" : guides_cut_chr_pos, \
                               "guide_cut_gene_nt_pos" : guides_cut_gene_dna_pos, \
                               "guide_gene_exon" : guides_exon_num})
        
        
        out_df[ ["CDS_len_nts", "guide0_chr_pos", "guide_cut_chr_pos", "guide_cut_gene_nt_pos", "guide_gene_exon"] ] = \
            out_df[ ["CDS_len_nts", "guide0_chr_pos", "guide_cut_chr_pos", "guide_cut_gene_nt_pos", "guide_gene_exon"] ].astype(int)
        
        gene_chr = self.genome_seq[self.chrom].seq
        
        #if (self.strand == '+'):
        #else
        
        
        #guides.revcomp =  [ str(gene_chr[(guide0-20):guide0+3])  for guide0  in guides0_chr_pos]
        
        #print out_df
        
        
        if (is_neg_strand):
            out_df['guide_strand'] = '-'
            #out_df['guide_id'] =  [ self.name + "_mS_" + str(guide_i)  for guide_i  in range(len(guides0_chr_pos))]
            
            #print [ gene_chr[(guide0-2  for guide0  in guides0_chr_pos]
            #print [ gene_chr[(guide0-2  for guide0  in guides0_chr_pos]
                   
            guide          =  [ str(gene_chr[(int(guide0)-2):(int(guide0)+21)].reverse_complement())  for guide0  in guides0_chr_pos]
            guide_noPAM    =  [ str(gene_chr[(int(guide0)-2+3):(int(guide0)+21)].reverse_complement())  for guide0  in guides0_chr_pos]
            guide_PAM_p7   =  [ str(gene_chr[(int(guide0)-2-7):(int(guide0)+21)].reverse_complement())  for guide0  in guides0_chr_pos]
            guide_PAM_m4p3 =  [ str(gene_chr[(int(guide0)-2-3):(int(guide0)+21+4)].reverse_complement())  for guide0  in guides0_chr_pos]
            guide_PAM_m4p3 =  [ str(gene_chr[(int(guide0)-2-3):(int(guide0)+21+4)].reverse_complement())  for guide0  in guides0_chr_pos]
            
            # the nucleotide left of the cut
            #guide_cut_chr_pos = guides0_chr_pos+2
            seq_around_cut =  [ str(gene_chr[(int(guide_cut)-seq_len_around_cut+1):(int(guide_cut)+seq_len_around_cut+1)].reverse_complement())  for guide_cut  in guides_cut_chr_pos]
            
            #what nt are covered by the guide and can be mutated
            if (self.is_neg_strand()):  # is gene negative strand
              out_df['guide_start_gene_nt_pos'] =  guides_cut_gene_dna_pos - 17
              out_df['guide_end_gene_nt_pos'] =  guides_cut_gene_dna_pos + 2
            else:
              out_df['guide_start_gene_nt_pos'] =  guides_cut_gene_dna_pos - 3
              out_df['guide_end_gene_nt_pos'] =  guides_cut_gene_dna_pos + 16
            
            
            
        else:
            out_df['guide_strand'] = '+'
            #out_df['guide_id'] =  [ self.name + "_pS_" + str(guide_i)  for guide_i  in range(len(guides0_chr_pos))]
            
            guide          =  [ str(gene_chr[(int(guide0)-20):(int(guide0)+3)])  for guide0  in guides0_chr_pos]
            guide_noPAM    =  [ str(gene_chr[(int(guide0)-20):(int(guide0))])  for guide0  in guides0_chr_pos]
            guide_PAM_p7   =  [ str(gene_chr[(int(guide0)-20):(int(guide0)+3+7)])  for guide0  in guides0_chr_pos]
            guide_PAM_m4p3 =  [ str(gene_chr[(int(guide0)-20-4):(int(guide0)+3+3)])  for guide0  in guides0_chr_pos]
            
            #guide_cut_chr_pos = guides0_chr_pos-3
            seq_around_cut =  [ str(gene_chr[(int(guide_cut)-seq_len_around_cut+1):(int(guide_cut)+seq_len_around_cut+1)])  for guide_cut  in guides_cut_chr_pos]
            
            #what nt are covered by the guide and can be mutated
            if (self.is_neg_strand()):  # is gene negative strand
              out_df['guide_start_gene_nt_pos'] =  guides_cut_gene_dna_pos - 3
              out_df['guide_end_gene_nt_pos'] =  guides_cut_gene_dna_pos + 16
            else:
              out_df['guide_start_gene_nt_pos'] =  guides_cut_gene_dna_pos -17
              out_df['guide_end_gene_nt_pos'] =  guides_cut_gene_dna_pos + 2
        
        
        #out_df['guide_cut_chr_pos'] = guide_cut_chr_pos
    
        out_df['guide']          = guide
        out_df['guide_noPAM']    = guide_noPAM
        out_df['guide_PAM_p7']   = guide_PAM_p7
        out_df['guide_PAM_m4p3'] = guide_PAM_m4p3
        out_df['guide_cut_gene_aa_pos'] = np.floor(guides_cut_gene_dna_pos/3) # zero based otherwise need to be ceil
        out_df['guide_cut_gene_aa_frame'] = np.remainder(guides_cut_gene_dna_pos,3)
        out_df['dna_around_cut'] = seq_around_cut
        out_df['dna_around_cut'] = out_df['dna_around_cut'].astype(str)
        
        num_lines_before_filtering = out_df.shape[0]
        
        # filtering guides that contain excluded sequences around them
        out_df = out_df[~out_df['dna_around_cut'].str.contains( '|'.join(excluded_seqs) )]
        
        #print 'Num lines that were filtered out because they contain excluded sequences %d (left %d lines)' % (num_lines_before_filtering-out_df.shape[0],out_df.shape[0])
       
        if (is_neg_strand):
            out_df['guide_id'] =  [ self.name + "_mS_" + str(guide_i)  for guide_i  in range(out_df.shape[0])]
        else:
            out_df['guide_id'] =  [ self.name + "_pS_" + str(guide_i)  for guide_i  in range(out_df.shape[0])]
            
    
 
        #print '--------------------'
        #print 'is negative strand: ' + str(is_neg_strand)
        
        #print guide[1:4]
        #print guide_noPAM[1:4]
        #print guide_PAM_p7[1:4]
        #print guide_PAM_m4p3[1:4]
        
        #print guides_cut_gene_dna_pos[1:4]
        #print out_df['guide_cut_gene_aa_pos'].iloc[1:4]
        #print out_df['guide_cut_gene_aa_frame'].iloc[1:4]
        
        #print out_df['guide_cut_chr_pos']
        
        #print '--------------------'
        
        # position in gene in percent
        out_df["guide_cut_gene_pos_frac"] = out_df["guide_cut_gene_nt_pos"] / out_df["CDS_len_nts"]

        
        return(out_df)
        
    ###############################################################################################################
    def get_all_guides_that_cut_in_cds(self,pam, seq_len_around_cut,
                                       min_mut_pos_in_guide, max_mut_pos_in_guide,
                                       excluded_seqs, mapping_cmd, sort_by = '5to3'):
        """
        Get all guides that cut within the coding sequence
        order from start to end
        
        the cut site is always before the cut nt. so for positive strand genes it is left and for negative it is right 
        """
        
        ordered_lefts = self.cds_lefts
        ordered_rights = self.cds_rights
        if (self.is_neg_strand()):
            ordered_lefts = ordered_lefts[::-1]
            ordered_rights = ordered_rights[::-1]
        
        ######
        # search positive strand for pam
        ######
        cur_cds_nt_start = 0
        exon_num = 0
        guides0_chr_pos = np.empty(0,dtype=int)
        guides_cut_chr_pos = np.empty(0,dtype=int)
        guides_cut_gene_dna_pos = np.empty(0,dtype=int)
        guides_exon_num = np.empty(0,dtype=int)

        for left,right in zip(ordered_lefts,ordered_rights):
            
            # cut is to the right of the nucleotide
            cur_left_for_pam  = left  + (self.CRISPR_CUT_INDEX + len(pam) - 1) + (1 * self.is_neg_strand()) 
            cur_right_for_pam = right + (self.CRISPR_CUT_INDEX + len(pam) - 1) + (1 * self.is_neg_strand()) 
              

            
            seq = self.genome_seq[self.chrom].seq[cur_left_for_pam:cur_right_for_pam]
            
            # returns a list of all the positions in that cut in cds
            cur_pam_dists = np.array([m.start() for m in re.finditer("(?="+pam+")", str(seq))])
            
            # removing guides that are not entirely in the CDS
            if ( (not np.isnan(min_mut_pos_in_guide))  and (not np.isnan(max_mut_pos_in_guide)) ):
                min_mut_pos_in_guide = int(min_mut_pos_in_guide)
                max_mut_pos_in_guide = int(max_mut_pos_in_guide)
                
                cur_pam_dists = cur_pam_dists[cur_pam_dists >= (-min_mut_pos_in_guide) - (self.CRISPR_CUT_INDEX + 1 * self.is_neg_strand() ) ]
                cur_pam_dists = cur_pam_dists[cur_pam_dists <= (len(seq) - 1  + len(pam) - 1 )  + ( (-max_mut_pos_in_guide) - (self.CRISPR_CUT_INDEX + 1 * self.is_neg_strand()) ) ]
                   
            
            cur_guides0_chr_pos = (cur_pam_dists-1) + cur_left_for_pam 
            
            if (self.is_neg_strand()): # negative
                cur_guides_cut_gene_dna_pos = (len(seq)-1-cur_pam_dists) + cur_cds_nt_start
                cur_guides_cut_chr_pos      = cur_guides0_chr_pos - (self.CRISPR_CUT_INDEX + 1)  # the cut is right of the nt
            else:
                cur_guides_cut_gene_dna_pos = cur_pam_dists + cur_cds_nt_start
                cur_guides_cut_chr_pos      = cur_guides0_chr_pos - self.CRISPR_CUT_INDEX  # the cut is left of the nt
            
            
            cur_guides_exon_num = np.full_like(cur_guides_cut_gene_dna_pos,exon_num)
            
            
            guides0_chr_pos = np.concatenate((guides0_chr_pos,cur_guides0_chr_pos))
            guides_cut_chr_pos  = np.concatenate((guides_cut_chr_pos,cur_guides_cut_chr_pos))
            guides_cut_gene_dna_pos = np.concatenate((guides_cut_gene_dna_pos,cur_guides_cut_gene_dna_pos))
            guides_exon_num = np.concatenate((guides_exon_num,cur_guides_exon_num))
            
            
            cur_cds_nt_start = cur_cds_nt_start + (right - left)
            exon_num = exon_num + 1
            
            
        pos_strand_guides_df = self.__guide_positions_to_df(pam, False, seq_len_around_cut, excluded_seqs, \
                                                            guides0_chr_pos, guides_cut_chr_pos, guides_cut_gene_dna_pos, guides_exon_num) 
        
        ######
        # search negative strand for pam
        ######
        cur_cds_nt_start = 0
        exon_num = 0
        guides0_chr_pos = np.empty(0,dtype=int)
        guides_cut_chr_pos = np.empty(0,dtype=int)
        guides0_gene_dna_pos = np.empty(0,dtype=int)
        guides_cut_gene_dna_pos = np.empty(0,dtype=int)
        guides_exon_num = np.empty(0,dtype=int)
        
        for left,right in zip(ordered_lefts,ordered_rights):
            
            
            cur_left_for_pam  = int(left)  - (self.CRISPR_CUT_INDEX + len(pam)+1) + (1 * self.is_neg_strand())
            cur_right_for_pam = int(right) - (self.CRISPR_CUT_INDEX + len(pam)+1) + (1 * self.is_neg_strand())
            
            seq = self.genome_seq[self.chrom].seq[cur_left_for_pam:cur_right_for_pam]
            
            revcomp_pam = Seq(pam,generic_dna).reverse_complement()
            
            # returns a list of all the positions in that cut in cds
            cur_pam_dists = np.array([m.start() for m in re.finditer("(?="+str(revcomp_pam)+")", str(seq))])
            
            
            if ( (not np.isnan(min_mut_pos_in_guide))  and (not np.isnan(max_mut_pos_in_guide)) ):
                min_mut_pos_in_guide = int(min_mut_pos_in_guide)
                max_mut_pos_in_guide = int(max_mut_pos_in_guide)
                
                cur_pam_dists = cur_pam_dists[cur_pam_dists >= (-min_mut_pos_in_guide) - (self.CRISPR_CUT_INDEX + 1 * self.is_neg_strand() ) ]
                cur_pam_dists = cur_pam_dists[cur_pam_dists <= (len(seq) - 1  + len(pam) - 1 )  + ( (-max_mut_pos_in_guide) - (self.CRISPR_CUT_INDEX + 1 * self.is_neg_strand()) ) ]
            
      
            
            cur_guides0_chr_pos = (cur_pam_dists+2) + cur_left_for_pam
            
            if (self.is_neg_strand()): # negative         
                cur_guides_cut_gene_dna_pos =  (len(seq)-1-cur_pam_dists) + cur_cds_nt_start
                cur_guides_cut_chr_pos      =  cur_guides0_chr_pos + self.CRISPR_CUT_INDEX  # the cut is right of the nt
            else: # positive
                cur_guides_cut_gene_dna_pos = cur_pam_dists + cur_cds_nt_start
                cur_guides_cut_chr_pos      = cur_guides0_chr_pos + self.CRISPR_CUT_INDEX + 1  # the cut is leftot the nt
            
            
            cur_guides_exon_num = np.full_like(cur_guides_cut_gene_dna_pos,exon_num)
            
            
            guides0_chr_pos = np.concatenate((guides0_chr_pos,cur_guides0_chr_pos))
            guides_cut_chr_pos  = np.concatenate((guides_cut_chr_pos,cur_guides_cut_chr_pos))
            guides_cut_gene_dna_pos = np.concatenate((guides_cut_gene_dna_pos,cur_guides_cut_gene_dna_pos))
            guides_exon_num = np.concatenate((guides_exon_num,cur_guides_exon_num))
            
            cur_cds_nt_start = cur_cds_nt_start + (right - left)
            exon_num = exon_num + 1
            
        
        neg_strand_guides_df = self.__guide_positions_to_df(pam, True, seq_len_around_cut, excluded_seqs, \
                                                            guides0_chr_pos, guides_cut_chr_pos, guides_cut_gene_dna_pos, guides_exon_num)
        
        
        # concating the positive and negative strands guides
        guides_df = pd.concat([pos_strand_guides_df, neg_strand_guides_df])
        
        # adding for each guide its location in the gene (5' -> 3'; fraction)
        guides_df["guide_cut_gene_pos_frac"] = guides_df["guide_cut_gene_nt_pos"] / guides_df["CDS_len_nts"]

        
        # if the 'sort' method is onlyStopCodon then leave only guide that cut the stop codon
        if sort_by == 'onlyStopCodon':
            guides_df = guides_df.ix[( ( (guides_df['CDS_len_nts']).values / 3) == ( (guides_df['guide_cut_gene_aa_pos']).values + 1) ) ,:]
        
        # calculating Azimuth score
        #print "Calculating Azimuth score"
        guides_df = cal_azimuth_score(guides_df, output_filename_GUIDE_withScores = "", guides_PAMm4p3_col_name="guide_PAM_m4p3")
        
        # calculating off targets
        #print "Testing off targets"
        guides_df = eval_guides_off_targets(guides_df, self.genome_seq, 'guide_id', 'guide_noPAM', pam, mapping_cmd)
        
        
        return (guides_df)
      
    ############################################################################################################### 
    def get_K_best_guide_ids(self, guides_df, K, off_targets_min_mismatch = 10, min_azimuth_score = 0, max_guide_pos_frac_in_gene = None, sort_by = '5to3'):
        """ 
        Form a data frame containing guides targeting the gene:
        Select K best guides (maximal number after filtering guides)   
        off_targets_min_mismatch - minimal number of mismatch allowed for off targets (> 4 means no off targets are allowed, 0 means no filtering)
        min_azimuth_score - min score
        sort_by: '5to3' - by position in the gene, 'Azimuth' - by Azimuth score
        """
        
        iloc_list =  np.array(range(guides_df.shape[0]))
        
        #print "In get_K_best_guide_ids, starting from %d guides" % (len(iloc_list))
        
        # filtering off tragets 
        if (off_targets_min_mismatch > 0):
            iloc_list = iloc_list[ (guides_df['guide_map_mismatch_0'].iloc[iloc_list] == 0).values ]
        if (off_targets_min_mismatch > 1):
            iloc_list = iloc_list[ (guides_df['guide_map_mismatch_1'].iloc[iloc_list] == 0).values ]
        if (off_targets_min_mismatch > 2):
            iloc_list = iloc_list[ (guides_df['guide_map_mismatch_2'].iloc[iloc_list] == 0).values ]
        if (off_targets_min_mismatch > 3):
            iloc_list = iloc_list[ (guides_df['guide_map_mismatch_3'].iloc[iloc_list] == 0).values ]
        
        
        # filtering low azimuth score
        iloc_list = iloc_list[ (guides_df['Azimuth'].iloc[iloc_list] >= min_azimuth_score).values ]
        
        
        # filtering guides too downstream in the gene
        if max_guide_pos_frac_in_gene is not None:
            iloc_list = iloc_list[ (guides_df['guide_cut_gene_pos_frac'].iloc[iloc_list] <= max_guide_pos_frac_in_gene).values ]
        
        
        #print "In get_K_best_guide_ids, after removing low Azimuth score: %d" % (len(iloc_list))
        
        if sort_by != 'onlyStopCodon':
            # filtering cut in the stop codon
            iloc_list = iloc_list[ ( (guides_df['CDS_len_nts'].iloc[iloc_list] / 3) > (guides_df['guide_cut_gene_aa_pos'].iloc[iloc_list]).values + 1 ).values]
            
            #print "In get_K_best_guide_ids, after removing guides that cut the stop codon: %d" % (len(iloc_list))
        #else:
            #print "In get_K_best_guide_ids, Not removing guides that cut the stop codon"
            
            
        iloc_list = iloc_list[ ( (guides_df['guide_cut_gene_aa_pos'].iloc[iloc_list]).values > 0)]
        
        
        #print "In get_K_best_guide_ids, after removing guides that cut the start codon: %d" % (len(iloc_list))
       
        
        # sorting
        if (sort_by == '5to3'):
            iloc_list = iloc_list[ np.argsort( (guides_df['guide_cut_gene_nt_pos'].iloc[iloc_list]).values ) ]
        elif (sort_by == 'Azimuth'):
            iloc_list = iloc_list[ np.argsort( (guides_df['Azimuth'].iloc[iloc_list]).values )[::-1] ]
        
        
        # get maximum K guides
        iloc_list = iloc_list[0:min(len(iloc_list), K)]
        
        #print (debug in iloc_list)
        #print '----+++'
        #print iloc_list
        
        
        #print "In get_K_best_guide_ids, after selecting top K: %d" % (len(iloc_list))
       
        return (guides_df['guide_id'].iloc[iloc_list].copy() )

    ###############################################################################################################
    def get_donor_mut_for_guides(self, guides_df, guides_ids, donor_mut_type, num_donor_variants = 1,
                                 mut_pos_in_guide = None, donor_length=100, donor_seq_offsets=[0],
                                 set_name = "", min_dist_cut_to_donor_edge = 30, excluded_seqs = [], do_revcomp_donor = False, scramble_guide_and_donor = False):
        """
        donor_mut_type - the type of mutations <'synonymous'|'nonsynonymous'|'nonsense'|'frameshift_<in|del>_<1|2>'>
        mut_pos_in_guide 1-20, default is 3 which is the cut site (between 3 and 4)
        num_donor_variants -  how many variants to produce from each mutation type for each guide
                            notice that frame shifts can have up to 4^insertion length 
                            deletions num of deletions
                            and other types depending on the number of options
        donor_length - total length of the sequence

        """
        out_guide_donor_df = pd.DataFrame(data=None)

        for donor_seq_offset in donor_seq_offsets:
            
            # TODO - this should be removed after implementing donors with offsets 
            if np.isnan(donor_seq_offset):
                raise ValueError("None/NaN offsets not implemented for gene guides - TODO...")
            
            for guide_id in guides_ids:
                
              #print "guide_id %s , donor_seq_offset %d" % (guide_id, donor_seq_offset) #DEBUG
                
              cur_guide_guide_donor_df = self.get_donor_mut_for_guide(guides_df, guide_id, donor_mut_type, num_donor_variants,
                                                                      mut_pos_in_guide, donor_length,int(donor_seq_offset), set_name,
                                                                      min_dist_cut_to_donor_edge,do_revcomp_donor, scramble_guide_and_donor)
              out_guide_donor_df = out_guide_donor_df.append(cur_guide_guide_donor_df,ignore_index=True)
              
              #print cur_guide_guide_donor_df.shape[0]
              #print out_guide_donor_df.shape[0]
        
        if out_guide_donor_df.shape[0] > 0:
            out_guide_donor_df['contain_excluded_sequences']  =  out_guide_donor_df['donor_seq'].str.contains( '|'.join(excluded_seqs) )

        return(out_guide_donor_df.reset_index(drop=True))

    ###############################################################################################################
    def get_donor_mut_for_guide(self, guides_df, guide_id, donor_mut_type, num_donor_variants = 1,
                                mut_pos_in_guide = None, donor_length=100,donor_seq_offset=0, set_name = "",
                                min_dist_cut_to_donor_edge = 30, do_revcomp_donor = False,
                                scramble_guide_and_donor = False):
        """
        donor_mut_type - the type of mutations <'synonymous'|'nonsynonymous'|'nonsense'|'frameshift_<in|del>_<1|2>'>
        mut_pos_in_guide 1-20 or None. default is 3.  cut is between 3 and 4 but the cut position depends on the gene orientation (cut before nt)
                                 None uses the cut site position
                                 ONLY None is IMPLEMENTED 
        num_donor_variants -  how many variants to produce from each mutation type for each guide
                            notice that frame shifts can have up to 4^insertion length 
                            deletions - only one
                            and other types depending on the number of options
        donor_length - total length of the sequence
        donor_seq_offset - allow to move donor left and right few bp to create more constructs
        
        Important: this function counts on the edited position and the cut site being in the CDS
        """
        
        #print '++++++++++++++++++++++++++1'
        #print guides_df[guides_df['guide_id'] == guide_id]
        #print '++++++++++++++++++++++++++2'
        
        ################################
        # extracting guide details
        ################################
        
        # single row of the guide
        guides_df = guides_df[guides_df['guide_id'] == guide_id]

        guide_gene = str(guides_df['Gene'].iloc[0])
        guide_cut_chr_pos = int(guides_df['guide_cut_chr_pos'].iloc[0])
        guide0_chr_pos = int(guides_df['guide0_chr_pos'].iloc[0])
        guide_is_negative_strand = (guides_df['guide_strand'].iloc[0] == '-')
        
        guide_cut_gene_aa_pos = int(guides_df['guide_cut_gene_aa_pos'].iloc[0])
        guide_cut_gene_aa_frame = int(guides_df['guide_cut_gene_aa_frame'].iloc[0])
   
        
        # TODO currently only mut_pos_in_guide = None is implemented which uses the cut site position and frame
        if mut_pos_in_guide is not None and not np.isnan(mut_pos_in_guide):
            #raise ValueError("get_donor_mut_for_guide: only None (cut site) is implemented for mut_pos_in_guide:" +  str(mut_pos_in_guide))
            # should take gene and guide orientations into account
            
            
            mut_pos_in_guide = int(mut_pos_in_guide)
            
            #DEBUG
            #print 'mutation postion in guide %d' % (mut_pos_in_guide)  
            
            cut_nt_pos_in_guide = -self.CRISPR_CUT_INDEX
            if (not self.is_neg_strand() and guide_is_negative_strand) or  (self.is_neg_strand() and not guide_is_negative_strand):
                cut_nt_pos_in_guide = int(-self.CRISPR_CUT_INDEX - 1)
                
            #DEBUG
            #print "cut pos in guide %d" % (cut_nt_pos_in_guide)
           
            edit_pos2cut_gene_nt_diff =  -mut_pos_in_guide +  cut_nt_pos_in_guide
            edit_pos2cut_chr_pos_diff =  -mut_pos_in_guide +  cut_nt_pos_in_guide

            if (not self.is_neg_strand() and not guide_is_negative_strand) or  (self.is_neg_strand() and guide_is_negative_strand):
               edit_pos2cut_gene_nt_diff = -edit_pos2cut_gene_nt_diff
            
            if (not guide_is_negative_strand):
                edit_pos2cut_chr_pos_diff = -edit_pos2cut_chr_pos_diff
                edit_chr_pos = guide0_chr_pos + mut_pos_in_guide
            else:
                edit_chr_pos = guide0_chr_pos - mut_pos_in_guide
            
            
            #DEBUG
            #print "edit_pos2cut_gene_nt_diff %d" % (edit_pos2cut_gene_nt_diff)
            #print "edit_pos2cut_chr_pos_diff %d" % (edit_pos2cut_chr_pos_diff)
            #print "guide0_chr_pos %d " % (guide0_chr_pos)
            #print "guide_cut_chr_pos %d " % (guide_cut_chr_pos)
            #print "edit_chr_pos %d" % (edit_chr_pos)
           
             
            edit_aa_pos = guide_cut_gene_aa_pos +  np.floor((guide_cut_gene_aa_frame+edit_pos2cut_gene_nt_diff)/3) 
            edit_aa_frame = np.remainder(guide_cut_gene_aa_frame+edit_pos2cut_gene_nt_diff,3)
        
        else:
            edit_chr_pos = guide_cut_chr_pos
            edit_aa_pos = guide_cut_gene_aa_pos
            edit_aa_frame = guide_cut_gene_aa_frame
            edit_pos2cut_chr_pos_diff = 0
            
            #DEBUG
            #print "-----------mut in cut"
        
        # parsing frameshift mutation type
        donor_mut_type_splt = donor_mut_type.split('_')
        donor_mut_name = donor_mut_type_splt[0]
        if (len(donor_mut_type_splt)>1):
            donor_mut_subtype = donor_mut_type_splt[1]
        if (len(donor_mut_type_splt)>2):
            if  (donor_mut_name == 'frameshift'):
                donor_mut_len = int(donor_mut_type_splt[2])
            else:
                donor_mut_infoStr = donor_mut_type_splt[2]
        
        
        ######################################
        # allocating out dataframe
        ######################################
        out_guide_donor_df = pd.DataFrame(data=None)
        
        # flag of the donor strand (should match the guide so we first create the donor and then fix the orientation)
        is_donor_revcomp = False
        donor_seq = Seq("", generic_dna)
        
        if  (donor_mut_name == 'frameshift'):
            
            if (donor_mut_subtype == "in"): # insertion
                # num of nt to add left and right (-1 for the cut nucleotide, -donor_mut_len for the added nucleotide)
                donor_nt_add_left  = int( np.floor((donor_length-donor_mut_len+edit_pos2cut_chr_pos_diff)/2) - donor_seq_offset)
                donor_nt_add_right = int((donor_length-1-donor_mut_len) - donor_nt_add_left)
                
                if (donor_nt_add_left < min_dist_cut_to_donor_edge or donor_nt_add_right <  min_dist_cut_to_donor_edge):
                     warnings.warn("Insert to big %d, there are not enough overlap with genome left: %d, right: %d - NOT returning " % (donor_mut_len, donor_nt_add_left, donor_nt_add_right))
                else:
                
                    # generate insertion sequences
                    insert_seqs, insert_seqs_num_requested_ok = get_K_random_dna_nts(donor_mut_len,num_donor_variants)

                    if not insert_seqs_num_requested_ok:
                        warnings.warn("can not create that num_donor_variants: " + str(num_donor_variants) + " of DNA length donor_mut_len:" +  str(donor_mut_len))
           
                    # DEBUG
                    #print str(int(np.floor((donor_length-1-donor_mut_len)/2) + donor_seq_offset + edit_pos2cut_chr_pos_diff))
            
                    #print edit_chr_pos
                    #print type(edit_chr_pos)
                    #print donor_nt_add_left
                    #print type(donor_nt_add_left)
                    #print donor_nt_add_right
                    #print type(donor_nt_add_right)
                    

                    for idx,cur_insert_seq in enumerate(insert_seqs):
                        # TODO ? make sure that does not create a stop codon? may be not because also cutting can create that
                        if (self.is_neg_strand()):
                            donor_seq = self.genome_seq[self.chrom].seq[(edit_chr_pos-donor_nt_add_left):(edit_chr_pos+1)] + \
                                        cur_insert_seq + \
                                        self.genome_seq[self.chrom].seq[(edit_chr_pos+1): (edit_chr_pos + donor_nt_add_right + 1)]
                        else:
                            donor_seq = self.genome_seq[self.chrom].seq[(edit_chr_pos-donor_nt_add_left):(edit_chr_pos)] + \
                                        cur_insert_seq + \
                                        self.genome_seq[self.chrom].seq[edit_chr_pos: (edit_chr_pos + donor_nt_add_right + 1)]
                        
                        # DEBUG
                        #print cur_insert_seq
                        #print donor_seq
                            
                        
                        donor_info_str = "in:" + str(cur_insert_seq)
                        
                        
                        # convert to match the guide orientation
                        if (guide_is_negative_strand):
                            donor_seq = donor_seq.reverse_complement()
                            
                        if do_revcomp_donor:
                            donor_seq = donor_seq.reverse_complement()
                        
                        if scramble_guide_and_donor:
                            donor_info_str = donor_info_str + ":" + "scramble"
                            tmp_donor_seq_lst = list(str(donor_seq))
                            random.shuffle(tmp_donor_seq_lst)
                            donor_seq = Seq(''.join(tmp_donor_seq_lst),alphabet=generic_dna)
                            
                        
                        
                        # appending to donor sequences matrix
                        cur_donor_line = pd.DataFrame({'Gene' : guide_gene, 'guide_id' : pd.Series(guide_id), 
                                          'donor_id' : pd.Series(guide_id + ':' + donor_mut_type + ':offset' + str(donor_seq_offset) + ':donorID' + str(idx) + ':EditPosInGuide' +  str(mut_pos_in_guide)), 
                                          'donor_seq': pd.Series(str(donor_seq)), 
                                          'donor_seq_shift' : pd.Series(int(donor_seq_offset)), 
                                          'donor_mut_pos_in_guide' : pd.Series(str(mut_pos_in_guide)), 
                                          'donor_info_str' : pd.Series( donor_info_str ),
                                          'set_name' : pd.Series(str(set_name)) })
                        
                        out_guide_donor_df = out_guide_donor_df.append(cur_donor_line)
                    
            elif (donor_mut_subtype == "del"): # deletion (currently only deletion after the cut (in term of the gene) are implemented)
                
                
                print " ------ In del --------------"
                
                if num_donor_variants != 1:
                    raise ValueError("Currently a deletion can produce only a single varient: " + str(num_donor_variants))
                else:
                    idx = 0
       
                if (self.is_neg_strand()):
                    left_side_deletion = int(np.floor(donor_mut_len/2))
                else:
                    left_side_deletion = int(np.ceil(donor_mut_len/2))
                right_side_deletion = int(donor_mut_len - left_side_deletion)
                
                # num of nt to add left and right (-1 for the cut nucleotide, +donor_mut_len for the deleted nucleotide)
                donor_nt_add_left = int(np.floor(donor_length/2)) # + edit_pos2cut_chr_pos_diff + donor_seq_offset
                donor_nt_add_right = int(donor_length - donor_nt_add_left)
                
                cur_donor_offset = int(np.floor(edit_pos2cut_chr_pos_diff/2 + donor_seq_offset))
                
                
                #print "donor_mut_len: %d" % (donor_mut_len)
                #print "right_side_deletion: %d" % (right_side_deletion)
                #print "donor_nt_add_left: %d" % (donor_nt_add_left)
                #print "donor_nt_add_right: %d" % (donor_nt_add_right)
                #print "edit_pos2cut_chr_pos_diff: %d" % (edit_pos2cut_chr_pos_diff)
                #print "donor_seq_offset: %d" % (donor_seq_offset)
                #print "cur_donor_offset: %d" % (cur_donor_offset)
                
                # TODO ? make sure that does not create a stop codon? may be not because also cutting can create that
                # deleting downstream to the cut
                #if (self.is_neg_strand()):
                #    donor_seq = self.genome_seq[self.chrom].seq[(edit_chr_pos-donor_nt_add_left):(edit_chr_pos-donor_mut_len)] + \
                #                self.genome_seq[self.chrom].seq[(edit_chr_pos+1): (edit_chr_pos + donor_nt_add_right + 1)]   
                
                #else:
                donor_seq = self.genome_seq[self.chrom].seq[(edit_chr_pos-donor_nt_add_left-left_side_deletion+cur_donor_offset):(edit_chr_pos-left_side_deletion+1)] + \
                            self.genome_seq[self.chrom].seq[(edit_chr_pos+1+right_side_deletion): (edit_chr_pos + 1 + right_side_deletion + donor_nt_add_right + cur_donor_offset)]

                #if (self.is_neg_strand()):
                #    donor_seq = self.genome_seq[self.chrom].seq[(edit_chr_pos-donor_nt_add_left):(edit_chr_pos+1-donor_mut_len)] + \
                #                self.genome_seq[self.chrom].seq[(edit_chr_pos+1): (edit_chr_pos + donor_nt_add_right + 1)]
                #else:
                #    donor_seq = self.genome_seq[self.chrom].seq[(edit_chr_pos-donor_nt_add_left):(edit_chr_pos)] + \
                #                self.genome_seq[self.chrom].seq[(edit_chr_pos+donor_mut_len): (edit_chr_pos + donor_nt_add_right + 1)]
                
                #print str(self.genome_seq[self.chrom].seq[(edit_chr_pos-donor_nt_add_left-left_side_deletion+cur_donor_offset+1):(edit_chr_pos-left_side_deletion+1)])
                #print str(self.genome_seq[self.chrom].seq[(edit_chr_pos+1+right_side_deletion): (edit_chr_pos + 1 + right_side_deletion + donor_nt_add_right + cur_donor_offset)])
                
                #print "donor seq:"
                #print str(donor_seq)



                # convert to match the guide orientation
                if (guide_is_negative_strand):
                    donor_seq = donor_seq.reverse_complement()
                
                if do_revcomp_donor:
                    donor_seq = donor_seq.reverse_complement()
                
                
                donor_info_str = "del:" + str(donor_mut_len)
                
                if scramble_guide_and_donor:
                    donor_info_str = donor_info_str + ":" + "scramble"
                    tmp_donor_seq_lst = list(str(donor_seq))
                    random.shuffle(tmp_donor_seq_lst)
                    donor_seq = Seq(''.join(tmp_donor_seq_lst),alphabet=generic_dna)

                # appending to donor sequences matrix
                cur_donor_line = pd.DataFrame({'Gene' : guide_gene, 'guide_id' : pd.Series(guide_id), 
                                  'donor_id' : pd.Series(guide_id + ':' + donor_mut_type + ':offset' + str(donor_seq_offset) + ':donorID' + str(idx) + ':EditPosInGuide' +  str(mut_pos_in_guide)), 
                                  'donor_seq': pd.Series(str(donor_seq)), 
                                  'donor_seq_shift' : pd.Series(int(donor_seq_offset)), 
                                  'donor_mut_pos_in_guide' : pd.Series(str(mut_pos_in_guide)), 
                                  'donor_info_str' : pd.Series(donor_info_str),
                                  'set_name' : pd.Series(str(set_name)) })
                out_guide_donor_df = out_guide_donor_df.append(cur_donor_line,ignore_index=True)
                
            else:
                raise ValueError('get_donor_mut_for_guide unknown donor_mut_subtype:' + donor_mut_subtype)
                
            
        elif (donor_mut_name == 'synonymous' or donor_mut_name == 'nonsynonymous' or
              donor_mut_name == 'nonsense' or donor_mut_name == 'stop2aa' or
              donor_mut_name == '1bp'):
            
            # num of nt to add left and right (-3 for the codon)
            donor_nt_add_left = int(np.floor((donor_length + edit_pos2cut_chr_pos_diff)/2) - donor_seq_offset)
            donor_nt_add_right = int((donor_length-1) - donor_nt_add_left)
            
            org_donor_seq = self.genome_seq[self.chrom].seq[(edit_chr_pos-donor_nt_add_left):(edit_chr_pos + donor_nt_add_right + 1)]
                        
            # num of nts before the cut 
            if (self.is_neg_strand()):
                nt_before_cut = donor_nt_add_right
                org_donor_seq = org_donor_seq.reverse_complement()
            else:
                nt_before_cut = donor_nt_add_left
                
            # the relevant codon is where the cut is or downstream to the cut
            i_first_nt_in_codon = (nt_before_cut-edit_aa_frame)
            cut_codon_dna_seq = org_donor_seq[(i_first_nt_in_codon):(i_first_nt_in_codon+3)]
            
            # the original codon
            #print "XXXXXXXX1"
            #print cut_codon_dna_seq
            #print nt_before_cut
            #print guide_cut_gene_aa_frame
            #print i_first_nt_in_codon
            #print edit_chr_pos
            #print donor_nt_add_left
            #print donor_seq_offset
            #print nt_before_cut
            #print self.chrom
            #print org_donor_seq
            #print (edit_chr_pos-donor_nt_add_left)
            #print (edit_chr_pos + donor_nt_add_right + 1)
            
            
            #DEBUG print
            #print type(cut_codon_dna_seq)
            #print str(cut_codon_dna_seq)
            
            org_codon = AminoAcidCodon(cut_codon_dna_seq)
            
            # getting synonymous mutations sequences
            if (donor_mut_name == 'synonymous'):
                if org_codon.is_start_codon():
                    codons = []
                    num_requested_ok = False
                    warnings.warn("can NOT do synonymous mutations in a start codon. codon: %s, edit_chr_pos: %d, nt_before_cut: %d, i_first_nt_in_codon: %d, org_donor_seq: %s" % 
                                  (cut_codon_dna_seq, edit_chr_pos, nt_before_cut, i_first_nt_in_codon, org_donor_seq))
                else:
                    codons, num_requested_ok = org_codon.get_K_synonymous_mut_codon_seqs(num_donor_variants)
            elif (donor_mut_name == 'nonsynonymous'):
                if org_codon.is_start_codon() or org_codon.is_stop_codon():
                    codons = []
                    num_requested_ok = False
                    warnings.warn("can NOT do nonsynonymous mutations in a start/stop codons. codon: %s, edit_chr_pos: %d, nt_before_cut: %d, i_first_nt_in_codon: %d, org_donor_seq: %s" % 
                                  (cut_codon_dna_seq, edit_chr_pos, nt_before_cut, i_first_nt_in_codon, org_donor_seq))
                else:
                    codons, num_requested_ok = org_codon.get_K_nonsynonymous_mut_codon_seqs(num_donor_variants)
            elif (donor_mut_name == 'nonsense'):
                if org_codon.is_stop_codon(): # this should happen only if early stop codon exist
                    codons = []
                    num_requested_ok = False
                    warnings.warn("can NOT do nonesense mutations in a stop codon. codon: %s, edit_chr_pos: %d, nt_before_cut: %d, i_first_nt_in_codon: %d, org_donor_seq: %s" % 
                                  (cut_codon_dna_seq, edit_chr_pos, nt_before_cut, i_first_nt_in_codon, org_donor_seq))
                else:
                    codons, num_requested_ok = org_codon.get_K_nonsense_mut_codon_seqs(num_donor_variants)
            elif (donor_mut_name == 'stop2aa'):
                if not org_codon.is_stop_codon():
                    codons = []
                    num_requested_ok = False
                    warnings.warn("can NOT do stop2aa mutations not in a stop codon. codon: %s, edit_chr_pos: %d, nt_before_cut: %d, i_first_nt_in_codon: %d, org_donor_seq: %s" % 
                                  (cut_codon_dna_seq, edit_chr_pos, nt_before_cut, i_first_nt_in_codon, org_donor_seq))
                else:
                    codons, num_requested_ok = org_codon.get_K_stop2aa_mut_codon_seqs(num_donor_variants)
            elif (donor_mut_name == '1bp'):
                
                # donor_mut_infoStr, donor_mut_subtype
                # nonsense_requiresOptionalSynonymous
                # 1bp_synonymous_requiresOptionalNonesense
                
                # this should happen only if early stop codon exist
                if (  ( (donor_mut_subtype == 'nonsense'   and donor_mut_infoStr == "requiresOptionalSynonymous") or
                        (donor_mut_subtype == 'synonymous' and donor_mut_infoStr == "requiresOptionalNonesense") or
                        (donor_mut_subtype == 'nonsynonymous' and donor_mut_infoStr == "requiresOptionalNonesenseAndSynonymous") ) and #requires that the codon can be mutated in the position to syn/nonsense
                        (not org_codon.is_mutable_2_synonymous_and_nonsese(edit_aa_frame) )  ):
                    codons = []
                    num_requested_ok = False
                    warnings.warn("can NOT do both synonymous and nonesense mutations in this codon. codon: %s, edit_chr_pos: %d, nt_before_cut: %d, i_first_nt_in_codon: %d, org_donor_seq: %s" % 
                                  (cut_codon_dna_seq, edit_chr_pos, nt_before_cut, i_first_nt_in_codon, org_donor_seq))
                else:
                    
                    if (donor_mut_subtype == 'synonymous'):
                        codons, num_requested_ok = org_codon.get_K_synonymous_mut_codon_seqs(num_donor_variants, edit_nt_i = edit_aa_frame)
                    elif (donor_mut_subtype == 'nonsynonymous'):
                        codons, num_requested_ok = org_codon.get_K_nonsynonymous_mut_codon_seqs(num_donor_variants, edit_nt_i = edit_aa_frame)
                    elif (donor_mut_subtype == 'nonsense'):
                        codons, num_requested_ok = org_codon.get_K_nonsense_mut_codon_seqs(num_donor_variants, edit_nt_i = edit_aa_frame)
                    else:
                        raise ValueError("get_donor_mut_for_guide - unknown mutation name:" + donor_mut_name)
            
            
            else:
                raise ValueError("get_donor_mut_for_guide - unknown 1bp mutation subtype name:" + donor_mut_subtype)
            
            if not num_requested_ok:
                warnings.warn("can not create that num_donor_variants: " + str(num_donor_variants) + " " + donor_mut_name + " mutations for codon:" +  str(cut_codon_dna_seq))
       
        
            for idx,alt_codon in enumerate(codons):
                alt_donor_seq = org_donor_seq[0:i_first_nt_in_codon] + \
                                alt_codon + \
                                org_donor_seq[(i_first_nt_in_codon+3):len(org_donor_seq)]
            

                # convert to match the guide orientation
                if ( (not guide_is_negative_strand and self.is_neg_strand()) or (guide_is_negative_strand and not self.is_neg_strand())  ):
                    alt_donor_seq = alt_donor_seq.reverse_complement()
                
                if do_revcomp_donor:
                    alt_donor_seq = alt_donor_seq.reverse_complement()
                
                donor_info_str = donor_mut_name + ":" + \
                                 str(cut_codon_dna_seq)   + ">" + str(alt_codon) + ":" + \
                                 org_codon.get_aa_seq() + ">" + org_codon.get_codon_aa_seq(alt_codon)
                
                if scramble_guide_and_donor:
                    donor_info_str = donor_info_str + ":" + "scramble"
                    tmp_donor_seq_lst = list(str(alt_donor_seq))
                    random.shuffle(tmp_donor_seq_lst)
                    alt_donor_seq = Seq(''.join(tmp_donor_seq_lst),alphabet=generic_dna)

                # appending to donor sequences matrix
                cur_donor_line = pd.DataFrame({'Gene' : guide_gene, 'guide_id' : pd.Series(guide_id), 
                                  'donor_id' : pd.Series(guide_id + ':' + donor_mut_type + ':offset' + str(donor_seq_offset) + ':donorID'  + str(idx) + ':EditPosInGuide' +  str(mut_pos_in_guide)), 
                                  'donor_seq': pd.Series(str(alt_donor_seq)), 
                                  'donor_seq_shift' : pd.Series(int(donor_seq_offset)), 
                                  'donor_mut_pos_in_guide' : pd.Series(str(mut_pos_in_guide)), 
                                  'donor_info_str' : pd.Series(donor_info_str),
                                  'set_name' : pd.Series(str(set_name)) })
                
 
                out_guide_donor_df = out_guide_donor_df.append(cur_donor_line,ignore_index=True)

        else:
            raise ValueError('get_donor_mut_for_guide unknown donor_mut_type:' + donor_mut_type +  " and donor_mut_name:" + donor_mut_name)
          

        return(out_guide_donor_df.reset_index(drop=True))
        



####################################################################################################################
def design_guide_and_donor_for_genes(input_design_table_filename, input_annot_gff_filename, input_genome_fasta_filename,
                                    PAM_seq = 'GA', excluded_seqs = ['TTTT'], BOWTIE_exe = '',
                                    off_targets_min_mismatch = 10,
                                    min_azimuth_score = 0.5,
                                    sort_by = '5to3', 
                                    min_dist_cut_to_donor_edge = 30,
                                    filter_out_donors_containing_excluded_seqs = False,
                                    max_guide_pos_frac_in_gene = 1.1,
                                    output_guides_df_filename = '', output_donors_df_filename = ''):
    """
    design_guide_and_donor_for_genes - A function for designing guides and donor sequences for genes according to the disgn table  
    input_design_table_filename = hand created table that defined the guides and donor seuences to be designed 
                                  contains the following columns (can contain other columns order does not matter):
                                  1. gene_id, example: YAL001C (gene id)
                                  2. guide_num, example: 5 (int)
                                  3. donor_mut_type, example:  nonsynonymous | synonymous | nonsense | frameshift_<in|out>_<1|2>
                                  4. num_donor_variants, example: 3 (will have a warning if can't desing that much)
                                  5. mut_pos_in_guide, example: currently only 'None' (which is represented by an empty string) that mean taht the mutation is in the cut site
                                  6. donor_length, example: 100 the full length of the donor sequence
                                  7. donor_seq_offsets, example: [-5,0,5] a list of shifts of the donor sequence center to create more constructs that do the same action 
    
    input_annot_gff_filename - sgd gene gff file (cleaned)
    input_genome_fasta_filename - sgd genome sequence in fasta file format, 
                                  NOTICE: the fasta files should have a bowtie index ready if bowtie is used for the mapping 
    
    # parameters for guide selection
    pam_seq - two letter PAM sequence ('GA')
    excluded_seqs - list of sequences that if in the donor or guide will make it nonfunctional
    BOWTIE_exe - bowtie 2 executable path. For mapping_cmd the command used for mapping the off tragets if empty string will use the default 
    
    off_targets_min_mismatch  -  min number of mismatch that are allowed in off targets (10, above 3 means no off targets allowed)
    excluded_seqs - used NOT as a filter just to add a column stating if contains excluded sequences
    min_azimuth_score - min of ok Azimuth score (0.5). Below it the guides are filtered  
    sort_by = how to sort the guides from best to worst. currently only '5to3' is implemented (select guides according to their position in the genes)
    
    # parameter for donor design
    min_dist_cut_to_donor_edge - parameter that control that donor with insertions have overlap with the genome. if not, returns empty data frame and warn
    
    # filtering guides and donors
    filter_out_donors_containing_excluded_seqs = filter donors that contain excluded sequences (True or false)
    max_guide_pos_frac_in_gene = max position of a guide in a gene. 5'->3', fraction of total gene length (between 0-1; larger than one or None means no filter)
    
    # output file names
    output_guides_df_filename - file name for designed guides dataframe (if empty, does not write)
    output_donors_df_filename - file name for designed donor sequences dataframe (if empty, does not write)
    
    # output
    Outputs the guides and donor sequences dataframes
    
    """
    
    ##################
    # processing input
    ##################

    # fixing the mapping command to be the default bowtie
    if not BOWTIE_exe or len(BOWTIE_exe) < 3:
        BOWTIE_exe = '/Users/eilon/software/bowtie2-2.2.8/bowtie2'
    
    # defining the bowtie mapping cmd
    bowtie_cmd = shlex.split(BOWTIE_exe + ' -x ' + os.path.splitext(input_genome_fasta_filename)[0] + ' -U - -f -D 20 -R 3 -N 1 -L 10 -i S,1,0.50 --gbar 3 --end-to-end -k 30 --no-head -t  --rdg 10,6 --rfg 10,6')
    mapping_cmd = bowtie_cmd

    
    # loading design matrix 
    # fields: gene_id	guide_num	donor_mut_type	num_donor_variants	mut_pos_in_guide	donor_length	donor_seq_offsets

    design_df = pd.read_table(input_design_table_filename, sep='\t', na_values = "")
    design_df['donor_seq_offsets'] = design_df['donor_seq_offsets'].apply(ast.literal_eval)
    design_df['donor_seq_offsets'] = design_df['donor_seq_offsets'].apply(np.array,dtype=np.float)
    
    design_df['gene_id'] = design_df['gene_id'].str.strip()
    design_df['donor_mut_type'] = design_df['donor_mut_type'].str.strip()
    
    
    # loading gene gff matrix
    print "loading genes gff file: " + input_annot_gff_filename
    genes_gff_df = sgd_gff2dataframe(input_annot_gff_filename, ['CDS'])
    
    # loading genome fasta file
    print "loading genome fasta file: " + input_genome_fasta_filename
    genome_seq = SeqIO.to_dict(SeqIO.parse(open(input_genome_fasta_filename),'fasta', alphabet=generic_dna))
 

    
    # init output dataframes
    out_guide_df = pd.DataFrame(data=None)
    out_guide_donor_df = pd.DataFrame(data=None)
    
    ####################
    # running gene by gene such that the donor sequences will use the same guides
    # if differnet number of guides are specified for each donor desing the guides will be selected from best to worst
    ####################
    
    # grouping by gene
    design_grouped = design_df.groupby('gene_id')
    

    # iterating over the genes
    for cur_gene_id, cur_gene_df in design_grouped:

        print "Designing guides for gene:" + cur_gene_id
        
        # the current gene object
        cur_gene = CRISPR_Gene(cur_gene_id, genes_gff_df, genome_seq)
        
        # if differnet number of guides are specified for each donor desing the guides will be selected from best to worst
        cur_gene_max_guide_num =  max(cur_gene_df['guide_num'])
        
        # if even one row requires the guide to be entirely in the CDS - find only guides that are entirely in the CDS
        #cur_require_entire_guide_to_be_in_CDS = np.any(~ np.isnan(cur_gene_df['mut_pos_in_guide'].values))
        # np.any(cur_gene_df['require_entire_guide_in_CDS'])
        #
        cur_min_mut_pos_in_guide = min(int(-cur_gene.CRISPR_CUT_INDEX - 1*(cur_gene.is_neg_strand())),np.nanmin(cur_gene_df['mut_pos_in_guide'].values))
        cur_max_mut_pos_in_guide = max(2,np.nanmax(cur_gene_df['mut_pos_in_guide'].values))
        
        
        # max half region length to look for excluded sequences (should be at least ~20 to test the guide)
        

        cur_max_donor_seq_offset = np.nanmax(np.abs(cur_gene_df['donor_seq_offsets'].apply(np.nanmax).values))
        if np.isnan(cur_max_donor_seq_offset):
            cur_max_donor_seq_offset = 0
        
        
        
        cur_max_seq_len_around_cut =  int( np.ceil(max(cur_gene_df['donor_length']) / 2) + \
                                           cur_max_donor_seq_offset)
        
        #print "Extracting all guides"
        
        # get all guides (computing filter guides that contain excluded sequences in the sequences around them)
        cur_all_gene_guides_df = cur_gene.get_all_guides_that_cut_in_cds(pam = PAM_seq, seq_len_around_cut = cur_max_seq_len_around_cut, 
                                                                         min_mut_pos_in_guide = cur_min_mut_pos_in_guide,
                                                                         max_mut_pos_in_guide = cur_max_mut_pos_in_guide,
                                                                         excluded_seqs = excluded_seqs, mapping_cmd = mapping_cmd,
                                                                         sort_by = sort_by)
        
        #DEBUG
        #print "before get_K_best_guide_ids" #DEBUG
        #print sum(cur_all_gene_guides_df['guide_id']== 'YAL001C_pS_33') #DEBUG
        
        #print cur_all_gene_guides_df[ cur_all_gene_guides_df['guide_id']== 'YAL001C_pS_33' ]
        #print cur_all_gene_guides_df[ cur_all_gene_guides_df['guide_id']== 'YAL001C_pS_20' ]
        #print cur_gene_max_guide_num
        #print off_targets_min_mismatch
        #print min_azimuth_score
        #print max_guide_pos_frac_in_gene

        # select top 'best' guides
        selected_guide_ids = cur_gene.get_K_best_guide_ids(guides_df = cur_all_gene_guides_df, K = cur_gene_max_guide_num, 
                                                           off_targets_min_mismatch = off_targets_min_mismatch,
                                                           min_azimuth_score = min_azimuth_score,
                                                           max_guide_pos_frac_in_gene = max_guide_pos_frac_in_gene,
                                                           sort_by = sort_by)
        #print "after get_K_best_guide_ids" #DEBUG
        #for gi in selected_guide_ids:
        #    print gi
        #print ('YAL001C_pS_33' in list(selected_guide_ids) ) 
        #print 'XXX-------------------' #DEBUG
        
        # selected guides (each guide should have at least one donor designed)
        cur_selected_guide_df = cur_all_gene_guides_df[ cur_all_gene_guides_df['guide_id'].isin( list(selected_guide_ids) ) ]
           
        print "--- Designing donor sequences for gene:" + cur_gene_id + ", for " + str(len(selected_guide_ids)) + " guides (out of " + str(cur_gene_max_guide_num) + " requested)"
        #print cur_selected_guide_df[ cur_selected_guide_df['guide_id'] == 'YAL001C_pS_33' ]
        
        # concating with the output dataframe
        out_guide_df = pd.concat([out_guide_df,cur_selected_guide_df],ignore_index=True)

        
        # desinging the donor sequences
        for idx,row in cur_gene_df.iterrows():
            
            if len(selected_guide_ids) <  int(row['guide_num']):
                warnings.warn("There are NOT enough guides in --- %s --- for the design. There are %d guides and the design if for %d" % (str(row["gene_id"]),len(selected_guide_ids), row['guide_num'] ))
            
            if len(selected_guide_ids) <= 0:
                continue
            
            cur_selected_guide_ids = selected_guide_ids.iloc[range( min(int(row['guide_num']),len(selected_guide_ids))  )]
            
            
            # do reverse complement for donor sequences
            if 'do_revcomp_donor' in cur_gene_df.columns:
                do_revcomp_donor = bool(row['do_revcomp_donor']==True)
            else:
                do_revcomp_donor = False
                
            
            
            
            # do reverse complement for donor sequences
            if 'do_scramble_guide_and_donor' in cur_gene_df.columns:
                scramble_guide_and_donor = bool(row['do_scramble_guide_and_donor']==True)
            else:
                scramble_guide_and_donor = False
                
            #print "do_revcomp_donor %d" % (do_revcomp_donor)
            #print "scramble_guide_and_donor %d" % (scramble_guide_and_donor)
            
            cur_all_gene_guides_df_fordonors = cur_all_gene_guides_df
            
            # permuting the guides and adding it to the guide df
            if scramble_guide_and_donor:
                print 'Scramble donor guides...'
                cur_all_gene_guides_df_scramble = cur_all_gene_guides_df[ cur_all_gene_guides_df['guide_id'].isin(cur_selected_guide_ids) ].copy()
                # scrambling the guides and updating their ids
                cur_all_gene_guides_df_scramble['guide_id'] = cur_all_gene_guides_df_scramble['guide_id'] + '_scramble'
                scramble_func = lambda x: ''.join(random.sample(str(x),len(str(x))))
                cur_all_gene_guides_df_scramble['guide'] = cur_all_gene_guides_df_scramble['guide'].apply(scramble_func)
                cur_all_gene_guides_df_scramble['guide_noPAM'] = cur_all_gene_guides_df_scramble['guide_noPAM'].apply(scramble_func)
                cur_all_gene_guides_df_scramble['guide_PAM_p7'] = cur_all_gene_guides_df_scramble['guide_PAM_p7'].apply(scramble_func)
                cur_all_gene_guides_df_scramble['guide_PAM_m4p3'] = cur_all_gene_guides_df_scramble['guide_PAM_m4p3'].apply(scramble_func)
                # adding the scrambeles guides to the guides table
                out_guide_df = pd.concat([out_guide_df,cur_all_gene_guides_df_scramble],ignore_index=True)

                
                # for donors design
                cur_all_gene_guides_df_fordonors = cur_all_gene_guides_df_scramble
                cur_selected_guide_ids = cur_selected_guide_ids + '_scramble'
            
            
            cur_gene_donor_df = cur_gene.get_donor_mut_for_guides(cur_all_gene_guides_df_fordonors, cur_selected_guide_ids, 
                                                                  donor_mut_type=str(row['donor_mut_type']), 
                                                                  num_donor_variants = int(row['num_donor_variants']), 
                                                                  mut_pos_in_guide = row['mut_pos_in_guide'], 
                                                                  donor_length=int(row['donor_length']), 
                                                                  donor_seq_offsets = row['donor_seq_offsets'],
                                                                  set_name = row['set_name'],
                                                                  min_dist_cut_to_donor_edge = min_dist_cut_to_donor_edge,
                                                                  excluded_seqs = excluded_seqs,
                                                                  do_revcomp_donor = do_revcomp_donor,
                                                                  scramble_guide_and_donor = scramble_guide_and_donor)
            
            out_guide_donor_df = pd.concat([out_guide_donor_df,cur_gene_donor_df],ignore_index=True)

        #print "---------------------------- Finished designing guide for the above gene"


    # adding for each guide tis location in the gene (now added in )
    #out_guide_df["guide_cut_gene_pos_frac"] = out_guide_df["guide_cut_gene_nt_pos"] / out_guide_df["CDS_len_nts"]

    # filtering out donor sequences with excluded sequences
    if filter_out_donors_containing_excluded_seqs:
        if out_guide_donor_df.shape[1] > 1: # if not null df
            out_guide_donor_df = out_guide_donor_df[ (~ out_guide_donor_df["contain_excluded_sequences"]).values ]
        
        
    if len(output_guides_df_filename) > 3:
        print "saving guides to: " +  output_guides_df_filename
        out_guide_df.to_csv(output_guides_df_filename, sep='\t', index = False)
        
        
    
    if len(output_donors_df_filename) > 3:
        print "saving donor sequences to: " +  output_donors_df_filename
        out_guide_donor_df.to_csv(output_donors_df_filename, sep='\t', index = False)
            
    
    return( (out_guide_df, out_guide_donor_df) )

