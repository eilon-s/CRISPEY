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

sys.path.append(os.path.expanduser('~/bin/python/crispr'))
from extract_guides_functions import *



class AminoAcidCodon(object):
    """ Holds a single amino acid codon 
        can mutate the DNA sequence to change the codon sequence
    """
    
    ###############################################################################################################
    def __init__(self, codon_dna_seq):
        """
        Initialize CRISPR_Gene object.
        
        Parameters:
        -----------
        codon_dna_seq : codon sequence (A/C/G/T string length 3)
      
        Attributes:
        -----------
        dna_seq : Seq(,generic_dna)
            codon DNA sequence
            gene name
        """
        
        # amino acid single letter code
        self.AA_dict = { \
                        "A":["GCA","GCC","GCG","GCT"], \
                        "B":["AAC","AAT","GAC","GAT"], \
                        "C":["TGC","TGT"], \
                        "D":["GAC","GAT"], \
                        "E":["GAA","GAG"], \
                        "F":["TTC","TTT"], \
                        "G":["GGA","GGC","GGG","GGT"], \
                        "H":["CAC","CAT"], \
                        "I":["ATA","ATC","ATT"], \
                        "K":["AAA","AAG"], \
                        "L":["CTA","CTC","CTG","CTT","TTA","TTG"], \
                        "M":["ATG"], \
                        "N":["AAC","AAT"], \
                        "P":["CCA","CCC","CCG","CCT"], \
                        "Q":["CAA","CAG"], \
                        "R":["AGA","AGG","CGA","CGC","CGG","CGT"], \
                        "S":["AGC","AGT","TCA","TCC","TCG","TCT"], \
                        "T":["ACA","ACC","ACG","ACT"], \
                        "V":["GTA","GTC","GTG","GTT"], \
                        "W":["TGG"], \
                        "Y":["TAC","TAT"], \
                        "Z":["CAA","CAG","GAA","GAG"], \
                        "*":["TAA","TAG","TGA"] \
                       }
        
        #self.AA_dict = { \
        #                "A":"GCA","A":"GCC","A":"GCG","A":"GCT", \
        #                "B":"AAC","B":"AAT","B":"GAC","B":"GAT", \
        #                "C":"TGC","C":"TGT", \
        #                "D":"GAC","D":"GAT", \
        #                "E":"GAA","E":"GAG", \
        #                "F":"TTC","F":"TTT", \
        #                "G":"GGA","G":"GGC","G":"GGG","G":"GGT", \
        #                "H":"CAC","H":"CAT", \
        #                "I":"ATA","I":"ATC","I":"ATT", \
        #                "K":"AAA","K":"AAG", \
        #                "L":"CTA","L":"CTC","L":"CTG","L":"CTT","L":"TTA","L":"TTG", \
        #                "M":"ATG", \
        #                "N":"AAC","N":"AAT", \
        #                "P":"CCA","P":"CCC","P":"CCG","P":"CCT", \
        #                "Q":"CAA","Q":"CAG", \
        #                "R":"AGA","R":"AGG","R":"CGA","R":"CGC","R":"CGG","R":"CGT", \
        #                "S":"AGC","S":"AGT","S":"TCA","S":"TCC","S":"TCG","S":"TCT", \
        #                "T":"ACA","T":"ACC","T":"ACG","T":"ACT", \
        #                "V":"GTA","V":"GTC","V":"GTG","V":"GTT", \
        #                "W":"TGG", \
        #                "Y":"TAC","Y":"TAT", \
        #                "Z":"CAA","Z":"CAG","Z":"GAA","Z":"GAG", \
        #                "*":"TAA","*":"TAG","*":"TGA" \
        #               }
        
        self.codon_table = self.__fill_codon_table()
        self.codon2aa    = self.__fill_codon2aa()
        
        # saving input parameters (pointers)
        self.dna_seq = Seq(str(codon_dna_seq),generic_dna)
        
        
    ###############################################################################################################
    def __repr__(self):
        return "AminoAcidCodon()"
    
    ###############################################################################################################
    def __str__(self):
        return "DNA codon: " + str(self.dna_seq)
    
    ###############################################################################################################
    def __fill_codon2aa(self):
        codon2aa = {}
        for aa, codons in self.AA_dict.iteritems():
            for codon in codons:
                codon2aa[codon] = aa
        return(codon2aa)
    
    ###############################################################################################################
    def __fill_codon_table(self):
        codon_table = np.chararray( (4,4,4) )
        for aa, codons in self.AA_dict.iteritems():
            #print aa
            for codon in codons:
                #print codon
                codon_table[self.__nt2int(codon[0]), \
                            self.__nt2int(codon[1]), \
                            self.__nt2int(codon[2])] = aa
        return(codon_table)

    ###############################################################################################################
    def __nt2int(self, nt_char):
        if nt_char == 'A':
            return(int(0))
        elif nt_char == 'C':
            return(int(1))
        elif nt_char == 'G':
            return(int(2))
        elif nt_char == 'T':
            return(int(3))
        else:
            raise ValueError(nt_char + " unknown nucleotide")
    
    ###############################################################################################################
    def __get_aa_seq(self):
        return(self.codon2aa[str(self.dna_seq)])
    
    ###############################################################################################################
    def get_aa_seq(self):
        return(self.codon2aa[str(self.dna_seq)])
    
    ###############################################################################################################
    def get_codon_aa_seq(self, codon):
        return(self.codon2aa[str(codon)])
  
    ###############################################################################################################
    def __int2nt(self, nt_int):
        if (nt_int < 0 or 3 < nt_int):
            raise ValueError(str(nt_int) + " not in 0-3 range")
        
        return('ACGT'[nt_int])
    
    ###############################################################################################################
    def is_stop_codon(self):
        return(self.__get_aa_seq() == "*")
    
    ###############################################################################################################
    def is_start_codon(self):
        return(self.__get_aa_seq() == "M")
    
    
    ###############################################################################################################
    def __keep_only_codons_differing_in_position(self,aa_codons, edit_nt_i = None):
        """
        if edit_nt_i is not None - returns only codon that differ by single bp mutations in the input position (0-2)
        """
        
        if (edit_nt_i is not None):
            org_codon_str = str(self.dna_seq)
            similarity_postios = [0,1,2]
            similarity_postios.remove(edit_nt_i)    
            # removing codons that are not similar in the other two positions 
            aa_codons  = [mut_codon for mut_codon in aa_codons if mut_codon[similarity_postios[0]] == org_codon_str[similarity_postios[0]] and mut_codon[similarity_postios[1]] == org_codon_str[similarity_postios[1]] ]

        return(aa_codons)
        
    
    
    ###############################################################################################################
    def get_all_synonymous_mut_codon_seqs(self, edit_nt_i = None):
        aa_codons = list(self.AA_dict[self.__get_aa_seq()])
        aa_codons.remove(self.dna_seq)
        
        aa_codons = self.__keep_only_codons_differing_in_position(aa_codons, edit_nt_i)
        
        return(aa_codons)
    
    ###############################################################################################################
    def get_a_synonymous_mut_codon_seq(self, edit_nt_i = None):
        
        aa_codons = self.get_all_synonymous_mut_codon_seqs(edit_nt_i)
        
        if len(aa_codons) == 0:
            return None
        else:
            return(random.choice(aa_codons))
        
    ###############################################################################################################
    def get_K_synonymous_mut_codon_seqs(self,K, edit_nt_i = None):
        """
        returns a tuple 
        the first element is minimum of K or the available number of mutations randomly selection
        the second element is whether K was not larger than the total number of possible mutations
        
        if edit_nt_i is not None - returns only single bp mutations that can in this position (0-2)
        """
        
        aa_codons = self.get_all_synonymous_mut_codon_seqs(edit_nt_i)

        if (K>len(aa_codons)):
            warnings.warn("get_K_synonymous_mut_codon_seqs: number of optional codon %d is smaller than the requested %d" % (len(aa_codons), K) )
       
        
        if len(aa_codons) == 0:
            return ( ([],K<=0) )
        else:
            return( (random.sample(aa_codons,min(K,len(aa_codons))), K<=len(aa_codons) ) )
    
    ###############################################################################################################
    def get_all_nonsynonymous_mut_codon_seqs(self, edit_nt_i = None):
        """ returns all codon sequences that are not of the same AA or stop codons 
            notice that differnt AA have different number of codons and therefore can be mutated to at different rates
        """
        
        if self.__get_aa_seq() == "*":
            raise ValueError("Can not do nonsynonymous mutation - the original codon is a stop codon")
        
        aa_codons = list()
        
        for aa, codons in self.AA_dict.iteritems():
            if aa != self.__get_aa_seq() and aa != "*":
                aa_codons = aa_codons + codons
        
        aa_codons = self.__keep_only_codons_differing_in_position(aa_codons, edit_nt_i)
        
        return(aa_codons)
    
    ###############################################################################################################
    def get_a_nonsynonymous_mut_codon_seq(self, edit_nt_i = None):
        aa_codons = self.get_all_nonsynonymous_mut_codon_seqs(edit_nt_i)
        
        if len(aa_codons) == 0:
            return None
        else:
            return(random.choice(aa_codons))
        
    ###############################################################################################################
    def get_K_nonsynonymous_mut_codon_seqs(self,K, edit_nt_i = None):
        """
        returns a tuple 
        the first element is minimum of K or the available number of mutations randomly selection
        the second element is whether K was not larger than the total number of possible mutations
        """
        aa_codons = self.get_all_nonsynonymous_mut_codon_seqs(edit_nt_i)
        
        if (K>len(aa_codons)):
            warnings.warn("get_K_nonsynonymous_mut_codon_seqs: number of optional codon %d is smaller than the requested %d" % (len(aa_codons), K) )
            
        if len(aa_codons) == 0:
            return ( ([],K<=0) )
        else:
            return( (random.sample(aa_codons,min(K,len(aa_codons))), K<=len(aa_codons) ) )
        
    ###############################################################################################################
    def get_all_nonsense_mut_codon_seqs(self, edit_nt_i = None):
        """ returns all stop codons 
        """
        
        if self.__get_aa_seq() == "*":
            raise ValueError("Can not do nonesense mutation - the original codon (%s) is a stop codon" % ( str(self.dna_seq) ))
        
        aa_codons = list(self.AA_dict['*'])
        
        aa_codons = self.__keep_only_codons_differing_in_position(aa_codons, edit_nt_i)
        
        return(aa_codons)
    
    ###############################################################################################################
    def get_a_nonsense_mut_codon_seq(self, edit_nt_i = None):
        aa_codons = self.get_all_nonsense_mut_codon_seqs(edit_nt_i)
        
        if len(aa_codons) == 0:
            return None
        else:
            return(random.choice(aa_codons))
        
    ###############################################################################################################
    def get_K_nonsense_mut_codon_seqs(self,K, edit_nt_i = None):
        """
        returns a tuple 
        the first element is minimum of K or the available number of mutations randomly selection
        the second element is whether K was not larger than the total number of possible mutations
        """
        aa_codons = self.get_all_nonsense_mut_codon_seqs(edit_nt_i)
        
        if (K>len(aa_codons)):
            warnings.warn("get_K_nonsense_mut_codon_seqs: number of optional codon %d is smaller than the requested %d" % (len(aa_codons), K) )
            
        if len(aa_codons) == 0:
            return ( ([], K<=0) )
        else:
            return( (random.sample(aa_codons,min(K,len(aa_codons))), K<=len(aa_codons) ) )
    
    ###############################################################################################################
    def is_mutable_2_synonymous_and_nonsese(self,edit_aa_nt):        
        """
        Is this codon can be mutated by single nt replacement in the position to both synonymous and stop codopn
        This will be true only for the third (=2) position of Cys and Tyr codons
        """
        return( (edit_aa_nt == 2) and (self.get_aa_seq() in ['Y','C']) )


    ###############################################################################################################
    def get_all_stop2aa_mut_codon_seqs(self, edit_nt_i = None):
        """
        returns all codon sequences that are some AA and not stop codon (the original codon is a stop codon)
        
        """
        
        if self.__get_aa_seq() != "*":
            raise ValueError("Can not do stop2aa mutation - the original codon has to be a stop codon %s" % (self.__get_aa_seq()))
        
        aa_codons = list()
        
        for aa, codons in self.AA_dict.iteritems():
            if aa != self.__get_aa_seq():
                aa_codons = aa_codons + codons
        
        aa_codons = self.__keep_only_codons_differing_in_position(aa_codons, edit_nt_i)
        
        return(aa_codons)
    
    ###############################################################################################################
    def get_a_stop2aa_mut_codon_seq(self, edit_nt_i = None):
        aa_codons = self.get_all_stop2aa_mut_codon_seqs(edit_nt_i)
        
        if len(aa_codons) == 0:
            return None
        else:
            return(random.choice(aa_codons))
        
    ###############################################################################################################
    def get_K_stop2aa_mut_codon_seqs(self,K, edit_nt_i = None):
        """
        returns a tuple 
        the first element is minimum of K or the available number of mutations randomly selection
        the second element is whether K was not larger than the total number of possible mutations
        """
        aa_codons = self.get_all_stop2aa_mut_codon_seqs(edit_nt_i)
        
        if (K>len(aa_codons)):
            warnings.warn("get_K_stop2aa_mut_codon_seqs: number of optional codon %d is smaller than the requested %d" % (len(aa_codons), K) )
            
        if len(aa_codons) == 0:
            return ( ([],K<=0) )
        else:
            return( (random.sample(aa_codons,min(K,len(aa_codons))), K<=len(aa_codons) ) )

