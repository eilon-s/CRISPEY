#!/usr/bin/env python3

import pandas as pd
import os
#import pysam
from Bio import SeqIO
import gzip

def crispr_generate_sim_reads(cur_sample,read_length, ref_fa, pool_with_counts_sim_filename, output_fastq_gz_r1, output_fastq_gz_r2):
    
    read_length = int(read_length)
    
    record_dict = SeqIO.to_dict(SeqIO.parse(ref_fa, "fasta"))

    pool_with_counts_sim_df = pd.read_table(pool_with_counts_sim_filename, sep='\t', na_values = "")
    

    with gzip.open(output_fastq_gz_r1, 'wb') as out_r1:
            with gzip.open(output_fastq_gz_r2, 'wb') as out_r2:

                for r in range(pool_with_counts_sim_df.shape[0]):


                    if r % 1000 == 0:
                        print(r)

                    seq_id = str(pool_with_counts_sim_df.oligo_agilent_id[r]) + '#' + pool_with_counts_sim_df.set_name[r]
                    seq = record_dict[seq_id]


                    cur_num_reads = int(pool_with_counts_sim_df.ix[r,cur_sample])

                    for i in range(cur_num_reads):
                        #print('@R1#' + seq_id + "#" + str(i) + '\n')
                        out_r1.write( bytes('@' + seq_id + ":" + str(i) + '\n', encoding='utf-8'))
                        out_r1.write( bytes( str(seq.seq[0:(0+read_length)]) + '\n', encoding='utf-8'))
                        out_r1.write( bytes('+\n', encoding='utf-8'))
                        out_r1.write( bytes('G'*read_length + '\n', encoding='utf-8'))


                        out_r2.write(bytes('@' + seq_id + ":" + str(i) + '\n', encoding='utf-8'))
                        out_r2.write(bytes('GCTACGATCGATCTGT' + str(seq.reverse_complement().seq[0:(read_length-16)]) + '\n', encoding='utf-8'))
                        out_r2.write(bytes('+\n', encoding='utf-8'))
                        out_r2.write(bytes('G'*read_length + '\n', encoding='utf-8'))

# for running as part of a snakemake pipeline
crispr_generate_sim_reads(snakemake.params["cur_smaple"], snakemake.params["read_length"], snakemake.input["ref_fa"], snakemake.input["sim_merged_count_table"], snakemake.output["sim_fastq_R1_file"], snakemake.output["sim_fastq_R2_file"])


