This workflow maps paired ends sequencing reads to the library in order to get guide/donor oligo counts for each sample.

## installing
Please install the dependencies in crispey_seq2counts_requirments.txt
using Miniconda 3. 

## Running
When running the Snakefile please provide a configuration file similar to the files under config_files. The config file points to a samples file (see examples). This file contains the required information about each sample (see below)

We recommend using raw reads mapped by STAR with zero mismatch parameter.

## Samples table columns are:

1. sample_index – unique sample indes
2. induction_time – time (h; float) from the start of the editing (competition sample should contain the details of the last editing sample that the competition started from) 
3. induction_time_index – time running index (int) from the start of the editing (competition sample should contain the details of the last editing sample that the competition started from)
4. induction_generation – generation (culture doubling; float) from the start of the editing (competition sample should contain the details of the last editing sample that the competition started from)
5. induction_replicate – editing replicate number (int) (competition sample should contain the details of the last editing sample that the competition started from)

6. competition_time - time (h; float) from the start of the growth competition
7. competition_time_index - time running index (int) from the start of the growth competition
8. competition_generation - generation (culture doubling; float) from the start of the growth competition
9. competition_replicate – growth competition replicate number (int)
10. competition_medium – name of growth competition media (str)
11. cell_num_start – number of cells in the previous time point, based on volume and OD measurements
12. cell_num_end - number of cells in the current time point, based on volume and OD measurements

13. growth_ratio-  cell_num_end /  cell_num_start
14. OD – measured OD
15. Sequencing_reads_R1 – path to sequencing read 1 fastq file
16. Sequencing_reads_R2 – path to sequencing read 2 fastq file

Column order is not flexible and the table can contain any additional columns
# Snakemake workflow rules graph
![sequences2counts](crispey_seq2counts_rule_graph.png)
