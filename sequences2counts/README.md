Mapping paired ends sequencing reads to get oligo counts

Snakefile requires a configuration file as input.
The config file points to samples file (see examples).
We found that raw reads mapped by STAR with zero mismatch parameter give good performances.

## The samples tables required columns are:

1. sample_index
2. induction_time
3. induction_time_index
4. induction_generation
5. induction_replicate
6. competition_time
7. competition_time_index
8. competition_generation
9. competition_replicate
10. competition_medium
11. cell_num_start
12. cell_num_end
13. growth_ratio
14. OD
15. Sequencing_reads_R1
16. Sequencing_reads_R2

Column order is not felxible and the table can contain any additional columns

# Snakemake workflow rules graph
![sequences2counts](crispey_seq2counts_rule_graph.png)

