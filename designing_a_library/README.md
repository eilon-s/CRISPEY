# crisprHDRlibrary

## Library design

The library was designed by running crispr_design_library.py.
Some paths of input files may need to be modified (input files are in the input directory)

## samples table

The samples tables required columns are:
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

