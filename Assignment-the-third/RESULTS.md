# Demultiplexing and Index Swapping - Results #
### Total number of records: 363246735 ###
### Matched Pairs: 331755033,  91.3% ###
### Unknown or Low Quality Pairs (if any base has qscore < 30): 30783962, 8.5% ###
### Index-hopped Pairs: 707740, 0.2% ###

### Percentage of Matched Reads from Each Sample ###
| sample number | group |	treatment |	percent |
| ------------- | ----- | --------- | ------- |
| 1	| 2A	| control	| 2.2 |
| 2	| 2B	| control	| 1.5 |
| 3	| 2B	| control	| 1.8 |
| 4	| 2C	| mbnl	| 2.4 |
| 6	| 2D	| mbnl	| 2.9 |
| 7	| 2E	| fox	| 1.4 |
| 8	| 2F	| fox	| 9.6 |
| 10	| 2G	| both	| 21.0 |
| 11	| 2H	| both	| 4.8 |
| 14	| 3B	| control	| 1.2 |
| 15	| 3C	| mbnl	| 2.0 |
| 16	| 3D	| mbnl	| 2.2 |
| 17	| 3E	| fox	| 3.1 |
| 19	| 3F	| fox	| 4.3 |
| 21	| 3G	| both	| 2.4 |
| 22	| 3H	| both	| 1.1 |
| 23	| 4A	| control	| 11.6 |
| 24	| 4A	| control	| 2.8 |
| 27	| 4C	| mbnl	| 1.9 |
| 28	| 4D	| mbnl	| 3.2 |
| 29	| 4E	| fox	| 1.3 |
| 31	| 4F	| fox	| 1.0 |
| 32	| 4G	| both	| 3.1 |
| 34	| 4H	| both	| 2.4 |


### Output Files ###
12619558_demult.txt (job number)
- printed contexts of dictionary containing index info and confirmed all files are closed at the end

summary_output_S1_L008.txt
- summarized results from the run
- Table with Percentage of Matched Reads from Each Sample

indexcombos_counts_S1_L008.tsv
- tab separated table with the counts for every possible combination of indexes, including matches

zipped output fastq files can be found on Talapas at 
```
/projects/bgmp/kroth/bioinfo/Bi622/demultiplexing/assignm_third/run3-success/
```
