# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 |
| 1294_S1_L008_R2_001.fastq.gz | index1 |
| 1294_S1_L008_R3_001.fastq.gz | index2 |
| 1294_S1_L008_R4_001.fastq.gz | read2 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. For quality scores of indexes, if a single base in the sequence has a qscore less than or equal to 12, then consider it low quality. Using a qscore of 12 as the cutoff would mean there must be a 94% - 95% chance that each and every base call is correct for it to be used for sample identification. For quality scores of the reads, the standard for average qscore cutoff seems to be 30, meaning the average probability that the bases are correct must be 99.9% for it to be analyzed downstream.
    3. For the index1 reads, 3976613 of the 363246735 (1.095%) total index1 sequences have at least 1 N. For index2 reads, 3328051 (0.9162%) have at least 1 N.
    ```
    $ /usr/bin/time -v zcat 1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | grep -E '[N]+' | wc -l
    $ /usr/bin/time -v zcat 1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | grep -E '[N]+' | wc -l
    ```

## Part 2
1. Define the problem
2. Describe output
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [4 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
