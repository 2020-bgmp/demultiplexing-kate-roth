#!/usr/bin/env python

import argparse
import gzip
import itertools
import re


def get_args():
    '''get names/paths of input files from the command line'''
    parser = argparse.ArgumentParser(description="A program to get names/paths of input files from the command line")
    parser.add_argument("-i", "--indexesfile", help="tab separated text file with index information", required=True, type=str)
    parser.add_argument("-q", "--qscore", help="index qscore cutoff (per individual base pair", default=30, type=int)
    parser.add_argument("-a", "--a", help="R1 (read1) fastq file", required=True, type=str)
    parser.add_argument("-b", "--b", help="R2 (index1) fastq file", required=True, type=str)
    parser.add_argument("-c", "--c", help="R3 (index2) fastq file", required=True, type=str)
    parser.add_argument("-d", "--d", help="R4 (index4) fastq file", required=True, type=str)
    return parser.parse_args()

# store args
args = get_args()
indexesf = args.indexesfile
threshold = args.qscore
read1f = args.a
read2f = args.b
read3f = args.c
read4f = args.d

# extract sample number and lane number from input file to use for naming other files
split_path = re.split('/', read1f)
split_name = re.split("\.", split_path[-1])
split_prefix = re.split("\_", split_name[0])

# use for naming output files
input_info = split_prefix[1] + '_' + split_prefix[2]


def get_revcomp(sequence):
    '''Takes a nucleotide sequence and returns the the reverse complement'''
    compl_seq = ''
    for bp in sequence:
        if bp == 'A':
            compl_seq+='T'
        elif bp == 'T':
            compl_seq+='A'
        elif bp == 'C':
            compl_seq+='G'
        elif bp == 'G':
            compl_seq+='C'
        else:
            compl_seq+=bp       # if there is an N or other random character
    revcomp = compl_seq[::-1]   # reverse the string
    return revcomp
# print(get_revcomp('GCT'))
# print(get_revcomp('NAC'))


def check_qscore(phred_seq):
    '''Takes a sequence of phred scores, converts to integer qscores (phred33 encoding only)
    and returns boolean for whether any one base is less than the specified cutoff.
    True for greater than or equal to cutoff, False for less than cutoff'''
    for phred in phred_seq:
        qscore = ord(phred) - 33
        if qscore < threshold:
            return False
    return True
# print(check_qscore('>JJJJJ'))
# print(check_qscore('?JJJJJ'))


# extract and store indexes from indexes file
# keys are seq of index, values are tuple of other info from indexes file
indexes_dict = {}
with open(indexesf) as indf:
    ind_ln = 0
    for ind_line in indf:
        if ind_ln < 1:
            ind_ln+=1
        else:
            ind_line = ind_line.strip()
            split = re.split('\t', ind_line)
            # {index sequence:(indexID, sample, group, treatment)}
            indexes_dict[split[4]] = (split[3], split[0], split[1], split[2])
# print(indexes_dict)


# make dictionary of all index combos using itertools
# keys are the permutations combo, values are the count
# {('AAA,ATC'):0}
permut_count_dict = {}
permutations = itertools.product(indexes_dict, repeat=2)    # product to get identity combos ('AAA','AAA')
for p in list(permutations):
    permut_count_dict[p] = 0
# print(permut_count_dict)


# generate filenames for matched filenames
matched_filenames = []
for index_info in indexes_dict.values():   
    # indexID_sample_group_treatment_illumina sample number_illumina lane number_read1/2.fq
    fname1 = index_info[0] + "_" + index_info[1] + "_" + index_info[2] + "_" + index_info[3] + "_" + input_info + "_read1.fq"
    fname2 = index_info[0] + "_" + index_info[1] + "_" + index_info[2] + "_" + index_info[3] + "_" + input_info + "_read2.fq"
    matched_filenames.append(fname1)
    matched_filenames.append(fname2)
# print(matched_filenames)


# make filename abbreviations for opening matched files
matched_seq_keys = []
for ind_seq,ind_info in indexes_dict.items():
    seqkey1 = ind_seq + 'read1'
    seqkey2 = ind_seq + 'read2'
    matched_seq_keys.append(seqkey1)
    matched_seq_keys.append(seqkey2)
# print(matched_seq_keys)
# print(len(matched_filenames)==len(matched_seq_keys))


# open matched files to write out to
# keys are index sequence + 'read1/2'; values are open file object
# {'AAAread1':<wrapper file object thing>}
open_files_dict = {}
for i in range(len(matched_filenames)):
    open_files_dict[matched_seq_keys[i]] = open(matched_filenames[i], 'w')

# check open_files_dict
# print(open_files_dict)
# for fn_open in open_files_dict.values():
#     print(fn_open, fn_open.closed)


# open index-hopped and uknown/low quality files to write out to
hopped_read1 = open('index-hopped_' + input_info + '_read1.fq', 'w')
hopped_read2 = open('index-hopped_' + input_info + '_read2.fq', 'w')
unkn_lowqual_read1 = open('unkn_lowqual_' + input_info + '_read1.fq', 'w')
unkn_lowqual_read2 = open('unkn_lowqual_' + input_info + '_read2.fq', 'w')

# open R1-4 files for reading in
R1file = gzip.open(read1f, 'rt')
R2file = gzip.open(read2f, 'rt')
R3file = gzip.open(read3f, 'rt')
R4file = gzip.open(read4f, 'rt')

# initialize counters and records
ln = 0
rn = 0
record_R1 = []
record_R2 = []
record_R3 = []
record_R4 = []

properly_matched_count = 0
hopped_general_count = 0
unkn_lowq_count = 0

# use zip to simultaneously read line by line
for line in zip(R1file,R2file,R3file,R4file):
    ln+=1
    line_r1 = line[0].strip()
    line_r2 = line[1].strip()
    line_r3 = line[2].strip()
    line_r4 = line[3].strip()
    if ln // 4 == rn:               # floor division to append first 3 lines of record
        record_R1.append(line_r1)
        record_R2.append(line_r2)
        record_R3.append(line_r3)
        record_R4.append(line_r4)
    else:                           # append last line of record
        record_R1.append(line_r1)
        record_R2.append(line_r2)
        record_R3.append(line_r3)
        record_R4.append(line_r4)
        rn+=1
        # print(record_R1, record_R2, record_R3, record_R4, '\n', sep='\n')

        # store index1 sequence and reverse complemented index2 sequence
        index1_seq = record_R2[1]
        revcompl_index2_seq = get_revcomp(record_R3[1])
        # print(index1_seq, revcompl_index2_seq, '\n', sep='\n')

        # concatenate index1 seq and rev complemented index2 seq to header lines of read1 and read2 (R4) records
        record_R1[0] = record_R1[0] + ";" + index1_seq + "-" + revcompl_index2_seq
        record_R4[0] = record_R4[0] + ";" + index1_seq + "-" + revcompl_index2_seq
        # print(record_R1, record_R2, record_R3, record_R4, '\n', sep='\n')

        if index1_seq and revcompl_index2_seq in indexes_dict:      # both need to match the given indexes to proceed
            good1 = check_qscore(index1_seq)
            good2 = check_qscore(revcompl_index2_seq)
            if good1 and good2:
                # if both indexes are known, good qual, and match each other, write to the appropriate matched files and incrememnt specific dictionary counter and general counter
                if index1_seq == revcompl_index2_seq:
                    # recreate keys to look up and write to appropriate file in the abbr
                    index_key1 = index1_seq + 'read1'
                    index_key2 = index1_seq + 'read2'
                    for l_r1 in record_R1:
                        open_files_dict[index_key1].write(l_r1 + '\n')
                    for l_r4 in record_R4:
                        open_files_dict[index_key2].write(l_r4 + '\n')
                    # recreate tuple key to look up and increment correct counter
                    search_key = (index1_seq, index1_seq)
                    permut_count_dict[search_key]+=1
                    properly_matched_count+=1

                # if both index seqs are known and good qual but are not the same, write to index hopped files and increment specific dictionary counter and general counter
                else:
                    for l_r1 in record_R1:
                        hopped_read1.write(l_r1 + '\n')
                    for l_r4 in record_R4:
                        hopped_read2.write(l_r4 + '\n')
                    # recreate tuple key to look up and increment correct counter
                    search_key = (index1_seq,revcompl_index2_seq)
                    permut_count_dict[search_key]+=1
                    hopped_general_count+=1

            # if a single base pair in either index read > 30, write to uknown/lowqual files and increment counter
            else:
                for l_r1 in record_R1:
                    unkn_lowqual_read1.write(l_r1 + '\n')
                for l_r4 in record_R4:
                    unkn_lowqual_read2.write(l_r4 + '\n')
                unkn_lowq_count+=1

        # if either index seq unknown, write to unknown/lowqual files and increment counter
        else:
            for l_r1 in record_R1:
                unkn_lowqual_read1.write(l_r1 + '\n')
            for l_r4 in record_R4:
                unkn_lowqual_read2.write(l_r4 + '\n')
            unkn_lowq_count+=1

        # clear records
        record_R1 = []
        record_R2 = []
        record_R3 = []
        record_R4 = []



print('Are Files Closed?')

# close matched files
for fn_close in open_files_dict.values():
    fn_close.close()
    print(fn_close, fn_close.closed)

# close all other files
hopped_read1.close()
hopped_read2.close()
unkn_lowqual_read1.close()
unkn_lowqual_read2.close()
R1file.close()
R2file.close()
R3file.close()
R4file.close()

# print whether other files are closed
print(hopped_read1, hopped_read1.closed)
print(hopped_read2, hopped_read2.closed)
print(unkn_lowqual_read1, unkn_lowqual_read1.closed)
print(unkn_lowqual_read2, unkn_lowqual_read2.closed)
print(R1file, R1file.closed)
print(R2file, R2file.closed)
print(R3file, R3file.closed)
print(R4file, R4file.closed)



### calulations and information for reporting ###

# keys are tuple of sample number, group number, treatment; values are percents
# {(sample,group,treatment):percent}
indiv_matched_counts = {}
for combo,counts in permut_count_dict.items():
    if combo[0] == combo [1]:
        sample_num = indexes_dict[combo[0]][1]
        group_num = indexes_dict[combo[0]][2]
        treatment = indexes_dict[combo[0]][3]
        indiv_matched_counts[(sample_num, group_num, treatment)] = ((counts / rn) * 100)
#print(indiv_matched_counts)

# calculate general percentages
matched_prcnt = properly_matched_count / rn * 100
hopped_prcnt = hopped_general_count / rn * 100
unkn_lowqual_prcnt = unkn_lowq_count / rn * 100


# write out summary info
with open('summary_output_' + input_info + '.txt', 'w') as outfile:
    outfile.write('Naming Convention for Matched Indexes Files:' + '\n')
    outfile.write('indexID_sample_group_treatment_IlluminaSampleNumber_IlluminaLaneNumber_read1/2.fq' + '\n' + '\n')
    outfile.write('Total number of records: ' + str(rn) + '\n')
    outfile.write('Matched Pairs: ' + str(properly_matched_count) + ', ' + str(format(matched_prcnt, '.1f')) + '%' + '\n')
    outfile.write('Uknown or Low Quality Pairs: ' + str(unkn_lowq_count) + ', ' + str(format(unkn_lowqual_prcnt, '.1f')) + '%' + '\n'+ '\n')
    outfile.write('Index-hopped Pairs: ' + str(hopped_general_count) + ', ' + str(format(hopped_prcnt, '.1f')) + '%' + '\n' + '\n')
    outfile.write('Percentage of Matched Reads from Each Sample' + '\n')
    outfile.write('sample number' + '\t' + 'group' + '\t' + 'treatment' + '\t' + 'percent' + '\n')
    for sample,s_prcnt in indiv_matched_counts.items():
        outfile.write(str(sample[0]) + '\t' + str(sample[1]) + '\t' + str(sample[2]) + '\t' + str(format(s_prcnt, '.1f')) + '\n')

# write out table of index combo counts
with open('indexcombos_counts_' + input_info + '.tsv', 'w') as table:
    table.write('index pair combination' + '\t' + 'number of pairs' + '\n')
    for combo,count in permut_count_dict.items():
        table.write(str(combo[0]) + ' ' + str(combo[1]) + '\t' + str(count) + '\n')



