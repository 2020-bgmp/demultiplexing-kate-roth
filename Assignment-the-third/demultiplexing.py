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
            compl_seq+=bp
    revcomp = compl_seq[::-1]
    return revcomp
# print(get_revcomp('GCT'))
# print(get_revcomp('NAC'))
# print(get_revcomp('CNG'))

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
permut_count_dict = {}
permutations = itertools.product(indexes_dict, repeat=2)
for p in list(permutations):
    permut_count_dict[p] = 0
# print(permut_count_dict)


# generate filenames for matched filenames
matched_filenames = []
for index_info in indexes_dict.values():
    # {index sequence:(indexID, sample, group, treatment)}
    fname1 = index_info[0] + "_" + index_info[1] + "_" + index_info[2] + "_" + index_info[3] + "_read1.fq"
    fname2 = index_info[0] + "_" + index_info[1] + "_" + index_info[2] + "_" + index_info[3] + "_read2.fq"
    matched_filenames.append(fname1)
    matched_filenames.append(fname2)
# print(matched_filenames)


# make filename abbreviations for opening matched files
matched_seq_keys = []
matched_filen_abbrev = []
for ind_seq,ind_info in indexes_dict.items():
    seqkey1 = ind_seq + 'read1'
    seqkey2 = ind_seq + 'read2'
    abbr1 = ind_info[0] + 'read1'
    abbr2 = ind_info[0] + 'read2'
    matched_seq_keys.append(seqkey1)
    matched_seq_keys.append(seqkey2)
    matched_filen_abbrev.append(abbr1)
    matched_filen_abbrev.append(abbr2)
# print(matched_seq_keys)
# print(matched_filen_abbrev)


# print(len(matched_filenames)==len(matched_filen_abbrev)==len(matched_seq_keys))

# open matched files to write out to
abbr_filenames_dict = {}
for i in range(len(matched_filenames)):
    abbr_filenames_dict[matched_seq_keys[i]] = open(matched_filenames[i], 'w')
# for fn_open in matched_filenames:
#     matched_filen_abbrev[abbr_count] = open(fn_open, 'w')
#     # print(matched_filen_abbrev[abbr_count], fn_open, matched_filen_abbrev[abbr_count].closed)
#     abbr_count+=1

# print(abbr_filenames_dict)
# for fn_open in abbr_filenames_dict.values():
#     print(fn_open, fn_open.closed)


# open index-hopped and uknown/low quality files to write out to
hopped_read1 = open('index-hopped_read1.fq', 'w')
hopped_read2 = open('index-hopped_read2.fq', 'w')
unkn_lowqual_read1 = open('unkn_lowqual_read1.fq', 'w')
unkn_lowqual_read2 = open('unkn_lowqual_read2.fq', 'w')

# open R1-4 files for reading in
R1file = gzip.open(read1f, 'rt')
R2file = gzip.open(read2f, 'rt')
R3file = gzip.open(read3f, 'rt')
R4file = gzip.open(read4f, 'rt')

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
    if ln // 4 == rn:
        record_R1.append(line_r1)
        record_R2.append(line_r2)
        record_R3.append(line_r3)
        record_R4.append(line_r4)
    else:
        record_R1.append(line_r1)
        record_R2.append(line_r2)
        record_R3.append(line_r3)
        record_R4.append(line_r4)
        rn+=1
        # print(record_R1, record_R2, record_R3, record_R4, '\n', sep='\n')

        index1_seq = record_R2[1]
        revcompl_index2_seq = get_revcomp(record_R3[1])
        # print(index1_seq, revcompl_index2_seq, '\n', sep='\n')

        record_R1[0] = record_R1[0] + ";" + index1_seq + "-" + revcompl_index2_seq
        record_R4[0] = record_R4[0] + ";" + index1_seq + "-" + revcompl_index2_seq
        # print(record_R1, record_R2, record_R3, record_R4, '\n', sep='\n')

        if index1_seq and revcompl_index2_seq in indexes_dict:
            good1 = check_qscore(index1_seq)
            good2 = check_qscore(revcompl_index2_seq)
            if good1 and good2:
                if index1_seq == revcompl_index2_seq:
                    index_key1 = index1_seq + 'read1'
                    index_key2 = index1_seq + 'read2'
                    for l_r1 in record_R1:
                        abbr_filenames_dict[index_key1].write(l_r1 + '\n')
                    for l_r4 in record_R4:
                        abbr_filenames_dict[index_key2].write(l_r4 + '\n')
                    search_key = (index1_seq, index1_seq)
                    permut_count_dict[search_key]+=1
                    properly_matched_count+=1

                else:
                    for l_r1 in record_R1:
                        hopped_read1.write(l_r1 + '\n')
                    for l_r4 in record_R4:
                        hopped_read2.write(l_r4 + '\n')
                    search_key = (index1_seq,revcompl_index2_seq)
                    permut_count_dict[search_key]+=1
                    hopped_general_count+=1

            else:
                for l_r1 in record_R1:
                    unkn_lowqual_read1.write(l_r1 + '\n')
                for l_r4 in record_R4:
                    unkn_lowqual_read2.write(l_r4 + '\n')
                unkn_lowq_count+=1


        else:
            for l_r1 in record_R1:
                unkn_lowqual_read1.write(l_r1 + '\n')
            for l_r4 in record_R4:
                unkn_lowqual_read2.write(l_r4 + '\n')
            unkn_lowq_count+=1

        record_R1 = []
        record_R2 = []
        record_R3 = []
        record_R4 = []


# close matched files
for fn_close in abbr_filenames_dict.values():
    fn_close.close()
    print(fn_close, fn_close.closed)

# close all other files
hopped_read1.close()
hopped_read2.close()
unkn_lowqual_read1.close()
unkn_lowqual_read1.close()


matched_prcnt = properly_matched_count / rn * 100
hopped_prcnt = hopped_general_count / rn * 100
unkn_lowqual_prcnt = unkn_lowq_count / rn * 100


print(permut_count_dict)
indiv_matched_counts = {}
for combo,counts in permut_count_dict.items():
    print(combo[0], combo[1])
    if combo[0] == combo [1]:
        sample_num = indexes_dict[combo[0]][1]
        indiv_matched_counts[sample_num] = counts

print(indiv_matched_counts)


# write out summary info
with open('summary_output.txt', 'w') as outfile:
    outfile.write('Total number of records: ' + str(rn) + '\n')
    outfile.write('Number of properly matched pairs: ' + str(properly_matched_count) + ', ' + str(format(matched_prcnt, '.1f')) + '%' + '\n')
    outfile.write('Number of index-hopped pairs: ' + str(hopped_general_count) + ', ' + str(format(hopped_prcnt, '.1f')) + '%' + '\n')
    outfile.write('Number of uknown or low quality pairs: ' + str(unkn_lowq_count) + ', ' + str(format(unkn_lowqual_prcnt, '.1f')) + '%' + '\n')

# print('index pair combination', 'number of pairs', sep='\t')
# for combo,counts in permut_count_dict.items():
#     print(combo[0],combo[1])
#     print(combo, counts, sep='\t')


