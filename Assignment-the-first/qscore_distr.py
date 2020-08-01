#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="A program to get filename of gzipped fastq with untrimmed reads, and the number of records contained, default read length is 101 base pairs")
    parser.add_argument("-f", "--file", help="zipped fastq file", required=True, type=str)
    parser.add_argument("-l", "--readlen", help="untrimmed read length", type=int, default=101)
    parser.add_argument("-p", "--pngname", help="what you want the mean plot png to be named", required=True, type=str)
    return parser.parse_args()

# store args
args = get_args()
qfile = args.file
read_length = args.readlen
plot_name = args.pngname

def init_list(value=0.0):
    '''Takes a value for populating a list of length 101. 
    If no value is passed, initializes list with 101 values of 0.0.'''
    ls = []
    for i in range(read_length):
        ls.append(value)
    return ls


def convert_phred(letter):
    '''Converts a single ASCII character into a numeric qscore. Phred33 encoding only'''
    qscore = ord(letter) - 33
    return qscore


def populate_list(file):
    '''Takes a gzipped fastq file with read lengths of 101 base pairs.
    Returns list of summed phred scores for each position 
    and the total record count of the input file'''
    tot_scores = init_list()
    with gzip.open(file, 'rt') as f:
        ln = 0
        for line in f:
            ln+=1
            line = line.strip()     # strip whitespace
            if ln % 4 == 0:             # isolate qscore lines
                for pos,qual in enumerate(line):
                    phred = convert_phred(qual)    # convert ASCII to qscore
                    tot_scores[pos]+=phred        # running sum of scores for each position
    rn = ln / 4
    return tot_scores, rn

# populate list and get number of records
mean_scores,num_rec = populate_list(qfile)

# print(mean_scores, num_rec, sep='\n')

# divide in place to get mean
for i in range(len(mean_scores)):
    mean_scores[i] = mean_scores[i] / num_rec

# get indices for plotting
indices = []
for i in range(len(mean_scores)):
    indices.append(i)

plt.figure(figsize=[12,6])
plt.bar(indices, mean_scores)
plt.xlabel("position (0-based)")
plt.ylabel("mean quality score")
plt.title("Mean Quality Score at Each Position")

plt.savefig(plot_name)
