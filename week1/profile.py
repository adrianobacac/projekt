#! /usr/bin/env python

'''
Created on Oct 8, 2014

@author: ppx10
'''

from os.path import os
import argparse
import sys

import matplotlib.pyplot as plt
from numpy import average
from Bio import SeqIO
from collections import OrderedDict

def size_human_readable(num):
    measure = ["B", "KB", "MB", "GB", "TB", "PB"]
    index = 0
    while(num > 1024 and index < len(measure)):
        num /= 1024.
        index += 1
    return str(round(num, 2)) + measure[index]

def make_boxplot(average_qualities, row_title, col_title, title, output_folder, file_name="avg_nuc_quality.png"):
    plt.figure()
    
    plt.xlabel(row_title)
    plt.ylabel(col_title)
    plt.title(title)
    
    plt.boxplot(average_qualities)
    
    
    save_path = os.path.join(output_folder, file_name)
    plt.savefig(save_path)
    
def make_histogram(lengths, row_title, col_title, title, output_folder, file_name="seq_len_dist.png"):
    
    ordered_lengths=OrderedDict(sorted(lengths.items(), key=lambda lengths: lengths[1]))
    
    y = ordered_lengths.keys()
    x = ordered_lengths.values()
    
    plt.figure()
    
    plt.xlabel(row_title)
    plt.ylabel(col_title)
    plt.title(title)
    
    # x = [i for i in range(1, len(y) + 1)]
    plt.bar(range(len(y)), y, width=1, align="center", color="g")
    plt.xticks(range(len(y)), x)
    # plt.axis([0.5, len(x) + 0.5, 0, max(y)*1.05])
    plt.grid(True)
    
    save_path = os.path.join(output_folder, file_name)
    plt.savefig(save_path)
    
def make_statistics_file(lengths, qualities, input_file, output_folder):
    
    size = size_human_readable((os.stat(input_file)).st_size)
    
    total_length = 0
    sequence_count = 0
    for key in lengths:
        total_length += int(key) * lengths[key]
        sequence_count += lengths[key]
    
    # average_length = average(lengths)
    average_length = round(float(total_length) / sequence_count, 2)
    average_read_quality = round(average([average(quality) for quality in qualities])) 
    
    
    
    with open(output_folder + "/profile_statistics.txt", "w") as txt:
        text_lines = ["File: " + input_file,
                      "\n",
                      "Size: " + size,
                      "\n",
                      "Average length: " + str(average_length),
                      "\n",
                      "Average read quality: " + str(average_read_quality),
                      "\n",
                      "Number of reads: " + str(sequence_count)
                      ]
        txt.writelines(text_lines)
        

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="A fastq file to be parsed")
    parser.add_argument("output_folder", help="Folder where the output will be stored")
    
        
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)
    
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)
    
    lengths = {}
    qualities = []
    # counter=0
    with open(args.input_file, "rU") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            # counter+=1
            try:
                lengths[len(record.seq)] += 1
            except:
                lengths.update({len(record.seq):1})
            # lengths.append(len(record.seq))
            qualities.append(record.letter_annotations["phred_quality"])

    qualities.sort(key=lambda row: len(row))
    
    qualities_per_position = []
    for i in range(0, max(lengths.keys())):
        qualities_per_position.append([])
        for quality in reversed(qualities):
            if i < len(quality):
                qualities_per_position[i].append(quality[i])
            else:
                break;
    make_histogram(lengths,
                    "sequences",
                    "length",
                     "sequence length distribution",
                     args.output_folder
                     )
    make_boxplot(qualities_per_position,
                 "position in read",
                 "average nucleotide quality",
                 "average nucleotide quality per position in the read",
                 args.output_folder
                 )
    
    make_statistics_file(lengths, qualities_per_position, args.input_file, args.output_folder)
        
if(__name__ == "__main__"):
    main()
