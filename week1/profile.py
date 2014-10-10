#! /usr/bin/env python

'''
Created on Oct 8, 2014

@author: ppx10
'''

from Bio import SeqIO
from os.path import os
import matplotlib.pyplot as plt
from numpy import average
import argparse

import sys

def get_size(path):
    return (os.stat(path)).st_size

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
    
def make_histogram(y, row_title, col_title, title, output_folder, file_name="seq_len_dist.png"):
    plt.figure()
    
    plt.xlabel(row_title)
    plt.ylabel(col_title)
    plt.title(title)
  
    x = [i for i in range(1, len(y) + 1)]
    plt.bar(x, y, width=1, align="center", color="g")
    plt.xticks(range(1, len(y) + 1), x)
    plt.axis([0.5, len(x) + 0.5, 0, max(y) + 10])
    plt.grid(True)
    
    save_path = os.path.join(output_folder, file_name)
    plt.savefig(save_path)
    
def create_output(lengths, qualities, input_file, output_folder):
    make_histogram(lengths,
                    "sequences",
                    "length",
                     "sequence length distribution",
                     output_folder
                     )
    make_boxplot(qualities,
                 "position in read",
                 "average nucleotide quality",
                 "average nucleotide quality per position in the read",
                 output_folder
                 )
    

    size = size_human_readable(get_size("shortExample.fastq"))
    average_length = average(lengths)
    # average_read_quality 
    read_count = len(lengths)

def main(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="A fastq file to be parsed")
    parser.add_argument("output_folder", help="Folder where the output will be stored")
    
    if len(argv) != 2:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)
    
    handle = open(args.input_file, "rU")
    lengths = []
    qualities = []
    for record in SeqIO.parse(handle, "fastq"):
        lengths.append(len(record.seq))
        qualities.append(record.letter_annotations["phred_quality"])
    """
    
    qualities.sort(key = lambda row: len(row))
    
    
    qualities_per_position = []
    for i in range(0, max(lengths)):
        qualities_per_position.append([])
        for quality in reversed(qualities):
            if i < len(quality):
                qualities_per_position[i].append(quality[i])
        
    """
    
    create_output(lengths, qualities, args.input_file, args.output_folder)
        
if(__name__ == "__main__"):
    main(sys.argv[1:])
