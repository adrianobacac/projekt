#! /usr/bin/env python

'''
Created on Oct 8, 2014

@author: ppx10
'''
import sys
import os
import argparse
import math
import time
from collections import Counter

import matplotlib.pyplot as plt
from Bio import SeqIO


def make_histogram(y, row_title, col_title, title, output_folder, file_name):
    plt.figure()
    plt.xlabel(row_title)
    plt.ylabel(col_title)
    plt.title(title)
    """
    changed plt.bar to plt.hist
    """
    plt.hist(y, color="b")
    plt.grid(True)
    
    save_path = os.path.join(output_folder, file_name)
    plt.savefig(save_path)
"""
##@profile
def make_histogram(y, row_title, col_title, title, output_folder, file_name):
    plt.figure()
    
    plt.xlabel(row_title)
    plt.ylabel(col_title)
    plt.title(title)
    
    x = [i for i in range(1, len(y) + 1)]
    plt.bar(x, y, width=1, align="center", color="b")
    plt.xticks(range(1, len(y) + 1), x)
    plt.axis([0.5, len(x) + 0.5, 0, max(y) * 1.1 ])
    plt.grid(True)
    
    save_path = os.path.join(output_folder, file_name)
    plt.savefig(save_path)
  """

def size_human_readable(num):
    measure = ["B", "KB", "MB", "GB", "TB", "PB"]
    index = 0
    while(num > 1024 and index < len(measure)):
        num /= 1024.
        index += 1
    return str(round(num, 2)) + measure[index]
  
def entropy(seq):
    nucleotid_occurances = Counter(seq)
    seq_len = len(seq)
    enth = 0
    for occurance in nucleotid_occurances.values():
        freq = float(occurance) / seq_len
        enth -= freq * math.log(freq, 2)
        
    return enth

def max_entropy(seq):
    return math.log(len(seq), 2)

def create_statistics_file(sequences, entropies, total_count, complex_count , input_file, output_folder):
    
    size = size_human_readable((os.stat(input_file)).st_size)
    creation_time = time.asctime()
    
    with open(output_folder + "/filter_statistics.txt", "w") as txt:
        text_lines = ["File: " + input_file,
                      "\n",
                      "Creation date: " + creation_time,
                      "\n",
                      "Size: " + size,
                      "\n",
                      "Reads: " + str(total_count),
                      "\n",
                      "Kept reads: " + str(complex_count),
                      "\n",
                      "Discarded reads: " + str(total_count - complex_count)
                      ]
        txt.writelines(text_lines)

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="A FASTA or FASTQ file to be filtered for low complexity sequences")
    parser.add_argument("output_folder", help="Folder where the output will be stored")
        
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)
    
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)
    
    if(args.input_file.endswith(".fasta")):
        file_type = "fasta"
    elif(args.input_file.endswith(".fastq")):
        file_type = "fastq"
    else:
        parser.error("File format is not supported")
    
    total_count = 0 
    complex_count = 0
    complex_records = []
    entropies = []
    
    entropy_threshold = 1.5
    
    with open(args.input_file, "rU") as handle:
        for record in SeqIO.parse(handle, file_type):
            total_count += 1
            ent = entropy(record.seq)
            entropies.append(ent)
            if ent >= entropy_threshold:
                complex_count += 1
                complex_records.append(record)
    
    make_histogram(entropies,
                   "entropy",
                    "sequence",
                     "sequence entropy",
                     args.output_folder,
                     "Entropy per sequence"
                     )
    
    
    with open(args.output_folder + "/filtered_" + args.input_file, "w") as output_handle:
        SeqIO.write(complex_records, output_handle, file_type)
     
    create_statistics_file(complex_records, entropies, total_count, complex_count, args.input_file, args.output_folder)
    
if(__name__ == "__main__"):
    main()
