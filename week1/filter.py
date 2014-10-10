#! /usr/bin/env python

'''
Created on Oct 8, 2014

@author: ppx10
'''
import sys
import os
import argparse
from Bio import SeqIO
from collections import Counter
import math
import matplotlib.pyplot as plt
import time

def make_histogram(y, row_title, col_title, title, output_folder, file_name):
    plt.figure()
    
    plt.xlabel(row_title)
    plt.ylabel(col_title)
    plt.title(title)
  
    x = [i for i in range(1, len(y) + 1)]
    plt.bar(x, y, width=1, align="center", color="b")
    plt.xticks(range(1, len(y) + 1), x)
    plt.axis([0.5, len(x) + 0.5, 0, max(y)*1.1 ])
    plt.grid(True)
    
    save_path = os.path.join(output_folder, file_name)
    plt.savefig(save_path)
  
def get_size(path):
    return (os.stat(path)).st_size

def size_human_readable(num):
    measure = ["B", "KB", "MB", "GB", "TB", "PB"]
    index = 0
    while(num > 1024 and index < len(measure)):
        num /= 1024.
        index += 1
    return str(round(num, 2)) + measure[index]
  
def entropy(seq):
    nucleotid_occurances = Counter(seq)
    seq_len=len(seq)
    enth=0
    for occurance in nucleotid_occurances.values():
        freq = float(occurance)/seq_len
        enth-= freq * math.log(freq)
    return enth

def max_entropy(seq):
    return math.log(len(seq))

def create_output(sequences, entropies, total_count, complex_count ,input_file, output_folder):
    make_histogram(entropies,
                    "sequence",
                    "entropy",
                     "sequence entropy",
                     output_folder,
                     "Entropy per sequence"
                     )
    
    output_handle = open(output_folder + "/filtered_" + input_file, "w")
    SeqIO.write(sequences, output_handle, "fasta")
    output_handle.close()
    
    size = size_human_readable(get_size(input_file))
    creation_time=time.asctime()
    
    with open(output_folder + "/filter_statistics.txt", "w") as txt:
        text_lines = ["File: "+input_file,
                      "\n", 
                      "Creation date: "+creation_time,
                      "\n", 
                      "Size: " + size, 
                      "\n",
                      "Reads: " + str(total_count), 
                      "\n",
                      "Kept reads:"+str(complex_count), 
                      "\n",
                      "Discarded reads:" + str(total_count-complex_count)
                      ]
        txt.writelines(text_lines)

def main(argv):
    
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="A FATSA file to be filtered for low complexity sequences")
    parser.add_argument("output_folder", help="Folder where the output will be stored")
    
    if len(argv) != 2:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)
    
    handle = open(args.input_file, "rU")
    
    total_count = 0 
    complex_count = 0
    complex_records = []
    entropies=[]
    
    entropy_threshold=1.1
    for record in SeqIO.parse(handle, "fasta"):
        
        total_count+=1
        ent=entropy(record.seq)
        entropies.append(ent)
        if ent >= entropy_threshold:
            complex_count += 1
            complex_records.append(record)
            
    create_output(complex_records, entropies, total_count, complex_count, args.input_file, args.output_folder)
    
if(__name__ == "__main__"):
    main(sys.argv[1:])
