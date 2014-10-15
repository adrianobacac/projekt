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
from contextlib import nested

import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np
from scipy import weave

BUFFER_SIZE = 100

# @profile
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
  

# @profile
def size_human_readable(num):
    measure = ["B", "KB", "MB", "GB", "TB", "PB"]
    index = 0
    while(num > 1024 and index < len(measure)):
        num /= 1024.
        index += 1
    return str(round(num, 2)) + measure[index]
  
# @profile
def entropy(seq):
    """
    whole function rewriten in C
    """
    code = r'''
        int counter[26]={0};
        int len = seq.length();
        
        for(int i = 0; i < len; ++i) {
            counter[seq[i]-'A']=int(counter[seq[i]-'A'])+1;
        }
        
        float percentage;
        float entropy = 0;
        
        for(int i = 0; i<26; ++i){
            if(counter[i]>0){
                percentage = (float)counter[i]/len;
                entropy -= percentage * log2(percentage);
            }
        }
        return_val = entropy;
        '''
    return  weave.inline(code, ['seq'])
    

# @profile
def max_entropy(seq):
    return math.log(len(seq), 2)

# @profile
def create_statistics_file(total_count, complex_count , input_file, output_folder):
    
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

@profile
def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="A FASTA file to be filtered for low complexity sequences")
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
    entropies = []
     
    entropy_threshold = 1.5
    
    with nested(open(args.input_file, "rU"), open(args.output_folder + "/filtered_" + args.input_file, "w")) as (fin, fout):
        buffer_size = BUFFER_SIZE
        buffered_records = []
        for record in SeqIO.parse(fin, file_type):
            if buffer_size > 0:
                buffer_size-=1
                total_count += 1
                ent = entropy(str(record.seq))
                entropies.append(ent)
                if ent >= entropy_threshold:
                    complex_count += 1
                    buffered_records.append(record)
            else:
                SeqIO.write(buffered_records, fout, file_type)
                buffered_records = []
                buffer_size = BUFFER_SIZE
    
    make_histogram(entropies,
                    "sequence",
                    "entropy",
                     "sequence entropy",
                     args.output_folder,
                     "Entropy per sequence"
                     )
    
    create_statistics_file(total_count, complex_count, args.input_file, args.output_folder)
    
if(__name__ == "__main__"):
    main()
