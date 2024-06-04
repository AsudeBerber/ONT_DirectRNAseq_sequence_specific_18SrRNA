__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

import argparse
import os
import pod5 as p5
import matplotlib.pyplot as plt
import numpy as np
import pysam as ps


sample = "../../resources/alignments/p2s_aligned_sorted.bam"
samfile = ps.AlignmentFile(f"{sample}")
read_ids = "62b6f170-c678-465c-aafe-b31af1e94f19"


for read in samfile.fetch():
        i = 0
        if read.query_name == read_ids:
            read_ID = read.query_name
            seq = read.query_alignment_sequence
            seq = seq[::-1] #sequence order is 5' -> 3', mv and signal are 3' -> 5'

            # Workaround in cases where two ts tags per read exists:
            # read.set_tag("ts", None) #first ts tag is transcript strand, has to be removed

            mv = read.get_tag("mv")
            ts = read.get_tag("ts")
            #print(read.get_tags())
            #print(f"ts:{ts}; {read_ID}") 
            aligned_pairs = read.get_aligned_pairs()
            for i, pair in enumerate(aligned_pairs):
                print(i, pair)
        else: 
            #removing the else part makes the code only 1s faster
            i = i+1
            k = i/500000
            if k.is_integer():
                # print("currently at " +read.query_name + "\n")
                print(f"reads compared: {i} of max. {max_reads}")
            continue