import argparse
import os
import pod5 as p5
import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np
import pysam as ps
import sys
try:
    import align_signal
except ImportError: 
    raise ImportError("Import of module align_signal failed, is align_signal.py in the same folder as this script?")

#sample = "resources/alignments/c05e1233-6e1a-4698-bb74-9b47df9507f2.bam"
#gets movetable (mv), ts and corresponding base sequence (seq) for given read id
def bam_aligned(sample, read_ids, region, pos):

    samfile = ps.AlignmentFile(f"{sample}")

    max_reads = samfile.mapped
    i = 0

    # if index file present, fetch(region = region)
    for read in samfile.fetch(region = region):
        if read.query_name == read_ids:
            read_ID = read.query_name
            seq = read.query_sequence
            

            # Workaround in cases where two ts tags per read exists:
            # read.set_tag("ts", None) #first ts tag is transcript strand(+|-), has to be removed

            mv = read.get_tag("mv")
            ts = read.get_tag("ts")

            ref_seq = read.get_reference_sequence()

            aln_pairs = read.get_aligned_pairs(with_seq = True)
            
            # creates pairs of base positions (query, reference) -> looks up position in alignment sequence for corresponding reference base position
            for pair in aln_pairs:
                if pair [1] == pos: 
                    pos_read = pair[0]
                
            # print(read.get_tags())
            # print(f"ts:{ts}; {read_ID}") 
            return(seq, mv, ts, aln_pairs, ref_seq, read)
        else: 
            #removing the else part makes the code only 1s faster
            i = i+1
            k = i/500000
            if k.is_integer():
                # print("currently at " +read.query_name + "\n")
                print(f"reads compared: {i} of max. {max_reads}")
            continue


def main(argv = sys.argv[1:]):
    args, fetch = cmd_parser(argv= argv)

    seq2mv, pos_read = seq_to_mv(reads_ids = args.readID, region = args.region, sample = args.sample,
                    seq = args.seq, mv = args.mv, ts = args.ts, fetch = fetch, pos = args.pos, range = args.range)  
    
def seq_to_mv(reads_ids, region, sample, seq=None, mv=None, ts=0, fetch = True, pos=1841, range = 8):
    if fetch == True:
        seq, mv, ts, aln_pairs, ref_seq, read = bam_aligned(sample, reads_ids, region, pos)
    if fetch == False:
        seq, mv, ts = args.seq, args.mv, args.ts
    seq2mv, rev_loci = align_signal.access_mv(signal = None, moves = mv, offset = ts,
                                                               loci = [int(pos)], motif_length = 1, extra_window = range,
                                                                 read = read, mode = "single_read")
    breakpoint()
    pass


def cmd_parser(argv):
    parser = argparse.ArgumentParser(description="plots electric current/timepoint with corresponding basecalled base for given Read_ID")
    parser.add_argument("--sequencer", help= "name of sequencer (p2i|p2s)")
    parser.add_argument("--sample", help= "Name of sample bam file w/o .bam ending", action="store")
    parser.add_argument("--readID", help= "Sample ID in bam and pod5 file", action="store")
    parser.add_argument("--pos", help= "position/index of base in middle, 1-based", type=int, action="store")
    parser.add_argument("--range", help= "bases to display the left/right of middle base", type=int, action="store")
    parser.add_argument("--no_fetch", help= "Disables fetching mv, ts, seq from readID - requieres --mv, --ts, --seq", action="store_true")
    parser.add_argument("--seq", help= "complete sequence of read", action="store_true")
    parser.add_argument("--mv", help= "movetable", action="store_true")
    parser.add_argument("--ts", help= "", action="store_true")
    parser.add_argument("--get_readids", help= "store name of all read ids in .txt", action= "store_true")
    parser.add_argument("--region", type= str, action= "store")
    parser.add_argument("--pod5_dir", type= str, action= "store")
    args = parser.parse_args()

    # starts at pos 0
    args.pos = args.pos - 1

    if args.get_readids == True:
        read_id_list_bam(args.sequencer, args.sample)

    if args.no_fetch == True:
        fetch = False
    else:
        fetch = True 

    return args, fetch


if __name__ == "__main__":
    exit(main())