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

def main(argv = sys.argv[1:]):
    args, fetch = cmd_parser(argv= argv)

    seq2mv, pos_read = seq_to_mv(reads_ids = args.readID, region = args.region, sample = args.sample,
                    seq = args.seq, mv = args.mv, ts = args.ts, fetch = fetch, pos = args.pos)  
    




def cmd_parser(argv):
    parser = argparse.ArgumentParser(description="plots electric current/timepoint with corresponding basecalled base for given Read_ID")
    parser.add_argument("--sequencer", help= "name of sequencer (p2i/p2s)")
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