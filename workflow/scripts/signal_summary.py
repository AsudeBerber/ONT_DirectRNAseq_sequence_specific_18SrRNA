__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

"""
calculates features for all given reads around positions given below (could be made into argument) and outputs it to .npz file
signal alignment and feature calculation are done by called module align_signal.py (has to be in same directory)
"""
# Code mostly written by Christoph Engelhardt

############### Run time for 4x10^6 reads: ~5 h; memory usage piles up over time up to ~50%/ 30 GB RAM in the end

import numpy as np
import pod5 as p5
import pysam
from tqdm import tqdm
import argparse
import sys
import os
import json
from pathlib import Path

try:
    import align_signal
except ImportError: 
    raise ImportError("Import of module align_signal failed, is align_signal.py in the same folder as this script?")


# pod5_file = "resources/pod5/p2s/"
# bam_file = f"resources/alignments/p2s_aligned_sorted.bam"
# window_size = 21
# npz_file = f"resources/results/p2s/CCG_window_{window_size}_{Path(bam_file).stem}.npz"

# different positions can be set here, ##### index is 0-based  ######
ref_ac1 = 1336
ref_ac2 = 1841
ref_no_ac = 429 #unacetylated CCG with similar sequence as 1336/1841 (TTCCG)
ref_pos = [ref_ac1] + [ref_ac2] + [ref_no_ac]
motif_length = 1


# handle sys args
def parse_args(argv):
    """Read arguments from command line."""

    parser = argparse.ArgumentParser()

    parser.add_argument("-j", "--json", type=str, required=True)
    parser.add_argument("-b", "--bam", type=str, required=True)
    parser.add_argument("-o", "--output", type=str, required=True)
    parser.add_argument("-w", "--window", type=int, default=21)

    args = parser.parse_args()

    return args

# for getting position in signal (goes 3' to 5')

def main(argv=sys.argv[1:]):

    args = parse_args(argv=argv)

    json_file = args.json
    bam_file = args.bam
    window_size = args.window
    npz_file = args.output

    extra_window = int((window_size - 1) / 2)

    with pysam.AlignmentFile(bam_file, mode = "rb", check_sq=False) as bam: 


        features, qual, query_seq, ref_seq, id = [], [], [], [], []
        with open(json_file, "r") as f:
            pod5_index= json.load(f)

        count_keyErr = 0 # counts skipped reads (s. below)

        for read in tqdm(bam):
            if read.is_unmapped:
                continue
            
            # get loci on the reference matching the motif
            aligned_pairs = read.get_aligned_pairs(with_seq=True, matches_only = False)
            ac_ccg= np.array(list(filter(lambda x: x[1] in ref_pos, aligned_pairs)), dtype= "object")
            # pairs_dict = dict((y, x) for x, y, z in ac_ccg if y is not None)
            loci, ref_loci, rev_loci = align_signal.get_loci(read, ac_ccg, extra_window, motif_length, ref_pos)
        
            if len(loci) == 0:
                continue
            
            # extract features from bam file
            try:
                per_site_qual = np.array([list(read.qual[locus-extra_window: locus+motif_length+extra_window]) for locus in loci])
                per_site_query_seq = np.array([list(read.query_sequence[locus-extra_window: locus+motif_length+extra_window]) for locus in loci])
                seq_dict = dict((x, z) for x, y, z in aligned_pairs)
                per_site_ref_seq = np.array([[seq_dict[key] for key in range(locus-extra_window, locus+motif_length+extra_window)] for locus in loci])
            except:
                raise IndexError("Something went wrong during feature extraction (check given range and indexing of that)")
        
            # dorado sometimes splits reads (https://github.com/nanoporetech/dorado/blob/release-v0.6/documentation/SAM.md#split-read-tags),
            # it doesn't seem possible to align these back to one original read id --> these reads are ignored
            try:
                pod5_file = pod5_index[read.query_name]
            except KeyError:
                count_keyErr += 1
                continue
                

            with p5.Reader(pod5_file) as pod5:
                # Read the selected read from the pod5 file
                # next() is required here as Reader.reads() returns a Generator
                try:
                    
                    # lookup is very fast (1E-4 s)
                    pod5_record = next(pod5.reads(selection=[read.query_name])) 

                    # events: 0.01 s = 100/s
                    #events is inverted as the signal goes from 3' -> 5', but sequence from 5' -> 3'
                    dict_events = align_signal.access_mv(pod5_record.signal, read.get_tag("mv"), read.get_tag("ts"), rev_loci, motif_length, extra_window, read=read, mode="signal_stats")
                    
                    # locus_rev is corresponding pos in signal, as this goes from 3' to 5'
                    
                    per_site_features = np.array([[dict_events[key] for key in reversed(range(locus - extra_window , locus + extra_window + motif_length))] for locus in rev_loci])
                    per_site_id = np.array([read.query_name + ':' + str(locus+1) for locus in ref_loci])
                    

                    #appending is very fast (1E-6 s)
                    features.append(per_site_features)
                    qual.append(per_site_qual)
                    query_seq.append(per_site_query_seq)
                    ref_seq.append(per_site_ref_seq)
                    id.append(per_site_id)
                except:
                    raise Exception("something went wrong while accessing pod5 data: \n 1) is the right pod5 path given?\n 2) does python try to append empty position?")
                    

    
    features = np.vstack(features)
    qual = np.vstack(qual)
    query_seq = np.vstack(query_seq)
    ref_seq = np.vstack(ref_seq)
    id = np.hstack(id)
    

    # checks if results folder exists, creates otherwise
    # check if plot dir exists, creates it otherwise
    npz_dir = os.path.dirname(npz_file) 
    if not os.path.isdir(npz_dir): os.makedirs(npz_file)
    
    # None's in query seq convert array to object type, which is not liked by np.load (loading of npz file), Nones are converted to string therfore
    ref_seq = np.array(ref_seq, dtype="U")
    np.savez_compressed(npz_file, feat = features, qual = qual, query = query_seq, ref = ref_seq, id = id)

    print(f"{count_keyErr} reads in bam but not in pod5 (this is normal, don't panic)")
    # it should be somehow possible to convert this to a txt file, however the dimensions of the array have to be reduced for this, maybe via for loop?
    # does not work so far
    # np.savetxt('resources/results/p2s/summary.txt', features, qual, query_seq, ref_seq, id))
    
    return 0


if __name__ == "__main__":
    exit(main())
