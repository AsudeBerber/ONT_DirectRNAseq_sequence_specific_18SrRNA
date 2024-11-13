__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

"""
calculates features for all given reads around positions given below (could be made into argument) and outputs it to .npz file
signal alignment and feature calculation are done by called module align_signal.py (has to be in same directory)
"""
# Code mostly written by Christoph Engelhardt

import numpy as np
import pod5 as p5
import pysam
import argparse
import sys
import os
import json
from pathlib import Path

try:
    import align_signal
except ImportError: 
    raise ImportError("Import of module align_signal failed, is align_signal.py in the same folder as this script?")

# Define constants
ref_ac1 = 1336
ref_ac2 = 1841
ref_no_ac = 429  # Unacetylated CCG position
ref_pos = [ref_ac1, ref_ac2, ref_no_ac]
motif_length = 1

# Parse command-line arguments
def parse_args(argv):
    """Read arguments from command line."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--json", type=str, required=True)
    parser.add_argument("-b", "--bam", type=str, required=True)
    parser.add_argument("-o", "--output", type=str, required=True)
    parser.add_argument("-w", "--window", type=int, default=21)
    return parser.parse_args(argv)

# Main function
def main(argv=sys.argv[1:]):
    args = parse_args(argv=argv)
    json_file = args.json
    bam_file = args.bam
    window_size = args.window
    npz_file = args.output
    extra_window = int((window_size - 1) / 2)

    print(f"Opening BAM file: {bam_file}")
    try:
        with pysam.AlignmentFile(bam_file, mode="rb", check_sq=False) as bam:
            print("BAM file opened successfully.")
            
            features, qual, query_seq, ref_seq, id = [], [], [], [], []
            
            print(f"Opening JSON file: {json_file}")
            with open(json_file, "r") as f:
                pod5_index = json.load(f)
            print("JSON file opened and loaded successfully.")
            
            count_keyErr = 0  # Counts skipped reads
            
            print("Starting to process reads...")
            for read in bam:
                print(f"Processing read {read.query_name}")
                if read.is_unmapped:
                    print(f"Skipped unmapped read {read.query_name}")
                    continue

                # Get loci on the reference matching the motif
                aligned_pairs = read.get_aligned_pairs(with_seq=True, matches_only=False)
                ac_ccg = np.array(list(filter(lambda x: x[1] in ref_pos, aligned_pairs)), dtype="object")
                loci, ref_loci, rev_loci = align_signal.get_loci(read, ac_ccg, extra_window, motif_length, ref_pos)
            
                if len(loci) == 0:
                    print(f"Skipped read {read.query_name}: No loci found.")
                    continue

                # Extract features from BAM file
                try:
                    per_site_qual = np.array([list(read.qual[locus-extra_window: locus+motif_length+extra_window]) for locus in loci])
                    per_site_query_seq = np.array([list(read.query_sequence[locus-extra_window: locus+motif_length+extra_window]) for locus in loci])
                    seq_dict = dict((x, z) for x, y, z in aligned_pairs)
                    per_site_ref_seq = np.array([[seq_dict[key] for key in range(locus-extra_window, locus+motif_length+extra_window)] for locus in loci])
                except IndexError:
                    print(f"Error during feature extraction for read {read.query_name}")
                    continue
                
                # Retrieve pod5 file path and record
                try:
                    pod5_file = pod5_index[read.query_name]
                except KeyError:
                    count_keyErr += 1
                    print(f"Skipped read {read.query_name}: KeyError in pod5 index.")
                    continue

                print(f"Processing read {read.query_name}: loci={loci}")

                try:
                    with p5.Reader(pod5_file) as pod5:
                        # Retrieve the read from pod5 file
                        pod5_record = next(pod5.reads(selection=[read.query_name])) 

                        # Get signal features from align_signal module
                        dict_events = align_signal.access_mv(
                            pod5_record.signal,
                            read.get_tag("mv"),
                            read.get_tag("ts"),
                            rev_loci,
                            motif_length,
                            extra_window,
                            read=read,
                            mode="signal_stats"
                        )
                        print(dict_events)
                        # Prepare per-site features and identifiers
                        per_site_features = np.array([[dict_events[key] for key in reversed(range(locus - extra_window, locus + extra_window + motif_length))] for locus in rev_loci])
                        per_site_id = np.array([read.query_name + ':' + str(locus + 1) for locus in ref_loci])

                        # Append results
                        features.append(per_site_features)
                        qual.append(per_site_qual)
                        query_seq.append(per_site_query_seq)
                        ref_seq.append(per_site_ref_seq)
                        id.append(per_site_id)
                
                except Exception as e:
                    print(f"Error accessing pod5 data for read {read.query_name}: {e}")
                    continue
            
            print("Finished processing reads.")
            print(features, id)
            # Convert lists to arrays for saving
            features = np.vstack(features)
            qual = np.vstack(qual)
            query_seq = np.vstack(query_seq)
            ref_seq = np.array(ref_seq, dtype="U")
            id = np.hstack(id)
            
            # Ensure output directory exists
            npz_dir = os.path.dirname(npz_file) 
            if not os.path.isdir(npz_dir):
                os.makedirs(npz_dir)
            
            # Save output in compressed format
            np.savez_compressed(npz_file, feat=features, qual=qual, query=query_seq, ref=ref_seq, id=id)
            print(f"Output saved to {npz_file}")
            print(f"{count_keyErr} reads in BAM but not in pod5 (this is normal, don't panic)")

    except FileNotFoundError as e:
        print(f"File not found error: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    exit(main())
