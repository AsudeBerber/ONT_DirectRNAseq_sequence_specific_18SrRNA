# Code mostly written by Christoph Engelhardt

############### Run time for 4x10^7 reads: ~5 h; memory usage piles up over time up to ~50%/ 30 GB RAM in the end

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

    parser.add_argument("-p", "--pod5", type=str)
    parser.add_argument("-b", "--bam", type=str)
    parser.add_argument("-o", "--output", type=str)
    parser.add_argument("-w", "--window", type=int, default=21)

    args = parser.parse_args()

    return args

# this part is from https://github.com/WGLab/DeepMod2/blob/main/src/detect.py
def get_events(signal, moves, offset, rev_loci, motif_length, extra_window):
    """
    Normalises and collapses the signal based on the moves table. Outputs an array with the
    following values for each called based:
    4: log10 signal length
    5: mean signal intensity
    6: standard deviation of signal intensity
    7: median signal intensity
    8: median absolute deviation of signal intensity
    0-3: mean signal intensity for each quartile
    """

    # normalise signal
    median = np.median(signal)
    mad = np.median(np.abs(signal-median))
    signal=(signal-median)/mad
    
    stride = moves.pop(0)
    move_index = np.where(moves)[0]
    rlen = len(move_index)
    
    
    dict_events = {}

    # last position is only to get sig_end
    pos_get_signal = [np.arange((locus - extra_window), locus + motif_length+extra_window) for locus in rev_loci]
    pos_get_signal = np.reshape(pos_get_signal, -1)
    
    for locus in pos_get_signal:

        data_tmp= np.zeros((9))
            
        prev = move_index[locus]*stride+offset
        sig_end = move_index[locus+1]*stride+offset
        
        
        sig_len = sig_end-prev
        data_tmp[4]=np.log10(sig_len)
        data_tmp[5]=np.mean(signal[prev:sig_end])
        data_tmp[6]=np.std(signal[prev:sig_end])
        data_tmp[7]=np.median(signal[prev:sig_end])
        data_tmp[8]=np.median(np.abs(signal[prev:sig_end]-data_tmp[4]))
        
        # get the mean signal for each quarter of the base signal
        for j in range(4):
            tmp_cnt=0
            for t in range(j*sig_len//4,min(sig_len, (j+1)*sig_len//4)):
                data_tmp[j]+=signal[t+prev]
                tmp_cnt+=1
            data_tmp[j]=data_tmp[j]/tmp_cnt
        dict_events.update({locus:data_tmp})
    return dict_events

def get_loci(read, pairs, wd, motif_length):
    """
    find positions that match motif
    """    
    ref_loci = []
    ref_loci_index = []
    # reversed loci, for operations from 3' end
    rev_loci = []

    #pairs[0]: query pos; [1]: ref pos; [2] ref base
    for i, pos in enumerate(ref_pos):
        # for empty array
        if pairs.shape == (0,):
            continue
        else:
            if (pos in pairs[:,1]): ref_loci.append(pos)
     
        # in cases where a read neither aligns to any of the ref positions, ref_loci_index would append [], which raises an index error --> solved by checking if array is empty
        index_pos = np.where(pairs[:,1] == pos)[0]
        if index_pos.shape == (0,):
            continue
        else:
            ref_loci_index.append(index_pos[0])
            
 
    loci = [pairs[locus, 0] for locus in ref_loci_index]
    # Remove loci that are not present on the query or too close to the ends of the alignment
    loci = [locus for locus in loci if locus is not None and locus > wd and locus < read.query_length - wd - (motif_length)]
    ref_loci = [pairs[index, 1] for index in ref_loci_index if pairs[index, 0] != None and pairs[index,0] in loci]
    if len(loci) != len(ref_loci):
        raise Exception ("length of reference and query sequence index not matching")
    
    for locus in loci:
        rev_loci.append(rev_locus(locus, read))
    return loci, ref_loci, rev_loci

# for getting position in signal (goes 3' to 5')
def rev_locus(locus, read):
    locus_rev = read.query_length -1 - locus
    return locus_rev


def main(argv=sys.argv[1:]):

    args = parse_args(argv=argv)

    pod5_file = args.pod5
    bam_file = args.bam
    window_size = args.window
    npz_file = args.output

    extra_window = int((window_size - 1) / 2)

    with pysam.AlignmentFile(bam_file, mode = "rb", check_sq=False) as bam: 


        features, qual, query_seq, ref_seq, id = [], [], [], [], []
        with open('resources/results/p2s/pod5.json', "r") as f:
            pod5_index= json.load(f)

        pod5_path = "resources/pod5/p2s/"
        count_keyErr = 0 # counts skipped reads (s. below)

        for read in tqdm(bam):
            if read.is_unmapped:
                continue
            
            # get loci on the reference matching the motif
            aligned_pairs = read.get_aligned_pairs(with_seq=True, matches_only = False)
            ac_ccg= np.array(list(filter(lambda x: x[1] in ref_pos, aligned_pairs)), dtype= "object")
            # pairs_dict = dict((y, x) for x, y, z in ac_ccg if y is not None)
            loci, ref_loci, rev_loci = get_loci(read, ac_ccg, extra_window, motif_length)
        
            if len(loci) == 0:
                continue
            
            # extract features from bam file
            try:
                per_site_qual = np.array([list(read.qual[locus-extra_window: locus+motif_length+extra_window]) for locus in loci])
                per_site_query_seq = np.array([list(read.query_sequence[locus-extra_window: locus+motif_length+extra_window]) for locus in loci])
                seq_dict = dict((x, z) for x, y, z in aligned_pairs)
                per_site_ref_seq = np.array([[seq_dict[key] for key in range(locus-extra_window, locus+motif_length+extra_window)] for locus in loci])
            except:
                breakpoint()
                pass
        
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
                    dict_events = align_signal.get_events(pod5_record.signal, read.get_tag("mv"), read.get_tag("ts"), rev_loci, motif_length, extra_window, signal_stats=True)
                    
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
                    breakpoint()
                    continue

    
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
