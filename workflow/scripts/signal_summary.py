# Code mostly written by Christoph Engelhardt

import numpy as np
import pod5 as p5
import pysam
from tqdm import tqdm
import re
import argparse
import sys
import os
import pdb
import time


pod5_file = "resources/pod5/p2s/"
bam_file = "resources/alignments/p2s_aligned_subsample_000001.bam"
motif = "CCG" # "HCG" is possible ("[ACT]CG"), highest specificity is "CCG"
window_size = 21
npz_file = f"resources/results/p2s/{motif}_window_{window_size}.txt"

# phred={}
# for x in range(0,94):
#     phred.update({(chr(x+33).encode('ascii').decode("utf-8")): x}) 


# handle sys args
def parse_args(argv):
    """Read arguments from command line."""

    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--pod5", type=str)
    parser.add_argument("-b", "--bam", type=str)
    parser.add_argument("-o", "--output", type=str)
    parser.add_argument("-m", "--motif", type=str, default="CCG")
    parser.add_argument("-w", "--window", type=int, default=21)

    args = parser.parse_args()

    return args


def get_events(signal, moves, offset):
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
    
    data = np.zeros((rlen,9))
    
    for i in range(rlen-1):
        prev = move_index[i]*stride+offset
        sig_end = move_index[i+1]*stride+offset
        
        sig_len = sig_end-prev
        data[i, 4]=np.log10(sig_len)
        data[i, 5]=np.mean(signal[prev:sig_end])
        data[i, 6]=np.std(signal[prev:sig_end])
        data[i, 7]=np.median(signal[prev:sig_end])
        data[i, 8]=np.median(np.abs(signal[prev:sig_end]-data[i, 4]))
        
        # get the mean signal for each quarter of the base signal
        for j in range(4):
            tmp_cnt=0
            for t in range(j*sig_len//4,min(sig_len, (j+1)*sig_len//4)):
                data[i, j]+=signal[t+prev]
                tmp_cnt+=1
            data[i, j]=data[i, j]/tmp_cnt

    return data


def get_loci(read, pairs, motif, wd, ml):
    """
    find positions that match motif
    """    
    ref_ac1 = range (1335,1338)
    ref_ac2 = range (1840,1843)
    ref_loci = []
    for m in ([*ref_ac1] + [*ref_ac2]):
        ref_loci.append(m)

    print(pairs)
    loci = [pairs[locus] for locus in ref_loci]
    # Remove loci that are not present on the query or too close to the ends of the alignment
    # loci = [locus for locus in loci if locus is not None and locus > wd-1 and locus < read.alen - wd - ml]
    loci = [locus for locus in loci if locus is not None and locus > wd-1 and locus < read.alen - wd - ml-1]

    return loci


def main(argv=sys.argv[1:]):

    args = parse_args(argv=argv)

    # pod5_file = args.pod5
    # bam_file = args.bam
    # motif = args.motif
    # window_size = args.window
    # npz_file = args.output

    compiled_motif = re.compile(motif)
    motif_length = len(motif)
    extra_window = int((window_size - motif_length) / 2)

    with pysam.AlignmentFile(bam_file, mode = "rb", check_sq=False) as bam: 

        features, qual, query_seq, ref_seq, id = [], [], [], [], []

        for read in tqdm(bam):
            if read.is_unmapped:
                continue
            
            # get loci on the reference matching the motif
            aligned_pairs = read.get_aligned_pairs(with_seq=True)
            ac_ccg= list(filter(lambda x: x[1] in range(1335,1338) or x[1] in range (1840,1843), aligned_pairs))

            pairs_dict = dict((y, x) for x, y, z in ac_ccg if y is not None)
            loci = get_loci(read, pairs_dict, compiled_motif, extra_window, motif_length)
        
            if len(loci) == 0:
                continue
            

            # extract features from bam file
            try:
                per_site_qual = [list(read.query_qualities[locus-extra_window: locus+motif_length+extra_window]) for locus in loci]
                per_site_query_seq = list(read.query_sequence[locus-extra_window: locus+motif_length+extra_window] for locus in loci)
                seq_dict = dict((x, z) for x, y, z in aligned_pairs)
                per_site_ref_seq = [[seq_dict[key] for key in range(locus-extra_window, locus+motif_length+extra_window)] for locus in loci]
            except:
                breakpoint()
                pass

            # extract features from pod5 file

            # with p5.DatasetReader(args.pod5) as dataset:
            time_st = time.process_time()
            #loops through all pod5 files in folder 
            for filename in os.listdir(args.pod5): #loops through all pod5 files in folder 
                pod5_file = os.path.join(args.pod5, filename)
                with p5.Reader(pod5_file) as pod5:
                    # Read the selected read from the pod5 file
                    # next() is required here as Reader.reads() returns a Generator
                    try:
                        pod5_record = next(pod5.reads(selection=[read.query_name])) 
                        events = get_events(pod5_record.signal, read.get_tag("mv"), read.get_tag("ts"))
                        per_site_features = np.array([events[locus-extra_window: locus+motif_length+extra_window] for locus in loci])
                        per_site_id = np.array([read.query_name + ':' + str(locus) for locus in loci])

                        features.append(per_site_features)
                        qual.append(per_site_qual)
                        query_seq.append(per_site_query_seq)
                        ref_seq.append(per_site_ref_seq)
                        id.append(per_site_id)
                    except:
                        continue
            time_pod = time.process_time() - time_st
            print (f"time pod5 loop: {time_pod}")

            # breakpoint()
            features = np.vstack(features)
            qual = np.vstack(qual)
            query_seq = np.vstack(query_seq)
            ref_seq = np.vstack(ref_seq)
            id = np.hstack(id)

    breakpoint()

    # checks if results folder exists, creates otherwise
    # check if plot dir exists, creates it otherwise
    if not os.path.isdir(npz_file): os.makedirs(npz_file)
        

    np.savez_compressed(npz_file, features, qual, query_seq, ref_seq, id)
    np.savetxt(npz_file, features, qual, query_seq, ref_seq, id)

    return 0


if __name__ == "__main__":
    exit(main())
