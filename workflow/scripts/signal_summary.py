# Code mostly written by Christoph Engelhardt

import numpy as np
import pod5 as p5
import pysam
from tqdm import tqdm
import re
import argparse
import sys


# pod5_file = "data/pod5/PAQ77977_pass_barcode01_362f656f_255b3204_0.pod5"
# bam_file = "data/mapped/PAQ77977_pass_barcode01_362f656f_255b3204_0.bam"
# motif = "CCG" # "HCG" is possible ("[ACT]CG"), highest specificity is "CCG"
# window_size = 21
# npz_file = ""

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
    ref_loci = []
    for m in motif.finditer(read.get_reference_sequence()):
        ref_loci.append(m.start())

    loci = [pairs[locus] for locus in ref_loci]
    # Remove loci that are not present on the query or too close to the ends of the alignment
    loci = [locus for locus in loci if locus is not None and locus > wd and locus < read.alen - wd - ml]

    return loci


def main(argv=sys.argv[1:]):

    args = parse_args(argv=argv)

    pod5_file = args.pod5
    bam_file = args.bam
    motif = args.motif
    window_size = args.window
    npz_file = args.output

    compiled_motif = re.compile(motif)
    motif_length = len(motif)
    extra_window = int((window_size - motif_length) / 2)

    with p5.Reader(pod5_file) as pod5, pysam.AlignmentFile(bam_file, mode = "rb", check_sq=False) as bam: 

        features, qual, query_seq, ref_seq, id = [], [], [], [], []

        for read in tqdm(bam):
            if read.is_unmapped:
                continue

            # get loci on the reference matching the motif
            aligned_pairs = read.get_aligned_pairs(with_seq=True)
            pairs_dict = dict((y-read.reference_start, x) for x, y, z in aligned_pairs if y is not None)
            loci = get_loci(read, pairs_dict, compiled_motif, extra_window, motif_length)
        
            if len(loci) == 0:
                continue

            # extract features from bam file
            per_site_qual = np.array([list(read.qual[locus-extra_window: locus+motif_length+extra_window]) for locus in loci])
            per_site_query_seq = np.array([list(read.query_sequence[locus-extra_window: locus+motif_length+extra_window]) for locus in loci])
            seq_dict = dict((x, z) for x, y, z in aligned_pairs)
            per_site_ref_seq = np.array([[seq_dict[key] for key in range(locus-extra_window, locus+motif_length+extra_window)] for locus in loci])

            # extract features from pod5 file
            pod5_record = next(pod5.reads(selection=[read.query_name])) 
            events = get_events(pod5_record.signal, read.get_tag("mv"), read.get_tag("ts"))
            per_site_features = np.array([events[locus-extra_window: locus+motif_length+extra_window] for locus in loci])
            per_site_id = np.array([read.query_name + ':' + str(locus) for locus in loci])

            features.append(per_site_features)
            qual.append(per_site_qual)
            query_seq.append(per_site_query_seq)
            ref_seq.append(per_site_ref_seq)
            id.append(per_site_id)


    features = np.vstack(features)
    qual = np.vstack(qual)
    query_seq = np.vstack(query_seq)
    ref_seq = np.vstack(ref_seq)
    id = np.hstack(id)

    # np.savez_compressed(npz_file, features, qual, query_seq, ref_seq, id)
    np.savetxt(npz_file, features, qual, query_seq, ref_seq, id)

    return 0


if __name__ == "__main__":
    exit(main())