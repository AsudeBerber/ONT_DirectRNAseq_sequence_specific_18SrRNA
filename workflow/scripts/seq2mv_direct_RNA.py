__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

"""  plots signal/squigle for one given read against basecalled sequence
bam file has to include following tags: mv, ts, ns (from basecaller), MD(minimap --MD) & also needs to be indexed!
runtime ~ 30s, faster when exact region is given, low memory usage """

# calls align_signal as module, has to be in same directory

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


##########################################################################
#  not executed, but might be useful

#creates textfile with all read_ids within bam file
def read_id_list_bam(sample=None):
    try:
        f = open("resources/read_id_list_bam.txt", "x")
    except:
        # if file exists, it is wiped
        f = open("resources/read_id_list_bam.txt", "w")
        f.write("")
        f.close()

    f = open("resources/read_id_list_bam.txt", "a")

    with ps.AlignmentFile(f"{sample}") as samfile:
        samfile.fetch()
        for reads in samfile.fetch():
            f.write(str(reads.query_name) +"\n")
            print(reads.query_name)
    f.close()


#writes generated array to txt file
def seq2mv_to_txt(seq2mv):
    try:
        a = open("resources/seq2mv.txt", "x")
    except: 
        # if file exists, it is wiped
        a = open("resources/seq2mv.txt", "w")
        a.write("")
        a.close()

    a = open("resources/seq2mv.txt", "a")
    for line in seq2mv:
        a.write(" ".join(line) + "\n")
    a.close()

##########################################################################
    

#assigns every base a color for annotation
def base_color(base):
    if base == "A": col = "#E9452C"
    elif base == "C": col = "#4596F7" 
    elif base == "T": col = "#55B73B"
    elif base == "U": col = "#F3B642" 
    elif base == "G": col = "#F3B642" 
    elif base == "X": col = "#00ffc3" #for acetylation?
    else: col = "black"
    return col    

#gets movetable (mv), ts and corresponding base sequence (seq) for given read id, as well as aligned pairs
def bam_aligned(sample, read_ids, region, pos):

    samfile = ps.AlignmentFile(f"{sample}")

    max_reads = samfile.mapped
    i = 0

    # loops through samfile until it finds matching read, rather slow 
    for read in samfile.fetch(region = region):
        if read.query_name == read_ids:
            read_ID = read.query_name
            seq = read.query_sequence
            
            # Workaround in cases where two ts tags per read exists:
            # read.set_tag("ts", None) #first ts tag is transcript strand(+|-), has to be removed

            mv = read.get_tag("mv")
            ts = read.get_tag("ts")

            ref_seq = read.get_reference_sequence()

            # list in form [query_pos, ref_pos, ref_base]; all 0-based
            aln_pairs = np.array(read.get_aligned_pairs(with_seq = True, matches_only = False))
            
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




# alignes reference & query sequence to matching positions in signal
# seq2mv: [start_signal, end_signal, position in query_sequence]; this is only done for the given window arround the required position
# rev_loci: position when going from 3' -> 5' instead of 5' -> 3'
def seq_to_mv(reads_ids, region, sample, seq=None, mv=None, ts=0, fetch = True, ref_pos=1841, range = 8):
    if fetch == True:
        seq, mv, ts, aln_pairs, ref_seq, read = bam_aligned(sample, reads_ids, region, ref_pos)
        pass
    if fetch == False:
        seq, mv, ts = args.seq, args.mv, args.ts

    q_pos, r_pos, rev_loci = align_signal.get_loci(read=read, pairs=aln_pairs, wd= range, motif_length=1, ref_pos=[ref_pos])                                                                                                                                                                            

    seq2mv, rev_loci = align_signal.access_mv(signal = None, moves = mv, offset = ts,
                                                               rev_loci = rev_loci, motif_length = 1, extra_window = range,
                                                                 read = read, mode = "single_read")
    return seq2mv, seq, aln_pairs


#plots array of [start, end, base] to position (start - end) in signal
def plot_signal_plus_seq(seq2mv, read_ids, pos, qseq, aln_pairs, range_bp, sequencer, full_read=False, range_var = "bases", pod5_dir = "resources/pod5/p2s"):
     
    if pod5_dir == None:
        pod5_dir = "resources/pod5/p2s"
    else: 
        pod5_dir = f"{pod5_dir}"

    # query sequence from 3' end on    
    rev_qseq = qseq[::-1]

    for filename in os.listdir(pod5_dir): #loops through all pod5 files in folder 
        pod5_file = os.path.join(pod5_dir, filename)
        with p5.Reader(pod5_file) as reader:
            # Read the selected read from the pod5 file
            # next() is required here as Reader.reads() returns a Generator
            try:
                read = next(reader.reads([read_ids]))
                print(f"read found in the following pod5_file: {filename}")
            except:
                continue
            
            signal = read.signal

            time = np.arange(len(signal)) #arbitrary time units

            # gets range for plotting signal
            start_pos = seq2mv[0][0]
            end_pos = seq2mv[-1][1]
            
            time_slice = time[start_pos:end_pos+1]
            signal_slice = signal[start_pos:end_pos+1]

        #constructs colormap for plot
        cmap_plot = list(range(range_bp*2 + 1))
        for i in  cmap_plot:
            if i == range_bp:
                cmap_plot[i] = "#b2182b"
            elif i %2 == 0:
                cmap_plot[i] = "#2166ac"
            else:
                cmap_plot[i] = "#67a9cf"

        # there is probably a better way to inverse this
        ref_seq_rev = aln_pairs[:,2]
        ref_seq_rev = ref_seq_rev[::-1]

        # Plot using matplotlib

        # for powerpoint slide:
        # fig, ax = plt.subplots(figsize=(18, 12))    
        # viridis = colormaps["viridis"].resampled(range_bp*2 +1)

        fig, ax = plt.subplots(figsize=(18, 4))
        xticks = []

        # plots signal
        ax.plot (time_slice, signal_slice,linewidth = 1, color = "#B9B9B9", zorder = 1)
        for i, [start, stop, rev_pos] in enumerate(seq2mv):
            signal_slice_base = signal[start:stop]
            time_slice_base = time [start:stop]
            
            #plots individual measure point on top of signal
            ax.scatter(time_slice_base, signal_slice_base,
            linewidth = 1, marker= "o", facecolor = cmap_plot[i], zorder = 2, alpha = 0.5, edgecolor = "none")
            x_coord = (start+stop)/2 

            # for ref_seq annotation
            ax.annotate("Reference sequence:", xy = (-0.05, -0.12), xycoords=("axes fraction", "axes fraction"), ha = "center", color = "grey")
            ax.annotate(ref_seq_rev[rev_pos], xy = (x_coord, -0.12), fontsize = 8, xycoords=("data", "axes fraction"), ha = "center", color = "grey")

            # read seq annotation, also gives 9-mer
            ax.annotate(rev_qseq[rev_pos], xy = (x_coord, 0.02), fontsize = 8, xycoords=("data", "axes fraction"), ha = "center", color = base_color(rev_qseq[rev_pos]))
            ax.annotate(rev_qseq[rev_pos:rev_pos-range_bp-1:-1], xy = (x_coord, 0.97), fontsize = 6, xycoords=("data", "axes fraction"), ha = "center", color = "black")
            ax.annotate(i+1, xy= (x_coord, -0.04), fontsize = 8, xycoords=("data", "axes fraction"), ha = "center")

            # plots vertical line between each base position + makes x axis tick there
            ax.axvline(start-0.5, linestyle = ":", linewidth = 0.5, color = "lightgrey")
            xticks.append(start-0.5)

            if i == len(seq2mv)-1: # plots axline after last base
                ax.axvline(stop+0.5, linestyle = ":", linewidth = 0.5, color = "lightgrey")
                xticks.append(stop+0.5)
    
        # titles, labeling
        ax.margins(0.05, 0.1)
        ax.set(xlabel = "base", ylabel = "signal (pA)")
        plt.title(str("Read ID: "+ read_ids))
        ax.annotate(f"18S rRNA transcript - position {pos+1} ± {range_bp} bp", xy= (0.1, 1.04), xycoords="axes fraction", ha = "center")

        # ugly way to set x axis ticks to start/end positions of each base
        ax.set_xticks(ticks = xticks)
        ax.set_xticklabels([])

        # check if plot dir exists, creates it otherwise
        if not os.path.isdir(f"resources/signal/{sequencer}/plots/{read_ids}"):
            os.makedirs(f"resources/signal/{sequencer}/plots/{read_ids}")

        # saves plot
        plt.savefig(f"resources/signal/{sequencer}/plots/{read_ids}/{read_ids}_{pos+1}-pm{range_bp}.svg", dpi = 300) 

# command line arguments
def cmd_parser(argv):
    parser = argparse.ArgumentParser(description="plots electric current/timepoint with corresponding basecalled base for given Read_ID")
    parser.add_argument("--sequencer", help= "name of sequencer (p2i|p2s)")
    parser.add_argument("--sample", help= "Name of sample bam file w/o .bam ending", action="store", required=True)
    parser.add_argument("--readID", help= "Sample ID in bam and pod5 file", action="store", required=True)
    parser.add_argument("--pos", help= "position/index of base in middle, 1-based", type=int, action="store", required=True)
    parser.add_argument("--range", help= "bases to display the left/right of middle base", type=int, action="store", required=True)
    parser.add_argument("--no_fetch", help= "Disables fetching mv, ts, seq from readID - requieres --mv, --ts, --seq", action="store_true")
    parser.add_argument("--seq", help= "complete sequence of read", action="store_true")
    parser.add_argument("--mv", help= "movetable", action="store_true")
    parser.add_argument("--ts", help= "", action="store_true")
    parser.add_argument("--get_readids", help= "store name of all read ids in .txt", action= "store_true")
    parser.add_argument("--region", type= str, action= "store")
    parser.add_argument("--pod5-dir", type= str, action= "store")

    args = parser.parse_args()

    if args.sequencer == None:
        args.sequencer = "p2s"
        print("assuming p2s as sequencer")

    # starts at pos 0
    args.pos = args.pos - 1

    # prints read id list
    if args.get_readids == True:
        read_id_list_bam(args.sequencer, args.sample)

    # when no-fetch is turned on, mv, ts, seq have to be given, for testing purposes
    if args.no_fetch == True:
        fetch = False
    else:
        fetch = True 

    return args, fetch


def main(argv = sys.argv[1:]):
    args, fetch = cmd_parser(argv= argv)

    # gets seq2mv array with positions in signal for given reference position +- range; query sequence and aligned query-reference pairs (query pos - reference position)
    seq2mv, qseq, aln_pairs = seq_to_mv(reads_ids = args.readID, region = args.region, sample = args.sample,
                    seq = args.seq, mv = args.mv, ts = args.ts, fetch = fetch, ref_pos = args.pos, range = args.range) 
    
    # plots requested position of given read
    plot_signal_plus_seq(seq2mv=seq2mv, read_ids = args.readID, pos=args.pos, qseq = qseq, aln_pairs = aln_pairs, range_bp = args.range,
                          sequencer = args.sequencer, full_read=False, range_var = "bases", pod5_dir = args.pod5_dir)
    

if __name__ == "__main__":
    exit(main())