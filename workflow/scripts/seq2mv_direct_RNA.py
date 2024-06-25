__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

#############  plots signal/squigle for one given read against basecalled sequence
###############  bam file has to include following tags: mv, ts, ns (from basecaller), MD(minimap --MD)
#################  runtime ~ 30s, faster when exact region is given, low memory usage

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
            return(seq, mv, ts, aln_pairs, ref_seq)
        else: 
            #removing the else part makes the code only 1s faster
            i = i+1
            k = i/500000
            if k.is_integer():
                # print("currently at " +read.query_name + "\n")
                print(f"reads compared: {i} of max. {max_reads}")
            continue



# assigns individual bases in sequence to corresponding part of signal through movetable information
def seq_to_mv(reads_ids, region, sample, seq=None, mv=None, ts=0, fetch = True, pos=42):
    if fetch == True:
        seq, mv, ts, aln_pairs, ref_seq = bam_aligned(sample, reads_ids, region, pos)
    
    seq = seq[::-1] #sequence order is 5' -> 3', mv and signal are 3' -> 5': therefore sequence is turned around
    s = mv[0] #stride length
    p = 1 #itinerates through movetable array
    x = 0 #number of additional strides (stride amount - 1)
    start = ts + 1

    seq2mv = np.array([[1, ts, "-", "-"]])
    for i, base in enumerate(seq):
        # print (base)
        while p < len(mv)-1 and mv[p + 1] == 0: #last movetable index(p): mv[p+1] doesn't exist, would cause error
            x = x + 1 #counts 0's
            p = p + 1 #0 found -> movetable index moves by one
            # print (f"{p},{x}\n")
        p = p + 1 
        # print(f"p={p}")
        # print(f"x={x}")
        end = start + (s-1) + (x*s) #stride 1: start + 4 (as start number is already first position in stride), all further strides: additional +5
        x = 0 #resets number of additional strides
        ref_base = ref_seq [i]
        seq2mv = np.append(values=[[start, end, base, ref_base]], arr=seq2mv, axis=0) # appends 
        start = end+1 #next base starts 1 after end of previous base

    print ("sequence-to-signal alignment finished,", len(seq), "bases, signal length =", end)    
    return seq2mv, pos_read


def reverse_seq_mv(seq2mv):
    k = int(seq2mv[-1][1]) #signal range
    rev_seq2mv = np.flip(seq2mv,0)
    for base_data in rev_seq2mv:
        start_old = int(base_data[0])
        end_old = int(base_data[1])
        base_data[0] = k - end_old #start and end switch when string is read in other direction, former last base is now first and vice versa
        base_data[1] = k - start_old
    return rev_seq2mv


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


#plots array of [start, end, base] to position (start - end) in signal
def plot_signal_plus_seq(seq2mv, read_ids, pos, pos_read, range_bp, sequencer, full_read=False, range_var = "bases", pod5_dir = "pod5"):
     
    if pod5_dir == None:
        pod5_dir = "resources/pod5"
    else: 
        pod5_dir = f"{pod5_dir}"

    #for output file naming

    start = pos_read - range_bp
    end = pos_read + range_bp

    start_b = start
    end_b = end

    for filename in os.listdir(pod5_dir): #loops through all pod5 files in folder 
        pod5_file = os.path.join(pod5_dir, filename)
        with p5.Reader(pod5_file) as reader:
            # Read the selected read from the pod5 file
            # next() is required here as Reader.reads() returns a Generator
            try:
                read = next(reader.reads([read_ids]))
            except:
                continue
            #
            signal = read.signal
            signal = signal[::-1] #reverses signal as RNA is sequenced 3' -> 5' but sequence is read 5' -> 3'
#
            # when range is set by bases, starting and ending time are pulled from seq2mv array --> base range is converted to time range
            time = np.arange(len(signal)) #arbitrary time units
            if range_var == "bases":
                print("plot for base range", start, "-", end)
                start = int(seq2mv[start][0])
                end = int(seq2mv[end][1])
                print("corresponding time range:", start,"-", end)
            elif range_var == "time":
                pass
            else:
                raise Exception ("range_var has to be bases or time")
            #
            if full_read == True:
                start = 1
                end = len(signal)
            #
            signal_slice = signal[start:end]
            time = np.arange(len(signal)) #arbitrary time units
            time_slice = time[start:end]

            
            #

        #constructs colormap for plot
        cmap_plot = list(range(range_bp*2 + 1))
        for i in  cmap_plot:
            if i == range_bp:
                cmap_plot[i] = "#b2182b"
            elif i %2 == 0:
                cmap_plot[i] = "#2166ac"
            else:
                cmap_plot[i] = "#67a9cf"
                

        # Plot using matplotlib

        # for powerpoint slide:
        # fig, ax = plt.subplots(figsize=(18, 12))    
        # viridis = colormaps["viridis"].resampled(range_bp*2 +1)

        fig, ax = plt.subplots(figsize=(18, 4))
        #
        ax.plot (time_slice, signal_slice,linewidth = 1, color = "#B9B9B9", zorder = 1)
        for i, base in enumerate(range(start_b, end_b +1)):
                    # like slice above, just for every base -> signal per base can be colored differently
                    base_start = int(seq2mv[base][0])
                    base_end = int(seq2mv[base][1])
                    signal_slice_base = signal[base_start:base_end]
                    time_slice_base = time [base_start:base_end]
                    ax.scatter(time_slice_base, signal_slice_base,
                                linewidth = 1, marker= "o", facecolor = cmap_plot[i], zorder = 2, alpha = 0.5, edgecolor = "none")
                                # linewidth = 1, marker= "o", facecolor = viridis.colors[i], zorder = 2, alpha = 0.5, edgecolor = "none", s = 600)
            
        
        
        
        ax.margins(0.05, 0.1)
        ax.set(xlabel = "base", ylabel = "signal (pA)")
        plt.title(str("Read ID: "+ read_ids))
        ax.annotate(f"18S RNA transcript - position {pos+1} ± {range_bp} bp", xy= (0.5, 0.95), xycoords="axes fraction", ha = "center")
#
        #annotation of bases to signal plot
        i = -range_bp
        xticks = []
        for base_data in seq2mv:
            x_coord = (int(base_data[0])+int(base_data[1]))/2 
            if x_coord < start:
                pass
            elif x_coord > start and x_coord < end: 
                # read seq
                ax.annotate(base_data[2], xy = (x_coord, 0.02), fontsize = 8, xycoords=("data", "axes fraction"), ha = "center", color = base_color(base_data[2]))
                # ref seq
                # ax.annotate(base_data[3], xy = (x_coord, 0.06), fontsize = 8, xycoords=("data", "axes fraction"), ha = "center", color = "grey")
                ax.annotate(i, xy= (x_coord, -0.04), fontsize = 8, xycoords=("data", "axes fraction"), ha = "center")
                ax.axvline(int(base_data[0])-0.5, linestyle = ":", linewidth = 0.5, color = "lightgrey")
                xticks.append(int(base_data[0])-0.5)
                i = i + 1
            else:
                ax.axvline(int(base_data[0])-0.5, linestyle = ":", linewidth = 0.5, color = "lightgrey")
                xticks.append(int(base_data[0])-0.5)
                break
        
        ax.set_xticks(ticks = xticks)
        ax.set_xticklabels([])

        # check if plot dir exists, creates it otherwise
        if not os.path.isdir(f"resources/signal/{sequencer}/plots/{read_ids}"):
            os.makedirs(f"resources/mapped/signal/{sequencer}/plots/{read_ids}")

        plt.savefig(f"resources/signal/{sequencer}/plots/{read_ids}/{read_ids}_{pos+1}-pm{range_bp}.svg", dpi = 300)           
        plt.show()

def main(argv = sys.argv[1:]):
    args, fetch = cmd_parser(argv= argv)

    seq2mv, pos_read = seq_to_mv(reads_ids = args.readID, region = args.region, sample = args.sample,
                    seq = args.seq, mv = args.mv, ts = args.ts, fetch = fetch, pos = args.pos)  


    # as RNA is sequenced 3' -> 5', but convention + basecalled sequence is 5' -> 3' the seq2mv has to be "turned around" 
    # for 5' -> 3' sequencing direction, only this part has to be removed
    rev_seq2mv=reverse_seq_mv(seq2mv)


    #plots array to signal
    plot_signal_plus_seq(rev_seq2mv, read_ids = args.readID, pos = args.pos, pos_read = pos_read,
                        range_bp = args.range, sequencer = args.sequencer, full_read=False, pod5_dir = args.pod5_dir)


##########################################################################################################################################################
# exemplary command line argument:
#
# python seq2mv_direct_RNA.py RNAseq_test_noisy_correction_sorted becda29e-a3f7-4647-a0b4-1dfd9134c3cc 300 350 --region SIRV6001:16-6000
#
# if known, specify region/reference span(=chromosome) (e.g. "ERCC-00076[:start-end]" ; []<- this part is optional"); time for execution will be much slower otherwise (minutes range)
# for region argument thousand seperators have to be removed, e.g. SIRV6001:16-6000 works, SIRV6001:16-6.000 (as given by IGV) does not
######################################################################################################################################################
        

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
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)


#saves seq2mv array (base aligned to signal position) to txt file in ../resources
# seq2mv_to_txt(seq2mv) 
# seq2mv_to_txt(rev_seq2mv)
       