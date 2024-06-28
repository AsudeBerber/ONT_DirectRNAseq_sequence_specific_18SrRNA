__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import logging
import os
import pathlib
import pdb

logger = logging.getLogger(__name__)
logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
# npz_file = f"resources/results/p2s/CCG_window_21_p2s_aligned_sorted.npz"
# frame=9


def load_npz(npz_file, no_mmap=False):
    logger.info("loading npz")
    ## mmap_mode writes in binary to disk to save RAM, can be deleted for smaller npz files
    if no_mmap == True:
        try:
            loaded = np.load(npz_file)
            print("loading in no mmap mode")
        except MemoryError:
            raise MemoryError("Not enough memory for loading without mmap-mode")
        
    
    else:
        loaded = np.load(npz_file, mmap_mode= "r")

    """
    FEATURES:
    3D array in shape (read,base,event[0:8])
    EVENTS:
    4: log10 signal length
    5: mean signal intensity
    6: standard deviation of signal intensity
    7: median signal intensity
    8: median absolute deviation of signal intensity
    0-3: mean signal intensity for each quartile
    """
    
    features = loaded["feat"]
    qual = loaded["qual"]
    query = loaded["query"]
    # ref contains "None" strings, in cases where deletions occured
    ref = loaded["ref"]
    id = loaded["id"]

    logger.info("npz loaded")

    window_size = ref.shape[1]

    pos = []
    for read in id:
        pos.append(read.rsplit(":")[1])

    return features, qual, query, ref, id, window_size, pos


def parse_args(argv):
    """Read arguments from command line."""

    parser = argparse.ArgumentParser()

    parser.add_argument("-w", "--window", type=int)
    parser.add_argument("-f", "--file", type=str)
    parser.add_argument("--no-mmap", action="store_true")
    parser.add_argument("--output-dir")

    args = parser.parse_args()

    return args

#frame: #number of bases to plot arround CCG (including CCG)
def slice_bases(event, query, features, ref):
    whole_window = len(query[0])
    middle = whole_window //2
    index_bases = np.arange(0,whole_window) [middle-frame:middle+frame+1]
    sliced_event = features[:,index_bases,event]
    sliced_refseq = ref[:,index_bases]
    return index_bases, sliced_event, sliced_refseq
    
def filter_by_pos(pos, df_event_pos):
    df = df_event_pos
    df_filtered = df[df["pos"] == str(pos)]
    df_plot = df_filtered.iloc[:,:frame*2+1]
    df_plot_flip = np.fliplr(df_plot)
    return df_plot_flip



def make_dfs(event, query, features, ref, pos): 
    index_bases, sliced_event, sliced_ref_seq = slice_bases(event=event, query=query, features=features, ref=ref)

    df_event = pd.DataFrame((sliced_event), columns = list(range(1,frame*2+2)))
    df_pos = pd.DataFrame(pos, columns = ["pos"])

    df_event_pos = pd.concat([df_event, df_pos], axis=1)

    return df_event_pos, sliced_event, sliced_ref_seq, index_bases

def boxplot_ann(refseq, axis):
    for i, base in enumerate(refseq.iloc[0]):
        if i == frame: 
            axis.annotate(base, xy = (i+1, 0.04), xycoords=("data", "axes fraction"), ha = "center", color = "darkred", fontsize = 12)
        else:
            axis.annotate(base, xy = (i+1, 0.04), xycoords=("data", "axes fraction"), ha = "center", fontsize = 12)


#1337 | 1842
def make_plot(event, window_size, ref, df_pos, index_bases, df_event_pos, file, output_dir):
    # df_event_filtered = filter_by_pos(pos)

    df_refseq = pd.DataFrame(np.fliplr(ref), columns = list(range(1,window_size+1)))
    df_refseq = pd.concat([df_refseq, df_pos], axis=1)

    df_refseq_430 = df_refseq[df_refseq["pos"] == "430"]
    df_refseq_1337 = df_refseq[df_refseq["pos"] == "1337"]
    df_refseq_1842 = df_refseq[df_refseq["pos"] == "1842"]
    df_refseq_430_sliced = df_refseq_430.iloc[:,index_bases]
    df_refseq_1337_sliced = df_refseq_1337.iloc[:,index_bases]
    df_refseq_1842_sliced = df_refseq_1842.iloc[:,index_bases]


    fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize = (20,4), sharey= True)
    ax1.violinplot(filter_by_pos(1337, df_event_pos), showmeans = False, showextrema = False)
    ax1.boxplot(filter_by_pos(1337, df_event_pos), showfliers = False)
    boxplot_ann(df_refseq_1337_sliced, ax1)
    ax1.set_yscale("symlog")
    ax1.set_xlabel("Steps")
    ax1.set_title(f"Pos 1337 ± {frame} bp")


    ax2.violinplot(filter_by_pos(1842, df_event_pos), showmeans = False, showextrema = False)
    ax2.boxplot(filter_by_pos(1842, df_event_pos), showfliers = False)
    boxplot_ann(df_refseq_1842_sliced, ax2)
    ax2.tick_params(labelleft = True)
    ax2.set_yscale("symlog")
    ax2.set_title(f"Pos 1842 ± {frame} bp")
    ax2.set_xlabel("Steps")

    ax3.violinplot(filter_by_pos(430, df_event_pos), showmeans = False, showextrema = False)
    ax3.boxplot(filter_by_pos(430, df_event_pos), showfliers = False)
    boxplot_ann(df_refseq_430_sliced, ax3)
    ax3.tick_params(labelleft = True)
    ax3.set_yscale("symlog")
    ax3.set_title(f"Pos 430 ± {frame} bp")
    ax3.set_xlabel("Steps")

    # creates directory if not present
    if not os.path.isdir(f"resources/signal/p2s/signal_summary/{output_dir}"): os.makedirs(f"resources/signal/p2s/signal_summary/{output_dir}")

    plt.savefig(f"resources/signal/p2s/signal_summary/{output_dir}/1337_1842_430_event_{event}.svg", dpi = 300)

def main(argv=sys.argv[1:]):

    args = parse_args(argv)

    global frame
    frame = args.window

    features, qual, query, ref, id, window_size, pos = load_npz(args.file, args.no_mmap)

    for event in range(0,9):
        logging.info(f"creating plot {event+1} of 9")
        df_event_pos, sliced_event, sliced_ref_seq, index_bases = make_dfs(event= event, query=query, features=features, ref=ref, pos=pos)
        make_plot(event = event, window_size = window_size, ref=ref, df_pos = df_event_pos, index_bases=index_bases,
                   df_event_pos=df_event_pos, file = args.file, output_dir= args.output_dir)

if __name__ == "__main__":
     exit(main())