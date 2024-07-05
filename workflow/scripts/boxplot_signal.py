__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

"""
plots boxplot for all features given by .npz file, currently for pos 403, 1337 and 1842
for future: no plotting, instead only loading of .npz file and production of small .csv file containing only data needed to construct boxplot (mean, quartils, limits)
--> could then be loaded in R (currently .npz is to big for that)

! FRAME is currently only the number of bases to plot left/right of the middle position (whole window = 2*FRAME + 1); could be changed to avoid confusion
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import logging
import os


# starts logger for giving status updates while program runs
logger = logging.getLogger(__name__)
logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))


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
    
    # splits .npz in original parts
    features = loaded["feat"]
    qual = loaded["qual"]
    query = loaded["query"]
    # ref contains "None" strings, in cases where deletions occured
    ref = loaded["ref"]
    id = loaded["id"]

    logger.info("npz loaded")

    window_size = ref.shape[1]

    # gets position (in refseq) from read id (at the end of read id)
    pos = []
    for read in id:
        pos.append(read.rsplit(":")[1])

    return features, qual, query, ref, id, window_size, pos


def parse_args(argv):
    """Read arguments from command line."""

    parser = argparse.ArgumentParser()

    parser.add_argument("-w", "--window", type=int, help= "positions to display left/right of middle position")
    parser.add_argument("-f", "--file", type=str)
    parser.add_argument("--no-mmap", action="store_true")
    parser.add_argument("--output-dir")

    args = parser.parse_args()

    return args


# restricts all arrays to relevant frame around middle (acetylated?) position
def slice_bases(event, query, features, ref):
    whole_window = len(query[0])

    # middle position is the investigated position (1842, 1337...)
    middle = whole_window //2
    
    #FRAME: number of bases to plot arround CCG (including CCG), global variable
    #as the beginnings and ends of the arrays are cut off, the index is needed to still assign the correct index in non-cut arrays
    index_bases = np.arange(0,whole_window) [middle-FRAME:middle+FRAME+1] 

    # events and reference sequence arrays are restricted to plotted window (start/end cut off, window given by index)
    sliced_event = features[:,index_bases,event]
    sliced_refseq = ref[:,index_bases]
    return index_bases, sliced_event, sliced_refseq
    

# creates dataframes from arrays for plotting
def filter_by_pos(pos, df_event_pos):
    df = df_event_pos
    df_filtered = df[df["pos"] == str(pos)]
    df_plot = df_filtered.iloc[:,:FRAME*2+1]
    df_plot_flip = np.fliplr(df_plot)
    return df_plot_flip


#makes dataframes to plot data
def make_dfs(event, query, features, ref, pos): 
    index_bases, sliced_event, sliced_ref_seq = slice_bases(event=event, query=query, features=features, ref=ref)

    # data frame with only events or position, events only for window that is plotted
    df_event = pd.DataFrame((sliced_event), columns = list(range(1,FRAME*2+2)))
    df_pos = pd.DataFrame(pos, columns = ["pos"])

    # dataframe with events and corresponding position in ref seq (1842, 1337...); e.g. [feature 0-8, 1842], [feature 0-8, 1337]...
    df_event_pos = pd.concat([df_event, df_pos], axis=1)

    return df_event_pos, sliced_event, sliced_ref_seq, index_bases


# annotates boxplot with corresponding bases from ref seq, this is done the same for all 3 subplots
def boxplot_ann(refseq, axis):
    for i, base in enumerate(refseq.iloc[0]):
        if i == FRAME: 
            axis.annotate(base, xy = (i+1, 0.04), xycoords=("data", "axes fraction"), ha = "center", color = "darkred", fontsize = 12) #center base is colored red
        else:
            axis.annotate(base, xy = (i+1, 0.04), xycoords=("data", "axes fraction"), ha = "center", fontsize = 12)


# creates plot (with 3 subplots for each position)
def make_plot(event, window_size, ref, df_pos, index_bases, df_event_pos, file, output_dir):

    # reference seq is flipped, as the plot goes from 3'-> 5'; corresponding position is again concatenated
    df_refseq = pd.DataFrame(np.fliplr(ref), columns = list(range(1,window_size+1)))
    df_refseq = pd.concat([df_refseq, df_pos], axis=1)

    # splits reversed refseq_df into 3 for individual positions
    df_refseq_430 = df_refseq[df_refseq["pos"] == "430"]
    df_refseq_1337 = df_refseq[df_refseq["pos"] == "1337"]
    df_refseq_1842 = df_refseq[df_refseq["pos"] == "1842"]

    # dfs are restricted to positions to plot (removes first and last positions)
    df_refseq_430_sliced = df_refseq_430.iloc[:,index_bases]
    df_refseq_1337_sliced = df_refseq_1337.iloc[:,index_bases]
    df_refseq_1842_sliced = df_refseq_1842.iloc[:,index_bases]

    # plots feature for position 1337 (ax1), 1842 (ax2) and 430 (ax3) next to each other
    # boxplot is added on top of violin plot, outliers are not shown, as there are far to many in the 4x10E6 reads
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize = (20,4), sharey= True) # plots the following 3 plots next to each other in 1 row, makes them share y-axis

    ## ax1
    # plots for 1337
    ax1.violinplot(filter_by_pos(1337, df_event_pos), showmeans = False, showextrema = False)
    ax1.boxplot(filter_by_pos(1337, df_event_pos), showfliers = False)
    
    # symlog scale = log scale in + and - direction, better for showing boxes
    ax1.set_yscale("symlog")

    # labeling
    boxplot_ann(df_refseq_1337_sliced, ax1)
    ax1.set_xlabel("Steps")
    ax1.set_title(f"Pos 1337 ± {FRAME} bp")

    ## ax2
    # plots for 1842
    ax2.violinplot(filter_by_pos(1842, df_event_pos), showmeans = False, showextrema = False)
    ax2.boxplot(filter_by_pos(1842, df_event_pos), showfliers = False)
    boxplot_ann(df_refseq_1842_sliced, ax2)

    # sharey = TRUE makes all figures the same y-scale, but removes ticks for ax2 + ax3 plots, this is solved with tick_params
    ax2.tick_params(labelleft = True)
    ax2.set_yscale("symlog")
    ax2.set_title(f"Pos 1842 ± {FRAME} bp")
    ax2.set_xlabel("Steps")

    ## ax3
    # plots for 430 pos
    ax3.violinplot(filter_by_pos(430, df_event_pos), showmeans = False, showextrema = False)
    ax3.boxplot(filter_by_pos(430, df_event_pos), showfliers = False)
    boxplot_ann(df_refseq_430_sliced, ax3)
    ax3.tick_params(labelleft = True)
    ax3.set_yscale("symlog")
    ax3.set_title(f"Pos 430 ± {FRAME} bp")
    ax3.set_xlabel("Steps")

    # creates directory if not present
    if not os.path.isdir(f"resources/signal/p2s/signal_summary/{output_dir}"): os.makedirs(f"resources/signal/p2s/signal_summary/{output_dir}")

    # saves figures in same directory with different endings for events
    plt.savefig(f"resources/signal/p2s/signal_summary/{output_dir}/1337_1842_430_event_{event}.svg", dpi = 300)

def main(argv=sys.argv[1:]):

    args = parse_args(argv)

    # frame is set as global variable, as it is needed in almost every function (!FRAME is only the number of bases to plot left/right, whole window is twice as big)
    global FRAME
    FRAME = args.window

    # loads data generated by signal_summary.py
    features, qual, query, ref, id, window_size, pos = load_npz(args.file, args.no_mmap)

    # makes 1 plot for each of the 9 events
    for event in range(0,9):
        logging.info(f"creating plot {event+1} of 9")

        # makes data frames for plotting
        df_event_pos, sliced_event, sliced_ref_seq, index_bases = make_dfs(event= event, query=query, features=features, ref=ref, pos=pos)

        # plotting
        make_plot(event = event, window_size = window_size, ref=ref, df_pos = df_event_pos, index_bases=index_bases,
                   df_event_pos=df_event_pos, file = args.file, output_dir= args.output_dir)


if __name__ == "__main__":
     exit(main())