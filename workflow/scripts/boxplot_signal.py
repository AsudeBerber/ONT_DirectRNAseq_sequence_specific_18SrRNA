__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import pdb



motif = "CCG"
window_size = 21
npz_file = f"resources/results/p2s/CCG_window_21_p2s_aligned.npz"
arround=9
event = 5
# event = args.feature

loaded = np.load(npz_file, mmap_mode= "w+")

def parse_args(argv):
    """Read arguments from command line."""

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--feature", type=int)

    args = parser.parse_args()

    return args

#arround: #bases to plot arround ac4C
def slice_bases(event, arround):
    whole_window = len(query[0])
    middle = whole_window //2
    index_bases = np.arange(0,whole_window) [middle-arround:middle+arround+1]
    sliced_event = features[:,index_bases,event]
    sliced_refseq = ref[:,index_bases]
    return index_bases, sliced_event, sliced_refseq
    
def filter_by_pos(pos):
    df = df_event_pos
    df_filtered = df[df["pos"] == str(pos)]
    df_plot = df_filtered.iloc[:,:arround*2+1]
    df_plot_flip = np.fliplr(df_plot)
    return df_plot_flip


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

pos = []
for read in id:
    pos.append(read.rsplit(":")[1])


index_bases, sliced_event, sliced_ref_seq = slice_bases(event=event, arround=arround)

df_event = pd.DataFrame((sliced_event), columns = list(range(1,arround*2+2)))
df_id = pd.DataFrame(id)
df_pos = pd.DataFrame(pos, columns = ["pos"])

df_event_pos = pd.concat([df_event, df_pos], axis=1)

def boxplot_ann(refseq, axis):
    for i, base in enumerate(refseq.iloc[0]):
        if i == arround: 
            axis.annotate(base, xy = (i+1, 0.96), xycoords=("data", "axes fraction"), ha = "center", color = "darkred")
        else:
            axis.annotate(base, xy = (i+1, 0.96), xycoords=("data", "axes fraction"), ha = "center")


#1337 | 1842
def make_plot(event):
    range_window_plot = ref.shape[1] 
    df_event_filtered = filter_by_pos(pos)

    df_refseq = pd.DataFrame(np.fliplr(ref), columns = list(range(1,ref.shape[1]+1)))
    df_refseq = pd.concat([df_refseq, df_pos], axis=1)

    df_refseq_430 = df_refseq[df_refseq["pos"] == "430"]
    df_refseq_1337 = df_refseq[df_refseq["pos"] == "1337"]
    df_refseq_1842 = df_refseq[df_refseq["pos"] == "1842"]
    df_refseq_430_sliced = df_refseq_430.iloc[:,index_bases]
    df_refseq_1337_sliced = df_refseq_1337.iloc[:,index_bases]
    df_refseq_1842_sliced = df_refseq_1842.iloc[:,index_bases]


    fig, (ax1, ax2) = plt.subplots(1,2)
    ax1.violinplot(filter_by_pos(1337), showmeans = False, showextrema = False)
    ax1.boxplot(filter_by_pos(1337), showfliers = False)
    boxplot_ann(df_refseq_1337_sliced, ax1)
    ax1.set_yscale("symlog")
    ax1.set_xlabel("time steps")
    ax1.set_title(f"Pos {1337} ± {arround} bp")


    ax2.violinplot(filter_by_pos(1842), showmeans = False, showextrema = False)
    ax2.boxplot(filter_by_pos(1842), showfliers = False)
    boxplot_ann(df_refseq_1842_sliced, ax2)
    ax2.set_yscale("symlog")
    ax2.set_title(f"Pos 1842 ± {arround} bp")
    ax2.set_xlabel("time steps")

    plt.savefig(f"resources/signal/signal_summary/1337_1842_event_{event}.svg", dpi = 300)

    fig, (ax1, ax2) = plt.subplots(1,2)

    ax1.violinplot(filter_by_pos(430), showmeans = False, showextrema = False)
    ax1.boxplot(filter_by_pos(430), showfliers = False)
    boxplot_ann(df_refseq_430_sliced, ax1)
    ax1.set_yscale("symlog")
    ax1.set_title(f"Pos 430 ± {arround} bp")
    ax1.set_xlabel("time steps")


    ax2.violinplot(filter_by_pos(1842), showmeans = False, showextrema = False)
    ax2.boxplot(filter_by_pos(1842), showfliers = False)
    boxplot_ann(df_refseq_1842_sliced,ax2)
    ax2.set_yscale("symlog")
    ax2.set_title(f"Pos 1842 ± {arround} bp")
    ax2.set_xlabel("time steps")

    plt.savefig(f"resources/signal/signal_summary/430_1842_event_{event}.svg", dpi = 300)

def main(argv=sys.argv[1:]):
    
    for event in range(0,9):
        make_plot(event=event)

if __name__ == "__main__":
     exit(main())