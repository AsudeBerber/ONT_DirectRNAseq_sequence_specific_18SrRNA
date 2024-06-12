__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

motif = "CCG"
window_size = 21
npz_file = f"../../resources/results/p2s/{motif}_window_21_subsample_00001.npz"
arround=3

loaded = np.load(npz_file)


#arround: #bases to plot arround ac4C
def slice_bases(event, arround):
    whole_window = len(query[0])
    middle = whole_window //2
    index_bases = np.arange(0,whole_window) [middle-arround:middle+arround+1]
    sliced_event = features[:,index_bases,event]
    sliced_refseq = ref[:,index_bases]
    return index_bases, sliced_event, sliced_refseq
    

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
    pos.append(read[-4:])

mean_signal_int = features[:,:,5]

event = 5 # mean intensity
index_bases, sliced_event, sliced_ref_seq = slice_bases(event=event, arround=arround)

df_event = pd.DataFrame((sliced_event), columns = list(range(1,8)))
df_id = pd.DataFrame(id)
df_pos = pd.DataFrame(pos, columns = ["pos"])

df_event_pos = pd.concat([df_event, df_pos], axis=1)

def filter_by_pos(pos):
    df = df_event_pos
    df_filtered = df[df["pos"] == str(pos)]
    df_plot = df_filtered.iloc[:,:7]
    return df_plot

#1337 | 1843
pos = 1843
df_event_filtered = filter_by_pos(pos)


fig, (ax1, ax2) = plt.subplots(1,2)
ax1.violinplot(filter_by_pos(1337), showmeans = False, showextrema = False)
ax1.boxplot(filter_by_pos(1337))

ax2.violinplot(filter_by_pos(1843), showmeans = False, showextrema = False)
ax2.boxplot(filter_by_pos(1843))