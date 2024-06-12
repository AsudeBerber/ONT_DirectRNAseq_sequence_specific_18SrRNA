__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

motif = "CCG"
window_size = 21
npz_file = f"../../resources/results/p2s/{motif}_window_21.npz"
arround=3

loaded = np.load(npz_file)


#arround: #bases to plot arround ac4C
def slice_bases(event, arround):
    whole_window = len(query[0])
    middle = whole_window //2
    index_bases = np.arange(0,whole_window) [middle-arround:middle+arround+1]
    sliced_event = features[:,index_bases,event]
    return index_bases, sliced_event
    

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

event = 6 # mean intensity
index_bases, sliced_event = slice_bases(event=event, arround=arround)

# df = pd.DataFrame((sliced_event, pos), columns = list(range(1,8)) + ["pos"])
df = pd.DataFrame()


fig, ax = plt.subplots()
ax.violinplot(sliced_event, showmeans = False, showextrema = False)
ax.boxplot(sliced_event)