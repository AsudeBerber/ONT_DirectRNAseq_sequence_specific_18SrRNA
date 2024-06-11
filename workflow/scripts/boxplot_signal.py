import argparse
import pod5 as p5
import matplotlib.pyplot as plt
import numpy as np
import os
# import pysam as ps

motif = "CCG"
window_size = 21
npz_file = f"resources/results/p2s/{motif}_window_{window_size}.npz"

loaded = np.load(npz_file)

features = loaded["feat"]
qual = loaded["qual"]
query = loaded["query"]
ref = loaded["ref"]
id = loaded["id"]

os.path.dirname(npz_file)
