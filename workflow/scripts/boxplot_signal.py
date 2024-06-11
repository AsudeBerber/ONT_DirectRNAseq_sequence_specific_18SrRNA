import argparse
import pod5 as p5
import matplotlib.pyplot as plt
import numpy as np
import pysam as ps

npz_file = f"resources/results/p2s/{motif}_window_{window_size}"

loaded = np.load(npz_file)