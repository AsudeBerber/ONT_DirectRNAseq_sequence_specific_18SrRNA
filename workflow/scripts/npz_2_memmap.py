__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

mean_q1, mean_q2, mean_q3, mean_q4 = features[:,:,0], features[:,:,1], features[:,:,2], features[:,:,3]
log10_len, mean, sd = features[:,:,4], features[:,:,5], features[:,:,6]
mdn, mdn_sd = features[:,:,7], features[:,:,8]

fp = np.memmap("../../resources/results/p2s/memmap/", mode = "w+")


fp.flush()