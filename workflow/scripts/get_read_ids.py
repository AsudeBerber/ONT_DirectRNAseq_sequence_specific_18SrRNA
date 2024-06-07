__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

import argparse
import os
import pod5 as p5
import numpy as np
import pysam as ps

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
