import pod5 as p5
import json
import argparse
from pathlib import Path
import os
import sys

# creates index/dict containg read + corresponding pod5 file

def parse_args(argv):
    """Read arguments from command line."""

    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--pod5", required=True)
    parser.add_argument("-o", "--output", required=True)

    args = parser.parse_args()

    return args


##### copied from https://github.com/WGLab/DeepMod2/blob/main/plot_utils/plot.py, Copyright (c) 2022 Wang Genomics Lab
def get_file_names(base_path):
    read_filename_dict={}
    
    if os.path.isdir(base_path):
        files=Path(base_path).rglob('*.pod5')
    else:
        files=[base_path]

    for read_path in files:
        read_path=str(read_path)
        with p5.Reader(read_path) as reader:
            for rname in reader.read_ids:
                read_filename_dict[rname]=read_path
                
    return read_filename_dict
                
try:
    argv=sys.argv[1:]
    args = parse_args(argv)
    pod5_path_dict = get_file_names(base_path=args.pod5)
    json_name = args.output
except:
    raise Exception("pod5 file not found in given path")

with open(f'resources/pod5/index/p2s/{json_name}.json', 'w') as fp:
    json.dump(pod5_path_dict, fp, sort_keys=True, separators=[",\n",":"], allow_nan=False, default=str)

