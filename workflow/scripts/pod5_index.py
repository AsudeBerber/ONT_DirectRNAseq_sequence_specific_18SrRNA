import pod5 as p5
import json
import os
import pdb

pod5_file = "resources/pod5/p2s/"


##### copied from https://github.com/WGLab/DeepMod2/blob/main/plot_utils/plot.py, Copyright (c) 2022 Wang Genomics Lab
def get_file_names(base_path):
    read_filename_dict={}
    
    if os.path.isdir(base_path):
        files=os.path(base_path).rglob('*.pod5')
    else:
        files=[base_path]
    for read_path in files:
        read_path=str(read_path)
        with p5.Reader(read_path) as reader:
            for rname in reader.read_ids:
                read_filename_dict[rname]=read_path
                
    return read_filename_dict

try:
    pod5_path_dict = get_file_names(base_path=pod5_file)
except:
    breakpoint()
    pass

# with open('resources/results/p2s/pod5.json', "r") as f:
#     test= json.load(f)

with open('resources/results/p2s/pod5.json', 'w') as fp:
    json.dump(pod5_path_dict, fp, sort_keys=True, separators=[",\n",":"], allow_nan=False, default=str)

# test = {"ID3": "jkalsdjfk", "ID5": "333", "ID1":"ksdfjlj"}
