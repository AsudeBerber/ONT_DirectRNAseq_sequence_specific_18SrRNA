import pod5
import json
import pdb

pod5_file = "resources/pod5/p2s/PAW35875_9fd38647_68d05f77_0.pod5"


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
    pod5_path_dict = get_file_names(base_path=base_path)
except:
    breakpoint()
    pass

# with open('resources/results/p2s/pod5.json', "r") as f:
#     test= json.load(f)

with open('resources/results/p2s/pod5.json', 'w') as fp:
    json.dump(pod5_path_dict, fp, sort_keys=True, separators=[",\n",":"], allow_nan=False, default=str)

# test = {"ID3": "jkalsdjfk", "ID5": "333", "ID1":"ksdfjlj"}
