import json
import pod5 as p5
import time
import pdb
import timeit
from sortedcontainers import SortedDict
import os

read_ID = "fffff88b-c3c0-4ed3-ae1a-afff6053366a"

pod5_path = "resources/pod5/p2s/"

time_st = time.process_time()
with open('resources/results/p2s/pod5.json', "r") as f:
    pod5_index= json.load(f)

pod5_index_sd = SortedDict(pod5_index)





time_st = time.time_ns()

pod5_file = pod5_index[read_ID]

time_index_unsrt = time.time_ns() - time_st


time_st = time.process_time_ns()

pod5_file = pod5_index_sd[read_ID]

time_index_srt = time.process_ns() - time_st


time_index_load = time.process_time() - time_st

breakpoint()
time_st = time.process_time()
#loops through all pod5 files in folder 
with p5.Reader(pod5_file) as pod5:
        # Read the selected read from the pod5 file
        # next() is required here as Reader.reads() returns a Generator
            try:
                pod5_record = next(pod5.reads(selection=[read_ID])) 
                idk = (pod5_record.signal)
                print(pod5_record.read_id, idk)
            except:
                breakpoint()
                pass
time_index = time.process_time() - time_st
print(time_index, time)


time_st = time.process_time()
for filename in os.listdir(pod5_path): #loops through all pod5 files in folder 
    pod5_file = os.path.join(pod5_path, filename)
    with p5.Reader(pod5_file) as pod5:
            read_ID = pod5.query_name

time_loop = time.process_time() - time_st

print (f"time index: {time_index_load}, {time_index} \n time loop: {time_loop}")