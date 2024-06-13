import json
import pod5 as p5
import time
import os

read_ID = "b0d44b01-953e-4020-9c91-350f0e16e138"

pod5_path = "resources/pod5/p2s/"

with open('resources/results/p2s/pod5.json', "r") as f:
    pod5_index= json.load(f)

pod5_file = pod5_index[read_ID]

time_st = time.process_time()
#loops through all pod5 files in folder 
with p5.Reader(os.path.join(pod5_path, pod5_file)) as pod5:
        # Read the selected read from the pod5 file
        # next() is required here as Reader.reads() returns a Generator
        try:
            pod5_record = next(pod5.reads(selection=[read.query_name])) 
            idk = (pod5_record.signal, read.get_tag("mv"), read.get_tag("ts"))
            print(pod5_record.read_id)
        except:
            continue
time_index = time.process_time() - time_st


print (f"time index: {time_index}")