import json
import pod5 as p5
import time
import os

read_ID = "b0d44b01-953e-4020-9c91-350f0e16e138"

pod5_path = "resources/pod5/p2s/"

time_st = time.process_time()
with open('resources/results/p2s/pod5.json', "r") as f:
    pod5_index= json.load(f)

pod5_file = pod5_index[read_ID]

time_index_load = time.process_time() - time_st

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


time_st = time.process_time()
for filename in os.listdir(pod5_path): #loops through all pod5 files in folder 
    pod5_file = os.path.join(pod5_path, filename)
    with p5.Reader(pod5_file) as pod5:
            # Read the selected read from the pod5 file
            # next() is required here as Reader.reads() returns a Generator
            try:
                pod5_record = next(pod5.reads(selection=[read_ID])) 
                idk = (pod5_record.signal)
                print(pod5_record.read_id, idk)
            except:
                continue

time_loop = time.process_time() - time_st

print (f"time index: {time_index_load}, {time_index} \n time loop: {time_loop}")