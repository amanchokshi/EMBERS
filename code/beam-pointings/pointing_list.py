import os
import json
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="""
        Reads data from downloaded pointing matadata
        json files. Extracts the start, stop time and
        grid-pointing number.
        """)

parser.add_argument('--meta_dir', metavar='\b', default='./../../outputs/beam-pointings/', help='Directory where json metadata files live. Default=./../../outputs/beam-pointings/')

args = parser.parse_args()
meta_dir = args.meta_dir

start = []
stop = []
pointings = []


for file in os.listdir(meta_dir):
    if file.endswith('.json'):
        f_path = os.path.join(meta_dir, file)
        
        with open(f_path) as table:
            data = json.load(table)
            for i in range(len(data)):
                start.append(data[i][0])
                stop.append(data[i][1])
                pointings.append(data[i][-1])

# sort the lists based on start_gps, incase files weren't read in correct order.

start, stop, pointings = zip(*sorted(zip(start, stop, pointings)))

start_gps = list(start)
stop_gps = list(stop)
pointings = list(pointings)

# Every morning, tests are run on the recievers.
# There are always 20 test observations of ~72 seconds
# These tests are assigned gridpointing == None
# The last of these is actually a zenith obs. 
# I manually set the 20th obs from None to 0
# because the time following this test, until the next obs is 
# usually a Looonnnnng time, and the telesocope remains at zenith

# if list of 20 Nones, 20th reset to 0!

for i in range(len(pointings)):
    if (i+19) < len(pointings):
        if (pointings[i] == None and pointings[i+19] == None):
            pointings[i+19] = 0

obs_length = list(np.diff(np.asarray(start_gps)))
pointings = pointings[:-1]
end_gps = start_gps[1:]
start_gps = start_gps[:-1]

for i in range(len(pointings)):
    print(start_gps[i], end_gps[i], obs_length[i], pointings[i])

#pointing_list = {}
#pointing_list['start_gps'] = []
#pointing_list['stop_gps'] = []
#pointing_list['pointing'] = []

