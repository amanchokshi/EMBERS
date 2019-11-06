import os
import json
import argparse

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
                start.extend(data[i][0])
                stop.extend(data[i][1])
                pointings.extend(data[i][-1])

# sort the lists based on start_gps, incase files weren't read in correct order.

start, stop, pointing = zip(*sorted(zip(start, stop, pointings)))

print(pointing)       


#pointing_list = {}
#pointing_list['start_gps'] = []
#pointing_list['stop_gps'] = []
#pointing_list['pointing'] = []

