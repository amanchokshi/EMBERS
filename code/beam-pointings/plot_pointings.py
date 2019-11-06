import json
import argparse
import numpy as np


parser = argparse.ArgumentParser(description="""
        Reads ultimate_pointing_times.json, and plots
        a histogram of total integration time at various
        pointings.
        """)

parser.add_argument('--meta_dir', metavar='\b', default='./../../outputs/beam-pointings/', help='Directory where json metadata files live. Default=./../../outputs/beam-pointings/')
parser.add_argument('--f_name', metavar='\b', default='ultimate_pointing_times.json', help='File name of json to be plotted. Default=ultimate_pointing_times.json')

args = parser.parse_args()
meta_dir = args.meta_dir
f_name = args.f_name

with open('{}{}'.format(meta_dir, f_name), 'r') as data:
    pointings = json.load(data)
    grid_pt = pointings['grid_pt']
    obs_length = pointings['obs_length']

# find unique pointings. 
# Replace None with -1, so sorting is possible.
# Then change -1 to None, again
unique_pointings = list(set(grid_pt))
unique_pointings[unique_pointings.index(None)] = -1
unique_pointings = sorted(unique_pointings)
unique_pointings[0] = None


pointings = []
integrations = []

for i in unique_pointings:
    pointings.append(i)
    time = 0
    for j in range(len(grid_pt)):
        if i == grid_pt[j]:
            time = time + obs_length[j]
    integrations.append(time)       
            
        


for i in range(len(pointings)):
    print(pointings[i], integrations[i])
            

