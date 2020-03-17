import json
import argparse
import numpy as np
from pathlib import Path

parser = argparse.ArgumentParser(description="""
        Reads data from downloaded pointing matadata
        json files. Extracts the start, stop time and
        grid-pointing number.
        """)

parser.add_argument('--meta_dir', metavar='\b', default='./../../outputs/beam_pointings/', help='Directory where json metadata files live. Default=./../../outputs/beam_pointings/')

args = parser.parse_args()
meta_dir = Path(args.meta_dir)

start_gps = []
stop_stop = []
obs_name = []
pointings = []

sweet_points = [0, 2, 4, None]

files = meta_dir.glob('pointings*.json')

for f in files:
    with open(f) as table:
        data = json.load(table)
        for i in range(len(data)):
            name = data[i][2]
            pointing = data[i][-1]
            start = data[i][0]
            if i != (len(data)-1):
                # Telescope pointing remains constant till it is changed
                # So, stop time is the start of next observation
                stop = data[i+1][1]
            else:
                stop = data[i][1]
            
            if pointing in sweet_points:
                print(f'{data[i][0]}: {data[i][1]}: {data[i][2]}: {data[i][-1]}')
