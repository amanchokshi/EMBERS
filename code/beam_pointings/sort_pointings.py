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

start_gps   = []
stop_gps    = []
obs_length  = []
pointings   = []

point_0 = ['All_0', 'Zenith_Test']
point_2 = ['EOR_Point_2', 'EOR_Point_2_Delays0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3_Ch100']
point_4 = ['EOR_Point_4', 'EOR_Point_4_Delays3,2,1,0,3,2,1,0,3,2,1,0,3,2,1,0_Ch100']

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
            
            # filter data by pointing
            if name in point_0:
                pointing = 0
            elif name in point_2:
                pointing = 2
            elif name in point_4:
                pointing = 4
            else:
                pass
            
            if pointing in [0, 2, 4]:
                start_gps.append(start)
                stop_gps.append(stop)
                obs_length.append(stop - start)
                pointings.append(pointing)
                #print(f'{stop - start}: {start}: {stop}: {name}: {pointing}')

for i in range(len(stop_gps)):
    print(f'{start_gps[i]}: {stop_gps[i]}: {pointings[i]}: {obs_length[i]}')
