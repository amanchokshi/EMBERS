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
            
                if i != (len(data)-1):
                # Telescope pointing remains constant till it is changed
                # So, stop time is the start of next observation
                    stop = data[i+1][0]
                else:
                    stop = data[i][1]

                length = stop - start

                if length >= 120:
                
                    start_gps.append(start)
                    stop_gps.append(stop)
                    obs_length.append(stop - start)
                    pointings.append(pointing)

## CAUTION
## Let's find where the pointings array changes values. Basically, I want to integrate all obs with consecutive pointings
## The code below works only if the first and last pointings aren't the same. If they are the same, the index of the first 
## element won't be returned
#p_changes = np.where(np.roll(pointings,1)!=pointings)[0]

#for i in range(len(start_gps)):
#    if i != (len(start_gps)-1):
#        if ((stop_gps[i] == start_gps[i+1]) and (pointings[i] == pointings[i+1])):

for i in range(len(start_gps)):
    print(f'{start_gps[i]}: {stop_gps[i]}: {pointings[i]}: {obs_length[i]}')

