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

for i in range(len(start_gps)):
    print(f'{pointings[i]}: {start_gps[i]}: {stop_gps[i]}: {obs_length[i]}')


#begin_gps = []
#end_gps = []
#obs_length = []
#obs_pointing = []
#
#for j in range(len(start_gps)):
#    if j != (len(start_gps)-1):
#        if ((stop_gps[j] != start_gps[j+1]) and (pointings[j] == pointings[j+1])):
#            begin_gps.append(start_gps[j])
#            end_gps.append(start_gps[j])
#            obs_length.append(stop_gps[j]-start_gps[j])
#            obs_pointing.append(pointings[j])



## if the stop time of one obs, and the start time of the next obs are same, they are consevutive.
## if they also have the pointing, let us combine them. To do this, create new start, stop list, 
## with corresponding pointing lists
#
## remove first element
#start_1 = start_gps[1:]
#p_1 = pointings[1:]
#
## remove last element
#stop_0 = stop_gps[:-1]
#p_0 = pointings[:-1]
#
##consecutive_obs = np.where(np.asarray(start_1) == np.asarray(stop_0))
##consecutive_poi =  np.where(np.asarray(p_1) == np.asarray(p_0)) 
##
##print(len(consecutive_obs[0]))
##print(len(consecutive_poi[0]))
#
#good_idx = np.where((np.asarray(start_1) == np.asarray(stop_0)) & (np.asarray(p_1) == np.asarray(p_0)))[0]
#
#
##for i in good_idx:
##    print(i)
#
#
#def ranges(nums):
#    nums = sorted(set(nums))
#    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s+1 < e]
#    edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
#    return list(zip(edges, edges))
#
#cont_ranges = ranges(good_idx)
#
#
#for c in cont_ranges:
#    if c[0] != c[1]:
#        t1 = start_gps[c[0]]
#        t2 = stop_gps[c[1]+1]
#        del_t = t2 - t1
#        p = pointings[c[0]]
#        begin_gps.append(t1)
#        end_gps.append(t2)
#        obs_length.append(del_t)
#        obs_pointing.append(p)
#    else:
#        t1 = start_gps[c[0]]
#        t2 = stop_gps[c[1]]
#        del_t = t2 - t1
#        p = pointings[c[0]]
#        begin_gps.append(t1)
#        end_gps.append(t2)
#        obs_length.append(del_t)
#        obs_pointing.append(p)
#
#begin_gps, end_gps, obs_pointing, obs_length = zip(*sorted(zip(begin_gps, end_gps, obs_pointing, obs_length )))
#
#for i in range(len(begin_gps)):
#    print(f'{obs_pointing[i]}: {begin_gps[i]}: {end_gps[i]}: {obs_length[i]}')
