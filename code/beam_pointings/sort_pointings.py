import json
import numpy as np
from pathlib import Path


def org_pointing_json(meta_dir):
    """Organize json files. Clean up quirks in data such 
    and Andrew's custom pointing names"""

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
                
                # Many obs are named None, but their title indicates their true pointing
                # These are mostly obs which Andrew Williams schduled for this experiment
                if name in point_0:
                    pointing = 0
                elif name in point_2:
                    pointing = 2
                elif name in point_4:
                    pointing = 4
                else:
                    pass
                
                ## Filter obs by passes
                #if pointing in [0, 2, 4]:
                if pointing in list(range(197)):
                
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
                        obs_length.append(length)
                        pointings.append(pointing)
    
    return (start_gps, stop_gps, obs_length, pointings)
    
    
def combine_obs(start_gps, stop_gps, obs_length, pointings, out_dir):
    # if consecutive obs have same pointing, combine them.
    for i in range(len(start_gps)):
        if i != len(start_gps)-1:
            if ((stop_gps[i] == start_gps[i+1]) and (pointings[i] == pointings[i+1])):
                start_gps[i+1] = start_gps[i]
                obs_length[i+1] += obs_length[i]
                
                pointings[i] = None
    
    # Apply mask, couldn't think of a better way to implement this
    good_idx    = np.where(np.asarray(pointings) != None)[0]
    start_gps   = np.asarray(start_gps)[good_idx].tolist()
    stop_gps    = np.asarray(stop_gps)[good_idx].tolist()
    obs_length  = np.asarray(obs_length)[good_idx].tolist()
    pointings   = np.asarray(pointings)[good_idx].tolist()
    
    
    # Create dictionary to be saved to json, for plotting
    pointing_list = {}
    pointing_list['grid_pt'] = pointings
    pointing_list['start_gps'] = start_gps
    pointing_list['stop_gps'] = stop_gps
    pointing_list['obs_length'] = obs_length


    with open(f'{out_dir}/ultimate_pointing_times.json', 'w') as outfile:
        json.dump(pointing_list, outfile, indent=4)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="""
            Reads data from downloaded pointing matadata
            json files. Extracts the start, stop time and
            grid-pointing number.
            """)
    
    parser.add_argument(
            '--meta_dir', metavar='\b', default='./../../outputs/beam_pointings/mwa_pointings/', 
            help='Directory where json metadata files live. Default=./../../outputs/beam_pointings/mwa_pointings/')
    
    parser.add_argument(
            '--out_dir', metavar='\b', default='./../../outputs/beam_pointings/', 
            help='Directory where json metadata files live. Default=./../../outputs/beam_pointings/')
    
    args = parser.parse_args()
    meta_dir    = Path(args.meta_dir)
    out_dir     = Path(args.out_dir)
    
    start_gps, stop_gps, obs_length, pointings = org_pointing_json(meta_dir)
    combine_obs(start_gps, stop_gps, obs_length, pointings, out_dir) 
