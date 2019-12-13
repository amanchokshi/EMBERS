import os
import json
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="""
        Reads data from downloaded pointing matadata
        json files. Extracts the start, stop time and
        grid-pointing number.
        """)

parser.add_argument('--meta_dir', metavar='\b', default='./../../outputs/beam_pointings/', help='Directory where json metadata files live. Default=./../../outputs/beam_pointings/')

args = parser.parse_args()
meta_dir = args.meta_dir

start = []
stop = []
obs_name = []
pointings = []


for file in os.listdir(meta_dir):
    if file.endswith('.json'):
        f_path = os.path.join(meta_dir, file)
        
        with open(f_path) as table:
            data = json.load(table)
            for i in range(len(data)):
                start.append(data[i][0])
                stop.append(data[i][1])
                obs_name.append(data[i][2])
                pointings.append(data[i][-1])

# sort the lists based on start_gps, incase files weren't read in correct order.
start, stop, pointings = zip(*sorted(zip(start, stop, pointings)))

start_gps = list(start)
stop_gps = list(stop)
obs_name = list(obs_name)
pointings = list(pointings)

# Andrew William's alternating sweetpoint 2,4 obs at 7:30 AM
# have gridpoint=None, but obs_name=EOR_Point_4 or EOR_Point_2
# Use this to modify the pointing list

for i in range(len(pointings)):
    if obs_name[i] == 'EOR_Point_2':
        pointings[i] = 2
    elif obs_name[i] == 'EOR_Point_4':
        pointings[i] = 4
    else:
        pass



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


# Find the index of pointings, where pointings change
p_change = [i for i in range(1,len(pointings)) if pointings[i]!=pointings[i-1]]

# Add 0 to the beginning and the lenght to the end, to have an inclusive list
p_change.insert(0, 0)
p_change.append(len(pointings))

# This section is a bit confusing.
# We have multiple consecutive observations at the same pointing.
# This code extracts the beginning of a pointing block == t_start
# The ending of pointing block == t_stop
# The total integration over the block == t_length
# The gridpointing number of the block
# all of these data points are added to lists, whose index corresponds to pointing blocks

t_start = []
t_stop = []
t_length = []
grid_pt = []
for i in range(len(p_change)):
    if i+1 < len(p_change):
        t_start.append(start_gps[p_change[i]])
        t_stop.append(stop_gps[p_change[i+1] - 1 ])
        grid_pt.append(pointings[p_change[i]])
        t_length.append(sum(obs_length[p_change[i]:p_change[i+1]]))

# Because of weird problem with int64 and json
t_length = (np.asarray(t_length)).tolist()

# Create dictionary to be saved to json
pointing_list = {}
pointing_list['grid_pt'] = grid_pt
pointing_list['start_gps'] = t_start
pointing_list['stop_gps'] = t_stop
pointing_list['obs_length'] = t_length


with open('{}/ultimate_pointing_times.json'.format(meta_dir), 'w') as outfile:
    json.dump(pointing_list, outfile, indent=4)
