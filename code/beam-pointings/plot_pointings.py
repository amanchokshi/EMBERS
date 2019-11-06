import json
import argparse
import numpy as np
#import matplotlib
## Force matplotlib to not use X-Server backend
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
    

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


int_hours = np.asarray(integrations)/(60*60)

# Only plot if integration for pointing is > threshold

time_point = []
point = []

for i in range(len(int_hours)):
    if int_hours[i] > 10:
        point.append(pointings[i])
        time_point.append(int_hours[i])

x = range(len(time_point))
leg = [int(i) for i in time_point]

fig, ax = plt.subplots(figsize=(12,7))
fig = plt.bar(x, time_point, color=sns.color_palette("rocket", len(time_point)))
if point[0] == None:
    point[0] = 'None'
plt.xticks(x, point,  rotation='vertical')
plt.ylabel('Hours')
plt.xlim(-1,len(time_point))
plt.xlabel('MWA Grid Poining Number')
#plt.legend(fig, leg, loc = "upper right", title = 'Hours',prop={'size':7})
plt.tight_layout()
plt.show()


