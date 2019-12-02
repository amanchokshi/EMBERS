import json
import argparse
import numpy as np
import matplotlib
# Force matplotlib to not use X-Server backend
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
    

parser = argparse.ArgumentParser(description="""
        Reads ultimate_pointing_times.json, and plots
        a histogram of total integration time at various
        pointings.
        """)

parser.add_argument('--meta_dir', metavar='\b', default='./../../outputs/beam_pointings/', help='Directory where json metadata files live. Default=./../../outputs/beam_pointings/')
parser.add_argument('--f_name', metavar='\b', default='ultimate_pointing_times.json', help='File name of json to be plotted. Default=ultimate_pointing_times.json')
parser.add_argument('--threshold', metavar='\b', default=15, help='Plot if integration at pointing is > threshold. Default=15')
parser.add_argument('--plt_name', metavar='\b', default='pointing_integration', help='Name of plot to be save, default=pointing_integration')

args = parser.parse_args()
meta_dir = args.meta_dir
f_name = args.f_name
threshold = args.threshold
plt_name = args.plt_name

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
    if int_hours[i] > threshold:
        point.append(pointings[i])
        time_point.append(int_hours[i])

x = range(len(time_point))
leg = [int(i) for i in time_point]

plt.style.use('seaborn')
fig, ax = plt.subplots(figsize=(12,6))
barplot = plt.bar(x, time_point, color=sns.color_palette("rocket", len(time_point)))

def autolabel(rects):
    for idx,rect in enumerate(barplot):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., height,
                leg[idx],
                ha='center', va='bottom', rotation=0)

autolabel(barplot)

if point[0] == None:
    point[0] = 'X'
plt.xticks(x, point,  rotation='vertical')
plt.ylabel('Hours')
plt.xlim(-0.7,len(time_point)-0.3)
plt.xlabel('MWA Grid Pointing Number')
plt.title('Integration at MWA Grid Pointings [Threshold: {} Hours]'.format(threshold))
plt.tight_layout()
plt.savefig('{}/{}.png'.format(meta_dir,plt_name))


