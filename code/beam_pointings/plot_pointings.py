import json
import matplotlib
import numpy as np
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
    


def point_integration(f_name, out_dir, int_thresh):
    """Bin data to calculate integration at each pointing"""
    
    with open(f'{out_dir}/{f_name}', 'r') as data:
        pointings = json.load(data)
        grid_pt = pointings['grid_pt']
        obs_length = pointings['obs_length']
    
    # find unique pointings. 
    unique_pointings = list(set(grid_pt))
    unique_pointings = sorted(unique_pointings)
    
    
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
    
    time_point = []
    point = []
    
    # time threshold at pointing in hours
    point_threshold = int_thresh
    
    for i in range(len(int_hours)):
        if int_hours[i] >= point_threshold:
            point.append(pointings[i])
            time_point.append(int_hours[i])

    return (time_point, point)


def pointing_hist(time_point, point, out_dir):
    """Pointing histogram"""

    x = range(len(time_point))
    leg = [int(i) for i in time_point]
    
    plt.style.use('seaborn')
    fig, ax = plt.subplots(figsize=(8,6))
    pal = sns.cubehelix_palette(len(time_point), start=.4, rot=-.5, dark=.4, reverse=True)
    barplot = plt.bar(x, time_point, color=sns.color_palette(pal))
    
    def autolabel(rects):
        for idx,rect in enumerate(barplot):
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., height,
                    leg[idx],
                    ha='center', va='bottom', rotation=0)
    
    autolabel(barplot)
    
    plt.xticks(x, point)
    plt.ylabel('Hours')
    plt.xlim(-0.7,len(time_point)-0.3)
    plt.xlabel('MWA Grid Pointing Number')
    plt.title('Integration at MWA Grid Pointings')
    plt.tight_layout()
    plt.savefig(f'{out_dir}/pointing_integration.png')


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="""
            Reads ultimate_pointing_times.json, and plots
            a histogram of total integration time at various
            pointings.
            """)
    
    parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/beam_pointings/', help='Directory where json metadata lives. Default=./../../outputs/beam_pointings/')
    parser.add_argument('--f_name', metavar='\b', default='ultimate_pointing_times.json', help='File name of json to be plotted. Default=ultimate_pointing_times.json')
    parser.add_argument('--int_thresh', metavar='\b', default=200, type=int,  help='Pointing integration time threshold. Default=200')
    
    args = parser.parse_args()
    out_dir     = args.out_dir
    f_name      = args.f_name
    int_thresh  = args.int_thresh

    time_point, point = point_integration(f_name, out_dir, int_thresh)
    pointing_hist(time_point, point, out_dir)

