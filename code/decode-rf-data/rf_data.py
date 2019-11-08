import numpy as np
from datetime import datetime
from astropy.time import Time, TimezoneInfo
import matplotlib.pyplot as plt


def read_data(filename):
    '''Converts raw rf data from the RF Explorer
    to a power and time array. The raw data is 
    saved in binary format. Time is in UNIX format.
    '''
    
    with open(filename, 'rb') as f:
        next(f)
        lines = f.readlines()
        
        times = []
        data_lines = []

        for line in lines:
            time, data = line.split('$Sp'.encode())
            times.append(time.decode())
            
            # List converts bytes to list of bytes
            # The last two charachters are excluded - Newline char
            data_lines.append(list(data[:-2]))
        
        # I have no idea where the (-1/2) factor comes from
        # Remenant of Jared Rasti's code. Seems to work
        # Must be something to do with how data is written from RF
        power = np.asarray(data_lines) * (-1 / 2)
        
        return (power, times)

def plot_waterfall(power, times, n_axes=True):
    
    power_median = np.median(power)
    image = power - power_median
    vmin = 0
    vmax = 30
    
    plt.style.use('dark_background')
    fig = plt.figure(figsize = (7,10))
    ax = fig.add_axes([0.12, 0.1, 0.72, 0.85])
    #ax = fig.add_subplot(111)
    im = ax.imshow(image, vmin=vmin, vmax=vmax, interpolation='none', cmap='Spectral')
    ax.set_aspect('auto')
    cax = fig.add_axes([0.88, 0.1, 0.03, 0.85])
    fig.colorbar(im, cax=cax)
    if n_axes == True:
        times_iso = []
        for i in range(len(times)):
            times_iso.append(Time(float(times[i]), format='unix').iso)
        print(times_iso)    
    else:
        ax.set_xlabel('Freq Channel')
        ax.set_ylabel('Time Step')
    
    plt.show()


power, times = read_data('./../../sample-data/rf0XX_2019-10-10-15:00.txt')
plot_waterfall(power, times)    



