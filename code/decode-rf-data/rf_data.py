import time
import numpy as np
from colormap import spectral
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


def plot_waterfall(power, times, name, n_axes=True):
    
    power_median = np.median(power)
    image = power - power_median
    vmin = 0
    vmax = 30

    # Custom spectral colormap
    cmap = spectral()
    
    plt.style.use('dark_background')
    fig = plt.figure(figsize = (7,10))
    ax = fig.add_axes([0.12, 0.1, 0.72, 0.85])
    im = ax.imshow(image, vmin=vmin, vmax=vmax, interpolation='none', cmap=cmap)
    cax = fig.add_axes([0.88, 0.1, 0.03, 0.85])
    fig.colorbar(im, cax=cax)
    ax.set_aspect('auto')
    ax.set_title('Waterfall Plot: {}'.format(name))

    if n_axes == True:
    
        # Number of time steps on y-axis
        number_t = 5
        t_step = int(len(times)/(number_t-1))
        times = times[::t_step]
        
        t_tz = []
        
        # Convert UNIX time to local HH:MM time
        for i in range(len(times)):
            
            perth_t = float(times[i])+28800 #28800=+8GMT @ PERTH
            hms = time.strftime('%H:%M', time.gmtime(perth_t))
            t_tz.append(hms)

        # Frequency: x-axis
        start_freq = 137.15
        stop_freq = 138.55

        # X-axis stuff
        x_ax = image.shape[1]
        freqs = np.arange(start_freq, stop_freq, 0.25)
        x_ticks = np.arange(0, x_ax, (0.25/0.0125)) #.0125MHz/ch
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(freqs)
        ax.set_xlabel('Freq [MHz]')

        # Y-axis stuff
        y_ax = image.shape[0]
        y_ticks = np.arange(0, y_ax, t_step)
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(t_tz)
        ax.set_ylabel('Time [HH:MM]')

    else:
        ax.set_xlabel('Freq Channel')
        ax.set_ylabel('Time Step')
    
    plt.show()


power, times = read_data('./../../sample-data/rf0XX_2019-10-10-15:00.txt')
plot_waterfall(power, times, 'Test File Name', False)    




