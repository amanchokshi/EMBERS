import time
import numpy as np
from colormap import spectral
import matplotlib
# Force matplotlib to not use X-Server backend
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_data(filename):
    '''Converts raw rf data from the RF Explorer
    to a power and time array. The raw data is 
    saved in binary format. Time is in UNIX format.
    '''
   
    try:
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
            # The (-1/2) converts an unsigned byte to a real value
            power = np.single(np.asarray(data_lines) * (-1 / 2))
            times = np.double(np.asarray(times))
            
            return (power, times)
    except Exception:
        print(f'Missing {filename}, moving on...')


def tile_names():
    '''A list of tile reference and MWA tile names '''
    tiles = [
            'rf0XX', 'rf0YY', 'rf1XX', 'rf1YY',
            'S06XX', 'S06YY', 'S07XX', 'S07YY',
            'S08XX', 'S08YY', 'S09XX', 'S09YY',
            'S10XX', 'S10YY', 'S12XX', 'S12YY',
            'S29XX', 'S29YY', 'S30XX', 'S30YY',
            'S31XX', 'S31YY', 'S32XX', 'S32YY',
            'S33XX', 'S33YY', 'S34XX', 'S34YY',
            'S35XX', 'S35YY', 'S36XX', 'S36YY'
            ]
    return tiles


def plot_waterfall(power, times, name, n_axes=True):
    '''Using power and times array, creates a waterfall
    plot depicting satellite passes in the raw rf data.
    By default this function will use n_axes = True to
    convert unix times to human readable HH:MM local time
    for the y-axis, and use the rf freqencies for the x-azis.
    
    If n_axes=False, will plot waterfall with time step
    vs freq channel (0-112)
    '''

    power_median = np.median(power)
    image = power - power_median
    # setting dynamic range of waterfall to be 30 dB above the median
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
    ax.set_title(f'Waterfall Plot: {name}')

    if n_axes == True:
    
        # Number of time steps on y-axis
        number_t = 5
        t_step = int(len(times)/(number_t-1))
        times = list(times)
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
    
    return plt


