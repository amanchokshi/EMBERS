import numpy as np
import matplotlib.pyplot as plt


def read_data(filename):
    '''Converts raw rf data from the RF Explorer
    to a power and time array. The raw data is 
    saved using unicode charachters.'''
    
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
        
        fig = plt.figure(figsize = (7,7))
        ax = fig.add_subplot(111)
        ax.imshow(power)
        ax.set_aspect('auto')
        plt.show()


read_data('./../../sample-data/rf0XX_2019-10-10-15:00.txt')



