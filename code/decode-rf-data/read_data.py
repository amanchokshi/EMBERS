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
            #data_lines.append(data[:-2].decode('latin-1'))
            data_lines.append(list(data[:-2]))
        #print(list(data_lines[0]))
        #power_lines = []
        print(len(data_lines[0]))
        water = np.array(data_lines)
        print(water.shape)
        #water.shape = (len(times), len(data_lines))
        fig = plt.figure(figsize = (7,7))
        ax = fig.add_subplot(111)
        ax.imshow(water)
        ax.set_aspect('auto')
        plt.show()

        #for d in data_lines:
        #    powers = []
        #    for char in d[:-2]:
        #        powers.append(ord(char))
        #    power_lines.append(powers)    


        #times = np.asarray(times).astype(np.float)
        #print(power_lines)    


read_data('./../../sample-data/rf0XX_2019-10-10-15:00.txt')



