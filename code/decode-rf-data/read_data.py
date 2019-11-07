import numpy as np


def read_data(filename):
    '''Converts raw rf data from the RF Explorer
    to a power and time array. The raw data is 
    saved using unicode charachters.'''

    with open(filename, 'rb') as table:
        next(table)
        table = table.readlines()
        for line in table:
            #data = line.read()
            time, dt = line.split('$Sp'.encode())
            print(time.decode())
#        for line in table:
            #print(line)
#        lines = table.read().split('\n')
#        
#        times = []
#        data_lines = []
#
#
#        for line in table:
#            time, data = line.split($Sp)
#            times.append(time)
#            data_lines.append(data)
#
#        times = np.asarray(times).astype(np.float)
#        print(times)


read_data('./../../sample-data/S08XX_2019-10-10-11:00.txt')


