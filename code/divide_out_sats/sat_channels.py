#TODO Import a dataset
#TODO Maybe use reference data to find channels. It's much cleaner
#TODO Find the median of the data. This will give us an approximate of the noise floor
#TODO Loop through freq channels and find any data above the median - potential satellite
#TODO Read chrono_ephem json file for that particular obs. Sort passes by lenght(time in sky)
#TODO If the potential sat occupies more than 80% of the sat pass length, classify it as a sat!
#TODO Or, come up with alternative thresholding scheme.
#TODO Exclude that channel from next loop
#TODO Find a way to include more channels - weather sats
#TODO Use reference data for sat, with corresponding chrono json for ephem

    
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../decode_rf_data')
from rf_data import read_data, plot_waterfall


ref_file = '../../../tiles_data/rf0XX/2019-10-07/rf0XX_2019-10-07-00:00.txt'

power, times = read_data(ref_file)

# To first order, let us consider the median to be the noise floor
noise_f = np.median(power)
print(noise_f)
print(power)


plot_waterfall(power, times, 'rf0XX')
plt.show()
    
