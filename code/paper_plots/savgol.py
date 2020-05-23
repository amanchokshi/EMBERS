import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('../align_data')
from align_data import savgol_interp

ch =8
ref_t, ref_p, tile_t, tile_p, ref_p_aligned, tile_p_aligned, time_array = savgol_interp(
        '../../../tiles_data/rf0XX/2019-09-15/rf0XX_2019-09-15-11:00.txt',
        '../../../tiles_data/S06XX/2019-09-15/S06XX_2019-09-15-11:00.txt',
        savgol_window_1 =11,
        savgol_window_2 =15,
        polyorder=2,
        interp_type='cubic',
        interp_freq=1
        )

plt.style.use('seaborn')

nice_fonts = {
        # Use LaTeX to write all text
        #"text.usetex": True,
        "font.family": "sans-serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 10,
        "font.size": 10,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 6,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        }

plt.rcParams.update(nice_fonts)


fig = plt.figure(figsize=(3.6,2.4))

tile_t = tile_t - tile_t[0]
time_array = time_array - time_array[0]

med = np.median(tile_p)
tile_p =tile_p - med
tile_p_aligned = tile_p_aligned - med

plt.plot(time_array,tile_p_aligned[::, ch], linewidth=1, color='#2c5d63', alpha=0.9, label='SavGol')
plt.scatter(tile_t, tile_p[::, ch], color='#7fa998',marker='.', s=2, alpha=0.6, label='AUT raw')

leg = plt.legend(loc="upper right", frameon=True, markerscale=4, handlelength=1)
leg.get_frame().set_facecolor('white')
for l in leg.legendHandles:
    l.set_alpha(1)

plt.ylabel('Power [dB]')
plt.xlabel('Time [s]')
plt.tight_layout()
#plt.savefig(f'savgol.png', dpi=600, bbox_inches='tight')
plt.savefig(f'./../../outputs/paper_plots/savgol.pdf', bbox_inches='tight')
