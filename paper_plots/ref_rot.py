import healpy as hp
import numpy as np
from embers.rf_tools.colormaps import jade, spectral
from embers.tile_maps.beam_utils import plot_healpix, rotate_map
from matplotlib import pyplot as plt

spec, _ = spectral()
jade, _ = jade()
nside = 32

ref_model = "../embers_out/tile_maps/ref_models/ref_dipole_models.npz"

# Load reference FEE model
# Rotate the fee models by -pi/2 to move model from spherical (E=0) to Alt/Az (N=0)
ref_fee_model = np.load(ref_model, allow_pickle=True)

ref_xx = rotated_fee = rotate_map(
    nside, angle=-(1 * np.pi) / 2.0, healpix_array=ref_fee_model["XX"]
)
ref_yy = rotated_fee = rotate_map(
    nside, angle=-(1 * np.pi) / 2.0, healpix_array=ref_fee_model["YY"]
)

plot_healpix(data_map=ref_xx, cmap=spec, vmin=-20, vmax=0)
plt.savefig("test_ref.png")
