import sys
import unittest
import numpy as np
import healpy as hp
from pathlib import Path
import matplotlib.pyplot as plt
sys.path.append('../../code/reproject_ref')
sys.path.append('../../code/decode_rf_data')
import project_beams2healpix as pb
from colormap import spectral


class Test_project_beams2healpix(unittest.TestCase):
        
    nside = 32
    len_empty_healpix = hp.nside2npix(nside)   #12288
    epsilon = 1.e-12
    refXX = '../../data/FEE_Reference_Models/MWA_reference_tile_FarField_XX.ffe'
    healpix_XX, _, _, _, _ = pb.create_model(len_empty_healpix, epsilon, nside, file_name=refXX)

    def test_create_model(self):    
        self.assertEqual(self.healpix_XX.shape[0], 12288)
        self.assertTrue(np.isclose(0, self.healpix_XX[0], atol=1e-2))

    def test_plot_healpix(self):
        pb.plot_healpix(data_map=self.healpix_XX, fig=1, sub=None,title=None,vmin=None,vmax=None,cmap=None, cbar=True)
        plt.savefig('healpix.png')
        heal_png = Path('healpix.png')
        self.assertTrue(heal_png.is_file())

        if heal_png.is_file():
            heal_png.unlink()


if __name__ == "__main__":
    unittest.main()


