import sys
import numpy as np
import healpy as hp
from scipy.stats import median_absolute_deviation as mad


def good_maps(raw_tile):
    '''Here, we extract "good satellie" data from the tile_maps_raw 
    generated by tile_maps.py. We also normalize the maps by
    the reference beam'''

    f_name, _ = raw_tile.name.split('.')
    tile, ref, _, _ = f_name.split('_')
    
    # Load reference FEE model
    ref_fee_model = np.load(ref_model, allow_pickle=True)
    if 'XX' in tile:
        ref_fee = ref_fee_model['XX']
    else:
        ref_fee = ref_fee_model['YY']
   
    # load data from map .npz file
    tile_data   = np.load(raw_tile, allow_pickle=True)
    tile_data   = {key:tile_data[key].item() for key in tile_data}
    ref_map     = tile_data['ref_map'] 
    tile_map    = tile_data['tile_map']  
    
    tile_maps_norm = {p:[] for p in pointings}
     
    for p in pointings:

        # This contains the data from all the good sat ref data
        ref_map_good = [[] for pixel in range(hp.nside2npix(nside))]
        
        # This contains the data from all the good sat ref data
        tile_map_good = [[] for pixel in range(hp.nside2npix(nside))]
        
        for sat in good_sats:

            for pix in range(hp.nside2npix(nside)):

                ref_map_good[pix].extend(ref_map[p][sat][pix])
                tile_map_good[pix].extend(tile_map[p][sat][pix])
        
        # Here, we divide tile power by ref power and multiply by the fee ref model in log space
        tile_map_norm = [[] for pixel in range(hp.nside2npix(nside))]
        
        for pixel in range(hp.nside2npix(nside)):
        
            ratio = np.asarray(np.subtract(tile_map_good[pixel], ref_map_good[pixel]))
            tile_map_norm[pixel].extend(ratio + ref_fee[pixel])

        tile_maps_norm[p].extend(tile_map_norm)


    # Save map arrays to npz file
    np.savez_compressed(f'{out_dir}/{tile}_{ref}_tile_maps.npz', **tile_maps_norm)


if __name__=='__main__':
    
    import argparse
    from pathlib import Path
    import concurrent.futures

    parser = argparse.ArgumentParser(description="""
        Extract the data from a list of 'good satellites' from the raw tile maps made my tile_maps.py. 
        Divide the tile/ref beams and apply the FEE ref model.
        """)
    
    parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/tile_maps/tile_maps_norm/',
            help='Output directory. Default=./../../outputs/tile_maps/tile_maps_norm/')
    parser.add_argument('--map_dir', metavar='\b', default='./../../outputs/tile_maps/tile_maps_raw',
            help='Raw Map directory. Default=./../../outputs/tile_maps/tile_maps_raw')
    parser.add_argument('--ref_model', metavar='\b', default='../../outputs/reproject_ref/ref_dipole_models.npz',
            help='Healpix reference FEE model file. default=../../outputs/reproject_ref/ref_dipole_models.npz')
    parser.add_argument('--nside', metavar='\b', type=int,  default=32,help='Healpix Nside. Default = 32')
    
    args = parser.parse_args()
    
    out_dir     = Path(args.out_dir)
    map_dir     = Path(args.map_dir)
    ref_model   = Path(args.ref_model)
    nside       = args.nside
    
    # make outpud dir
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Good sats from which to make plots
    good_sats = [
            25338, 25984, 25985,
            28654, 40086, 40087,
            40091, 41179, 41180,
            41182, 41183, 41184,
            41185, 41187, 41188,
            41189, 44387
            ]
    
    # list of beam pointings
    pointings = ['0','2','4']
   
    # list of all raw tile maps
    map_files = [item for item in map_dir.glob('*.npz')]

    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(good_maps, map_files)

    #good_maps(map_files[0])






