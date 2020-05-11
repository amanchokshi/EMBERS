import json
import wget
import time
import argparse
import numpy as np
from pathlib import Path


parser = argparse.ArgumentParser(description="""
        Download metafits files to find
        dipole flagging of tiles
        """)

parser.add_argument(
        '--obs_ids', metavar='\b', default='./../../outputs/beam_pointings/ultimate_pointing_times.json', 
        help='Directory where json metadata files live. Default=./../../outputs/beam_pointings/ultimate_pointing_times.json')

parser.add_argument(
        '--out_dir', metavar='\b', default='./../../outputs/beam_pointings/metafits/', 
        help='Directory where json metadata files live. Default=./../../outputs/beam_pointings/metafits/')

args = parser.parse_args()
obs_ids     = args.obs_ids
out_dir     = Path(args.out_dir)

out_dir.mkdir(parents=True, exist_ok=True)

with open(obs_ids) as gps:
    gps_times = json.load(gps)['start_gps']

cerberus_url = 'http://ws.mwatelescope.org/metadata/fits?obs_id='

for i in gps_times[10::10]:
    print(f'Downloading {i} metafits')
    wget.download(f'{cerberus_url}{i}', f'{out_dir}/{i}.metafits')
    print('\n')
    time.sleep(10)   
    





