import json
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
        '--out_dir', metavar='\b', default='./../../outputs/metafits/', 
        help='Directory where json metadata files live. Default=./../../outputs/metafits/')

args = parser.parse_args()
obs_ids     = args.obs_ids
out_dir     = Path(args.out_dir)

out_dir.mkdir(parents=True, exist_ok=True)

with open(obs_ids) as gps:
    gps_times = json.load(gps)['start_gps']

for i in gps_times[::10]:
    print(i)





