import os
import json
import argparse

parser = argparse.ArgumentParser(description="""
        Collates satellite pass data from all ephem json files.
        Sorts the data based on rise time and saves it to a new
        file - ultimate_ephem_list.json
        """)

parser.add_argument('--json_dir', metavar='\b', default='./../../outputs/sat_ephemeris/ephem_json/', help='Directory where ephem json files live. Default=./../../outputs/sat_ephemeris/ephem_json/')

args = parser.parse_args()
json_dir = args.json_dir


# Open ultimate_ephem_list and extract lists from dict
with open(f'{json_dir}/ultimate_ephem_list.json', 'r') as ephem:
    sat_ephem = json.load(ephem)
    
    sat_id  =   sat_ephem['sat_id']
    t_rise  =   sat_ephem['t_rise']
    t_set   =   sat_ephem['t_set']
    sat_alt =   sat_ephem['sat_alt']
    sat_az  =   sat_ephem['sat_az']

print((t_set[0]-t_rise[0]) / 20)
print(len(sat_az[0]))
print((t_set[1]-t_rise[1])/20)
print(len(sat_az[1]))
print((t_set[2]-t_rise[2])/20)
print(len(sat_az[2]))


