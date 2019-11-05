import os
import json
import argparse

parser = argparse.ArgumentParser(description="""
        Collates satellite pass data from all ephem json files.
        Sorts the data based on rise time and saves it to a new
        file - ultimate_ephem_list.json
        """)

parser.add_argument('--json_dir', metavar='\b', default='./../../outputs/sat-ephemeris/ephem_json/', help='Directory where ephem json files live. Default=./../../outputs/sat-ephemeris/ephem_json/')

args = parser.parse_args()
json_dir = args.json_dir

sat_id = []
t_rise = []
t_set = []
alt = []
az = []

for file in os.listdir(json_dir):
    if file.endswith('.json'):
        f_path = os.path.join(json_dir, file)

        with open(f_path) as ephem:
            sat_ephem = json.load(ephem)
            az.extend(sat_ephem['sat_az'])
            alt.extend(sat_ephem['sat_alt'])
            t_set.extend(sat_ephem['t_set'])
            t_rise.extend(sat_ephem['t_rise'])
            sat_id.extend([sat_ephem['sat_id'][0] for i in range(len(sat_ephem['t_rise']))])


# sort all the lists based on t_rise
t_rise, t_set, alt, az, sat_id = zip(*sorted(zip(t_rise, t_set, alt, az, sat_id)))

chrono_ephem = {}
chrono_ephem['sat_id'] = list(sat_id)
chrono_ephem['t_rise'] = list(t_rise)
chrono_ephem['t_set'] = list(t_set)
chrono_ephem['sat_alt'] = list(alt)
chrono_ephem['sat_az'] = list(az)


with open('{}/ultimate_ephem_list.json'.format(json_dir), 'w') as outfile:
    json.dump(chrono_ephem, outfile)


