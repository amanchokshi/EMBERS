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

time_array  = []
sat_id      = []
alt         = []
az          = []

for file in os.listdir(json_dir):
    if file.endswith('.json'):
        f_path = os.path.join(json_dir, file)

        with open(f_path) as ephem:
            sat_ephem = json.load(ephem)
            az.extend(sat_ephem['sat_az'])
            alt.extend(sat_ephem['sat_alt'])
            time_array.extend(sat_ephem['time_array'])
            sat_id.extend([sat_ephem['sat_id'][0] for i in range(len(sat_ephem['time_array']))])


# sort all the lists based on time_array
time_array, alt, az, sat_id = zip(*sorted(zip(time_array, alt, az, sat_id)))

chrono_ephem = {}
chrono_ephem['sat_id']      = list(sat_id)
chrono_ephem['time_array']  = list(time_array)
chrono_ephem['sat_alt']     = list(alt)
chrono_ephem['sat_az']      = list(az)


with open(f'{json_dir}/ultimate_ephem_list.json', 'w') as outfile:
    json.dump(chrono_ephem, outfile, indent=4)



