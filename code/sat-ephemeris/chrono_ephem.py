import os
import json
import numpy as np
from sat_ephemeris import sat_plot

sat_id = []
t_rise = []
t_set = []
alt = []
az = []

for file in os.listdir('./../../outputs/ephem_json'):
    if file.endswith('.json'):
        f_path = os.path.join('./../../outputs/ephem_json', file)

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


with open('./../../outputs/ephem_json/ultimate_ephem_list.json', 'w') as outfile:
    json.dump(chrono_ephem, outfile)


