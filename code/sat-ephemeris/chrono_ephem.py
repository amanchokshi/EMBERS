import os
import json

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
            sat_id.extend([sat_ephem['sat_id'][0] for i in range(len(t_rise))])
    
print(len(sat_id))    
print(len(t_set))
print(len(az))

