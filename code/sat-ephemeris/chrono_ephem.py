import os
import json

sat_id = []
t_rise = []
t_set = []
alt = []
az = []



with open('./../../outputs/ephem_json/21576.json') as ephem:
    sat_ephem = json.load(ephem)
    az.extend(sat_ephem['sat_az'])
    alt.extend(sat_ephem['sat_alt'])
    t_set.extend(sat_ephem['t_set'])
    t_rise.extend(sat_ephem['t_rise'])
    sat_id.extend([sat_ephem['sat_id'][0] for i in range(len(t_rise))])
    
print(sat_id)    
print(t_set)

