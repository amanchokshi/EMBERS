import json

pointing_list = {}
pointing_list['start_gps'] = []
pointing_list['stop_gps'] = []
pointing_list['pointing'] = []


with open('./../../outputs/beam-pointings/pointings_01.json', 'r') as table:
    data = json.load(table)
    for i in range(len(data)):
        pointing_list['start_gps'].append(data[i][0])
        pointing_list['stop_gps'].append(data[i][1])
        pointing_list['pointing'].append(data[i][-1])

print(pointing_list[0])        
