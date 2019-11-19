# A dictionary of satellites which emit between 137.150 and 138.550 MHz
# ORBCOMM Communication satellites, NOAA & METEOR Weather satellites
# Satellite Name : NORAD CATALOGUE ID

norad_ids = {
        'ORBCOMM-X': 21576,
        'ORBCOMM FM 1': 23545,
        'ORBCOMM FM 2': 23546,
        'ORBCOMM FM 8': 25112,
        'ORBCOMM FM 10': 25113,
        'ORBCOMM FM 11': 25114,
        'ORBCOMM FM 12': 25115,
        'ORBCOMM FM 9': 25116,
        'ORBCOMM FM 5': 25117,
        'ORBCOMM FM 6': 25118,
        'ORBCOMM FM 7': 25119,
        'ORBCOMM FM 3': 25158,
        'ORBCOMM FM 4': 25159,
        'ORBCOMM FM 17': 25413,
        'ORBCOMM FM 18': 25414,
        'ORBCOMM FM 19': 25415,
        'ORBCOMM FM 20': 25416,
        'ORBCOMM FM 16': 25417,
        'ORBCOMM FM 15': 25418,
        'ORBCOMM FM 14': 25419,
        'ORBCOMM FM 13': 25420,
        'ORBCOMM FM 21': 25475,
        'ORBCOMM FM 22': 25476,
        'ORBCOMM FM 23': 25477,
        'ORBCOMM FM 24': 25478,
        'ORBCOMM FM 25': 25479,
        'ORBCOMM FM 26': 25480,
        'ORBCOMM FM 27': 25481,
        'ORBCOMM FM 28': 25482,
        'ORBCOMM FM 30': 25980,
        'ORBCOMM FM 31': 25981,
        'ORBCOMM FM 32': 25982,
        'ORBCOMM FM 33': 25983,
        'ORBCOMM FM 36': 25984,
        'ORBCOMM FM 35': 25985,
        'ORBCOMM FM 34': 25986,
            'ORBCOMM FM 38': 33060,
            'ORBCOMM FM 41': 33061,
            'ORBCOMM FM 29': 33062,
            'ORBCOMM FM 39': 33063,
            'ORBCOMM FM 37': 33064,
            'ORBCOMM FM 40': 33065,
            'ORBCOMM FM 109': 40086,
            'ORBCOMM FM 107': 40087,
            'ORBCOMM FM 106': 40088,
            'ORBCOMM FM 111': 40089,
            'ORBCOMM FM 104': 40090,
            'ORBCOMM FM 103': 40091,
            'ORBCOMM FM 114': 41179,
            'ORBCOMM FM 119': 41180,
            'ORBCOMM FM 105': 41181,
            'ORBCOMM FM 110': 41182,
            'ORBCOMM FM 118': 41183,
            'ORBCOMM FM 112': 41184,
            'ORBCOMM FM 113': 41185,
            'ORBCOMM FM 115': 41186,
            'ORBCOMM FM 108': 41187,
            'ORBCOMM FM 117': 41188,
            'ORBCOMM FM 116': 41189,
            'ORBCOMM FM 16 DEB 1': 44018,
            'ORBCOMM FM 16 DEB 2': 44019,
            'ORBCOMM FM 16 DEB 3': 44020,
            'ORBCOMM FM 16 DEB 4': 44022,
            'ORBCOMM FM 16 DEB 5': 44023,
            'ORBCOMM FM 16 DEB 6': 44024,
            # 'ORBCOMM FM 16 DEB 7': 44025, [NO CURRENT TLEs FOUND!]
            'ORBCOMM FM 16 DEB 8': 44026,
            'ORBCOMM FM 16 DEB 9': 44027,
            'ORBCOMM FM 16 DEB 10': 44028,
            'NOAA 18': 28654,
            'NOAA 15': 25338,
            'Meteor M2': 40069,
            'Meteor M2-2': 44387
           }
    
if __name__ == '__main__':
    for ids in norad_ids.values():
        print(ids)