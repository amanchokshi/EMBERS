import argparse
from astropy.time import Time

parser = argparse.ArgumentParser(description="""
        Convert date-time [iso] to gps time
        """)

def convert_gps(date):
    return int(Time(date, format='iso').gps)

parser.add_argument('--date', metavar='\b', help='Date-time to be converted to gps-seconds. Ex:2019-01-01 00:00:00')

args = parser.parse_args()
date = args.date

gps = convert_gps(date)
print(f'gps: {gps}')



