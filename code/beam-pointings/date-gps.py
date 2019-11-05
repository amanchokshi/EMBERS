import argparse
from astropy.time import Time

parser = argparse.ArgumentParser(description="""
        Convert date-time [iso] to gps time
        """)

parser.add_argument('--date', metavar='\b', help='Date-time to be converted to gps-seconds. Ex:2019-01-01 00:00:00')

args = parser.parse_args()
date = args.date

print(Time(date, format='iso').gps)


