#! /bin/bash

python sat_ids.py >> sat_list.txt

while read p; do
  echo "python sat_ephemeris.py --sat=$p" >> parallel_ephem.txt
done <sat_list.txt

# 8 Jobs in Parallel
parallel --jobs 35 < parallel_ephem.txt

rm parallel_ephem.txt
rm sat_list.txt
