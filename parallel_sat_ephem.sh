#! /bin/bash

while read p; do
  echo "python sat_ephemeris.py --sat=$p" >> parallel_ephem.txt
done <sat_list.txt

# GNU Parallel. 8 Jobs at a time
parallel --jobs 8 < parallel_ephem.txt 

rm parallel_ephem.txt
