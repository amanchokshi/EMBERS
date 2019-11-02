#! /bin/bash

while read p; do
  echo "python sat_ephemeris.py --sat=$p" >> parallel_ephem.txt
done <sat_list.txt

# 8 Jobs in Parallel
parallel --jobs 7 < parallel_ephem.txt 

rm parallel_ephem.txt
