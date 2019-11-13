python sat_ids.py >> sat_list.txt
while read p; do python plot_ephemeris.py --sat=$p; done <sat_list.txt
rm sat_list.txt
