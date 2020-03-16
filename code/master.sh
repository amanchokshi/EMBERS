#! /bin/bash

## ONE SCRIPT TO RULE THEM ALL
## MUHAHAHAHA HA HA ...

# Batch process all data :)

# Plot waterfalls for one month of data
python batch_waterfall.py --data_dir=../../../tiles_data --start_date=2019-10-01 --stop_date=2019-11-01

# align all data
python batch_align.py --data_dir=./../../../tiles_data --start_date=2019-09-12 --stop_date=2020-03-16


