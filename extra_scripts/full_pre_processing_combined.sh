#!/usr/bin/env bash


python match_to_target_combined_df.py -d first_plate -t first_plate_targets.csv -i A01-B23 -p 50 

python filter_data.py -i "A01-B23_matched_new.csv"