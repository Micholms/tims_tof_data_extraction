#!/usr/bin/env bash


for dir in first_plate/Exported_data/*/     # list directories in the form "/tmp/dirname/"
do
    dir=${dir%*/}    # remove the trailing "/"
    dir=${dir##*/}   # print everything after the final "/"
    
    python match_to_target.py -t first_plate/first_plate_targets.csv  -d first_plate/Exported_data/ -w ${dir} -p 10
done
list=()
for i in *10.0_matched.csv
do
    
    list+=($i)
done

python merge_pre_processed.py ${list[@]}

python filter_data.py -i "merged_A01_A10_matchedonly.csv"