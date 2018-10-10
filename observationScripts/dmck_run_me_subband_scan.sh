#!/bin/bash

# Generate the base path of the observation
basepath='/data/home/user1/data/'`date +"%Y"`'/'`date +"%m"`'/'`date +"%d"`'/XST/mode5/'

# Get the start time for the first level subfolder
start_timestrp="`date +"%H:%M:%S"`"

# Ensure the path for the log exists
mkdir -p $basepath
bash ./dmck_run_subband_scan.sh $1 $start_timestrp | tee -a $basepath$start_timestrp'/allSky.log' 

# Remove the temporary folder when we finish
rm -r $basepath'/temp_00/'