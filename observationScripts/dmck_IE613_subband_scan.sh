#!/bin/bash

# Generate the base path of the observation
basepath='/data/home/user1/data/'`date +"%Y"`'/'`date +"%m"`'/'`date +"%d"`'/XST/'

# Get the start time for the first level subfolder
start_timestrp="`date +"%H:%M:%S"`"

# Ensure the path for the log exists
mkdir -p $basepath

# Scan subband 51 to 461, sample every subband for 60 seconds, given we are on IE613
bash ./dmck_generic_subband_scan.sh 51 462 1 60 $start_timestrp 1 | tee -a $basepath$start_timestrp'/allSky.log' 

# Remove the temporary folder when we finish
rm -r $basepath'/temp_00/'