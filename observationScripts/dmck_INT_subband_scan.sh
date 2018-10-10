#!/bin/bash

# Generate the base path of the observation
basepath='/data/home/user1/data/'`date +"%Y"`'/'`date +"%m"`'/'`date +"%d"`'/XST/mode5/'

# Get the start time for the first level subfolder
start_timestrp="`date +"%H%M%S"`"

# Ensure the path for the log exists
mkdir -p $basepath

# Total integration time: ~28 minutes.
# Total observation time: depends on how quickly HBA activations are applied.

# Scan subband 51 to 120, sample every 5th subband for 60 seconds
bash ./dmck_generic_subband_scan.sh 51 121 5 60 $start_timestrp 1 | tee -a $basepath$start_timestrp'/allSky.log' 

# Scan subband 120 to 461, sample every 25th subband for 60 seconds
bash ./dmck_generic_subband_scan.sh 120 462 25 60 $start_timestrp 1 | tee -a $basepath$start_timestrp'/allSky.log' 

# Remove the temporary folder when we finish
rm -r $basepath'/temp_00/'