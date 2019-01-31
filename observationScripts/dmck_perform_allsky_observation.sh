#!/bin/bash

if test $# != 8; then
	echo "Incorrect inputs: ./dmck_perform_allsky_observation.sh <minutes> <duty cycle> <subband mode 3 start> <subband mode 3 end> <subband mode 5 start> <subband mode 5 end> <subband mode 7 start> <subband mode 7 end>"
	exit;
fi

# Generate the base path of the observation
basepath='/data/home/user1/data/'`date +"%Y"`'/'`date +"%m"`'/'`date +"%d"`'/xst/'

# Get the start time for the first level subfolder
start_timestrp="`date +"%H:%M:%S"`"

# Ensure the path for the log exists
mkdir -p $basepath
bash ./dmck_allsky_backend.sh $1 $basepath $2 $3 $4 $5 $6 $7 $8 | tee -a $basepath$start_timestrp'/allSky.log' 

# Remove the temporary folder when we finish
rm -r $basepath'/temp_00/'