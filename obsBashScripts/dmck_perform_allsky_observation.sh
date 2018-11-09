#!/bin/bash

if test $# != 5; then
	echo "Incorrect inputs: ./dmck_perform_allsky_observation.sh <minutes> <duty cycle> <subband mode 3> <subband mode 5> <subband mode 7>"
	exit;
fi

# Generate the base path of the observation
basepath='/data/home/user1/data/'`date +"%Y"`'/'`date +"%m"`'/'`date +"%d"`'/xst/'

# Get the start time for the first level subfolder
start_timestrp="`date +"%H:%M:%S"`"

# Ensure the path for the log exists
mkdir -p $basepath
bash ./dmck_allsky_backend.sh $1 $basepath $2 $3 $4 $5 | tee -a $basepath$start_timestrp'/allSky.log' 

# Remove the temporary folder when we finish
rm -r $basepath'/temp_00/'