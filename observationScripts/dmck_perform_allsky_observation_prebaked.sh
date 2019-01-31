#!/bin/bash

if test $# != 2; then
	echo "Incorrect inputs: ./dmck_perform_allsky_observation_prebaked.sh <minutes> <duty cycle>"
	exit;
fi

# Generate the base path of the observation
basepath='/data/home/user1/data/'`date +"%Y"`'/'`date +"%m"`'/'`date +"%d"`'/xst/'

# Get the start time for the first level subfolder
start_timestrp="`date +"%H:%M:%S"`"

# Ensure the path for the log exists
mkdir -p $basepath
bash ./dmck_allsky_backend_prebaked.sh $basepath $1 $2 | tee -a $basepath$start_timestrp'/allSky.log' 

# Remove the temporary folder when we finish
rm -r $basepath'/temp_00/'