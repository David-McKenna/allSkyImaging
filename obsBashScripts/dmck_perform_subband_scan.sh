#!/bin/bash

echo `date`

if test $# != 1; then
	echo "Incorrect inputs: ./dmck_perform_subband_scan.sh <rcu mode>"
	exit;
fi

# Generate the base path of the observation
basepath='/data/home/user1/data/'`date +"%Y"`'/'`date +"%m"`'/'`date +"%d"`'/xst/mode$1/'

# Get the start time for the first level subfolder
start_timestrp="`date +"%H%M%S"`"

# Ensure the path for the log exists
mkdir -p $basepath$start_timestrp

basepath=$basepath$start_timestrp'/'

# Scan subband 51 to 451, sample every 5th subband for 10 seconds
bash ./dmck_generic_subband_scan.sh 51 452 5 10 $basepath $1 | tee -a $basepath'allSky.log' 
