#!/bin/bash

echo `date`
if test $# != 1; then
	echo "Incorrect inputs: ./dmck_INT_subband_scan.sh <Station Code>"
	exit;
fi
# Generate the base path of the observation
basepath='./'$1'/'`date +"%Y"`'/'`date +"%m"`'/'`date +"%d"`'/XST/mode5/'

# Get the start time for the first level subfolder
start_timestrp="`date +"%H%M%S"`"

# Ensure the path for the log exists
mkdir -p $basepath$start_timestrp

basepath=$basepath$start_timestrp'/'

# Scan subband 51 to 120, sample every 5th subband for 60 seconds
bash ./dmck_generic_subband_scan.sh 51 121 5 60 $basepath $1 | tee -a $basepath'allSky.log' 
echo `date`
