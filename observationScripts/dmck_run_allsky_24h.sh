﻿
#!/bin/bash
clear
swlevel 2


if test $# != 5; then
	echo "Incorrect inputs: ./dmck_run_allsky.sh <minutes> <duty cycle> <subband mode 3> <subband mode 5> <subband mode 7>"
	exit;
fi

# Initialise seconds (auto counts), a second cache and a run count (duty cycles completed)
SECONDS=0
CURRSEC=0
RUNCOUNT=0

# Determine the duty cycles for the given observation, time of the next cycle begining
let MAXRUNS=$1*60/$2 # 24 hour cycle
NEXTRUNOFFSET=$2
NEXTRUNOFFSET=$((NEXTRUNOFFSET*60))
NEXTRUN=$NEXTRUNOFFSET

echo "Running for $1 minutes ($MAXRUNS cycles) with an intended duty cycle of $2 minutes."

# Add an AllSky suffix to the input folder
folder_name='/data/home/user1/data/'`date +"%Y"`'/'`date +"%m"`'/'`date +"%d"`'/xst/mode'

mkdir -p $folder_name'3/'
mkdir -p $folder_name'5/'
mkdir -p $folder_name'7/'


while true; do
	# Call our generic observations script: <RCUMODE> <SUBBAND> <INT TIME (s)> <OUTPUT FOLDER SUFFIX>
        rspctl --mode=3
	bash ./dmck_xst_generic_mode.sh 3 $3 5 $folder_name'3/'



        rspctl --mode=5
        echo 'Activating HBA tiles by the default (Generic 2015) scheme.'
        python dmck_set_HBA_single_element_pattern_IE613.py -v 2>&1 | tee -a "$folder_name"'last_hba_activation.log'
        sleep 5

        bash ./dmck_xst_generic_mode.sh 5 $4 10 $folder_name'5/'


        rspctl --mode=7
        echo 'Activating HBA tiles by the default (Generic 2015) scheme.'
        python dmck_set_HBA_single_element_pattern_IE613.py -v 2>&1 | tee -a "$folder_name"'last_hba_activation.log'
        sleep 5

        bash ./dmck_xst_generic_mode.sh 7 $5 10 $folder_name'7/'

        # Check the current time to see if our duty cycle length has been met.
        # If it has, sleep until the next start time.
        CURRSEC=$SECONDS
        SLEEPLEN=$(($NEXTRUN-$CURRSEC))
	echo "Performing Sleep Test with $SLEEPLEN s sleep expected."
        if [ $SLEEPLEN -gt 0 ]; then
                echo "Sleeping for $SLEEPLEN s on run $RUNCOUNT"
                sleep $SLEEPLEN
        fi

        # Once the 24 hours are up, exit. Assumings all observations are completed within
        #       one duty cycle length.
        if [ $RUNCOUNT -gt $MAXRUNS ]; then
                exit;
        fi

        let RUNCOUNT++

        echo $SLEEPLEN
        NEXTRUN=$(($NEXTRUN+$NEXTRUNOFFSET))

        echo $CURRSEC, $NEXTRUN


done

