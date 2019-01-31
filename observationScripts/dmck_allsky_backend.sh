
#!/bin/bash
swlevel 2


if test $# != 9; then
	echo "Incorrect inputs: ./dmck_allsky_backend.sh <output folder> <minutes> <duty cycle> <subband mode 3 start> <subband mode 3 end> <subband mode 5 start> <subband mode 5 end> <subband mode 7 start> <subband mode 7 end>"
	exit;
fi

# Initialise seconds (auto counts), a second cache and a run count (duty cycles completed)
SECONDS=0
CURRSEC=0
RUNCOUNT=0

# Determine the duty cycles for the given observation, time of the next cycle begining
let MAXRUNS=$2*60/$3
let MAXTIME=$2*60
NEXTRUNOFFSET=$3
NEXTRUNOFFSET=$((NEXTRUNOFFSET*60))
NEXTRUN=$NEXTRUNOFFSET

echo "Running for $2 minutes ($MAXRUNS cycles) with an intended duty cycle of $3 minutes."

# Add an AllSky suffix to the input folder
folder_name=$1'mode'

while true; do
	# Call our generic observations script: <RCUMODE> <START SUBBAND> <END SUBBAND> <SKIP N SUBBANDS> <INT TIME (s)> <OUTPUT FOLDER SUFFIX> ('skip' if skip applying activation)
        rspctl --mode=3
        bash ./dmck_generic_subband_scan.sh 3 $4 $5 10 2 4 $folder_name'3/'


        rspctl --mode=5
        bash ./dmck_generic_subband_scan.sh 5 $6 $7 10 2 4 $folder_name'5/'


        rspctl --mode=7
        bash ./dmck_generic_subband_scan.sh 7 $8 $9 10 2 4 $folder_name'5/'

        # Check the current time to see if our duty cycle length has been met.
        # If it has, sleep until the next start time.
        CURRSEC=$SECONDS
        SLEEPLEN=$(($NEXTRUN-$CURRSEC))
	echo "Performing Sleep Test with $SLEEPLEN s sleep expected."
        if [ $SLEEPLEN -gt 0 ]; then
                echo "Sleeping for $SLEEPLEN s on run $RUNCOUNT"
                sleep $SLEEPLEN
        fi

        # Once the 24 hours are up, exit. 
        if [ $RUNCOUNT -gt $MAXRUNS ]; then
                exit;
        fi
        if [ $SECONDS -gt $MAXTIME ]; then
                exit;
        fi

        let RUNCOUNT++

        echo $SLEEPLEN
        NEXTRUN=$(($NEXTRUN+$NEXTRUNOFFSET))

        echo $CURRSEC, $NEXTRUN


done

