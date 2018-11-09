
#!/bin/bash
swlevel 2


if test $# != 6; then
	echo "Incorrect inputs: ./dmck_allsky_backend.sh <output folder> <minutes> <duty cycle> <subband mode 3> <subband mode 5> <subband mode 7>"
	exit;
fi

# Initialise seconds (auto counts), a second cache and a run count (duty cycles completed)
SECONDS=0
CURRSEC=0
RUNCOUNT=0

# Determine the duty cycles for the given observation, time of the next cycle begining
let MAXRUNS=$2*60/$3
NEXTRUNOFFSET=$3
NEXTRUNOFFSET=$((NEXTRUNOFFSET*60))
NEXTRUN=$NEXTRUNOFFSET

echo "Running for $2 minutes ($MAXRUNS cycles) with an intended duty cycle of $3 minutes."

# Add an AllSky suffix to the input folder
folder_name=$1'mode'

while true; do
	# Call our generic observations script: <RCUMODE> <SUBBAND> <INT TIME (s)> <OUTPUT FOLDER SUFFIX>
        rspctl --mode=3
	bash ./dmck_xst_generic_mode.sh 3 $4 5 $folder_name'3/'



        rspctl --mode=5
        echo 'Activating HBA tiles by the default (Generic 2015) scheme.'
        python dmck_set_HBA_single_element_pattern_IE613.py -v 2>&1 | tee -a "$folder_name"'5/last_hba_activation.log'
        sleep 5

        bash ./dmck_xst_generic_mode.sh 5 $5 10 $folder_name'5/'


        rspctl --mode=7
        echo 'Activating HBA tiles by the default (Generic 2015) scheme.'
        python dmck_set_HBA_single_element_pattern_IE613.py -v 2>&1 | tee -a "$folder_name"'7/last_hba_activation.log'
        sleep 5

        bash ./dmck_xst_generic_mode.sh 7 $6 10 $folder_name'7/'

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

