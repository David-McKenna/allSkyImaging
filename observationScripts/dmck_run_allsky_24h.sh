
#!/bin/bash
clear
swlevel 3


if test $# != 2; then
	echo "Incorrect inputs: ./script.sh <minutes> <start time/head folder name>"
	exit;
fi

# Initialise seconds (auto counts), a second cache and a run count (duty cycles completed)
SECONDS=0
CURRSEC=0
RUNCOUNT=0

# Determine the duty cycles for the given observation, time of the next cycle begining
let MAXRUNS=24*60/$1 # 24 hour cycle
NEXTRUNOFFSET=$1
NEXTRUNOFFSET=$((NEXTRUNOFFSET*60))
NEXTRUN=$NEXTRUNOFFSET

echo "Running for 24 hours ($MAXRUNS cycles) with an intended duty cycle of $1 minutes."

# Add an AllSky suffix to the input folder
folder_name=$2'/AllSky'

while true; do
	# Call our generic observations script: <RCUMODE> <SUBBAND> <INT TIME (s)> <OUTPUT FOLDER SUFFIX>
        # Currently takes ~ 7.45 to complete one duty cycle on IE613 with the HBA activations
        bash ./dmck_xst_generic_mode.sh 3 210 10 $folder_name
	bash ./dmck_xst_generic_mode.sh 3 210 5 $folder_name


        # Run 2x mode5/7 as we never got a purely clean subband; sampling multiple should give us a better chance.
        bash ./dmck_xst_generic_mode.sh 5 281 60 $folder_name
        bash ./dmck_xst_generic_mode.sh 5 359 60 $folder_name

        bash ./dmck_xst_generic_mode.sh 7 100 60 $folder_name
        bash ./dmck_xst_generic_mode.sh 7 130 60 $folder_name

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

