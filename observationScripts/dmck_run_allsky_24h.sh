
#!/bin/bash
clear
swlevel 3


if test $# != 1; then
	echo "Incorrect inputs: ./script.sh <minutes>"
	exit;
fi

SECONDS=0
CURRSEC=0
RUNCOUNT=0
let MAXRUNS=24*60/$1 # 24 hour cycle
NEXTRUNOFFSET=$1
NEXTRUNOFFSET=$((NEXTRUNOFFSET*60))
NEXTRUN=$NEXTRUNOFFSET

echo "Running for 24 hours ($MAXRUNS cycles) with a duty cycle of $1."
start_time="`date +"%H%M%S"`"

folder_name=$start_time'/AllSky'
while true; do
	#bash ./dmck_xst_generic_mode.sh 3 210 10 $folder_name
	#bash ./dmck_xst_generic_mode.sh 3 210 5 $folder_name


        # Run 2x mode5/7 as we never got a purely clean subband; sampling multiple should give us a better chance.
        bash ./dmck_xst_generic_mode.sh 5 281 60 $folder_name
        bash ./dmck_xst_generic_mode.sh 5 359 60 $folder_name

        #bash ./dmck_xst_generic_mode.sh 7 100 60 $folder_name
        #bash ./dmck_xst_generic_mode.sh 7 130 60 $folder_name

        CURRSEC=$SECONDS
        SLEEPLEN=$(($NEXTRUN-$CURRSEC))
	echo "Performing Sleep Test with $SLEEPLEN s sleep expected."
        if [ $SLEEPLEN -gt 0 ]; then
                echo "Sleeping for $SLEEPLEN s on run $RUNCOUNT"
                sleep $SLEEPLEN
        fi

        if [ $RUNCOUNT -gt $MAXRUNS ]; then
                # 24 hours with 5 minutes per run
                exit;
        fi

        let RUNCOUNT++

        echo $SLEEPLEN
        NEXTRUN=$(($NEXTRUN+$NEXTRUNOFFSET))

        echo $CURRSEC, $NEXTRUN


done

