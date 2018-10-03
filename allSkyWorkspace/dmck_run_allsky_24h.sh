
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
while true; do
        bash ./dmck_xst_generic_mode.sh 3 200 5 AllSky

        # Run 2x mode5/7 as we never got a purely clean subband; sampling multiple should give us a better chance.
        bash ./dmck_xst_generic_mode.sh 5 100 60 AllSky
        bash ./dmck_xst_generic_mode.sh 5 130 60 AllSky

        bash ./dmck_xst_generic_mode.sh 7 100 60 AllSky
        bash ./dmck_xst_generic_mode.sh 7 130 60 AllSky

        CURRSEC=$SECONDS
        SLEEPLEN=$(($NEXTRUN-$CURRSEC))
	echo "Performing Sleep Test with $SLEEPLEN s sleep expected."
        if [ '$SLEEPLEN' > '0' ]; then
                echo "Sleeping for $SLEEPLEN s on run $RUNCOUNT"
                sleep $SLEEPLEN
        fi

        if [ '$RUNCOUNT' > '$MAXRUNS']; then
                # 24 hours with 5 minutes per run
                exit;
        fi

        let RUNCOUNT++

        echo $SLEEPLEN
        NEXTRUN=$(($NEXTRUN+$NEXTRUNOFFSET))

        echo $CURRSEC, $NEXTRUN


done

