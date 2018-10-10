
#!/bin/bash
clear
swlevel 3


if test $# != 5; then
	echo "Incorrect inputs: ./script.sh <start subband> <end subband> <sample every nth subband> <integration time> <output folder suffix> <IE613 Bool>"
	exit;
fi

SUBBAND=$1
END_SUBBAND=$2
SPARSE_VAR=$3
INT=$4
# Initialise seconds (auto counts), a second cache and a run count (duty cycles completed)

# Add an AllSky suffix to the input folder
folder_name=$5'/AllSky'

# Generate the data output path, if we're on IE613 we will be copying data back.
basepath='/data/home/user1/data/'`date +"%Y"`'/'`date +"%m"`'/'`date +"%d"`'/XST'

while [ $SUBBAND -lt $END_SUBBAND]; do
        echo "Integrating subband $SUBBAND for $INT seconds."
        bash ./dmck_xst_generic_mode.sh 5 $SUBBAND $INT $folder_name

        # If we are on IE613, copy the data base if our subband is divisible by 5.
        if [ $6 -eq 1 ]; do
                if [ $(($SUBBAND % 5)) -eq 0 ]; do
                        rsync -rzu $basepath USERNEEDED@LGCIPNEEDED:/OUTFOLDERFOLDERNEEDED
                fi
        fi

        # Update the subband for the next run
        SUBBAND=$((SUBBAND + SPARSE_VAR))

done

