
#!/bin/bash
clear
swlevel 3


if test $# != 6; then
	echo "Incorrect inputs: ./script.sh <start subband> <end subband> <sample every nth subband> <integration time> <output folder> <Station Name>"
	exit;
fi

SUBBAND=$1
END_SUBBAND=$2
SPARSE_VAR=$3
INT=$4
# Initialise seconds (auto counts), a second cache and a run count (duty cycles completed)

# Add a time and AllSky suffix to the input folder
folder_name=$5'AllSky/'

# Generate the data output path, if we're on IE613 we will be copying data back and need this.
basepathIE='/data/home/user1/data/'`date +"%Y"`'/'`date +"%m"`'/'`date +"%d"`'/XST/'

while [ $SUBBAND -lt $END_SUBBAND ]; do
        echo "Integrating subband $SUBBAND for $INT seconds."
        cp  '/run/media/MCKENND2/16a416df-eb7b-42cc-be13-400117b8d5c1/offload_allsky/24hData/mode3/sb210/5s/20181001_111303_sb210_xst.dat' $folder_name'temp_00/'
        bash ./dmck_xst_generic_mode.sh 5 $SUBBAND $INT $folder_name

        # If we are on IE613, copy the data base if our subband is divisible by 5.
        if [ $6 = 'IE613' ]; then
                if [ $(($SUBBAND % 5)) -eq 0 ]; then
                        rsync -rzu $basepathIE USERNEEDED@LGCIPNEEDED:/OUTFOLDERFOLDERNEEDED
                fi
        fi

        # Update the subband for the next run
        SUBBAND=$((SUBBAND + SPARSE_VAR))

done

