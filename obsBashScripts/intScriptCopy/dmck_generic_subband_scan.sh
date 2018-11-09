
#!/bin/bash
clear
swlevel 2


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

rcumode=5
# If we are observing with HBAs, apply an activation pattern so we don't have tiled sidelobes
rspctl --mode=$rcumode

if [ $rcumode -gt 4 ]; then
	echo 'Activating HBA tiles by the default (Generic 2015) scheme.'
	python dmck_set_HBA_single_element_pattern_IE613.py -v 2>&1 | tee -a "$folder_name"'last_hba_activation.log'
fi

rspctl --realdelays  | tee 'realdelays1.txt'

while [ $SUBBAND -lt $END_SUBBAND ]; do
        echo "Integrating subband $SUBBAND for $INT seconds."

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
rspctl --realdelays  | tee 'realdelays2.txt'
