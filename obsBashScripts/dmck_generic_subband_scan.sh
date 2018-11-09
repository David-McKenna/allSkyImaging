
#!/bin/bash
swlevel 2


if test $# != 6; then
	echo "Incorrect inputs: ./dmck_generic_subband_scan.sh <rcu mode> <start subband> <end subband> <sample every nth subband> <integration time> <output folder>"
	exit;
fi

SUBBAND=$2
END_SUBBAND=$3
SPARSE_VAR=$4
INT=$5
# Initialise seconds (auto counts), a second cache and a run count (duty cycles completed)

# Add a time and AllSky suffix to the input folder
folder_name=$6'AllSky/'

rcumode=$1
# If we are observing with HBAs, apply an activation pattern so we don't have tiled sidelobes
rspctl --mode=$rcumode

if [ $rcumode -gt 4 ]; then
	echo 'Activating HBA tiles by the default (Generic 2015) scheme.'
	python dmck_set_HBA_single_element_pattern_IE613.py -v 2>&1 | tee -a "$folder_name"'last_hba_activation.log'
        sleep 5
fi

rspctl --realdelays  | tee 'realdelays1.txt'

while [ $SUBBAND -lt $END_SUBBAND ]; do
        echo "Integrating subband $SUBBAND for $INT seconds."

        bash ./dmck_xst_generic_mode.sh 5 $SUBBAND $INT $folder_name

        # Update the subband for the next run
        SUBBAND=$((SUBBAND + SPARSE_VAR))

done
rspctl --realdelays  | tee 'realdelays_end.txt'
