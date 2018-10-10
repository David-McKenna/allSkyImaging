#!/bin/bash

# Ensure we've been passed the right amount of parameters.
if test $# != 4; then
   echo "Usage: dmck_xst_generic_mode.sh RCUMODE SUBBAND INTTIME FOLDERNAME"
   exit;
fi

# Setup an output data structure in /data/home/user1/<YY>/<MM>/<DD>/XST
# We will use a temp folder for the given day to ingest observations before 
# 	adding metadata into their names and moving them to .../mode_X/<FOLDERNAME>/sb<SB>/
basepath='/data/home/user1/data/'`date +"%Y"`'/'`date +"%m"`'/'`date +"%d"`'/XST'
datapath=$basepath'/mode'$1'/'$4'/sb'$2'/'
datapathnosubband=$basepath'/mode'$1'/'$4'/'
temppath=$basepath'/temp_00/'

# Make the temporary folder if it doesn't exist already (-p creates parent folders too)
mkdir -p $temppath

# Rebind $1 to a human readable name and set the bitmode to 16
rcumode=$1
bits=16

# Log the start time of the observation
start_time="`date +"%Y/%m/%d@%H:%M:%S"`"
start_timestrp="`date +"%H%M%S"`"

# Log our intended actions
echo "*** RCU Mode $1, Subband $2, Integration $3 Data being saved to $datapath"

# Ensure the data directory exists and create a copy of our scripts inside it for future reference
mkdir -p $datapath
cp -n ./dmck_* $datapathnosubband

# Rebind $2 to a human readable name
subband="$2"


# If we are observing with HBAs, apply an activation pattern so we don't have tiled sidelobes
rspctl --mode=$rcumode
if [ $rcumode -gt 4 ]; then
	echo 'Activating HBA tiles by the default (Generic 2015) scheme.'
	python dmck_set_HBA_single_element_pattern_IE613.py -v 2>&1 | tee -a "$datapath"'last_hba_activation.log'
fi

# Switch over to 16 bit and the intended observation subband
rspctl --bitmode=$bits
rspctl --xcsubband=$subband

# Start the XST observation, log the intent
echo "Begining $3 second integration in Mode $1 Subband $2"
rspctl --xcstatistics --integration=$3 --duration=$3 --directory=$temppath

# Generate the output filename from the observation, move the file to the intended location
filename=$(ls -1 $temppath*dat)
newname=$(echo $filename | sed s/_xst/_sb"$subband"_xst/g)
newname=${datapath}$(echo $newname | sed -e 's,'$temppath',,g')

end_time="`date +"%Y/%m/%d@%H:%M:%S"`"
echo 'Moving mode '$1' subband '$2' data to ' $newname
echo "XST mode $1 subband $2 observation ended at $end_time"
mv $filename $newname

# Log everything for debug reasons.
printf 'Mode '$1'\nSubband '$2'\nIntegrationTime '$3's\nStartTime '$start_time'\nEndTime '$end_time >> $newname'.log'
cat $datapath'last_hba_activation.log' >> $newname'.log'


# Add some spacing before the next observation starts.
printf "\n\n\n\n\n"
exit;
