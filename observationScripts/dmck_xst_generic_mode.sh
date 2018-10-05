#!/bin/bash



#if test $# == 0; then
#
#fi
if test $# != 4; then
   echo "Usage: dmck_xst_generic_mode.sh RCUMODE SUBBAND INTTIME FOLDERNAME"
   exit;
fi

basepath='/data/home/user1/data/'`date +"%Y"`'/'`date +"%m"`'/'`date +"%d"`'/XST'
datapath=$basepath'/mode'$1'/'$4'/sb'$2'/'
temppath=$basepath'/temp_00/'

#rm -r $temppath
mkdir -p $temppath

#record until the following date/time
#stopdate=$1
#stoptime=$2

rcumode=$1
bits=16


start_time="`date +"%Y/%m/%d@%H:%M:%S"`"
start_timestrp="`date +"%H:%M:%S"`"

#write data datapath on LCU
echo "*** Mode $1 Subband $2 Data being saved to $datapath"

#create the data directory & copy this script to it as a record of what we've done.
mkdir -p $datapath
cp -n ./dmck_xst_generic_mode.sh $datapath

subband="$2"

rspctl --mode=$rcumode
if [ $rcumode -gt 4 ]; then
	echo 'Activating HBA tiles by the default (Generic 2015) scheme.'
	python set_HBA_single_element_pattern_IE613.py -v 2>&1 | tee -a $basepath$datapath$start_timestrp'last_hba_activation.log'
fi

rspctl --bitmode=$bits

rspctl --xcsubband=$subband

echo "Begining $3 second integration in Mode $1 Subband $2"
rspctl --xcstatistics --integration=$3 --duration=$3 --directory=$temppath

filename=$(ls -1 $temppath*dat)
newname=$(echo $filename | sed s/_xst/_sb"$subband"_xst/g)
newname=${datapath}$(echo $newname | sed -e 's,'$temppath',,g')

end_time="`date +"%Y/%m/%d@%H:%M:%S"`"

printf 'Mode '$1'\nSubband '$2'\nIntegrationTime '$3's\nStartTime '$start_time'\nEndTime '$end_time >> $newname'.log'
cat $basepath$datapath$start_timestrp'last_hba_activation.log' >> $newname'.log'
echo 'Moving mode '$1' subband '$2' data to ' $newname
mv $filename $newname
echo "XST mode $1 subband $2 observation ended at $end_time"
printf "\n\n\n\n\n"
exit;
