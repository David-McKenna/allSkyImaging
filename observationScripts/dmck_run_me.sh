#!/bin/bash

basepath='/data/home/user1/data/'`date +"%Y"`'/'`date +"%m"`'/'`date +"%d"`'/XST/'

start_timestrp="`date +"%H:%M:%S"`"

mkdir -p $basepath
bash ./dmck_run_allsky_24h.sh $1 | tee -a $basepath$start_timestrp'allSky.log' 

rm -r $basepath'/temp_00/'