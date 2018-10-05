#!/usr/bin/python
# set_HBA_single_element_pattern.py
# set the HBA tile single elements for all-sky imaging
# 2009 Nov 10  James M Anderson  --MPIfR  start
# 2015 Nov 26  Menno Norden - Astron upgrade core & remote
# 2015 Dec 03  Michiel Brentjens - Optimum elements determined
import os
import time
import sys
import optparse
import subprocess
NUM_HBA_ANTENNAS = 16


Effelsberg_elements_20091110 = [1,4,13,15,11,9,14,1,15,0,8,2,11,3,14,0,2,4,3,0,0,2,12,12,12,12,15,11,14,15,7,5,1,0,3,10,1,11,0,12 \
,12,1,6,7,0,10,9,6,15,14,11,7,2,0,7,12,15,8,13,3,7,6,3,15,11,1,4,11,8,1,8,15,4,0,5,6,12,0,12,15,3,7,14,8,3,12,12,2,9,8,14,2,5,6,12,0]

Generic_International_Station_20091110 = [15,0,15,3,9,15,14,2,0,3,4,14,10,8,5,15,12,0,2,11,3,12,12,1,5,4,4,8,6,3,0,5,3,11,3,2,8 \
,15,13,8,3,2,9,1,14,8,8,0,12,13,0,11,15,3,12,3,13,3,10,5,0,10,1,6,4,10,3,15,3,14,0,12,0,7,0,12,7,3,13,0,7,3,15,4,14,4,3,8,4,9,12,0,14,9,3,11]

# Optimum element calculation done by M.Brentjes (Dec 2015)
Generic_Int_201512 = [0,5,3,1,8,3,12,15,10,13,11,5,12,12,5,2,10,8,0,3,5,1,4,0,11,6,2,4,9,14,15,3,7,5,13,15,5,6,5,12,15,7,1,1,14 \
,9,4,9,3,9,3,13,7,14,7,14,2,8,8,0,1,4,2,2,12,15,5,7,6,10,12,3,3,12,7,4,6,0,5,9,1,10,10,11,5,11,7,9,7,6,4,4,15,4,1,15]
Generic_Core_201512 = [0,10,4,3,14,0,5,5,3,13,10,3,12,2,7,15,6,14,7,5,7,9,0,15,0,10,4,3,14,0,5,5,3,13,10,3,12,2,7,15,6,14,7,5,7,9,0,15]
Generic_Remote_201512 = [0,13,12,4,11,11,7,8,2,7,11,2,10,2,6,3,8,3,1,7,1,15,13,1,11,1,12,7,10,15,8,2,12,13,9,13,4,5,5,12,5,5,9,11,15,12,2,15]

def set_hba_antennas(elements, sleeptime, verbose=False, process=True):
    """Goes through antenna by antenna to select which RCUs to set to this
antenna.  For most stations, this results in fewer rspctl commands being sent
than going through tile by tile.
"""

    if(verbose): # Added by JmCC for analysing program
        print 'Using ' + options.element_list
    NUM_TILES = len(elements)
    for antenna in range(NUM_HBA_ANTENNAS): #for each of the 16 elements in a tile select which one will be on
        rcu_set = False
        rcu_ctrl_string = '='
        for tile in range(NUM_TILES): # for each of the 96 tiles
            if(elements[tile] == antenna):
                rcu_set = True
                # convert tile number to RCU number
                rcu_ctrl_string += "%d,%d,"%(2*tile,2*tile+1)
        if(not rcu_set):
            continue
        rcu_ctrl_string = rcu_ctrl_string[0:len(rcu_ctrl_string)-1]
        ant_ctrl_string='='
        for index in range(NUM_HBA_ANTENNAS):
            if index == antenna:
                ant_ctrl_string=ant_ctrl_string + '128,' # 128 means LNA on, RF ON, 0ns Delay
            else:	
                ant_ctrl_string=ant_ctrl_string + '2,' # 2 means LNA off, RF OFF (30dB supression of signal), 0ns Delay

        ant_ctrl_string=ant_ctrl_string[0:len(ant_ctrl_string)-1]
        cmd_str='rspctl --hbadelay' + ant_ctrl_string + ' --select' + rcu_ctrl_string # put them together
        if(verbose): # print the generated command
            sys.stderr.write("Running '%s'\n"%cmd_str)
        if(process): # run the generated command
            subprocess.call(cmd_str.split(' '))
            time.sleep(sleeptime)
	time.sleep(sleeptime)
	processObj = subprocess.Popen(['/opt/lofar/bin/rspctl', '--realdelays'], stdout = subprocess.PIPE)


#############################################################
if __name__ == "__main__":
    p = optparse.OptionParser()
    p.add_option("--verbose", "-v", action="store_true")
    p.add_option("--echo_only", action="store_true")
    p.add_option("--sleeptime",default=2.0,type='float',help="Time to delay between individual RCU elements being turned on, in seconds.  May be floating point")
    p.add_option("--element_list", default='International_2015', choices=['Effelsberg','International_2009', 'International_2015', 'Core', 'Remote'], help="The station antenna element list to use.\nOptions are:\nEffelsberg\nInternational_2009\nInternational_2015\nCopre\nRemote")
    options, arguments = p.parse_args()
    if(options.element_list == 'Effelsberg'):
        set_hba_antennas(Effelsberg_elements_20091110,
                         options.sleeptime,
                         options.verbose,not options.echo_only)
    elif(options.element_list == 'International_2009'):
        set_hba_antennas(Generic_International_Station_20091110,
                         options.sleeptime,
                         options.verbose,not options.echo_only)
    elif(options.element_list == 'International_2015'):
        set_hba_antennas(Generic_Int_201512,
                         options.sleeptime,
                         options.verbose,not options.echo_only)
    elif(options.station == 'Core'):
        set_hba_antennas(Generic_Core_201512,
                         options.sleeptime,
                         options.verbose,not options.echo_only)
    elif(options.station == 'Remote'):
        set_hba_antennas(Generic_Remote_201512,
                         options.sleeptime,
                         options.verbose,not options.echo_only)
    else:
        raise RuntimeError("Unknown element list '%s'"%options.element_list)
