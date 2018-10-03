#!/usr/bin/env python

import os,time,sys
import optparse

NUM_HBA_ANTENNAS = 16

#import random
#el=range(NUM_HBA_ANTENNAS)
#a=[]
#for i in range(96):
#    a.append(random.choice(el))
Effelsberg_elements_20091110 = [1,4,13,15,11,9,14,1,15,0,8,2,11,3,14,0,2,4,3,0,0,2,12,12,12,12,15,11,14,15,7,5,1,0,3,10,1,11,0,12,12,1,6,7,0,10,9,6,15,14,11,7,2,0,7,12,15,8,13,3,7,6,3,15,11,1,4,11,8,1,8,15,4,0,5,6,12,0,12,15,3,7,14,8,3,12,12,2,9,8,14,2,5,6,12,0]

Generic_International_Station_20091110 = [15,0,15,3,9,15,14,2,0,3,4,14,10,8,5,15,12,0,2,11,3,12,12,1,5,4,4,8,6,3,0,5,3,11,3,2,8,15,13,8,3,2,9,1,14,8,8,0,12,13,0,11,15,3,12,3,13,3,10,5,0,10,1,6,4,10,3,15,3,14,0,12,0,7,0,12,7,3,13,0,7,3,15,4,14,4,3,8,4,9,12,0,14,9,3,11]

def set_hba_antennas(elements,sleeptime,verbose=False,process=True):
    """Goes through antenna by antenna to select which RCUs to set to this
antenna.  For most stations, this results in fewer rspctl commands being sent
than going through tile by tile.
"""
    NUM_TILES = len(elements)
    for antenna in range(NUM_HBA_ANTENNAS):
        rcu_set = False
        rcu_ctrl_string = '='
        for tile in range(NUM_TILES):
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
                ant_ctrl_string=ant_ctrl_string + '128,'
            else:	
                ant_ctrl_string=ant_ctrl_string + '2,'

        ant_ctrl_string=ant_ctrl_string[0:len(ant_ctrl_string)-1]
        cmd_str='rspctl --hbadelay' + ant_ctrl_string + ' --select' + rcu_ctrl_string
        if(verbose):
            sys.stderr.write("Running '%s'\n"%cmd_str)
        if(process):
            os.popen(cmd_str)
            time.sleep(sleeptime)

if __name__ == "__main__":
    p = optparse.OptionParser()
    p.add_option("--verbose", "-v", action="store_true")
    p.add_option("--echo_only", action="store_true")
    p.add_option("--sleeptime",default=2.0,type='float',help="Time to delay between individual RCU elements being turned on, in seconds.  May be floating point")
    p.add_option("--element_list", default='International', choices=['Effelsberg','International'], help="The station antenna element list to use.\nOptions are:\nEffelsberg\nInternational")
    options, arguments = p.parse_args()
    if(options.element_list == 'Effelsberg'):
        set_hba_antennas(Effelsberg_elements_20091110,
                         options.sleeptime,
                         options.verbose,not options.echo_only)
    elif(options.element_list == 'International'):
        set_hba_antennas(Generic_International_Station_20091110,
                         options.sleeptime,
                         options.verbose,not options.echo_only)
    else:
        raise RuntimeError("Unknown element list '%s'"%options.element_list)

