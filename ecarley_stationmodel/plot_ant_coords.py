#!/usr/bin/python

########################################################################
##
## File:
##    plot_ant_coords.py
##
## Description:
##    Plot the positions of a station's antennas in ECEF coordinates
##
## Disclaimer:
##
## Notes:
##
## Examples:
##   plot_ant_coords.py /Users/eoincarley/LOFAR/data/IE613/AntennaField.conf
##
## Testing:
##
## Author:
##    Eoin Carley, School of Physics, TCD
##
## Created/Updated:
##    20-Nov-2017 -- v0.1, Eoin Carley
##
########################################################################

import numpy
import numpy as np
import ephem
import pymap3d as pm
import pylab
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys,os
import struct
import time
from time import sleep
import pdb

#########################################################   
#   Functions
#
def read_ant_xyz(ant_field_file, rcuInfo, rcumode, station):
    fh=open(ant_field_file)
    readXYZ=False
    readArrayHDR=False
    readArray=False
    linecnt=0
    nelem=0
    ants=[]

    for line in fh:
        if line.lower().startswith(rcuInfo[rcumode]['array_type'].lower()):
            readXYZ=True
            continue
        if readXYZ:
            arr_x=line.split(' ')[2]
            arr_y=line.split(' ')[3]
            arr_z=line.split(' ')[4]
            readXYZ=False
            readArrayHDR=True
            continue
        if readArrayHDR:
            if station == 'IE': 
                nelem=eval(line.split()[0])[1]+1
            else:    
                nelem=int(line.split(' ')[0])    #This code expects '96' here. But the IE613 file has '(0,95)'. 
            readArrayHDR=False
            readArray=True
            continue
        if readArray and linecnt < nelem:
            cl=' '.join(line.split())
            ants.append(map(float,cl.split(' ')))
            linecnt+=1   
    fh.close()        

    xyz=[]
    for a in ants: xyz.append([a[0]+float(arr_x),a[1]+float(arr_y),a[2]+float(arr_z)])
    xyz=numpy.array(xyz) # Pretty sure this is Earth Centered, Earth Fixed (ECEF) Coords.

    return xyz

#
#########################################################   


if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('Plot XYZ of IE613 LOFAR Antennas')
    o.set_description(__doc__)
    opts, args = o.parse_args(sys.argv[1:])  

#
#########################################################   

ant_field_file=args[0]

rcuInfo = [ {'mode':'OFF', 'rcuID':0, 'array_type':'LBA', 'bw':100000000.},            #0
            {'mode':'LBL_HPF10MHZ', 'rcuID':1, 'array_type':'LBA', 'bw':100000000.},   #1
            {'mode':'LBL_HPF30MHZ', 'rcuID':2, 'array_type':'LBA', 'bw':100000000.},   #2
            {'mode':'LBH_HPF10MHZ', 'rcuID':3, 'array_type':'LBA', 'bw':100000000.},   #3
            {'mode':'LBH_HPF30MHZ', 'rcuID':4, 'array_type':'LBA', 'bw':100000000.},   #4
            {'mode':'HBA_110_190MHZ', 'rcuID':5, 'array_type':'HBA', 'bw':100000000.}, #5
            {'mode':'HBA_170_230MHZ', 'rcuID':6, 'array_type':'HBA', 'bw':100000000.}, #6
            {'mode':'HBA_210_290MHZ', 'rcuID':7, 'array_type':'HBA', 'bw':100000000.}] #7

fig = plt.figure()
ax = fig.gca(projection='3d')

lba_xyz = read_ant_xyz(ant_field_file, rcuInfo, 3, 'IE')
LBAX  = [ coords[0]/1e3 for coords in lba_xyz]   
LBAY  = [ coords[1]/1e3 for coords in lba_xyz]   
LBAZ  = [ coords[2]/1e3 for coords in lba_xyz]

sun_ecef = pm.aer2ecef(45.0, 45.0, 0.0, 53.09472, -7.9213880, 300.0, deg=True) #pm.geodetic2ecef(53.09472, -7.921388, 75.0, deg=True)
meanx = np.mean(LBAX)
meany = np.mean(LBAY)
meanz = np.mean(LBAZ)
lba_to_sunx = [ meanx, sun_ecef[0]/1000.0 ]
lba_to_suny = [ meany, sun_ecef[1]/1000.0 ]
lba_to_sunz = [ meanz, sun_ecef[2]/1000.0 ]

#LBAX.append(sun_ecef[0]/1000.)
#LBAY.append(sun_ecef[1]/1000.)
#LBAZ.append(sun_ecef[2]/1000.)

hba_xyz = read_ant_xyz(ant_field_file, rcuInfo, 5, 'IE')
HBAX  = [ coords[0]/1e3 for coords in hba_xyz]   
HBAY  = [ coords[1]/1e3 for coords in hba_xyz]   
HBAZ  = [ coords[2]/1e3 for coords in hba_xyz]    

ax.plot(LBAX, LBAY, LBAZ, 'bo', label='LBAs')
ax.plot(HBAX, HBAY, HBAZ, 'go', label='HBAs')
ax.plot(lba_to_sunx, lba_to_suny, lba_to_sunz, 'r')
ax.plot(lba_to_sunx, lba_to_suny, lba_to_sunz, 'ro')


ax.set_xlim3d(3.8016e3, 3.8016e3+0.15)
ax.set_ylim3d(-5.2892e2, -5.2892e2-0.15)
ax.set_zlim3d(5.0769e3, 5.0769e3+0.15)


plt.title('IE613 antenna positions in ECEF coordinates')
ax.legend()
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
ax.legend()

"""
ant_field_file = '/Users/eoincarley/LOFAR/data/chilbolton/AntennaField.conf'  
lba_xyz = read_ant_xyz(ant_field_file, rcuInfo, 3, 'UK')
LBAX  = [ coords[0]/1e3 for coords in lba_xyz]   
LBAY  = [ coords[1]/1e3 for coords in lba_xyz]   
LBAZ  = [ coords[2]/1e3 for coords in lba_xyz]
hba_xyz = read_ant_xyz(ant_field_file, rcuInfo, 5, 'UK')
HBAX  = [ coords[0]/1e3 for coords in hba_xyz]   
HBAY  = [ coords[1]/1e3 for coords in hba_xyz]   
HBAZ  = [ coords[2]/1e3 for coords in hba_xyz]    
ax.plot(LBAX, LBAY, LBAZ, 'ro', label='UK608 LBAs')
ax.plot(HBAX, HBAY, HBAZ, 'ro', label='UK608 HBAs')
ax.legend()

ant_field_file = '/Users/eoincarley/LOFAR/data/CS006/AntennaField.conf'  
lba_xyz = read_ant_xyz(ant_field_file, rcuInfo, 3, 'CS')
LBAX  = [ coords[0]/1e3 for coords in lba_xyz]   
LBAY  = [ coords[1]/1e3 for coords in lba_xyz]   
LBAZ  = [ coords[2]/1e3 for coords in lba_xyz]
hba_xyz = read_ant_xyz(ant_field_file, rcuInfo, 5, 'CS')
HBAX  = [ coords[0]/1e3 for coords in hba_xyz]   
HBAY  = [ coords[1]/1e3 for coords in hba_xyz]   
HBAZ  = [ coords[2]/1e3 for coords in hba_xyz]    
ax.plot(LBAX, LBAY, LBAZ, 'go', label='CS006 LBAs')
ax.plot(HBAX, HBAY, HBAZ, 'go', label='CS006 HBAs')
"""


plt.show()