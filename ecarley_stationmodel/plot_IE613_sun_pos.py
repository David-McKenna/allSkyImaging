#!/usr/bin/python

########################################################################
##
## File:
##    plot_IE613_sun_pos.py
##
## Description:
##    Plot the positions of a station's antennas in ECEF coordinates along with 
##    the position of the sun at a particular time.
##
## Disclaimer:
##
## Notes:
##
## Examples:
##   plot_IE613_sun_pos.py /Users/eoincarley/LOFAR/data/IE613/AntennaField.conf
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

import numpy as np
import math
import ephem
import pysolar
from pysolar.solar import *
import pymap3d as pm
import pylab
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys,os
import struct
import time
import datetime
import calendar
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
rcuInfo = [ {'mode':'OFF', 'rcuID':0, 'array_type':'LBA', 'bw':100000000.},            #0
            {'mode':'LBL_HPF10MHZ', 'rcuID':1, 'array_type':'LBA', 'bw':100000000.},   #1
            {'mode':'LBL_HPF30MHZ', 'rcuID':2, 'array_type':'LBA', 'bw':100000000.},   #2
            {'mode':'LBH_HPF10MHZ', 'rcuID':3, 'array_type':'LBA', 'bw':100000000.},   #3
            {'mode':'LBH_HPF30MHZ', 'rcuID':4, 'array_type':'LBA', 'bw':100000000.},   #4
            {'mode':'HBA_110_190MHZ', 'rcuID':5, 'array_type':'HBA', 'bw':100000000.}, #5
            {'mode':'HBA_170_230MHZ', 'rcuID':6, 'array_type':'HBA', 'bw':100000000.}, #6
            {'mode':'HBA_210_290MHZ', 'rcuID':7, 'array_type':'HBA', 'bw':100000000.}] #7

ant_field_file=args[0]
ant_arr_file = "/Users/eoincarley/LOFAR/data/IE613/AntennaArrays.conf"
fh=open(ant_arr_file)
rcumode=3
readLatLon=False
for line in fh:
    if line.lower().startswith(rcuInfo[rcumode]['array_type'].lower()):
        readLatLon=True
        continue
    if readLatLon:
        lon=line.split(' ')[2]
        lat=line.split(' ')[3]
        elev=line.split(' ')[4]
        readLatLon=False
        continue       
fh.close()

#------------------------------------------------
#     Define local geographic coords and time
#
obs = ephem.Observer()
obs.lon = math.radians(float(lon))
obs.lat = math.radians(float(lat))
obs.elevation=float(elev)
time_0 = calendar.timegm( time.strptime("2017-06-21 05:00:00", "%Y-%m-%d %H:%M:%S") )
time_1 = calendar.timegm( time.strptime("2017-06-21 20:00:00", "%Y-%m-%d %H:%M:%S") )
image_num = 0
plt.ion()

#-----------------------------------------
#      Get IE613 antenna positions
lba_xyz = read_ant_xyz(ant_field_file, rcuInfo, 3, 'IE')
LBAX  = [ coords[0]/1e3 for coords in lba_xyz]   
LBAY  = [ coords[1]/1e3 for coords in lba_xyz]   
LBAZ  = [ coords[2]/1e3 for coords in lba_xyz]

hba_xyz = read_ant_xyz(ant_field_file, rcuInfo, 5, 'IE')
HBAX  = [ coords[0]/1e3 for coords in hba_xyz]   
HBAY  = [ coords[1]/1e3 for coords in hba_xyz]   
HBAZ  = [ coords[2]/1e3 for coords in hba_xyz]  

meanx = np.mean(LBAX)
meany = np.mean(LBAY)
meanz = np.mean(LBAZ)

u = np.linspace(0, 2.0*np.pi, 200)
v = np.linspace(0, np.pi, 200)
xsphere0 = 7e-3*np.outer(np.cos(u), np.sin(v)) 
ysphere0 = 7e-3*np.outer(np.sin(u), np.sin(v)) 
zsphere0 = 7e-3*np.outer(np.ones(np.size(u)), np.cos(v)) 

while time_0 < time_1:
    time_utc = time.strftime("%Y/%m/%d %H:%M", time.gmtime(time_0)) #UTC (Irish time in Winter)
    time_ist = time.strftime("%Y/%m/%d %H:%M", time.gmtime(time_0 + 60.0*60.0)) # Irish Standard Time (Irish time in Summer)
    obs.date = time_ist

    #-----------------------------------------
    #           Figure setup
    fig = plt.figure(figsize=(15, 15))
    ax = fig.gca(projection='3d')

    #----------------------------------------------------
    #    Define solar ephem object to get local alt-az     
    #
    sun = ephem.Sun()
    sun.compute(obs)
    alt = math.degrees(sun.alt)
    az = math.degrees(sun.az) 
    #print(alt) 
    #print(az)


    sun_ecef = pm.aer2ecef(az, alt, 150.0, 53.09472, -7.9213880, 75.0, deg=True) #pm.geodetic2ecef(53.09472, -7.921388, 75.0, deg=True)

    lba_to_sunx = [ meanx, sun_ecef[0]/1000.0 ]
    lba_to_suny = [ meany, sun_ecef[1]/1000.0 ]
    lba_to_sunz = [ meanz, sun_ecef[2]/1000.0 ]

    #-----------------------------------------  
    #        Make sphere for the sun
    xsphere = xsphere0 + lba_to_sunx[1]
    ysphere = ysphere0 + lba_to_suny[1]
    zsphere = zsphere0 + lba_to_sunz[1]


    #----------------------------------------- 
    #        Plot and format
    ax.plot(LBAX, LBAY, LBAZ, 'bo', label='LBAs (ECEF coords)')
    ax.plot(HBAX, HBAY, HBAZ, 'o', color='salmon', label='HBAs (ECEF coords)')
    ax.plot(lba_to_sunx, lba_to_suny, lba_to_sunz, 'g')
    ax.plot_surface(xsphere, ysphere, zsphere, color='y', linewidth=0.01)

    #----------------------------------------- 
    #        Plot projections
    #
    LBAZproj = [5.0769e3 for i in LBAZ]
    HBAZproj = [5.0769e3  for i in HBAZ]
    lba_to_sunzproj = [5.0769e3  for i in lba_to_sunz]
    zsphereproj = [5.0769e3  for i in zsphere]

    LBAXproj = [3.80155e3 for i in LBAX]
    HBAXproj = [3.80155e3  for i in HBAX]
    lba_to_sunxproj = [3.80155e3  for i in lba_to_sunx]
    xsphereproj = [3.80155e3  for i in xsphere]

    ax.plot(LBAXproj, LBAY, LBAZ, '.', color='lightgray')
    ax.plot(HBAXproj, HBAY, HBAZ, '.', color='lightgray')
    ax.plot(lba_to_sunxproj, lba_to_suny, lba_to_sunz, color='lightgray')
    ax.plot_surface(xsphereproj, ysphere, zsphere, color='lightgray', linewidth=0.01)

    ax.plot(LBAX, LBAY, LBAZproj, '.', color='lightgray')
    ax.plot(HBAX, HBAY, HBAZproj, '.', color='lightgray')
    ax.plot(lba_to_sunx, lba_to_suny, lba_to_sunzproj, color='lightgray')
    ax.plot_surface(xsphere, ysphere, zsphereproj, color='lightgray', linewidth=0.01)


    #----------------------------------------- 
    #        Plot format
    #
    ax.set_xlim3d(3.80155e3, 3.80155e3+0.25)
    ax.set_ylim3d(-5.2889e2 -0.25, -5.2889e2)
    ax.set_zlim3d(5.0769e3, 5.0769e3+0.25)
    plt.title('IE613 solar alt: %s$^{\circ}$, az: %s$^{\circ}$ @ %s IST' % (round(alt), round(az), time_ist))
    ax.legend()
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.legend()
    ax.view_init(elev=16, azim=-55.0)

    plt.show()
    #pdb.set_trace()
    plt.savefig('/Users/eoincarley/LOFAR/data/IE613/image_'+str(format(image_num, '03'))+'.png')   # save the figure to file
    image_num+=1
    
    time_0 = time_0 + 10.0*60.0



