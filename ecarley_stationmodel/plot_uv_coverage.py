#!/usr/bin/etc python

########################################################################
##
## File:
##    plot_uv_coverage.py
##
## Description:
##    Plot the uv coverage for a given set of antenna positions
##
## Disclaimer:
##
## Notes:
##
## Examples:
##   python2.7  plot_uv_coverage.py /Users/eoincarley/LOFAR/data/2017.08.21/AntennaField.conf
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
import calendar
import pdb

######################
#  
#   Functions
#
def eq2top_m(ha, dec):
    """Return the 3x3 matrix converting equatorial coordinates to topocentric
    at the given hour angle (ha) and declination (dec)."""
    sin_H, cos_H = numpy.sin(ha), numpy.cos(ha)
    sin_d, cos_d = numpy.sin(dec), numpy.cos(dec)
    zero = numpy.zeros_like(ha)
    map =  numpy.array([[    sin_H    ,       cos_H  ,       zero  ],
                        [ -sin_d*cos_H,   sin_d*sin_H,      cos_d  ],
                        [  cos_d*cos_H,  -cos_d*sin_H,      sin_d  ]])
    if len(map.shape) == 3: map = map.transpose([2, 0, 1])
    return map
def get_baseline(i, j, src, obs):
    """Return the baseline corresponding to i,j""" 
    bl = j - i
    try:
        if src.alt < 0:
            raise PointingError('Phase center below horizon')
        m=src.map
    except(AttributeError):
        ra,dec = src._ra,src._dec
        #opposite HA since the we want to move the source at zenith away to phase to the original zenith source
        m = eq2top_m(ra-obs.sidereal_time(), dec)
        #normal HA
        #m = eq2top_m(obs.sidereal_time() - ra, dec)
    return numpy.dot(m, bl).transpose()
def gen_uvw(i, j, src, obs, f):
    """Compute uvw coordinates of baseline relative to provided FixedBody"""
    x,y,z = get_baseline(i,j,src,obs)
    afreqs = numpy.reshape(f, (1,f.size))
    afreqs = afreqs/ephem.c #1/wavelength
    if len(x.shape) == 0: return numpy.array([x*afreqs, y*afreqs, z*afreqs])
    x.shape += (1,); y.shape += (1,); z.shape += (1,)
    return numpy.array([numpy.dot(x,afreqs), numpy.dot(y,afreqs), numpy.dot(z,afreqs)])
def xyz2uvw(xyz, src, obs, f):
    """Return an array of UVW values"""
    uvw=numpy.zeros((xyz.shape[0],xyz.shape[0],3))
    for i in range(xyz.shape[0]):
        for j in range(xyz.shape[0]):
            if i==j: continue
            uvw[i,j]=gen_uvw(xyz[i], xyz[j], src, obs, f)[:,0,0]
    return uvw
def dft2(d,k,l,u,v):
    """compute the 2d DFT for position (k,l) based on (d,uvw)"""
    return numpy.sum(d*numpy.exp(-2.*numpy.pi*1j*((u*k) + (v*l))))
    #psf:
    return numpy.sum(numpy.exp(-2.*numpy.pi*1j*((u*k) + (v*l))))
def dftImage(d,uvw,px,res,mask=False):
    """return a DFT image"""
    nants=uvw.shape[0]
    im=numpy.zeros((px[0],px[1]),dtype=complex)
    mid_k=int(px[0]/2.)
    mid_l=int(px[1]/2.)
    u=uvw[:,:,0]
    v=uvw[:,:,1]
    w=uvw[:,:,2]
    u/=mid_k
    v/=mid_l
    start_time=time.time()
    for k in range(px[0]):
        for l in range(px[1]):
            im[k,l]=dft2(d,(k-mid_k),(l-mid_l),u,v)
            if mask:        #mask out region beyond field of view
                rad=(((k-mid_k)*res)**2 + ((l-mid_l)*res)**2)**.5
                if rad > mid_k*res: im[k,l]=0
                #else: im[k,l]=dft2(d,(k-mid_k),(l-mid_l),u,v)
    print time.time()-start_time
    return im
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
ant_arr_file = "/Users/eoincarley/LOFAR/data/IE613/AntennaArrays.conf"
rcuInfo = [ {'mode':'OFF', 'rcuID':0, 'array_type':'LBA', 'bw':100000000.},            #0
            {'mode':'LBL_HPF10MHZ', 'rcuID':1, 'array_type':'LBA', 'bw':100000000.},   #1
            {'mode':'LBL_HPF30MHZ', 'rcuID':2, 'array_type':'LBA', 'bw':100000000.},   #2
            {'mode':'LBH_HPF10MHZ', 'rcuID':3, 'array_type':'LBA', 'bw':100000000.},   #3
            {'mode':'LBH_HPF30MHZ', 'rcuID':4, 'array_type':'LBA', 'bw':100000000.},   #4
            {'mode':'HBA_110_190MHZ', 'rcuID':5, 'array_type':'HBA', 'bw':100000000.}, #5
            {'mode':'HBA_170_230MHZ', 'rcuID':6, 'array_type':'HBA', 'bw':100000000.}, #6
            {'mode':'HBA_210_290MHZ', 'rcuID':7, 'array_type':'HBA', 'bw':100000000.}] #7

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

xyz = read_ant_xyz(ant_field_file, rcuInfo, rcumode, 'IE')

obs=ephem.Observer()
obs.long=lon
obs.lat=lat
obs.elevation=float(elev)
obs.epoch=2000.0
freq = numpy.array([50e6])
time_0 = calendar.timegm( time.strptime("2017-12-21 09:00:00", "%Y-%m-%d %H:%M:%S") )
time_1 = calendar.timegm( time.strptime("2017-12-21 20:00:00", "%Y-%m-%d %H:%M:%S") )
image_num = 0
plt.ion()
while time_0 < time_1:
    time_utc = time.strftime("%Y/%m/%d %H:%M", time.gmtime(time_0)) #UTC
    time_ist = time.strftime("%Y/%m/%d %H:%M", time.gmtime(time_0+ 60.0*60.0)) #Irish Standard Time
    obs.date = time_utc
    sun = ephem.Sun()
    sun.compute(obs)
    src=sun
    src._ra=sun.ra
    src._dec=sun.dec
    src.compute(obs)
    
    ####################################
    #       UVW calc function
    uvw=xyz2uvw(xyz, src, obs, freq[0])
    U=uvw[:,:,0]
    V=uvw[:,:,1]

    """
    pixels=opts.pixels
    px=[pixels,pixels]
    fov=numpy.pi    #Field of View in radians
    res=fov/px[0]   #pixel resolution
    """
    fig = plt.figure(figsize=(10,10))
    plt.plot(U, V, '.')
    plt.xlabel('u ($\lambda$)')
    plt.ylabel('v ($\lambda$)')
    plt.title('UV-coverage at %s MHz for IE613 LBAs. Source: Sun @ %s IST' % (round(freq[0]/1e6, 1), time_utc))
    plt.axis([-10, 10, -10, 10])
    plt.show( )
    plt.savefig('/Users/eoincarley/LOFAR/data/IE613/image_'+str(format(image_num, '03'))+'.png')   # save the figure to file
    plt.close(fig)
    image_num+=1
    time_0 = time_0+15.0*60.0

#ffmpeg -y -r 20 -i image_%03d.png -vb 50M IE613_uv_coverage_sun.mpg    


