#!/usr/bin/python
"""
########################################################################

    File:
        model_sun_pos_uv_beam.py

    Description:
        Plot the IE613 station antenna positions, position of the sun in the sky
        at a particular time, the uv-coverage for the sun direction and the resulting
        station beam for that direction.

    Disclaimer:

    Notes:
        Adapted in parts from dft_acc.py for the production of all sky images.

    Examples:
         model_sun_pos_uv_beam.py /Users/eoincarley/LOFAR/data/IE613/2017.08.21/AntennaField.conf

    Testing:

    Author:
        Eoin Carley, School of Physics, TCD

    Created/Updated:
        20-Nov-2017 -- v0.1, Eoin Carley

########################################################################
"""

import numpy
import math
import numpy as np
import ephem
import pymap3d as pm
import pylab
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
import sys,os
import struct
import time
import datetime
import calendar
import pdb

######################
#  
#     Functions
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
    #return numpy.sum(d*numpy.exp(-2.*numpy.pi*1j*((u*k) + (v*l))))
    #psf:
    return numpy.sum(numpy.exp(-2.*numpy.pi*1j*((u*k) + (v*l))))
def dftImage(d,uvw,pxPlot,res,mask=False):
    """return a DFT image"""
    start_time = time.time()
    nants=uvw.shape[0]
    im=numpy.zeros((pxPlot[0],pxPlot[1]),dtype=complex)

    u=uvw[:,:,0].reshape(-1)
    v=uvw[:,:,1].reshape(-1)

    # u,v = x,y / lambda implies theta = 1 / u,v

    uvMax = np.max([u, v])
    u /= uvMax
    v /= uvMax
    #w=uvw[:,:,2]

    res = 1. / (2. * uvMax) # 2x for Nyquist sampling
    px = int(2. * uvMax / res)
    
    # Temp debug
    pxPlot = [px, px]
    imSize = px + 11 # 5 padding on each side + 1 for rounding
    centralRef = imSize / 2

    uvGrid = np.zeros([imSize, imSize])

    sampleGrid = np.mgrid[-4:5, -4:5].reshape(1, 2, -1).astype(float) * res
    sampleCache = sampleGrid.copy()

    gaussCoeff = 1. / (2. * np.square(res))
    guassMatrix = np.array([[gaussCoeff, 0], [0, gaussCoeff]])
    #sampledPoints: (1,2,-1)
    gaussLambda = lambda sampledPoints: np.sum(np.exp(-1. * np.dot(guassMatrix, np.square(sampledPoints))), axis = (0,1))
    
    for sampleCoord in np.vstack([u, v]).T[:, np.newaxis, :]:
        offsets = res - (sampleCoord % res)[..., np.newaxis]
        sampleCache += offsets
        
        # Gaussian
        pointSamples = gaussLambda(sampleCache)
        
        # Rectangular
        #print(np.argwhere(sampleCache.reshape(1, 2, -1) == 0))
        #pointSamples = np.zeros([81])
        #pointSamples[30] = 1.

        sampleIndex = ((sampleCoord[..., np.newaxis] + offsets + sampleCache.reshape(1, 2, -1)) / (res) + centralRef).astype(int)

        uvGrid[sampleIndex[0, 0, :], sampleIndex[0, 1, :]] += pointSamples

        sampleCache -= offsets
    
    uvGrid = uvGrid[5:-6, 5:-6] # Remove padding
    im = np.fft.fft2(uvGrid)
    im = np.fft.fftshift(np.abs(im.real)) # Abs of FFT == DFT?

    print(res, imSize, res * imSize)
    
    if mask:
        sampleGrid = np.array(np.meshgrid(np.arange(im.shape[0]), np.arange(im.shape[1]))) - np.max(im.shape) / 2
        maskRad = np.sqrt(np.square(sampleGrid[0] * 1.22 * res) + np.square(sampleGrid[1] * 1.22 * res)) # 1.22 factor to convert to radians

        print(np.max(maskRad), np.max(sampleGrid), res, uvMax)
        maskRad = maskRad > (2. * np.pi)
        #print(np.argwhere(maskRad))
        im[maskRad] = 0.

        im = im[~np.all(im == 0, axis = 0)]
        im = im[:, ~np.all(im.T == 0, axis = 1)]

    print("DFT Image ended after {0}".format(time.time() - start_time))
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

print("Initialisation of functions Complete")

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('Plot XYZ of IE613 LOFAR Antennas')
    o.set_description(__doc__)
    opts, args = o.parse_args(sys.argv[1:])  

#
#########################################################   

if args:
    ant_field_file=args[0]
else:
    if not os.path.exists("./ant_field_file.conf"):
        import urllib
        urllib.urlretrieve("https://raw.githubusercontent.com/griffinfoster/SWHT/master/SWHT/data/LOFAR/StaticMetaData/IE613-AntennaField.conf", "./ant_field_file.conf")
    ant_field_file = "./ant_field_file.conf"


ant_arr_file = "/Users/eoincarley/LOFAR/data/IE613/AntennaArrays.conf"

if not os.path.exists(ant_arr_file):
    import urllib
    urllib.urlretrieve("https://raw.githubusercontent.com/griffinfoster/SWHT/master/SWHT/data/LOFAR/StaticMetaData/AntennaArrays/AntennaArrays_Int.conf", "AntennaArrays.conf")
    ant_arr_file = "./AntennaArrays.conf"

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

print("Parameter Initialisation Complete")

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

print("Antenna Array Read")
#------------------------------------------------
#     Define local geographic coords and time
#
obs = ephem.Observer()
obs.lon = math.radians(float(lon))
obs.lat = math.radians(float(lat))
obs.elevation=float(elev)
obs.epoch=2000.0
freq = numpy.array([50e6])
time_0 = calendar.timegm( time.strptime("2017-12-21 09:00:00", "%Y-%m-%d %H:%M:%S") )
time_1 = calendar.timegm( time.strptime("2017-12-21 20:00:00", "%Y-%m-%d %H:%M:%S") )
image_num = 0
plt.ion()

print("Coords / time processed")
#----------------------------------------------------
#      Get IE613 antenna positions and projections
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

zplane = 5.0769e3
LBAZproj = [zplane for i in LBAZ]
HBAZproj = [zplane  for i in HBAZ]
xplane = 3.80155e3
LBAXproj = [xplane for i in LBAX]
HBAXproj = [xplane  for i in HBAX]
yplane = -5.288e2
LBAYproj = [yplane for i in LBAY]
HBAYproj = [yplane  for i in HBAY]

u = np.linspace(0, 2.0*np.pi, 200)
v = np.linspace(0, np.pi, 200)
xsphere0 = 7e-3*np.outer(np.cos(u), np.sin(v)) 
ysphere0 = 7e-3*np.outer(np.sin(u), np.sin(v)) 
zsphere0 = 7e-3*np.outer(np.ones(np.size(u)), np.cos(v)) 

print("Projections Processed")
#------------------------------------------
#       For the instrument beam image
pixels=128       #opts.pixels
px=[pixels,pixels]
fov=numpy.pi    #Field of View in radians
res=fov/px[0]   #pixel resolution
dummy = 0.0   

while time_0 < time_1:

    print("Starting Loop")
    time_utc = time.strftime("%Y/%m/%d %H:%M", time.gmtime(time_0)) #UTC (Irish time in Winter)
    time_ist = time.strftime("%Y/%m/%d %H:%M", time.gmtime(time_0 + 60.0*60.0)) # Irish Standard Time (Irish time in Summer)
    obs.date = time_utc
    sun = ephem.Sun()
    sun.compute(obs)
    print("Sun Obtained")

    fig = plt.figure(figsize=(18, 10))

    ####################################
    #       Plot IE613 and sun
    src=sun
    src._ra=sun.ra
    src._dec=sun.dec
    src.compute(obs)
    ax = plt.subplot2grid((2, 2), (0, 0), rowspan=2, projection='3d')
    #ax = fig.add_subplot(3, 2, 1, projection='3d', rowspan=2, colspan=2)
    print("Plotted for IE613")
    #----------------------------------------------------
    #    Define solar ephem object to get local alt-az     
    #
    alt = math.degrees(sun.alt)
    az = math.degrees(sun.az) 

    sun_ecef = pm.aer2ecef(az, alt, 150.0, 53.09472, -7.9213880, 75.0, deg=True) #pm.geodetic2ecef(53.09472, -7.921388, 75.0, deg=True)

    lba_to_sunx = [ meanx, sun_ecef[0]/1000.0 ]
    lba_to_suny = [ meany, sun_ecef[1]/1000.0 ]
    lba_to_sunz = [ meanz, sun_ecef[2]/1000.0 ]
    print("Sun object created")
    #-----------------------------------------  
    #        Make sphere for the sun
    xsphere = xsphere0 + lba_to_sunx[1]
    ysphere = ysphere0 + lba_to_suny[1]
    zsphere = zsphere0 + lba_to_sunz[1]


    #----------------------------------------- 
    #        Plot and format
    ax.plot(LBAX, LBAY, LBAZ, 'bo', label='LBAs (ECEF coords)', zorder=-1)
    ax.plot(HBAX, HBAY, HBAZ, 'o', color='salmon', label='HBAs (ECEF coords)', zorder=-1)
    ax.plot_surface(xsphere, ysphere, zsphere, color='y', linewidth=0.01, zorder=1)
    ax.plot(lba_to_sunx, lba_to_suny, lba_to_sunz, 'g', zorder=2)
    print("Station Plotted")
    #----------------------------------------- 
    #    Plot projections on the XYZ planes
    #
    lba_to_sunzproj = np.full_like(lba_to_sunz, zplane)
    lba_to_sunxproj = np.full_like(lba_to_sunx, xplane)
    lba_to_sunyproj = np.full_like(lba_to_suny, yplane)

    zsphereproj = np.full_like(zsphere, zplane)
    xsphereproj = np.full_like(xsphere, xplane)
    ysphereproj = np.full_like(ysphere, yplane)

    ax.plot(LBAXproj, LBAY, LBAZ, '.', color='lightgray', zorder=-2)
    ax.plot(HBAXproj, HBAY, HBAZ, '.', color='lightgray', zorder=-2)
    ax.plot(lba_to_sunxproj, lba_to_suny, lba_to_sunz, color='lightgray', zorder=-2)
    ax.plot_surface(xsphereproj, ysphere, zsphere, color='lightgray', linewidth=0.01, zorder=-2)

    ax.plot(LBAX, LBAY, LBAZproj, '.', color='lightgray', zorder=-2)
    ax.plot(HBAX, HBAY, HBAZproj, '.', color='lightgray', zorder=-2)
    ax.plot(lba_to_sunx, lba_to_suny, lba_to_sunzproj, color='lightgray', zorder=-2)
    ax.plot_surface(xsphere, ysphere, zsphereproj, color='lightgray', linewidth=0.01, zorder=-2)

    ax.plot(LBAX, LBAYproj, LBAZ, '.', color='lightgray', zorder=-2)
    ax.plot(HBAX, HBAYproj, HBAZ, '.', color='lightgray', zorder=-2)
    ax.plot(lba_to_sunx, lba_to_sunyproj, lba_to_sunz, color='lightgray', zorder=-2)
    ax.plot_surface(xsphere, ysphereproj, zsphere, color='lightgray', linewidth=0.01, zorder=-2)
    print("Projecions plotted")

    #----------------------------------------- 
    #        Plot format
    #
    box_width = 0.2
    xbox_cen = 3.80155e3
    ybox_cen = -5.2889e2
    zbox_cen = 5.0769e3
    ax.set_xlim3d(xbox_cen, xbox_cen+box_width)
    ax.set_ylim3d(ybox_cen -box_width, ybox_cen)
    ax.set_zlim3d(zbox_cen, zbox_cen+box_width)
    plt.title('IE613 solar alt: %s$^{\circ}$, az: %s$^{\circ}$ @ %s IST' % (round(alt), round(az), time_utc))
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    plt.gca().set_position([0, 0.0, 0.52, 1.0])
    #ax.legend()
    ax.view_init(elev=16, azim=-55.0)
    print("Plot reoritentatied")

    ####################################
    #       Plot UV coverage
    uvw=xyz2uvw(lba_xyz, src, obs, freq[0])
    U=uvw[:,:,0]
    V=uvw[:,:,1]
    ax0 = plt.subplot2grid((2, 2), (0, 1))
    ax0.set_aspect('equal')
    plt.plot(U, V, '.')
    plt.xlabel('u ($\lambda$)')
    plt.ylabel('v ($\lambda$)')
    plt.title('UV-coverage at %s MHz for IE613 LBAs. Source: Sun' % (round(freq[0]/1e6, 1)))
    plt.axis([-10, 10, -10, 10])
    print("UVW Plotted")

    ####################################
    #       Plot instrument beam
    ax1 = plt.subplot2grid((2, 2), (1, 1))
    img = dftImage(dummy, uvw, px, res, mask=True)
    img=img.real
    img=numpy.log(img)

    ax1 = plt.subplot2grid((2, 2), (1, 1))
    im = plt.imshow(img, extent = (-1.0*pixels, pixels, pixels, -1.0*pixels), vmin = 0.51, vmax= 10.5)
    plt.title('IE613 %s MHz LBA station beam. Source: Sun' % (round(freq[0]/1e6, 1)))
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.gca().xaxis.set_major_locator(plt.NullLocator())

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label('log$_{10}$(I$_{beam}$)')
    print("Beam Plotted")
    plt.show( )
    plt.savefig('./sun_pos_uv_projection_'+str(format(image_num, '03'))+'.png')   # save the figure to file
    #pdb.set_trace()    
    plt.close(fig)
    print(image_num)
    image_num+=1
    time_0 = time_0+5.0*60.0
    raw_input()


    #ffmpeg -y -r 20 -i image_%03d.png -vb 50M IE613_uv_coverage_sun.mpg    




