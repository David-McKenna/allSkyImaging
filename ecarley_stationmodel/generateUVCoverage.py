import matplotlib

matplotlib.use('Agg')

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
import pytz
import skimage.measure as skim
import copy
import multiprocessing as mp

import subprocess

ffmpegLoc = "/cphys/ugrad/2015-16/JF/MCKENND2/Downloads/ffmpeg-4.0.2-64bit-static/ffmpeg"

global rcuInfo
rcuInfo = [ {'mode':'OFF', 'rcuID':0, 'array_type':'LBA', 'bw':100000000., 'refFreq': 0},            #0
			{'mode':'LBL_HPF10MHZ', 'rcuID':1, 'array_type':'LBA', 'bw':100000000., 'refFreq': 50e6},   #1
			{'mode':'LBL_HPF30MHZ', 'rcuID':2, 'array_type':'LBA', 'bw':100000000., 'refFreq': 50e6},   #2
			{'mode':'LBH_HPF10MHZ', 'rcuID':3, 'array_type':'LBA', 'bw':100000000., 'refFreq': 50e6},   #3
			{'mode':'LBH_HPF30MHZ', 'rcuID':4, 'array_type':'LBA', 'bw':100000000., 'refFreq': 50e6},   #4
			{'mode':'HBA_110_190MHZ', 'rcuID':5, 'array_type':'HBA', 'bw':100000000., 'refFreq': 150e6}, #5
			{'mode':'HBA_170_230MHZ', 'rcuID':6, 'array_type':'HBA', 'bw':100000000., 'refFreq': 200e6}, #6
			{'mode':'HBA_210_290MHZ', 'rcuID':7, 'array_type':'HBA', 'bw':100000000., 'refFreq': 225e6}] #7
			
randWalk = [ 0,  2,  7,  5,  7,  9,  5, 15,  3, 11, 15, 15,  2,  1,  3,  6,  6,
        6, 15, 11, 15,  2,  0,  0,  6, 11, 14, 14,  7,  7, 15, 12,  7,  4,
        1,  2,  0, 10,  9,  5, 10, 13, 14,  1,  0,  2,  4,  5,  2,  0, 11,
        8, 13,  2,  2,  1,  2,  9,  3, 15,  8, 14, 12, 15,  1,  8,  8,  0,
        5,  6,  9, 14,  5,  8, 12,  0,  1,  4,  0, 15,  1, 12, 14,  9,  0,
        0,  0, 12, 14, 12, 13,  4,  8,  8, 13,  8]


Effelsberg_elements_20091110 = [1,4,13,15,11,9,14,1,15,0,8,2,11,3,14,0,2,4,3,0,0,2,12,12,12,12,15,11,14,15,7,5,1,0,3,10,1,11,0,12,12,1,6,7,0,10,9,6,15,14,11,7,2,0,7,12,15,8,13,3,7,6,3,15,11,1,4,11,8,1,8,15,4,0,5,6,12,0,12,15,3,7,14,8,3,12,12,2,9,8,14,2,5,6,12,0]
Generic_International_Station_20091110 = [15,0,15,3,9,15,14,2,0,3,4,14,10,8,5,15,12,0,2,11,3,12,12,1,5,4,4,8,6,3,0,5,3,11,3,2,8,15,13,8,3,2,9,1,14,8,8,0,12,13,0,11,15,3,12,3,13,3,10,5,0,10,1,6,4,10,3,15,3,14,0,12,0,7,0,12,7,3,13,0,7,3,15,4,14,4,3,8,4,9,12,0,14,9,3,11]
custom_debug = (range(16) * 8)[:len(Effelsberg_elements_20091110)]

global activationSchemes
activationSchemes = [Effelsberg_elements_20091110, Generic_International_Station_20091110, custom_debug, randWalk]

def eq2top_m(ha, dec):
	"""Return the 3x3 matrix converting equatorial coordinates to topocentric
	at the given hour angle (ha) and declination (dec)."""
	sin_H, cos_H = np.sin(ha), np.cos(ha)
	sin_d, cos_d = np.sin(dec), np.cos(dec)
	
	zero = np.zeros_like(ha)

	mapMatrix =  np.array([[    sin_H    ,       cos_H  ,       zero  ],
						[ -sin_d*cos_H,   sin_d*sin_H,      cos_d  ],
						[  cos_d*cos_H,  -cos_d*sin_H,      sin_d  ]])

	if len(mapMatrix.shape) == 3: 
		mapMatrix = mapMatrix.transpose([2, 0, 1])
	return mapMatrix

def get_baseline(i, j, src, obs):
	"""Return the baseline corresponding to i,j""" 
	
	bl = j - i
	
	try:
		if src.alt < 0:
			print(src.date, src.alt, obs.date)
			raise PointingError('Phase center below horizon')
		m=src.map
	except(AttributeError):
		# Attribute error is consistently reached -- why?
		ra,dec = src._ra,src._dec
		#opposite HA since the we want to move the source at zenith away to phase to the original zenith source
		m = eq2top_m(ra-obs.sidereal_time(), dec)
		#normal HA
		#m = eq2top_m(obs.sidereal_time() - ra, dec)
	return np.dot(m, bl).transpose()

def gen_uvw(i, j, src, obs, f):
	"""Compute uvw coordinates of baseline relative to provided FixedBody"""
	x,y,z = get_baseline(i,j,src,obs)

	afreqs = np.reshape(f, (1,f.size))
	afreqs = afreqs/ephem.c #1/wavelength

	if len(x.shape) == 0: 
		return np.array([x*afreqs, y*afreqs, z*afreqs])
	
	x.shape += (1,); y.shape += (1,); z.shape += (1,)

	return np.array([np.dot(x,afreqs), np.dot(y,afreqs), np.dot(z,afreqs)])

def xyz2uvw(xyz, src, obs, f):
	"""Return an array of UVW values"""
#	uvw=np.zeros((xyz.shape[0],xyz.shape[0],3))
	uvw=np.zeros((xyz.shape[0],xyz.shape[0] - 1,3))

	for i in range(xyz.shape[0]):
		refIdx = 0
		for j in range(xyz.shape[0]):
			if i==j: continue
#			uvw[i,j]=gen_uvw(xyz[i], xyz[j], src, obs, f)[:,0,0]
			uvw[i,refIdx]=gen_uvw(xyz[i], xyz[j], src, obs, f)[:,0,0]
			refIdx += 1

	return uvw

def dft2(d,k,l,u,v):
	"""compute the 2d DFT for position (k,l) based on (d,uvw)"""
	#return np.sum(d*np.exp(-2.*np.pi*1j*((u*k) + (v*l))))
	#psf:

	return np.sum(np.exp(-2.*np.pi*1j*((u*k) + (v*l))))

def dftImage(d,uvw,px,res,mask=False):
	"""return a DFT image"""

	nants=uvw.shape[0]
	
	im=np.zeros((px[0],px[1]),dtype=complex)

	mid_k=int(px[0]/2.)
	mid_l=int(px[1]/2.)

	u=uvw[:,:,0]
	v=uvw[:,:,1]
	w=uvw[:,:,2]

	u/=mid_k
	v/=mid_l

	start_time=time.time()

	# Speedup with numpy arrays / numba port in future
	for k in range(px[0]):
		for l in range(px[1]):
			im[k,l]=dft2(d,(k-mid_k),(l-mid_l),u,v)
			if mask:        #mask out region beyond field of view
				rad=(((k-mid_k)*res)**2 + ((l-mid_l)*res)**2)**.5
				if rad > mid_k*res: im[k,l]=0

	print("DFT Image ended after {0}".format(time.time() - start_time))

	# Debug after removing (0,0) UV elements
	im = np.abs(im)
	return im

def fftImage(d,uvw,pxPlot,res,mask=False, useDVar = False, method = 'gaus'):
	"""return a FFT image"""

	print("Init FFT")
	start_time = time.time()
	nants=uvw.shape[0]
	im=np.zeros((pxPlot[0],pxPlot[1]),dtype=complex)

	u=uvw[:,:,0].reshape(-1)
	v=uvw[:,:,1].reshape(-1)

	print(u.size, v.size)
	print("Reshaped Arrays")
	# u,v = x,y / lambda implies theta = 1 / u,v

	uvMax = np.max([u, v])

	u /= uvMax
	v /= uvMax

	res = 1. / (2. * uvMax) # 2x for Nyquist sampling
	px = int(2. * uvMax / res)
	
	print(np.max(u), np.max(v), uvMax, res)
	# Temp debug
	pxPlot = [px, px]
	imSize = px + 11 # 5 padding on each side + 1 for rounding
	centralRef = imSize / 2

	uvGrid = np.zeros([imSize, imSize])
	sampleCache = np.mgrid[-4:5, -4:5].reshape(1, 2, -1).astype(float) * res

	gaussCoeff = 1. / (2. * np.square(res))
	guassMatrix = np.array([[gaussCoeff, 0], [0, gaussCoeff]])
	# Sample points expects shape (1,2,-1), where -1 denotes the amount of points to be sampled.
	gaussLambda = lambda sampledPoints: np.sum(np.exp(-1. * np.dot(guassMatrix, np.square(sampledPoints))), axis = (0,1))

	for sampleCoord in np.vstack([u, v]).T[:, np.newaxis, :]:
		offsets = res - (sampleCoord % res)[..., np.newaxis]
		sampleCache += offsets
		
		if method == 'gaus':
			pointSamples = gaussLambda(sampleCache)
		
		elif method == 'rect':
			pointSamples = np.zeros([81])
			pointSamples[30] = 1. # Why is it so low?

		elif method == 'rectgaus':
			pointSamples = np.zeros([81])
			pointSamples[30] = gaussLambda(offsets)

		else:
			print("Unknown method! Falling back to Gaussian.")
			pointSamples = gaussLambda(sampleCache)

		sampleIndex = ((sampleCoord[..., np.newaxis] + offsets + sampleCache.reshape(1, 2, -1)) / (res) + centralRef).astype(int)
		uvGrid[sampleIndex[0, 0, :], sampleIndex[0, 1, :]] += pointSamples

		sampleCache -= offsets

	print("Gridding Complete, size ", uvGrid.size, res, px, im.shape)
	uvGrid = uvGrid[5:-6, 5:-6] # Remove padding

	# Sometimes we get absolutely huge arrays (Thanks HBAs.). Let's resample these down to a more reasonable value of ~512 ** 2
	if uvGrid.shape[0] > 728:
		print("Attempting to shrink size from {0}".format(uvGrid.shape[0]))
		imDiv = min(1, int(uvGrid.shape[0] / 512) - 1)
		uvGrid = skim.block_reduce(uvGrid.real, (imDiv, imDiv), np.sum)	

	oddOffset = 1 - (im.shape[0] % 2)
	uvGrid = uvGrid[5:-6 + oddOffset, 5:-6 + oddOffset] # Remove padding but ensure an odd result
	im = np.fft.fft2(uvGrid)
	im = np.fft.fftshift(np.abs(im.real)) # Abs of FFT == DFT?

	if useDVar:
		im *= d

	# Needs retesting after image resizing implemented
	if mask:
		sampleGrid = np.array(np.meshgrid(np.arange(im.shape[0]), np.arange(im.shape[1]))) - np.max(im.shape) / 2
		maskRad = np.sqrt(np.square(sampleGrid[0] * 1.22 * res) + np.square(sampleGrid[1] * 1.22 * res)) # 1.22 factor to convert to radians

		print(np.max(maskRad), np.max(sampleGrid), res, uvMax)
		maskRad = maskRad > (2. * np.pi)
		#print(np.argwhere(maskRad))
		im[maskRad] = 0.

		im = im[~np.all(im == 0, axis = 0)]
		im = im[:, ~np.all(im.T == 0, axis = 1)]

	print("FFT Image ended after {0} for size {1}".format(time.time() - start_time, uvGrid.shape))
	return im

def read_ant_xyz(antFieldFile, rcumode, station, activationDeltas = [None, None]):

	# Grab HBA offsets iff needed, else initialise arrays of 0s
	if activationDeltas[0] is not None:
		activation, hbaDeltas = activationDeltas
		activationArr = activationSchemes[activation]
		deltaArr = []
		print(activationArr, hbaDeltas)

		with open(hbaDeltas) as deltaRef:
			readNext = False
			for line in deltaRef:
				if line.startswith(']'):
					readNext = False

				if readNext:
					lineCull = line[:-2].split(' ')
					deltaArr.append([float(delta) for delta in line[:-2].split(' ')[4:] if delta != ''])

				if line.startswith("(0,15)"):
					readNext = True
		deltaArr = np.vstack(deltaArr).astype(float)

	else:
		activationArr = np.zeros([96], dtype = int)
		deltaArr = np.zeros([16, 3])

	
	with open(antFieldFile) as fh:	
		readXYZ=False
		readArrayHDR=False # These two are never set it True?
		readArray=False # These two are never set to True?
		
		linecnt=0
		nelem=0
		ants=[]
	
		for line in fh:
			if line.lower().startswith(rcuInfo[rcumode]['array_type'].lower()):
				readXYZ=True
				continue

			if readXYZ:
				arr_x=float(line.split(' ')[2])
				arr_y=float(line.split(' ')[3])
				arr_z=float(line.split(' ')[4])
				readXYZ=False
				readArrayHDR=True
				continue

			# Next two if statements are never used?
			if readArrayHDR:
				if station == 'IE': 
					nelem=eval(line.split()[0])[1]+1
				else:    
					nelem=int(line.split(' ')[0])    # This code expects '96' here. But the IE613 file has '(0,95)'. 
				readArrayHDR=False
				readArray=True
				continue

			if readArray and linecnt < nelem:
				cl=' '.join(line.split())
				ants.append(map(float,cl.split(' ')))
	
				linecnt+=1         

	xyz = [[a[0] + arr_x + deltaArr[activationArr[idx]][0],
			a[1] + arr_y + deltaArr[activationArr[idx]][1],
			a[2] + arr_z + deltaArr[activationArr[idx]][2]]

				for idx, a in enumerate(ants)]

	xyz=np.array(xyz) # Earth Centered, Earth Fixed (ECEF) Coords.

	return xyz

def initialiseLocalFiles(args):
	if args:
		antFieldFile=args[0]
	else:
		antFieldFile = "./antFieldFile.conf"
		if not os.path.exists(antFieldFile):
			import urllib
			urllib.urlretrieve("https://raw.githubusercontent.com/griffinfoster/SWHT/master/SWHT/data/LOFAR/StaticMetaData/IE613-AntennaField.conf", antFieldFile)
	
	
	antArrFile = "./stations_etrs.txt" # TODO: Find public replacement for station lat/lon
	
	hbaDeltas = "./IE613-iHBADeltas.conf"
	
	if not os.path.exists(hbaDeltas):
		import urllib
		urllib.urlretrieve("https://raw.githubusercontent.com/griffinfoster/SWHT/master/SWHT/data/LOFAR/StaticMetaData/IE613-iHBADeltas.conf", "IE613-iHBADeltas.conf")
	
	return antFieldFile, antArrFile, hbaDeltas

def getLatLongElevation(antArrFile, arrayName):
	with open(antArrFile) as fh:
		for line in fh:
			if line.lower().startswith(arrayName.lower()):
				lon=float(line.split('\t')[4])
				lat=float(line.split('\t')[5])
				elev=float(line.split('\t')[6])
				print(lon, lat, elev)
				break

	return [lon, lat, elev]

def antProc(antFieldFile, antArrFile, rcuMode, obsTime, hbaActivation = [None, None]):
	# TODO: Generalise to any source, not just the sun. Maybe other stations.
	# May be harder than expected, not very familiar with pyehepm, but porting the entire script over to
	# 	astropy may be a challenge.

	antType = rcuInfo[rcuMode]['array_type']
	antLatLongEle = np.array(getLatLongElevation(antArrFile, 'IE613' + antType), dtype = float)

	freq = np.array([rcuInfo[rcuMode]['refFreq']])

	obs = ephem.Observer()
	sunObj = ephem.Sun()

	obs.lon, obs.lat = math.radians(antLatLongEle[1]), math.radians(antLatLongEle[0])
	obs.elevation = antLatLongEle[2]
	obs.pressure = 0

	obs.date = obsTime

	sunObj.compute(obs)

	initTime = str(obs.previous_rising(sunObj))
	endTime = str(obs.next_setting(sunObj))
	print(initTime, endTime)

	initTime = initTime[:-6] + initTime[-6:-4] + '0:00'
	endTime = endTime[:-6] + endTime[-6:-4] + '0:00'

	print(initTime, endTime)
	initTime = calendar.timegm( time.strptime(initTime, "%Y/%m/%d %H:%M:%S") )
	endTime = calendar.timegm( time.strptime(endTime, "%Y/%m/%d %H:%M:%S") )

	initTime = datetime.datetime.utcfromtimestamp(initTime) + datetime.timedelta(minutes = 30) # Have a 30 minute buffer after/before reaching the horizon
	endTime = datetime.datetime.utcfromtimestamp(endTime) - datetime.timedelta(minutes = 30)

	sunObj.compute(obs)
	sourceObj = sunObj
	sourceObj._ra = sunObj.ra
	sourceObj._dec = sunObj.dec
	sourceObj.compute(obs)

	xyz = read_ant_xyz(antFieldFile, rcuMode, 'IE', hbaActivation)

	uvw = xyz2uvw(xyz, sourceObj, obs, freq[0])

	return uvw, freq, initTime, endTime, xyz, obs, sunObj, sourceObj

def plotConsts(xyz, planeConsts):
	x = np.array([coords[0] / 1e3 for coords in xyz])
	y = np.array([coords[1] / 1e3 for coords in xyz])
	z = np.array([coords[2] / 1e3 for coords in xyz])

	meanX = np.mean(x)
	meanY = np.mean(y)
	meanZ = np.mean(z)

	xproj = np.full_like(x, planeConsts[0])
	yproj = np.full_like(y, planeConsts[1])
	zproj = np.full_like(z, planeConsts[2])

	u = np.linspace(0, 2. * np.pi, 200)
	v = np.linspace(0, np.pi, 200)

	xSphereInit = 7e-3*np.outer(np.cos(u), np.sin(v)) 
	ySphereInit = 7e-3*np.outer(np.sin(u), np.sin(v)) 
	zSphereInit = 7e-3*np.outer(np.ones(np.size(u)), np.cos(v))

	return x, y, z, [[meanX, meanY, meanZ], [xproj, yproj, zproj], [xSphereInit, ySphereInit, zSphereInit], planeConsts]

def processPlot(xyz, obs, sunObj, refTimeUtc, plotConstants, sun_ecef, ax, antType = 'LBA', plotSuns = True):
	x, y, z = xyz
	meanArr, projArr, sphereInit, planeConsts = plotConstants

	meanX, meanY, meanZ = meanArr
	xProj, yProj, zProj = projArr
	xSphereInit, ySphereInit, zSphereInit = sphereInit
	xplane, yplane, zplane = planeConsts

	obs.date = refTimeUtc
	sunObj.compute(obs)

	to_sunx = [ meanX, sun_ecef[0]/1000.0 ]
	to_suny = [ meanY, sun_ecef[1]/1000.0 ]
	to_sunz = [ meanZ, sun_ecef[2]/1000.0 ]

	xSphere = xSphereInit + to_sunx[1]
	ySphere = ySphereInit + to_suny[1]
	zSphere = zSphereInit + to_sunz[1]

	to_sunzproj = np.full_like(to_sunz, zplane)
	to_sunxproj = np.full_like(to_sunx, xplane)
	to_sunyproj = np.full_like(to_suny, yplane)
	
	zSphereProj = np.full_like(zSphere, zplane)
	xSphereProj = np.full_like(xSphere, xplane)
	ySphereProj = np.full_like(ySphere, yplane)

	if plotSuns:
		ax.plot_surface(xSphere, ySphere, zSphere, color='y', linewidth=0.01, zorder=-2)
		ax.plot_surface(xSphereProj, ySphere, zSphere, color='lightgray', linewidth=0.01, zorder=-2)
		ax.plot_surface(xSphere, ySphere, zSphereProj, color='lightgray', linewidth=0.01, zorder=-2)
		ax.plot_surface(xSphere, ySphereProj, zSphere, color='lightgray', linewidth=0.01, zorder=-2)

	if antType == 'LBA':
		ax.plot(x, y, z, 'o', color = 'blue', label='LBAs (ECEF coords)', zorder=-1)
	else:
		ax.plot(x, y, z, 'o', color = 'salmon', label='HBAs (ECEF coords)', zorder=-1)

	ax.plot(to_sunx, to_suny, to_sunz, color = 'g', zorder=2)

	ax.plot(xProj, y, z, '.', color='lightgray', zorder=-2)
	ax.plot(to_sunxproj, to_suny, to_sunz, color='lightgray', zorder=-2)

	ax.plot(x, y, zProj, '.', color='lightgray', zorder=-2)
	ax.plot(to_sunx, to_suny, to_sunzproj, color='lightgray', zorder=-2)

	ax.plot(x, yProj, z, '.', color='lightgray', zorder=-2) # Determine why this projection appears to go below the array. Side effect of the solar position approximation?
	ax.plot(to_sunx, to_sunyproj, to_sunz, color='lightgray', zorder=-2)

	ax.set_aspect('equal')
	ax.set_anchor('SW')

	return ax

def processUVCoverage(idx, sunObj, obs, xyz, sourceObj, obsArr, freq, colIdx, subplotTuple, antType):

	sourceObj=sunObj
	sourceObj._ra=sunObj.ra
	sourceObj._dec=sunObj.dec
	sourceObj.compute(obs)
	uvw=xyz2uvw(xyz, sourceObj, obs, freq[0])
	U=uvw[:,:,0]
	V=uvw[:,:,1]

	print(subplotTuple, colIdx)
	ax0 = plt.subplot2grid(subplotTuple, (0, colIdx), fig = plt.figure(idx))
	ax0.plot(U, V, '.')
	ax0.set_xlabel('u ($\lambda$)')
	ax0.set_ylabel('v ($\lambda$)')
	ax0.set_title('UV-coverage at {0} MHz for IE613 {1}s. Source: Sun'.format(round(freq[0]/1e6, 1), antType))

	if antType == 'LBA':
		limVar = 10
	if antType == 'HBA':
		limVar = 30
	limVar = [-1.1 * limVar, 1.1 * limVar]

	ax0.set_xlim(limVar)
	ax0.set_ylim(limVar)
	ax0.set_aspect('equal')

	print("UVW Plotted")

	return uvw

def beamGenerator(idx, dummy, uvw, px, res, fftType, maskVar):
	if fftType is None:
		img = dftImage(dummy, uvw, px, res, mask = maskVar)
	else:
		img = fftImage(dummy, uvw, px, res, mask = maskVar, method = fftType)

	print("Exiting process for figure {0}".format(idx))
	return img

def beamPlot(idx, img, dummy, pixels, freq, outputFolder, fftType, colIdx, subplotTuple, antType):
	fig = plt.figure(idx)

	img=img.real
	#img=np.log(img)
	ax1 = plt.subplot2grid(subplotTuple, (1, colIdx), fig = fig, projection = '3d')

	ax1 = fig.gca()
	ax1.w_xaxis.set_pane_color((0.75,0.75,0.75,1))
	ax1.w_yaxis.set_pane_color((0.75,0.75,0.75,1))
	ax1.zaxis.pane.fill = False
#	im = ax1.imshow(img, extent = (-1.0*pixels, pixels, pixels, -1.0*pixels), vmin = 0.51, vmax= 10.5)
	X = np.arange(-img.shape[0] / 2, img.shape[0] / 2)
	Y = np.arange(-img.shape[1] / 2, img.shape[1] / 2)

	X = np.repeat(X[:, np.newaxis], img.shape[1], axis = 1)
	Y = np.repeat(Y[:, np.newaxis], img.shape[0], axis = 1).T

	img = np.log10(img)

	plotLim = np.max(img)
	plotLim = [plotLim * 0.1, plotLim]
	limitPlots = img < plotLim[0]

	#im = ax1.plot_surface(X, Y, plotBelowAxis, cmap = 'viridis', vmax = plotLim[1], vmin = plotLim[0], alpha = 1.)
	midEle = [int((X.shape[0] / 2) * 1.1) + 1, int((X.shape[1] / 2) * 1.1) + 1]

	selection = (X + Y) >= 0
	imgCpy = img.copy()
	imgCpy[selection] = np.nan

	im = ax1.plot_surface(X, Y, imgCpy, cmap = 'viridis', vmax = plotLim[1], vmin = plotLim[0], alpha = .5)
	ax1.set_zlim(plotLim[0], plotLim[1])
	ax1.contour(X, Y, img, zdir = 'z', offset = ax1.get_zlim()[0],  cmap = 'viridis', vmax = plotLim[1], vmin = plotLim[0], linewidths = 1.)

	plotAboveAxis = img.copy()
	plotAboveAxis[limitPlots] = np.nan

	approx10Lines = plotAboveAxis.shape[0] / 10
	stridedPlotX = skim.block_reduce(plotAboveAxis, (approx10Lines, 1), np.nanmax)
	stridedPlotY = skim.block_reduce(plotAboveAxis, (1, approx10Lines), np.nanmax)

	ax1.contour(X[:midEle[0]:approx10Lines ,:], Y[:midEle[0]:approx10Lines ,:], stridedPlotX[:midEle[0]/approx10Lines + 1,:], zdir = 'x', offset = ax1.get_xlim()[0], cmap = 'viridis')
	ax1.contour(X[:, :midEle[1]:approx10Lines], Y[:, :midEle[1]:approx10Lines], stridedPlotY[:, :midEle[1]/approx10Lines + 1], zdir = 'y', offset = ax1.get_ylim()[0], cmap = 'viridis')
	cbar = plt.colorbar(im, shrink = 0.6)

	ax1.view_init(elev=25, azim=45.0)
	ax1.set_title('IE613 {0} MHz {1} station beam. Source: Sun'.format((round(freq[0]/1e6, 1)), antType))
	ax1.xaxis.set_major_locator(plt.NullLocator())
	ax1.yaxis.set_major_locator(plt.NullLocator())
	
	#divider = make_axes_locatable(ax1)
	#cax = divider.append_axes("right", size="2%", pad=0.05)
	#cbar = fig.colorbar(im)
	#cbar.set_label('log$_{10}$(I$_{beam}$)')
	print("{0} Beam Plotted for figure {1}".format(antType, idx))

	if colIdx == subplotTuple[1] - 1:
		plt.show()
		plt.savefig(outputFolder + '/station_uv_beam_'+str(format(idx, '03'))+'.png')   # save the figure to file 
		plt.close(fig)
	
def mpCheckAndCallFT(idx, dummy, uvw, px, res, fftType, maskVar, mpPool = None):
	if mpPool is not None:
		returnVar = mpPool.apply_async(beamGenerator, args = ([idx, dummy, uvw, px, res, fftType, maskVar]))
	else:
		returnVar = beamGenerator(idx, dummy, uvw, px, res, fftType, maskVar)

	return returnVar

def mainCall(opts, args):
	# Cleanup any past plots
	plt.close('all')

	# Extract parameters
	fftType = opts.fftType
	assert(fftType in [None, 'gaus', 'rect', 'rectgaus'])

	pltLBA = opts.LBA
	pltHBA = opts.HBA
	hbaActivation = opts.HBA_act
	obsTime = opts.time
	pixels = opts.pixels
	timeStep = float(opts.ts) * 60.
	outputFolder = opts.output_folder
	enableMp = opts.mp
	maskVar = opts.mask

	antFieldFile, antArrFile, hbaDeltas = initialiseLocalFiles(args)
	print("Files detected or acquired")

	lbaRcuMode = opts.lba_rcu_mode
	lbaLatLongEle = getLatLongElevation(antArrFile, 'IE613LBA')

	hbaRcuMode = opts.hba_rcu_mode
	hbaLatLongEle = getLatLongElevation(antArrFile, 'IE613HBA')

	if len(hbaActivation) == 1:
		hbaActivation = int(hbaActivation)
	if hbaActivation is not None:
		if hbaActivation in [0, '0', 'effelsberg', 'eff', 'effels']:
			hbaActivation = 0
			hbaActStr = '_Effelsberg'

		elif hbaActivation in [1, '1', 'generic', 'gen', 'orig']:
			hbaActivation = 1
			hbaActStr = '_Generic'

		elif hbaActivation in [2, '2', 'debug', 'db', 'dbg']:
			hbaActivation = 2
			hbaActStr = '_Debug'
		elif isinstance(hbaActivation, int):
			hbaActStr = '_unknown'
		else:
			raise LookupError("Unknown Activation Scheme {0}!".format(hbaActivation))
	else:
		hbaActStr = ''

	if outputFolder is None:
		outputFolder = 'IE613_station_uv_coverage'
		if pltLBA:
			outputFolder += '_LBA_{0}'.format(lbaRcuMode)
		if pltHBA:
			outputFolder += '_HBA_{0}{1}'.format(hbaRcuMode, hbaActStr)
				
		outputFolder += '_' + obsTime

		if fftType is not None:
			outputFolder += '_FFT_{0}'.format(fftType)
		else:
			outputFolder += '_DFT'


	if obsTime in ['summer', 'winter']:
		if obsTime == 'summer': obsTime = '06/21'
		elif obsTime == 'winter': obsTime = '12/21'
		else: print("Broken date statement, this message should be unreachable.")

	obsTime = "2018/{0} 15:00:00".format(obsTime) # Used to determine previous sunrise / next sunset

	uvwLBA, freqLBA, refTime, endTime, xyzLBA, obsLBA, sunLBAObj, sourceObjLBA = antProc(antFieldFile, antArrFile, lbaRcuMode, obsTime)
	obs = obsLBA
	sunObj = sunLBAObj

	uvwHBA, freqHBA, __, __, xyzHBA, obsHBA, sunHBAObj, sourceObjHBA = antProc(antFieldFile, antArrFile, hbaRcuMode, obsTime, [hbaActivation, hbaDeltas])

	image_num = 0
	plt.ion()


	zplane = 5.0769e3
	xplane = 3.80155e3
	yplane = -5.28875e2
	planeConsts = [xplane, yplane, zplane]


	xLBA, yLBA, zLBA, plotConstantsLBA = plotConsts(xyzLBA, planeConsts)
	xHBA, yHBA, zHBA, plotConstantsHBA = plotConsts(xyzHBA, planeConsts)


	px=[pixels,pixels]
	fov=np.pi    #Field of View in radians
	res=fov/px[0]   #pixel resolution
	dummy = 0.0   
	
	if not os.path.exists(outputFolder):
		os.mkdir(outputFolder)

	if enableMp:
		processes = mp.cpu_count() - 1
		mpPool = mp.Pool(processes = processes)
	else:
		mpPool = None

	idx = 0

	lbaRet = []
	hbaRet = []
	while refTime <= endTime:
		print("Starting Loop")
		print(obs.date)
		refTimeUtc = refTime #UTC (Use for Winter)
		refTimeIst = pytz.timezone('UTC').localize(refTimeUtc).astimezone(pytz.timezone('Europe/Dublin'))
		obs.date = refTimeUtc
		sunObj.compute(obs)
		print(refTimeUtc, refTimeIst)
			
		#TODO: Progamatically adjust
		uvPlots = sum([pltHBA, pltLBA])
		print(10 + uvPlots * 8, 'uvPlots')
		fig = plt.figure(idx, figsize=(10 + uvPlots * 8, 10))
		
		####################################
		#       Plot IE613 and sun
		subplotTuple = (2, 2 + uvPlots)
		ax = plt.subplot2grid(subplotTuple, (0, 0), rowspan= 2, colspan = 2, projection='3d', fig = fig)
		ax = fig.gca()

		alt = math.degrees(sunObj.alt)
		az = math.degrees(sunObj.az) 
		
		sun_ecef = pm.aer2ecef(az, alt, 100.0, 53.09472, -7.9213880, 75.0, deg=True) #pm.geodetic2ecef(53.09472, -7.921388, 75.0, deg=True)

		ax = processPlot([xLBA, yLBA, zLBA], obsLBA, sunLBAObj, refTimeUtc, plotConstantsLBA, sun_ecef, ax, antType = 'LBA', plotSuns = True) # Plot suns first to prevent z-buffer issues
		ax = processPlot([xHBA, yHBA, zHBA], obsHBA, sunHBAObj, refTimeUtc, plotConstantsHBA, sun_ecef, ax, antType = 'HBA', plotSuns = False) 
	
		print("Projecions plotted")

		box_width = 0.2
		xbox_cen = 3.80155e3
		ybox_cen = -5.2889e2
		zbox_cen = 5.0769e3
		ax.set_xlim3d(xbox_cen, xbox_cen+box_width)
		ax.set_ylim3d(ybox_cen -box_width, ybox_cen)
		ax.set_zlim3d(zbox_cen, zbox_cen+box_width)

		plt.title('IE613 solar alt: %s$^{\circ}$, az: %s$^{\circ}$ @ %s IST' % (round(alt), round(az), refTimeUtc))

		ax.set_xlabel('X (km)')
		ax.set_ylabel('Y (km)')
		ax.set_zlabel('Z (km)')

		print(2. / (2. + uvPlots) + 0.15)
		plt.gca().set_position([0., 0., 2. / (2. + uvPlots), 1.])
		ax.view_init(elev=16, azim=-55.0)
		

		ax.tick_params(axis = 'z', direction = 'out', pad = 18.0)
		ax.tick_params(axis = 'y', direction = 'out', pad = 14.0)
		print("Plot reoritentatied")

		####################################
		#       Plot UV coverage, beam
		#		If they aren't called separately, we'll randomly lose graphs as matplotlib hates multithreading.
		colIdx = 2
		if pltLBA:
			uvwLBA = processUVCoverage(idx, sunLBAObj, obs, xyzLBA, sourceObjLBA, obsLBA, freqLBA, colIdx, subplotTuple, 'LBA')
			colIdx += 1
		if pltHBA:
			uvwHBA = processUVCoverage(idx, sunHBAObj, obs, xyzHBA, sourceObjHBA, obsHBA, freqHBA, colIdx, subplotTuple, 'HBA')

		colIdx = 2
		if pltLBA:
			lbaRet.append([idx, mpCheckAndCallFT(idx, dummy, uvwLBA, px, res, fftType, maskVar, mpPool)])
			colIdx +=1
		if pltHBA:
			hbaRet.append([idx, mpCheckAndCallFT(idx, dummy, uvwHBA, px, res, fftType, maskVar, mpPool)])
		
		colIdx = 2
		if not enableMp:
			if pltLBA:
				beamPlot(idx, lbaRet[-1][1], dummy, pixels, freqLBA, outputFolder, fftType, colIdx, subplotTuple, 'LBA')
				colIdx += 1
			if pltHBA:
				beamPlot(idx, hbaRet[-1][1], dummy, pixels, freqHBA, outputFolder, fftType, colIdx, subplotTuple, 'HBA')
		
		if not (pltLBA or pltHBA):
			fig.savefig(outputFolder + '/station_uv_beam_'+str(format(idx, '03'))+'.png')   # save the figure to file
			plt.close(idx)

		refTime += datetime.timedelta(seconds = timeStep)
		idx += 1

	# Finish plots, cleanup if we used mp
	if enableMp:
		mpPool.close()
		mpPool.join()

		coldIdx = 2
		if pltLBA:
			for idx, asyncCallback in lbaRet:
				lbaImg = asyncCallback.get()
				beamPlot(idx, lbaImg, dummy, pixels, freqLBA, outputFolder, fftType, colIdx, subplotTuple, 'LBA')

		colIdx = 2 + pltLBA
		if pltHBA:
			for idx, asyncCallback in hbaRet:
				hbaImg = asyncCallback.get()
				beamPlot(idx, hbaImg, dummy, pixels, freqHBA, outputFolder, fftType, colIdx, subplotTuple, 'HBA')

	plt.close('all')

	# Generate the output video
	subprocess.call([ffmpegLoc, '-y',  '-r',  '20',  '-i',  '{0}/station_uv_beam_%03d.png'.format(outputFolder), '-vb', '50M', "{0}/{0}.mpg".format(outputFolder)])

	
#ffmpeg -y -r 20 -i image_%03d.png -vb 50M IE613_uv_coverage_sun.mpg    
	
#########################################################   
def mainCallBreakThings(opts, args):

	# Extract parameters
	fftType = opts.fftType
	assert(fftType in [None, 'gaus', 'rect', 'rectgaus'])

	pltLBA = opts.LBA
	pltHBA = opts.HBA
	hbaActivation = opts.HBA_act
	obsTime = opts.time
	pixels = opts.pixels
	timeStep = float(opts.ts) * 60.
	outputFolder = opts.output_folder
	enableMp = opts.mp
	maskVar = opts.mask

	antFieldFile, antArrFile, hbaDeltas = initialiseLocalFiles(args)
	print("Files detected or acquired")

	lbaRcuMode = opts.lba_rcu_mode
	lbaLatLongEle = getLatLongElevation(antArrFile, 'IE613LBA')

	hbaRcuMode = opts.hba_rcu_mode
	hbaLatLongEle = getLatLongElevation(antArrFile, 'IE613HBA')

	obsTime = '06/21'

	obsTime = "2018/{0} 15:00:00".format(obsTime) # Used to determine previous sunrise / next sunset

	uvwLBA, freqLBA, refTime, endTime, xyzLBA, obsLBA, sunLBAObj, sourceObjLBA = antProc(antFieldFile, antArrFile, lbaRcuMode, obsTime)
	obs = obsLBA
	sunObj = sunLBAObj

	px=[pixels,pixels]
	fov=np.pi    #Field of View in radians
	res=fov/px[0]   #pixel resolution
	dummy = 0.0   
	processes = mp.cpu_count() - 1
	mpPool = mp.Pool(processes = processes)
	obs.date = refTimeUtc + datetime.timedelta(hours = 6)
	sunObj.compute(obs)
	colIdx = 2
	while True:
		for i in range(200):
			uvwHBA, freqHBA, __, __, xyzHBA, obsHBA, sunHBAObj, sourceObjHBA = antProc(antFieldFile, antArrFile, hbaRcuMode, obsTime, ['random', hbaDeltas])

		uvwHBA = processUVCoverage(idx, sunHBAObj, obs, xyzHBA, sourceObjHBA, obsHBA, freqHBA, colIdx, subplotTuple, 'HBA')

		hbaRet.append([idx, mpCheckAndCallFT(idx, dummy, uvwHBA, px, res, fftType, maskVar, mpPool)])



		mpPool.close()
		mpPool.join()


		if pltHBA:
			for idx, asyncCallback in hbaRet:
				hbaImg = asyncCallback.get()
				beamPlot(idx, hbaImg, dummy, pixels, freqHBA, outputFolder, fftType, colIdx, subplotTuple, 'HBA')

	plt.close('all')


print("Initialisation of functions Complete")

if __name__ == '__main__':
	from optparse import OptionParser
	o = OptionParser()
	o.set_usage('Plot of IE613 LOFAR Antennas and Expected Responces')
	o.set_description(__doc__)

	o.add_option('--lba_rcu_mode', dest = 'lba_rcu_mode', help = "LBA RCU Mode considered for processing", default = 3)
	o.add_option('--hba_rcu_mode', dest = 'hba_rcu_mode', help = "HBA RCU Mode considered for processing", default = 5)
	
	o.add_option('-f', '--fft', dest = 'fftType', help = "Perform beam synthesis by FFT/gridding with the (rect)angular or (gaus)sian method for windowing.", default = None)
	
	o.add_option('-l', '--lba', action = 'store_false', dest = 'LBA', help = "Plot LBA UV/beam in output image.", default = True)
	o.add_option('-c', '--hba', action = 'store_true', dest = 'HBA', help = "Plot HBA UV/beam in output image.", default = False)
	
	o.add_option('-a', '--hba_activation', dest = 'HBA_act', help = "HBA activation pattern ('effelsberg', 'generic', 'debug', ...)", default = None)
	o.add_option('-t', '--time', dest = 'time', help = "Observation date (2018, provide date in mm/dd syntax, or 'winter' / 'summer' for solstices)", default = 'winter')
	o.add_option('-p', '--pixels', dest = 'pixels', help = "Size of the synthesized beam image in pixels (DFT METHOD ONLY)", default = 128)
	o.add_option('-o', '--output_folder', dest = 'output_folder', help = "Where to put the output files (can be auto generated)", default = None)
	o.add_option('-s', '--time_step', dest = 'ts', help = "Time (in minutes) between sampling steps.", default = 5.)
	o.add_option('-m', '--multi-thread', dest = 'mp', action = 'store_true', help = "Using this flag will enable multithreading with n-1 of available CPU cores being used.", default = False)
	o.add_option('-r', '--mask_radius', dest = 'mask', action = 'store_true', help = "Limit the output beam shape to a circle (TODO: Inputable radius)", default = False)

	opts, args = o.parse_args(sys.argv[1:])  
	mainCall(opts, args)
#
#########################################################   