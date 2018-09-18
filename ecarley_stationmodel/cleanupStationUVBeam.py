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

global rcuInfo
rcuInfo = [ {'mode':'OFF', 'rcuID':0, 'array_type':'LBA', 'bw':100000000., 'refFreq': 50e6},            #0
			{'mode':'LBL_HPF10MHZ', 'rcuID':1, 'array_type':'LBA', 'bw':100000000., 'refFreq': 50e6},   #1
			{'mode':'LBL_HPF30MHZ', 'rcuID':2, 'array_type':'LBA', 'bw':100000000., 'refFreq': 50e6},   #2
			{'mode':'LBH_HPF10MHZ', 'rcuID':3, 'array_type':'LBA', 'bw':100000000., 'refFreq': 50e6},   #3
			{'mode':'LBH_HPF30MHZ', 'rcuID':4, 'array_type':'LBA', 'bw':100000000., 'refFreq': 50e6},   #4
			{'mode':'HBA_110_190MHZ', 'rcuID':5, 'array_type':'HBA', 'bw':100000000., 'refFreq': 150e6}, #5
			{'mode':'HBA_170_230MHZ', 'rcuID':6, 'array_type':'HBA', 'bw':100000000., 'refFreq': 200e6}, #6
			{'mode':'HBA_210_290MHZ', 'rcuID':7, 'array_type':'HBA', 'bw':100000000., 'refFreq': 225e6}] #7

def eq2top_m(ha, dec):
	"""Return the 3x3 matrix converting equatorial coordinates to topocentric
	at the given hour angle (ha) and declination (dec)."""
	sin_H, cos_H = np.sin(ha), np.cos(ha)
	sin_d, cos_d = np.sin(dec), np.cos(dec)
	zero = np.zeros_like(ha)

	mapMatrix =  np.array([[    sin_H    ,       cos_H  ,       zero  ],
						[ -sin_d*cos_H,   sin_d*sin_H,      cos_d  ],
						[  cos_d*cos_H,  -cos_d*sin_H,      sin_d  ]])
	if len(mapMatrix.shape) == 3: mapMatrix = mapMatrix.transpose([2, 0, 1])
	return mapMatrix

def get_baseline(i, j, src, obs, includeSelfRef = True):
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
	return np.dot(m, bl).transpose()

def gen_uvw(i, j, src, obs, f):
	"""Compute uvw coordinates of baseline relative to provided FixedBody"""
	x,y,z = get_baseline(i,j,src,obs)

	afreqs = np.reshape(f, (1,f.size))
	afreqs = afreqs/ephem.c #1/wavelength

	if len(x.shape) == 0: return np.array([x*afreqs, y*afreqs, z*afreqs])
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
	for k in range(px[0]):
		for l in range(px[1]):
			im[k,l]=dft2(d,(k-mid_k),(l-mid_l),u,v)
			if mask:        #mask out region beyond field of view
				rad=(((k-mid_k)*res)**2 + ((l-mid_l)*res)**2)**.5
				if rad > mid_k*res: im[k,l]=0
				#else: im[k,l]=dft2(d,(k-mid_k),(l-mid_l),u,v)
	print time.time()-start_time

	# Debug after removing (0,0) UV elements
	im = np.abs(im)
	return im

def fftImage(d,uvw,pxPlot,res,mask=False, useDVar = False, method = 'gaus'):
	"""return a FFT image"""
	start_time = time.time()
	nants=uvw.shape[0]
	im=np.zeros((pxPlot[0],pxPlot[1]),dtype=complex)

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
		
		if method == 'gaus':
			pointSamples = gaussLambda(sampleCache)
		
		elif method == 'rect':
			pointSamples = np.zeros([81])
			pointSamples[30] = 1.

		elif method == 'rectgaus':
			pointSamples = np.zeros([81])
			pointSamples[30] = gaussLambda(offsets)

		else:
			print("Unknown method! Falling back to Gaussian.")
			pointSamples = gaussLambda(sampleCache)

		sampleIndex = ((sampleCoord[..., np.newaxis] + offsets + sampleCache.reshape(1, 2, -1)) / (res) + centralRef).astype(int)

		uvGrid[sampleIndex[0, 0, :], sampleIndex[0, 1, :]] += pointSamples

		sampleCache -= offsets
	
	uvGrid = uvGrid[5:-6, 5:-6] # Remove padding


	im = np.fft.fft2(uvGrid)
	im = np.fft.fftshift(np.abs(im.real)) # Abs of FFT == DFT?

	if useDVar:
		im *= d
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

def read_ant_xyz(antFieldFile, rcumode, station):
	fh=open(antFieldFile)
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
	xyz=np.array(xyz) # Pretty sure this is Earth Centered, Earth Fixed (ECEF) Coords.

	return xyz

def initialiseLocalFiles(args):
	if args:
		antFieldFile=args[0]
	else:
		antFieldFile = "./antFieldFile.conf"
		if not os.path.exists(antFieldFile):
			import urllib
			urllib.urlretrieve("https://raw.githubusercontent.com/griffinfoster/SWHT/master/SWHT/data/LOFAR/StaticMetaData/IE613-AntennaField.conf", antFieldFile)
	
	
	antArrFile = "./AntennaArrays.conf"
	
	if not os.path.exists(antArrFile):
		import urllib
		urllib.urlretrieve("https://raw.githubusercontent.com/griffinfoster/SWHT/master/SWHT/data/LOFAR/StaticMetaData/AntennaArrays/AntennaArrays_Int.conf", antArrFile)
	
	return antFieldFile, antArrFile

def getLatLongElevation(antArrFile, arrayType):
	with open(antArrFile) as fh:
		readLatLon=False

		for line in fh:
			if line.lower().startswith(arrayType.lower()):
				readLatLon=True
				continue
			if readLatLon:
				lon=line.split(' ')[2]
				lat=line.split(' ')[3]
				elev=line.split(' ')[4]
				readLatLon=False
				break

	return [lon, lat, elev]

def antProc(antFieldFile, antArrFile, rcuMode, obsTime, hbaActivation = None):
	antType = rcuInfo[rcuMode]['array_type']
	antLatLongEle = np.array(getLatLongElevation(antArrFile, antType), dtype = float)

	freq = np.array([rcuInfo[rcuMode]['refFreq']])

	obs = ephem.Observer()
	sunObj = ephem.Sun()

	obs.lon, obs.lat = math.radians(antLatLongEle[0]), math.radians(antLatLongEle[1])
	obs.elevation = antLatLongEle[2]
	obs.epoch=2000.0

	obs.date = obsTime

	initTime = str(obs.previous_rising(sunObj))[:-5] + '00:00'

	initTime = initTime[:-7] + str(int(initTime[-7]) + 1) + initTime[-6:]
	endTime = str(obs.next_setting(ephem.Sun()))[:-5] + '00:00'

	initTime = calendar.timegm( time.strptime(initTime, "%Y/%m/%d %H:%M:%S") )
	endTime = calendar.timegm( time.strptime(endTime, "%Y/%m/%d %H:%M:%S") )
	initTime = datetime.datetime.utcfromtimestamp(initTime)
	endTime = datetime.datetime.utcfromtimestamp(endTime)


	sunObj.compute(obs)
	sourceObj = sunObj
	sourceObj._ra = sunObj.ra
	sourceObj._dec = sunObj.dec
	sourceObj.compute(obs)

	xyz = read_ant_xyz(antFieldFile, rcuMode, 'IE')

	if hbaActivation is None:
		uvw = xyz2uvw(xyz, sourceObj, obs, freq[0])
	else:
		return #TODO

	return uvw, freq, initTime, endTime, xyz, obs, sunObj

def plotConsts(xyz, planeConsts):
	x = np.array([coords[0] / 1e3 for coords in xyz])
	y = np.array([coords[1] / 1e3 for coords in xyz])
	z = np.array([coords[2] / 1e3 for coords in xyz])

	meanX = np.mean(x)
	meanY = np.mean(y)
	meanZ = np.mean(z)

	zproj = np.full_like(z, planeConsts[0])
	xproj = np.full_like(x, planeConsts[1])
	yproj = np.full_like(y, planeConsts[2])

	u = np.linspace(0, 2. * np.pi, 200)
	v = np.linspace(0, np.pi, 200)

	xSphereInit = 7e-3*np.outer(np.cos(u), np.sin(v)) 
	ySphereInit = 7e-3*np.outer(np.sin(u), np.sin(v)) 
	zSphereInit = 7e-3*np.outer(np.ones(np.size(u)), np.cos(v))

	return x, y, z, [meanX, meanY, meanZ], zproj, xproj, yproj, [xSphereInit, ySphereInit, zSphereInit]


def mainCall(opts, args):

	fftType = opts.fftType
	print(fftType)
	assert(fftType in [None, 'gaus', 'rect', 'rectgaus'])
	pltLba = opts.LBA
	pltHba = opts.HBA
	obsTime = opts.time
	pixels = opts.pixels

	if obsTime in ['summer', 'winter']:
		if obsTime == 'summer': obsTime = '06/21'
		elif obsTime == 'winter': obsTime = '12/21'
		else: print("Broken date statement, this message should be unreachable.")
	obsTime = "2017/{0} 13:00:00".format(obsTime) # Apparently ephem adds a year when you get the next sunrise.

	antFieldFile, antArrFile = initialiseLocalFiles(args)
	print("Files detected or acquired")

	# TODO: move more logic into here, after cleanup
	lbaRcuMode = opts.lba_rcu_mode
	lbaLatLongEle = getLatLongElevation(antArrFile, 'LBA')

	hbaAct = opts.HBA_act
	hbaRcuMode = opts.hba_rcu_mode
	hbaLatLongEle = getLatLongElevation(antArrFile, 'HBA')

	print("Parameter Initialisation Complete")
	#------------------------------------------------
	#     Define local geographic coords and time
	#
	uvwLBA, freqLBA, refTime, endTime, xyzLBA, obsLba, sunLbaObj = antProc(antFieldFile, antArrFile, lbaRcuMode, obsTime)
	uvwHBA, freqHBA, __, __, xyzHBA, obsHba, sunHbaObj = antProc(antFieldFile, antArrFile, hbaRcuMode, obsTime, hbaAct)

	image_num = 0
	plt.ion()
	
	print("Coords / time processed")
	#----------------------------------------------------
	#      Get IE613 antenna positions and projections	
	zplane = 5.0769e3
	xplane = 3.80155e3
	yplane = -5.28875e2
	planeConsts = [zplane, xplane, yplane]

	xLBA, yLBA, zLBA, meanLbaArr, zLbaProj, xLbaProj, yLbaProj, sphereInit = plotConsts(xyzLBA, planeConsts)
	xHBA, yHBA, zHBA, meanHbaArr, zHbaProj, xHbaProj, yHbaProj, __ = plotConsts(xyzHBA, planeConsts)

	meanLbaX, meanLbaY, meanLbaZ = meanLbaArr
	meanHbaX, meanHbaY, meanHbaZ = meanHbaArr

	xSphereInit, ySphereInit, zSphereInit = sphereInit
	
	print("Projections Processed")
	#------------------------------------------
	#       For the instrument beam image
	px=[pixels,pixels]
	fov=np.pi    #Field of View in radians
	res=fov/px[0]   #pixel resolution
	dummy = 0.0   
	
	if not os.path.exists("./station_plots"):
		os.mkdir("./station_plots")

	while refTime < endTime:
	
		print("Starting Loop")

		refTimeUtc = refTime #UTC (Use for Winter)
		refTimeIst = pytz.timezone('UTC').localize(refTimeUtc).astimezone(pytz.timezone('Europe/Dublin'))
		
		print(refTimeUtc, refTimeIst)
		#if pltLba:
		#	 = processLba()

		#if pltHba:
		#	 = processHba()

		obsLba.date = refTimeUtc
		obsHba.date = refTimeUtc

		sunLbaObj.compute(obsLba)
		sunHbaObj.compute(obsHba)
		print("Sun Obtained")
		
		#TODO: Progamatically adjust
		fig = plt.figure(figsize=(18, 10))
	
		####################################
		#       Plot IE613 and sun
		sourceObjLBA=sunLbaObj
		sourceObjLBA._ra=sunLbaObj.ra
		sourceObjLBA._dec=sunLbaObj.dec
		sourceObjLBA.compute(obsLba)

		sourceObjHBA=sunHbaObj
		sourceObjHBA._ra=sunHbaObj.ra
		sourceObjHBA._dec=sunHbaObj.dec
		sourceObjHBA.compute(obsHba)


		#TODO: Programatically adjust
		ax = plt.subplot2grid((2, 2), (0, 0), rowspan=2, projection='3d')

		print("Plotted for IE613")

		#----------------------------------------------------
		#    Define solar ephem object to get local alt-az     
		#
		alt = math.degrees(sunLbaObj.alt)
		az = math.degrees(sunLbaObj.az) 
	
		sun_ecef = pm.aer2ecef(az, alt, 150.0, 53.09472, -7.9213880, 75.0, deg=True) #pm.geodetic2ecef(53.09472, -7.921388, 75.0, deg=True)
	
		#TODO: Abstract
		lba_to_sunx = [ meanLbaX, sun_ecef[0]/1000.0 ]
		lba_to_suny = [ meanLbaY, sun_ecef[1]/1000.0 ]
		lba_to_sunz = [ meanLbaZ, sun_ecef[2]/1000.0 ]

		hba_to_sunx = [ meanHbaX, sun_ecef[0]/1000.0 ]
		hba_to_suny = [ meanHbaY, sun_ecef[1]/1000.0 ]
		hba_to_sunz = [ meanHbaZ, sun_ecef[2]/1000.0 ]

		xSphereLBA = xSphereInit + lba_to_sunx[1]
		ySphereLBA = ySphereInit + lba_to_suny[1]
		zSphereLBA = zSphereInit + lba_to_sunz[1]

		xSphereHBA = xSphereInit + hba_to_sunx[1]
		ySphereHBA = ySphereInit + hba_to_suny[1]
		zSphereHBA = zSphereInit + hba_to_sunz[1]
	
	
		#----------------------------------------- 
		#        Plot and format
		#TODO: Abstract
		ax.plot(xLBA, yLBA, zLBA, 'bo', label='LBAs (ECEF coords)', zorder=-1)
		ax.plot(xHBA, yHBA, zHBA, 'o', color='salmon', label='HBAs (ECEF coords)', zorder=-1)

		ax.plot_surface(xSphereLBA, ySphereLBA, zSphereLBA, color='y', linewidth=0.01, zorder=1)
		ax.plot(lba_to_sunx, lba_to_suny, lba_to_sunz, 'g', zorder=2)

		ax.plot_surface(xSphereHBA, ySphereHBA, zSphereHBA, color='y', linewidth=0.01, zorder=1)
		ax.plot(hba_to_sunx, hba_to_suny, hba_to_sunz, 'g', zorder=2)
		print("Station Plotted")
		#----------------------------------------- 
		#    Plot projections on the XYZ planes
		#
		lba_to_sunzproj = np.full_like(lba_to_sunz, zplane)
		lba_to_sunxproj = np.full_like(lba_to_sunx, xplane)
		lba_to_sunyproj = np.full_like(lba_to_suny, yplane)

		hba_to_sunzproj = np.full_like(hba_to_sunz, zplane)
		hba_to_sunxproj = np.full_like(hba_to_sunx, xplane)
		hba_to_sunyproj = np.full_like(hba_to_suny, yplane)
	
		zSphereProjLBA = np.full_like(zSphereLBA, zplane)
		xSphereProjLBA = np.full_like(xSphereLBA, xplane)
		ySphereProjLBA = np.full_like(ySphereLBA, yplane)

		zSphereProjHBA = np.full_like(zSphereHBA, zplane)
		xSphereProjHBA = np.full_like(xSphereHBA, xplane)
		ySphereProjHBA = np.full_like(ySphereHBA, yplane)


		ax.plot(xLbaProj, yLBA, zLBA, '.', color='lightgray', zorder=-3)
		ax.plot(xHbaProj, yHBA, zHBA, '.', color='lightgray', zorder=-3)
		ax.plot(lba_to_sunxproj, lba_to_suny, lba_to_sunz, color='lightgray', zorder=-3)
		ax.plot_surface(xSphereProjLBA, ySphereLBA, zSphereLBA, color='lightgray', linewidth=0.01, zorder=-3)
		ax.plot(hba_to_sunxproj, hba_to_suny, hba_to_sunz, color='lightgray', zorder=-3)
		ax.plot_surface(xSphereProjHBA, ySphereHBA, zSphereHBA, color='lightgray', linewidth=0.01, zorder=-3)

		ax.plot(xLBA, yLBA, zLbaProj, '.', color='lightgray', zorder=-3)
		ax.plot(xHBA, yHBA, zHbaProj, '.', color='lightgray', zorder=-3)
		ax.plot(lba_to_sunx, lba_to_suny, lba_to_sunzproj, color='lightgray', zorder=-3)
		ax.plot_surface(xSphereLBA, ySphereLBA, zSphereProjLBA, color='lightgray', linewidth=0.01, zorder=-3)
		ax.plot(hba_to_sunx, hba_to_suny, hba_to_sunzproj, color='lightgray', zorder=-3)
		ax.plot_surface(xSphereHBA, ySphereHBA, zSphereProjHBA, color='lightgray', linewidth=0.01, zorder=-3)
	
		ax.plot(xLBA, yLbaProj, zLBA, '.', color='lightgray', zorder=-3)
		ax.plot(xHBA, yHbaProj, zHBA, '.', color='lightgray', zorder=-3)
		ax.plot(lba_to_sunx, lba_to_sunyproj, lba_to_sunz, color='lightgray', zorder=-3)
		ax.plot_surface(xSphereLBA, ySphereProjLBA, zSphereLBA, color='lightgray', linewidth=0.01, zorder=-3)
		ax.plot(hba_to_sunx, hba_to_sunyproj, hba_to_sunz, color='lightgray', zorder=-3)
		ax.plot_surface(xSphereHBA, ySphereProjHBA, zSphereHBA, color='lightgray', linewidth=0.01, zorder=-3)

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
		plt.title('IE613 solar alt: %s$^{\circ}$, az: %s$^{\circ}$ @ %s IST' % (round(alt), round(az), refTimeUtc))
		ax.set_xlabel('X (km)')
		ax.set_ylabel('Y (km)')
		ax.set_zlabel('Z (km)')
		plt.gca().set_position([0, 0.0, 0.52, 1.0])
		#ax.legend()
		ax.view_init(elev=16, azim=-55.0)
		print("Plot reoritentatied")
	
		####################################
		#       Plot UV coverage
		uvw=xyz2uvw(xyzLBA, sourceObjLBA, obsLba, freqLBA[0])
		U=uvw[:,:,0]
		V=uvw[:,:,1]
		ax0 = plt.subplot2grid((2, 2), (0, 1))
		ax0.set_aspect('equal')
		plt.plot(U, V, '.')
		plt.xlabel('u ($\lambda$)')
		plt.ylabel('v ($\lambda$)')
		plt.title('UV-coverage at %s MHz for IE613 LBAs. Source: Sun' % (round(freqLBA[0]/1e6, 1)))
		plt.axis([-10, 10, -10, 10])
		print("UVW Plotted")
	
		####################################
		#       Plot instrument beam
		ax1 = plt.subplot2grid((2, 2), (1, 1))

		if fftType is None:
			img = dftImage(dummy, uvw, px, res, mask=True)
		else:
			img = fftImage(dummy, uvw, px, res, mask = True, method = fftType)
		img=img.real
		img=np.log(img)
	
		ax1 = plt.subplot2grid((2, 2), (1, 1))
		im = plt.imshow(img, extent = (-1.0*pixels, pixels, pixels, -1.0*pixels), vmin = 0.51, vmax= 10.5)
		plt.title('IE613 %s MHz LBA station beam. Source: Sun' % (round(freqLBA[0]/1e6, 1)))
		plt.gca().yaxis.set_major_locator(plt.NullLocator())
		plt.gca().xaxis.set_major_locator(plt.NullLocator())
	
		divider = make_axes_locatable(ax1)
		cax = divider.append_axes("right", size="2%", pad=0.05)
		cbar = plt.colorbar(im, cax=cax)
		cbar.set_label('log$_{10}$(I$_{beam}$)')
		print("Beam Plotted")
		plt.show( )
		plt.savefig('./station_plots/station_uv_beam_'+str(format(image_num, '03'))+'.png')   # save the figure to file
		#pdb.set_trace()    
		plt.close(fig)
		print(image_num)
		image_num+=1
		refTime = refTime + datetime.timedelta(minutes = 5.)
	
	
		#ffmpeg -y -r 20 -i image_%03d.png -vb 50M IE613_uv_coverage_sun.mpg    
	
#########################################################   

print("Initialisation of functions Complete")

if __name__ == '__main__':
	from optparse import OptionParser
	o = OptionParser()
	o.set_usage('Plot XYZ of IE613 LOFAR Antennas')
	o.set_description(__doc__)

	o.add_option('--lba_rcu_mode', dest = 'lba_rcu_mode', help = "LBA RCU Mode considered for processing", default = 3)
	o.add_option('--hba_rcu_mode', dest = 'hba_rcu_mode', help = "HBA RCU Mode considered for processing", default = 5)
	
	o.add_option('-f', '--fft', dest = 'fftType', help = "Perform beam synthesis by FFT/gridding with the (rect)angular or (gaus)sian method for windowing.", default = None)
	
	o.add_option('-l', '--lba', action = 'store_true', dest = 'LBA', help = "Plot LBA UV/beam in output image.", default = True)
	o.add_option('-c', '--hba', action = 'store_true', dest = 'HBA', help = "Plot HBA UV/beam in output image.", default = False)
	
	o.add_option('-a', '--hba_activation', dest = 'HBA_act', help = "HBA activation pattern ('effelsberg', 'generic', 'debug', ...)", default = None)
	o.add_option('-t', '--time', dest = 'time', help = "Observation date (2018, provide date in mm/dd syntax, or 'winter' / 'summer' for solstaces)", default = 'winter')
	o.add_option('-p', '--pixels', dest = 'pixels', help = "Size of the synthesized beam image in pixels", default = 128)
	opts, args = o.parse_args(sys.argv[1:])  
	mainCall(opts, args)
#
#########################################################   