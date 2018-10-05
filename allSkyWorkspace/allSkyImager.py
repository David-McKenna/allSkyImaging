

import ast
import scipy.constants
import scipy.ndimage
import astropy.coordinates
import astropy.time
import astropy.units as u
import datetime
import subprocess
import os
import numpy as np
import matplotlib as mpl
import matplotlib.patheffects as mplPe
import matplotlib.pyplot as plt

from functools import wraps
import multiprocessing as mp


ffmpegLoc = "/cphys/ugrad/2015-16/JF/MCKENND2/Downloads/ffmpeg-4.0.2-64bit-static/ffmpeg"

def cache(function):
	cachedLookup = {}
	@wraps(function)
	def trueFunc(*args):
		if args in cachedLookup:
			return cachedLookup[args]

		result = function(*args)
		cachedLookup[args] = result
		return result
	return trueFunc

def calcFreq(rcuMode, subband):
	baseFreq = 100.0

	rcuMode = int(rcuMode)
	subband = int(subband)
	if rcuMode == 5:
		freqOff = 100e6

	elif rcuMode == 6:
		freqOff = 160e6
		baseFreq = 80.0

	elif rcuMode == 7:
		freqOff = 200e6

	else:
		freqOff = 0

	frequency = ((baseFreq / 512 ) * (subband) + (freqOff/1e6))

	return frequency


def parseiHBAField(afilename, hbaDeltasFile, activeElems, rcuMode, EU):
	# [posxpol, posypol, lon, lat, refpos, rotmat] =
	#     parseiHBAField(filename, hbaDeltasFile, selection, rcuMode, EU)
	#
	# Parser for AntennaFields.conf and iHBADeltas.conf files.
	#
	# arguments
	# afilename : filename (including path if necessary) of the AntennaFields.conf
	#            file
	# hbaDeltasFile: filename (including path if necessary) of the iHBADeltas.conf
	#            file
	# selection: element number per tile [0..15]
	# rcuMode  : RCU mode number (determines array configuration)
	# EU       : true for European stations
	#
	# return values
	# posxpol : nElem x 3 matrix with (x, y, z)-positions of x-dipoles (in m)
	#           w.r.t. the reference position in ITRF coordinates
	# posypol : nElem x 3 matrix with (x, y, z)-positions of y-dipoles (in m)
	#           w.r.t. the reference position in ITRF coordinates
	# lon     : geographic longitude of the station (in degrees)
	# lat     : geographic latitude of the station (in degrees)
	# refpos  : nElem x 3 vector with (x, y, z)-coordinate of ITRF reference
	#           position of the specified antenna array (in m)
	# rotmat  : 3 x 3 rotation matrix to convert ITRF to local coordinates
	#
	# SJW, January 2011
	# modified by SJW, January 2012
	# including iHBADeltas.conf by MJN, November 2015
	# update new format Antenna Field size by MJN, June 2017

	# added by Pearse Murphy, August 2018 to account for sparse HBA antennas 
	elemsDict={'Effelsberg_elements_20091110' : [1,4,13,15,11,9,14,1,15,0,8,2,11,3,14,0,2,4,3,0,0,2,12,12,12,12,15,11,14,15,7,5,1,0,3,10,1,11,0,12,12,1,6,7,0,10,9,6,15,14,11,7,2,0,7,12,15,8,13,3,7,6,3,15,11,1,4,11,8,1,8,15,4,0,5,6,12,0,12,15,3,7,14,8,3,12,12,2,9,8,14,2,5,6,12,0], 
	'Generic_International_Station_20091110' : [15,0,15,3,9,15,14,2,0,3,4,14,10,8,5,15,12,0,2,11,3,12,12,1,5,4,4,8,6,3,0,5,3,11,3,2,8,15,13,8,3,2,9,1,14,8,8,0,12,13,0,11,15,3,12,3,13,3,10,5,0,10,1,6,4,10,3,15,3,14,0,12,0,7,0,12,7,3,13,0,7,3,15,4,14,4,3,8,4,9,12,0,14,9,3,11],
	'Generic_Int_201512' : [0,5,3,1,8,3,12,15,10,13,11,5,12,12,5,2,10,8,0,3,5,1,4,0,11,6,2,4,9,14,15,3,7,5,13,15,5,6,5,12,15,7,1,1,14,9,4,9,3,9,3,13,7,14,7,14,2,8,8,0,1,4,2,2,12,15,5,7,6,10,12,3,3,12,7,4,6,0,5,9,1,10,10,11,5,11,7,9,7,6,4,4,15,4,1,15],
	'Generic_Core_201512' : [0,10,4,3,14,0,5,5,3,13,10,3,12,2,7,15,6,14,7,5,7,9,0,15,0,10,4,3,14,0,5,5,3,13,10,3,12,2,7,15,6,14,7,5,7,9,0,15],
	'Generic_Remote_201512' : [0,13,12,4,11,11,7,8,2,7,11,2,10,2,6,3,8,3,1,7,1,15,13,1,11,1,12,7,10,15,8,2,12,13,9,13,4,5,5,12,5,5,9,11,15,12,2,15]}
	
	if activeElems:
		activeElements = elemsDict[activeElems]

		with open(hbaDeltasFile, 'r') as hbaDeltasRef:
			deltaLines = [line for line in hbaDeltasRef]

		arrayDeltas = parseBlitzFile(deltaLines, 'HBADeltas', False)
		np.save('./hbadeltas.npy', arrayDeltas)

	
	if rcuMode in [1, 2, 3, 4]:
		arrayName = 'LBA'
	else:
		arrayName = 'HBA'

	with open(afilename, 'r') as arrayRef:
		arrayLines = [line for line in arrayRef]
		
	arrayLoc = parseBlitzFile(arrayLines, arrayName, False)
	antLocs = parseBlitzFile(arrayLines, arrayName, True)

	print(antLocs[0], antLocs[0, 0, :], antLocs[0, :, 0])
	posX = antLocs[:, 0, :]
	posY = antLocs[:, 1, :]

	np.save('antLoc.npy', np.stack([posX, posY], axis = -1))

	if activeElems:
		posX += arrayDeltas[activeElements]
		posY += arrayDeltas[activeElements]
	
	np.save('antLoc2.npy', np.stack([posX, posY], axis = -1))
	# select the right set of antennas
	if not EU:
		nElem = posX.shape[0]
		if rcuMode in [1, 2]:
			sel = range(nElem/2,nElem) #sel=25:48
		elif rcuMode in [3, 4]:
			sel = range(0,nElem/2)         #sel=1:24
		else:
			sel = range(0,nElem)
		
		posX = posX[sel, :]
		posY = posY[sel, :]
	
#    # obtain rotation matrix
	if 'HBA0\n' in arrayLines:
		####Matlab code. Needs porting for completeness
		#        while ~((strcmp(token, 'HBA1') && strcmp(prevtoken, 'ROTATION_MATRIX')) || isempty(token))
		#            prevtoken = token;
		#            token = fscanf(fid, '#s', 1);
		#        end
		#        fscanf(fid, '#s', 4);
		#        rotmat2 = zeros(3, 3);
		#        rotmat2(:) = fscanf(fid, '#f', 9);
		#        posxpol(1:nElem/2, :) = posxpol(1:nElem/2, :) / rotmat;
		#        posxpol(nElem/2+1:nElem, :) = posxpol(nElem/2+1:nElem, :) / rotmat2;
		#        posypol(1:nElem/2, :) = posypol(1:nElem/2, :) / rotmat;
		#        posypol(nElem/2+1:nElem, :) = posypol(nElem/2+1:nElem, :) / rotmat2;
		####Matlab code end

		# David McKenna: Port attempted based on python implementation of else statement by Joe McCauley (don't have access to orignal source). Untested as I don't have matlab to test the actual intended output and haven't pointed as a core/remote station yet.
		rotationMatrixOne = parseBlitzFile(arrayLines, 'ROTATION_MATRIX ' + arrayName + '0', False).T
		rotationMatrixTwo = parseBlitzFile(arrayLines, 'ROTATION_MATRIX ' + arrayName + '0', False).T

		posX = np.matrix(posX)
		posY = np.matrix(posY)

		posX[:, :nElem/2+1] = posX[:, :nElem/2+1] * np.linalg.pinv(rotationMatrixOne)
		posY[:, :nElem/2+1] = posY[:, :nElem/2+1] * np.linalg.pinv(rotationMatrixOne)

		posX[:, nElem/2+1:] = posX[:, nElem/2+1:] * np.linalg.pinv(rotationMatrixTwo)
		posY[:, nElem/2+1:] = posY[:, nElem/2+1:] * np.linalg.pinv(rotationMatrixTwo)

		roationMatrix = np.dstack([rotationMatrixOne, rotationMatrixTwo])

	else:
		rotationMatrix = parseBlitzFile(arrayLines, 'ROTATION_MATRIX ' + arrayName, False).T #has to be transposed to get into correct format for further operations below
		rotationMatrix = np.matrix(rotationMatrix)

		posX = np.matrix(posX)
		posY = np.matrix(posY)
		posX = posX * np.linalg.pinv( rotationMatrix ) #right matrix division
		posY = posY * np.linalg.pinv( rotationMatrix )

	posX = np.array(posX)
	posY = np.array(posY)
	# obtain longitude and latitude of the station
	# David McKenna: I don't want to touch this wizardry.
	wgs84_f = 1 / 298.257223563
	wgs84_a = 6378137
	wgs84_e2 = wgs84_f * (2 - wgs84_f)
	lon = np.arctan2(arrayLoc[1], arrayLoc[0])
	lon=lon* 180 / np.pi
	r = np.sqrt(arrayLoc[0]**2 + arrayLoc[1]**2)
	prev_lat = 100000
	lat = np.arctan2(arrayLoc[2], r)
	while (abs(lat - prev_lat) >= 1e-12):
		prev_lat = lat
		normalized_earth_radius = 1 / np.sqrt((1-wgs84_f)**2 * np.sin(lat)**2 + np.cos(lat)**2)
		lat = np.arctan2(wgs84_e2 * wgs84_a * normalized_earth_radius * np.sin(lat) + arrayLoc[2], r)
	lat = lat * 180 /np.pi

	return posX, posY, lon, lat, arrayLoc, rotationMatrix

def parseBlitzFile(linesArray, keyword, refLoc = False):
	linesArray = [line.strip('\n') for line in linesArray]

	triggerLine = [idx for idx, line in enumerate(linesArray) if line == keyword][0]
	triggerLine += refLoc + 1 # refLoc = add a buffer line to account for a location refernece before processing.


	if ']' in linesArray[triggerLine]:
		line = linesArray[triggerLine]
		startIdx = line.find('[') + 1 # Don't include it as a character
		endIdx = line.find(']')

		dataElements = __processLine(line[startIdx:endIdx])
		return np.array(dataElements)

	endLine = [idx + triggerLine for idx, line in enumerate(linesArray[triggerLine:]) if ']' in line][0]
	iterateLines = range(triggerLine + 1, endLine)

	arrayShapeParts = linesArray[triggerLine].count('x') + 1

	startArrayLoc = linesArray[triggerLine].index('[') 
	splitTuples =  linesArray[triggerLine][:startArrayLoc].split('x')

	arrayShape = filter(None, splitTuples)[:arrayShapeParts]
	arrayShape = [ast.literal_eval(strEle.strip(' ')) for strEle in arrayShape]
	arrayShape = [tuplePair[1] - tuplePair[0] + 1 for tuplePair in arrayShape]

	arrayData = []
	for lineIdx in iterateLines:
		procLine = __processLine(linesArray[lineIdx])
		arrayData.append(procLine)

	arrayData = np.vstack(arrayData).reshape(arrayShape)

	return arrayData

def __processLine(line):
	line = filter(None, line.split(' '))
	return [float(element) for element in line]


def corrToSkyImage(correlationMatrix, posX, posY, obsFreq, lVec, mVec):
	# corrMat: (96,96,nChan)
	nElem = np.shape(correlationMatrix)[0]
	frameCount = correlationMatrix.shape[2]
	c = scipy.constants.c

	# Temp for lazy debugging
	if len(lVec.shape) != 2:
		lVec = lVec[np.newaxis, :] # (1, 96)

	if len(mVec.shape) != 2:
		mVec = mVec[np.newaxis, :] # (1, 96)

	if len(posX.shape) != 2:
		posX = posX[:, np.newaxis] # (96, 1)

	if len(posY.shape) != 2:
		posY = posY[:, np.newaxis] # (96, 1)

	skyView = np.zeros([lVec.size, mVec.size, frameCount]) # (l,m,nChan)
	
	wavelength = c / obsFreq
	k = (2 * np.pi) / wavelength

	wx = np.exp(-1j * k * posX * lVec) # (l, 96)
	wy = np.exp(-1j * k * posY * mVec) # (m, 96)
	weight = np.multiply(wx[:, np.newaxis, :], wy[:, :, np.newaxis]).transpose((1,2,0))[..., np.newaxis] # (l,m,96,1)
	conjWeight = np.conj(weight).transpose((0,1,3,2)) # (l,m,1,96)
	# Should be able to fully vectorise this over all channels, should give  a nicespeedup if we pass frames as alternative channels... for another day.
	for frame in range(frameCount):
		correlationMatrixChan = correlationMatrix[..., frame] # (96,96)
		print("Processing Frame {0} of {1} at Frequency {2:.2F}E+06".format(frame + 1, frameCount, obsFreq / 1e6))

		tempProd = np.dot(conjWeight, correlationMatrixChan) # w.H * corr # (l,m, 1, 96)

		skyView[..., frame] = np.sum(np.multiply(tempProd.transpose((0,1,3,2)), weight).real, axis = (2,3)) # (w.H * corr) * w, faster than loop below by about 50ms (running at 620ms, down from 2s in original implementation)

		#for l in range(mVec.size):
		#	for m in range(lVec.size):
		#		skyView[l, m, frame] = np.dot(tempProd[l, m], weight[l, m]).real[0, 0]
		#skyView[..., frame] = np.tensordot(tempProd, weight, ([3],[2])) # Dot products producing 201x201x1x201x201x1 arrays, thanks numpy.
	return skyView.transpose((1,0,2))

# -11.9
def generatePlots(inputCorrelations, antPos, plotOptions, dateArr, rcuMode, subband, multiprocessing = True, stationRotation = -11.9, plotX = True, plotY = True, mask = True, lVec = None, mVec = None, 
calibrationX = None, calibrationY = None, baselineLimits = None):
	inputCorrelationsX = inputCorrelations[..., 0]
	inputCorrelationsY = inputCorrelations[..., 1]

	posX = antPos[..., 0]
	posY = antPos[..., 1]

	if baselineLimits:
		print('Test baseline limits.')
		if not baselineLimits[1]:
			print('Test no 0th correlations')
			selections = np.arange(inputCorrelationsX.shape[0])
			inputCorrelationsX[selections, selections] = 1.
			inputCorrelationsY[selections, selections] = 1.
		else:
			minBaseline, maxBaseline = baselineLimits

			baselines = posX - posY.T
			print(posX.shape, posY.shape, np.median(np.abs(baselines)), np.percentile(np.abs(baselines), 66),  np.percentile(np.abs(baselines), 75),  np.percentile(np.abs(baselines), 80),  np.percentile(np.abs(baselines), 90))
			baselines = np.logical_or(np.abs(baselines) > maxBaseline, np.abs(baselines) < minBaseline)

			inputCorrelationsX[baselines] = 0.
			inputCorrelationsY[baselines] = 0.

	frequency = calcFreq(rcuMode, subband) * 1e6

	labelOptions = [dateArr, rcuMode, subband, frequency]

	if not lVec:
		lVec = np.arange(-1., 1. + 0.01, 0.008 ) # Default to 200 pixels
	if not mVec:
		mVec = lVec.copy()

	pixels = np.max([lVec.size, mVec.size])

	if mask:
		mask = np.zeros([lVec.size, mVec.size]).astype(bool)
		dist = np.sqrt(np.square(np.meshgrid(lVec)) + np.square(np.meshgrid(mVec)).T)

		mask[dist >= 1] = True
	else:
		mask = np.zeros([lVec.size, mVec.size]).astype(bool)

	### Multithread? Might not need it if it's fast enough.
	### Either way, this should be abtracted to a higher degree.
	figNum = 0
	returnVar = {}


	if plotX:
		if calibrationX is not None:
			print('Calibrating X for subband {0} from shape {1}'.format(subband, calibrationX.shape))
			calSubband = calibrationX[:, subband]
			print(calSubband.shape)
			calMatrixX = np.outer(calSubband, np.conj(calSubband).T)[..., np.newaxis]
			print(inputCorrelationsX.shape, calMatrixX.shape)
			inputCorrelationsX = np.multiply(np.conj(calMatrixX), inputCorrelationsX)

		labelOptionsX = labelOptions + ['X', 0]
		allSkyImX, __, __ = __processAllSkyIm(inputCorrelationsX, posX, posY, frequency, lVec, mVec, stationRotation, mask, plotOptions, labelOptionsX)
		print("X Polarisation Processed, begining plotting")

		xFileLoc = []
		for i in range(allSkyImX.shape[-1]):
			labelOptionsX[0] = dateArr[i]
			labelOptionsX[-1] = figNum
			figNum += 1
			xFileLoc.append(plotAllSkyImage(allSkyImX[..., i], plotOptions, labelOptionsX, pixels))

		if plotOptions[5] and allSkyImX.shape[-1] > 20:
			filePrefix = xFileLoc[0].split(' ')[0]
			fileSuffix = '_'.join(xFileLoc[0].split('/')[-1].split('_')[1:])[:-4]
			print("Exporting frames to video at " + "./{0}{1}.mpg".format(filePrefix, fileSuffix))
			subprocess.call([ffmpegLoc, '-y',  '-r',  '20', '-pattern_type', 'glob',  '-i',  '{0}*{1}.png'.format(filePrefix, fileSuffix), '-vb', '50M', "{0}{1}.mpg".format(filePrefix, fileSuffix)])

		returnVar['X'] = allSkyImX
		print(returnVar)

	if plotY:
		if calibrationY is not None:
			print('Calibrating Y for subband {0} from shape {1}'.format(subband, calibrationX.shape))
			calSubband = calibrationY[:, subband]
			print(calSubband.shape)
			calMatrixY = np.outer(calSubband, np.conj(calSubband).T)[..., np.newaxis]
			print(inputCorrelationsX.shape, calMatrixY.shape)
			inputCorrelationsY = np.multiply(np.conj(calMatrixY), inputCorrelationsY)

		labelOptionsY = labelOptions + ['Y', 0]
		allSkyImY, __, __ = __processAllSkyIm(inputCorrelationsY, posX, posY, frequency, lVec, mVec, stationRotation, mask, plotOptions, labelOptionsY)
		print("Y Polarisation Processed, begining plotting")

		yFileLoc = []
		for i in range(allSkyImY.shape[-1]):
			labelOptionsY[0] = dateArr[i]
			labelOptionsX[-1] = figNum
			figNum += 1
			yFileLoc.append(plotAllSkyImage(allSkyImY[..., i], plotOptions, labelOptionsY, pixels))

		if plotOptions[5] and allSkyImY.shape[-1] > 20:
			filePrefix = yFileLoc[0].split(' ')[0]
			fileSuffix = '_'.join(yFileLoc[0].split('/')[-1].split('_')[1:])[:-4]
			print("Exporting frames to video at " + "./{0}{1}.mpg".format(filePrefix, fileSuffix))
			subprocess.call([ffmpegLoc, '-y',  '-r',  '20',  '-pattern_type', 'glob', '-i',  '{0}*{1}.png'.format(filePrefix, fileSuffix), '-vb', '50M', "{0}{1}.mpg".format(filePrefix, fileSuffix)])

		returnVar['Y'] = allSkyImY
		print(returnVar)

	print(returnVar)

	return returnVar

def __processCorrPlot():
	if calibrationArr is not None:
		print('Calibrating {0} for subband {1} from shape {2}'.format(polChar, subband, calibrationX.shape))
		calSubband = calibrationArr[:, subband]
		calMatrixArr = np.outer(calSubband, np.conj(calSubband).T)[..., np.newaxis]
		inputCorrelations = np.multiply(np.conj(calMatrixArr), inputCorrelations)

	labelOptions = labelOptions + [polChar, 0]
	allSkyIm, __, __ = __processAllSkyIm(inputCorrelationsArr, posX, posY, frequency, lVec, mVec, stationRotation, mask, plotOptions, labelOptions)
	print("{0} Polarisation Processed, begining plotting".format(polChar))

	fileLoc = []
	for i in range(allSkyIm.shape[-1]):
		labelOptions[0] = dateArr[i]
		labelOptions[-1] = figNum
		figNum += 1
		fileLoc.append(plotAllSkyImage(allSkyIm[..., i], plotOptions, labelOptions, pixels))

	if plotOptions[5] and allSkyIm.shape[-1] > 20:
		filePrefix = fileLoc[0].split(' ')[0]
		fileSuffix = '_'.join(fileLoc[0].split('/')[-1].split('_')[1:])[:-4]
		print("Exporting frames to video at " + "./{0}{1}.mpg".format(filePrefix, fileSuffix))
		subprocess.call([ffmpegLoc, '-y',  '-r',  '20', '-pattern_type', 'glob',  '-i',  '{0}*{1}.png'.format(filePrefix, fileSuffix), '-vb', '50M', "{0}{1}.mpg".format(filePrefix, fileSuffix)])

	#return __, __, __...

def __processAllSkyIm(inputCorrelations, posX, posY, frequency, lVec, mVec, stationRotation, mask, plotOptions, labelOptions):
	allSkyIm = corrToSkyImage(inputCorrelations, posX, posY, frequency, lVec, mVec)
	allSkyIm = np.rot90(allSkyIm, 3)
	allSkyIm = np.flipud(allSkyIm)
	allSkyIm = scipy.ndimage.rotate(allSkyIm, stationRotation, mode = 'constant', cval = 100., reshape = False)
	
	allSkyIm[mask] = np.nan

	# Extra returns for if we start multithreading
	return allSkyIm, plotOptions, labelOptions

def plotAllSkyImage(allSkyImage, plotOptions, labelOptions, pixels):
	logPlot, skyObjColor, gridThickness, backgroundColor, foregroundColor, saveImage, radialLabelAngle, colorBar, obsSite, outputFolder = plotOptions
	dateTime, rcuMode, subband, frequency, polarity, figNum = labelOptions

	if obsSite in ['Birr', 'IE613', 'IE', 'EIRE']:
		telescopeLoc = astropy.coordinates.EarthLocation( lat = 53.095 * u.deg, lon = -7.9218 * u.deg, height = 100 * u.m )
	else:
		telescopeLoc = astropy.coordinates.EarthLocation.of_site(obsSite)
	latLon = telescopeLoc.to_geodetic()

	if len(dateTime) == 15:
		dateTime = str(datetime.datetime.strptime(dateTime, '%Y%m%d-%H%M%S'))

	obsTime = astropy.time.Time(dateTime)

	try:
		knownSources = ['Polaris', 'Cas A', 'Cyg A', 'Sgr A', 'Tau A', 'Vir A', 'Cen A', 'Vela']
		referenceObject = []
		referenceObject = [cachedSkyCoords(name) for name in knownSources]
	
		knownBodies = ['Sun', 'Jupiter', 'Moon', 'Uranus', 'Neptune']
		referenceObject.extend([astropy.coordinates.get_body(sourceName, obsTime, telescopeLoc) for sourceName in knownBodies])
		knownSources[0] = 'NCP (Polaris)'
		knownSources.extend(knownBodies)

	except astropy.coordinates.name_resolve.NameResolveError as timeoutException:
		print("Unable to resolve all sources (likely timed out on database lookup), skipping plotting some sources.")
		print(timeoutException)

	altAzRef = astropy.coordinates.AltAz(obstime = obsTime, location = telescopeLoc)
	altAzObjects = [skyObj.transform_to(altAzRef) for skyObj in referenceObject]

	gs = mpl.gridspec.GridSpec(1, 2, width_ratios = [18, 1])
	gs2 = mpl.gridspec.GridSpec(1, 2, width_ratios = [18, 1])

	fig = plt.figure(figNum, figsize = (18, 14))
	fig.patch.set_facecolor(backgroundColor)
	plt.suptitle( 'LOFAR mode {0}{1} all sky plot at {2}MHz (sb{3}) for {4}\n'.format(rcuMode, polarity, round(frequency / 1e6, 2), subband, obsSite), fontsize = 28, color=foregroundColor )#, va = 'top') 
	plt.rcParams["text.color"] = foregroundColor
	plt.rcParams["axes.labelcolor"] = foregroundColor
	plt.rcParams["xtick.color"] =  foregroundColor
	plt.rcParams["ytick.color"] = foregroundColor
	plt.rcParams['axes.edgecolor'] = foregroundColor
	plt.rcParams['axes.linewidth'] = gridThickness

	axImage = fig.add_subplot(gs[0], label = 'ax_image')
	axImage.axis('off')
	if logPlot:
		allSkyImageLog = np.log(allSkyImage)
		vminVar = np.nanpercentile(allSkyImageLog, 99)
		vmaxVar = np.nanpercentile(allSkyImageLog, 5)
		pltIm = axImage.imshow(allSkyImageLog, alpha = 1, cmap='jet', label = 'ax_image', vmin = vminVar, vmax = vmaxVar)
	else:
		vminVar = np.nanpercentile(allSkyImage, 99)
		vmaxVar = np.nanpercentile(allSkyImage, 5)
		pltIm = axImage.imshow(allSkyImage, alpha = 1, cmap='jet', label = 'ax_image', vmin = vminVar, vmax = vmaxVar)
	axImage.axis('off')
	if colorBar:
		axColorBar = plt.subplot(gs[1])
		colorBarObj = plt.colorbar(pltIm, axColorBar)

		axColorBar.tick_params(which = 'minor', length = 2)
		axColorBar.tick_params(which = 'major', length = 4, width = 1)      
		axColorBar.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(10))
	   
	pltObj = fig.add_subplot(gs2[0], label = 'ax', polar = True)

	axImage.set_xlim((0, pixels))
	axImage.set_ylim((pixels, 0))
	pltObj.set_theta_zero_location("N")
	pltObj.set_theta_direction(1)

	plotStatus = [[True, __plotSkyObject(axImage, skyObj, pixels, skyObjColor, knownSources[idx])] if skyObj.alt.deg > 20. else [False, __plotSkyObject(axImage, skyObj, pixels, skyObjColor, knownSources[idx], offset = True)]  for idx, skyObj in enumerate(altAzObjects)]
	
	plt.rcParams["text.color"] = foregroundColor
	legend = axImage.legend( loc = 8, bbox_to_anchor = ( 0.5, -0.128 ), ncol = 4, framealpha = 0.0, fontsize = 14, title = str(obsTime)[:-4])
	legend.get_title().set_fontsize('22')

	for idx, skyText in enumerate(legend.get_texts()):
		if not plotStatus[idx][0]:
			plt.setp(skyText, color = 'red')
	radii = []

	for r in range(0, 90, 15): # r grid at 15 degree intervals
		radii.append(180 * np.cos(r * np.pi/180)) # plot the radii so as to display as an orthographic grid
		pltObj.set_rgrids(radii)
	if radialLabelAngle: # you would not want to put y ticks on 0 anyhow as it would be messy
		yLabel = [ '', '15' + u'\xb0', '30' + u'\xb0', '45' + u'\xb0', '60' + u'\xb0', '75' + u'\xb0' ]
		pltObj.set_yticklabels(yLabel, color = skyObjColor)
		pltObj.set_rlabel_position(radialLabelAngle)
	else:
		yLabel = []
		pltObj.set_yticklabels(yLabel)

	thetaticks = np.arange(0, 360, 45)
	pltObj.set_thetagrids(thetaticks, weight = 'bold', color = skyObjColor, fontsize = 18)
	pltObj.tick_params('x', pad = -35, rotation = 'auto')

	pltObj.grid(False, 'both', color = skyObjColor, linewidth = gridThickness)
	pltObj.patch.set(alpha = 0.0)
	plt.sca(axImage)

	if saveImage:
		if not os.path.exists(outputFolder):
			os.makedirs(outputFolder)

		plotFilename = "{6}{0}_{1}_sb{2}_mode{3}{4}_{5}MHz.png".format(dateTime, obsSite, subband, rcuMode, polarity, int(frequency/1e6), outputFolder)
		print("Saving output to {0}".format(plotFilename))
	 
		fig.savefig(plotFilename, facecolor=fig.get_facecolor(), edgecolor='none')
		plt.close(figNum)
		return plotFilename
	else:
		fig.show()
		plt.close(figNum)
		return

@cache
def cachedSkyCoords(name):
	return astropy.coordinates.SkyCoord.from_name(name)

def __plotSkyObject(axIm, skyObj, pixels, skyObjColor, sourceName, offset = False):
	if not offset:
		rho = np.sin(np.pi / 2. - skyObj.alt.rad )
		phi = skyObj.az.rad

		y, x = rho * np.cos(phi), rho * np.sin(phi)
		x, y = (pixels/2) - (pixels/2) * x, (pixels/2) - (pixels/2) * y
	else:
		x = 0. 
		y = 0.

	skyPlt = axIm.scatter(x, y, color = skyObjColor, marker = 'D', s = 50, label = u"{0} - Az={1}\xb0, El={2}\xb0".format(sourceName, round(skyObj.az.deg, 1), round(skyObj.alt.deg, 1)), alpha = 1)
	if not offset:
		textObj = axIm.annotate(sourceName, xy = (x,y), xytext = (x+2,y+2), color = skyObjColor, fontsize = 18)
		textObj.set_path_effects([mplPe.withStroke(linewidth = 5, foreground = 'w')])

#logPlot, skyObjColor, gridThickness, backgroundColor, foregroundColor, saveImage, radialLabelAngle, colorBar, obsSite = plotOptions
plotOptions = [False, 'black', .5, 'black', 'white', True, 0, True, 'Birr', './']
