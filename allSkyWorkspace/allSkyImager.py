

import ast
import scipy.constants
import scipy.ndimage
import astropy.coordinates
import astropy.time
import astropy.units as u
import datetime
import subprocess
import os
import collections
import numpy as np
import matplotlib as mpl
import matplotlib.patheffects as mplPe
import matplotlib.pyplot as plt

from functools import wraps
import multiprocessing as mp

from healpy.sphtfunc import Alm as almMap
import healpy

global vmaxCache
global vminCache

vmaxCache = collections.deque([], 20)
vminCache = collections.deque([], 20)

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
	# David McKenna: I don't want to touch this wizardry. Looks a lot like the code from Michiel Brentjens' lofar-antenna-positions though, so I supplemented it with his height code.
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
	height = r * np.cos(lat * np.pi / 180.) + arrayLoc[2] * np.sin(lat * np.pi / 180.) - wgs84_a * np.sqrt(1. - wgs84_e2 * np.sin(lat * np.pi / 180.) ** 2)


	return posX, posY, lon, lat, height, arrayLoc, antLocs

def parseBlitzFile(linesArray, keyword, refLoc = False):
	linesArray = [line.strip('\n') for line in linesArray]
	#print(linesArray)
	print(linesArray)
	triggerLine = [idx for idx, line in enumerate(linesArray) if line.startswith(keyword)][0]
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
	#print(splitTuples)
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
		lVec = lVec[np.newaxis, :] # (1, l)

	if len(mVec.shape) != 2:
		mVec = mVec[np.newaxis, :] # (1, m)

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
	for frame in np.arange(frameCount):
		correlationMatrixChan = correlationMatrix[..., frame] # (96,96)
		print("Processing Frame {0} of {1} at Frequency {2:.2F}E+06".format(frame + 1, frameCount, obsFreq / 1e6))

		tempProd = np.dot(conjWeight, correlationMatrixChan) # w.H * corr # (l,m, 1, 96)

		prodProd = np.multiply(tempProd.transpose((0,1,3,2)), weight).real
		skyView[..., frame] = np.sum(prodProd, axis = (2,3)) # (w.H * corr) * w, faster than loop below by about 50ms (running at 620ms, down from 2s in original implementation)

		#for l in range(mVec.size):
		#	for m in range(lVec.size):
		#		skyView[l, m, frame] = np.dot(tempProd[l, m], weight[l, m]).real[0, 0]
		#skyView[..., frame] = np.tensordot(tempProd, weight, ([3],[2])) # Dot products producing 201x201x1x201x201x1 arrays, thanks numpy.
	return skyView.transpose((1,0,2))

# -11.9
def generatePlots(inputCorrelations, antPos, plotOptions, dateArr, rcuMode, subband, multiprocessing = True, stationRotation = -11.9, plotX = True, plotY = True, mask = True, lVec = None, mVec = None, 
calibrationX = None, calibrationY = None, baselineLimits = None, stationLocation = None):
	inputCorrelationsX = inputCorrelations[..., 0]
	inputCorrelationsY = inputCorrelations[..., 1]

	antPos, rawAnts = antPos
	rawAnts = rawAnts[:, 0, :][:, np.newaxis]
	#raw_input([antPos.shape, rawAnts.shape, 'hi'])
	posX = antPos[..., 0]
	posY = antPos[..., 1]

	# Reference
	stackedArr = np.zeros(np.array(inputCorrelationsX.shape[:2]) * 2, dtype = 'complex')
	stackedArr[::2, ::2] = inputCorrelationsX[..., 0]
	stackedArr[1::2, 1::2] = inputCorrelationsY[..., 0]

	xxCorr = inputCorrelationsX
	yyCorr = inputCorrelationsY
	xyCorr = stackedArr[1::2, ::2]
	yxCorr = stackedArr[::2, 1::2]

	if baselineLimits:
		print('Test baseline limits.')
		if not baselineLimits[0]:
			baselineLimits[0] = 0.
			print('By not using autocorrelations the sky is a lot quieter: we are disabling log plotting as a result.')
			plotOptions[0] = False
		elif not baselineLimits[1]:
			baselineLimits[1] = 1e9
		minBaseline, maxBaseline = baselineLimits

		baselines = posX - posY.T
		print(posX.shape, posY.shape, np.median(np.abs(baselines)), np.percentile(np.abs(baselines), 66),  np.percentile(np.abs(baselines), 75),  np.percentile(np.abs(baselines), 80),  np.percentile(np.abs(baselines), 90))
		baselines = np.logical_or(np.abs(baselines) > maxBaseline, np.abs(baselines) < minBaseline)

		inputCorrelationsX[baselines] = 0.
		inputCorrelationsY[baselines] = 0.

	frequency = calcFreq(rcuMode, subband) * 1e6

	#raw_input((xxCorr + yyCorr).shape)
	test = swhtSkyImage(xxCorr + yyCorr, rawAnts, stationLocation, dateArr, 15, frequency)

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
			labelOptionsX[0] = str(dateArr[i])
			labelOptionsX[0] += '.{0:02d}'.format(i)
			labelOptionsX[-1] = figNum
			figNum += 1
			xFileLoc.append(plotAllSkyImage(allSkyImX[..., i], plotOptions, labelOptionsX, pixels, stationLocation, lVec, mVec))

		if plotOptions[5] and allSkyImX.shape[-1] > 20:
			filePrefix = xFileLoc[0].split(' ')[0][:-3]
			fileSuffix = '_'.join(xFileLoc[0].split('/')[-1].split('_')[1:])[:-4]
			print("Exporting frames to video at " + "./{0}{1}.mpg".format(filePrefix, fileSuffix))
			subprocess.call([ffmpegLoc, '-y',  '-r',  '20', '-pattern_type', 'glob',  '-i',  '{0}*{1}.png'.format(filePrefix, fileSuffix), '-vb', '50M', "{0}{1}.mpg".format(filePrefix, fileSuffix)])

		global vminCache
		global vmaxCache

		vminCache.clear()
		vmaxCache.clear()
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
			labelOptionsY[0] = str(dateArr[i])
			labelOptionsY[0] += '.{0:02d}'.format(i)
			labelOptionsY[-1] = figNum
			figNum += 1
			yFileLoc.append(plotAllSkyImage(allSkyImY[..., i], plotOptions, labelOptionsY, pixels, stationLocation, lVec, mVec))

		if plotOptions[5] and allSkyImY.shape[-1] > 20:
			filePrefix = yFileLoc[0].split(' ')[0][:-3]
			fileSuffix = '_'.join(yFileLoc[0].split('/')[-1].split('_')[1:])[:-4]
			print("Exporting frames to video at " + "./{0}{1}.mpg".format(filePrefix, fileSuffix))
			subprocess.call([ffmpegLoc, '-y',  '-r',  '20',  '-pattern_type', 'glob', '-i',  '{0}*{1}.png'.format(filePrefix, fileSuffix), '-vb', '50M', "{0}{1}.mpg".format(filePrefix, fileSuffix)])
		
		global vminCache
		global vmaxCache

		vminCache.clear()
		vmaxCache.clear()

		returnVar['Y'] = allSkyImY
		print(returnVar)

	if plotX and plotY:
		stokes_I = np.abs(allSkyImX) + np.abs(allSkyImY)
		labelOptionsI = labelOptions + ['I', 0]
		stokesILoc = []
		for i in range(stokes_I.shape[-1]):
			labelOptionsI[0] = str(dateArr[i])
			labelOptionsI[0] += '.{0:02d}'.format(i)
			labelOptionsI[-1] = figNum
			figNum += 1

			stokesILoc.append(plotAllSkyImage(stokes_I[..., i], plotOptions, labelOptionsI, pixels, stationLocation, lVec, mVec))


		if plotOptions[5] and stokes_I.shape[-1] > 20:
			filePrefix = stokesILoc[0].split(' ')[0][:-3]
			fileSuffix = '_'.join(stokesILoc[0].split('/')[-1].split('_')[1:])[:-4]
			print("Exporting frames to video at " + "./{0}{1}.mpg".format(filePrefix, fileSuffix))
			subprocess.call([ffmpegLoc, '-y',  '-r',  '20',  '-pattern_type', 'glob', '-i',  '{0}*{1}.png'.format(filePrefix, fileSuffix), '-vb', '50M', "{0}{1}.mpg".format(filePrefix, fileSuffix)])


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
		labelOptions[0] = dateArr[i] + 1
		labelOptions[-1] = figNum
		figNum += 1
		fileLoc.append(plotAllSkyImage(allSkyIm[..., i], plotOptions, labelOptions, pixels, stationLocation, lVec, mVec))

	if plotOptions[5] and allSkyIm.shape[-1] > 20:
		filePrefix = fileLoc[0].split(' ')[0]
		fileSuffix = '_'.join(fileLoc[0].split('/')[-1].split('_')[1:])[:-4]
		print("Exporting frames to video at " + "./{0}{1}.mpg".format(filePrefix, fileSuffix))
		subprocess.call([ffmpegLoc, '-y',  '-r',  '20', '-pattern_type', 'glob',  '-i',  '{0}*{1}.png'.format(filePrefix, fileSuffix), '-vb', '50M', "{0}{1}.mpg".format(filePrefix, fileSuffix)])

	#return __, __, __...

def __processAllSkyIm(inputCorrelations, posX, posY, frequency, lVec, mVec, stationRotation, mask, plotOptions, labelOptions):
	allSkyIm = corrToSkyImage(inputCorrelations, posX, posY, frequency, lVec, mVec)
	#allSkyIm = swhtSkyImage(inputCorrelations, posX, posY, frequency, lVec, mVec)
	allSkyIm = np.rot90(allSkyIm, 3)
	allSkyIm = np.flipud(allSkyIm)
	allSkyIm = scipy.ndimage.rotate(allSkyIm, stationRotation, mode = 'constant', cval = 100., reshape = False)
	
	allSkyIm[mask] = np.nan
	allSkyIm[allSkyIm == 100.] = np.nan

	# Extra returns for if we start multithreading
	return allSkyIm, plotOptions, labelOptions

def swhtSkyImage(inputCorrelations, xyz, stationLocation, timesArr, lMax, frequency):
	#Currently non functional. Stencilling in code at the moment, will flesh out input parameters later.
	xyz_raw = xyz.copy()
	xyz = xyz - xyz.transpose(1,0,2)
	xyz_triu = np.triu_indices(xyz.shape[0])

	inputCorrelations[np.eye(inputCorrelations.shape[0], dtype = bool)] = 0.
	inputDropped = inputCorrelations[xyz_triu]
	xyz = xyz[xyz_triu]

	print(stationLocation)
	trueTelescopeLoc = astropy.coordinates.EarthLocation( lat = stationLocation[0] * u.deg, lon = stationLocation[1] * u.deg, height = stationLocation[2] * u.m )
	telescopeLoc = astropy.coordinates.EarthLocation(lat = 0 * u.deg, lon = -90 * u.deg, height = 0 * u.m) # Not sure why, used by Griffin as a baseline for calculations and it works.
	uvw = np.zeros([0] + list(xyz.shape))
	print(uvw.shape)

	zenithArr = []
	altAzArr = []
	for timestamp in timesArr:
		print(timestamp)
		obsTime = astropy.time.Time(timestamp, scale = 'utc', location = telescopeLoc)
		utcObsTime = astropy.time.Time(timestamp, scale = 'utc')

		altAz = astropy.coordinates.AltAz(az = 0. * u.deg, alt = 90. * u.deg, location = trueTelescopeLoc, obstime = utcObsTime)
		zenith = altAz.transform_to(astropy.coordinates.Galactic)

		print(zenith.l.deg, zenith.b.deg)
		sidAngle = float(obsTime.sidereal_time('mean').radian)
		zenithArr.append(zenith)
		altAzArr.append(altAz)

		sinSA = np.sin(sidAngle)
		cosSA = np.cos(sidAngle)
		rotationMatrix = np.array([[sinSA, cosSA, 0.],
						[-1. * cosSA, sinSA, 0.],
						[0., 0., 1.]])

		#raw_input([sidAngle, obsTime, sinSA, cosSA])
		uvw = np.concatenate([uvw, np.dot(rotationMatrix, xyz.T).T[np.newaxis]], axis = 0)

	r, theta, phi = cartToSpherical(uvw)
	r = r[0] # All samples have the same r values

	k_0 = 2. * np.pi * frequency / scipy.constants.c

	#k_0 = frequency / scipy.constants.c
	preFac = 2. * (k_0 ** 2) / np.pi
	kRVec = k_0 * r

	arraySize = almMap.getsize(lMax - 1)
	results = np.empty(arraySize, dtype = complex)
	ores = np.zeros((lMax, 2* lMax + 1), dtype = complex)
	jv = []
	sph = []

	kZeroBool = kRVec == 0.


#	mGrid = xyz_raw.shape[0] - 1 - np.mgrid[0:xyz_raw.shape[0], 0:xyz_raw.shape[1]].transpose(1,2,0)[xyz_triu]
#	index = [almMap.getidx(lMax, lm[0], lm[1]) for lm in mGrid]
#	_cython_loop(kRVec, kZeroBool, phi, theta, lMax, results, preFac, index)
	'''
	try:
		knownBodies = ['Sun', 'Jupiter', 'Moon', 'Uranus', 'Neptune']
		plotLists = [[]] * len(knownBodies)
		bodyList = [[]] * len(timesArr)

		altAzRef = astropy.coordinates.AltAz(obstime = obsTime, location = telescopeLoc)
		altAzObjArr = []

		for timeIdx, obsTimeStr in enumerate(timesArr):
			print(timeIdx)
			#altAzRef.obstime = obsTimeStr # If only they allowed us to just change the time...
			bodyList[timeIdx] = [astropy.coordinates.get_body(sourceName, obsTime, telescopeLoc) for sourceName in knownBodies]

			altAzObjArr.append([skyObj.transform_to(altAzArr[timeIdx]) for skyObj in bodyList[timeIdx]])
		
		for bodyIdx, timeArr in enumerate(altAzObjArr):
			print(timeArr)
			for bodyObj in timeArr:
				print(plotLists[bodyIdx], bodyIdx, bodyObj.alt.deg)
				if bodyObj.alt.deg > 20.:
					plotLists[bodyIdx].extend((bodySkyObj.galactic.l.deg, bodySkyObj.galactic.b.deg))
				else:
					plotLists[bodyIdx].extend((np.nan, np.nan))

	except NameError:
		plotLists = []
	'''


	for l in range(lMax):
		j_l = scipy.special.spherical_jn(l, kRVec) # (rLen,)
		j_l[kZeroBool] = 0. # nan otherwise
		
		lRev = (4 * np.pi * (-1j)**l) # (1,)
		print(l)
		jv.append(np.sum(j_l))
		visBessel = inputDropped.flatten() * np.repeat(j_l, uvw.shape[0], axis = 0) 
		for m in range(l + 1):
			y_lm_star = np.conj(scipy.special.sph_harm(m, l, phi.T.flatten(), theta.T.flatten()))# (timeSamples * nants ** 2)
			print(phi.shape, theta.shape)
			#y_lm_star = np.conj( SWHT.Ylm.Ylm( l, m, phi.T.flatten(), theta.T.flatten())) # Same output as scipy...
			sph.append(np.sum(y_lm_star))
			#print(j_l, y_lm_star)

			resultsIndex = almMap.getidx(lMax - 1, l, m)
			#rI.append(resultsIndex)
			#print(inputCorrelations.shape, inputCorrelations[xyz_triu].shape)
			#print(inputCorrelations[xyz_triu].shape, j_l[:, np.newaxis].shape, y_lm_star.T.shape)

			#print(l, m, resultsIndex, lMax)
			results[resultsIndex] = preFac * np.sum(visBessel * y_lm_star) / lRev
			#print(inputDropped.shape, j_l.shape, y_lm_star.shape)
			#results[resultsIndex] += preFac * np.sum(np.conj(inputDropped) * j_l[:, np.newaxis] * y_lm_star.T) / lRev
			#ores[l, l + m] = results[resultsIndex]
			'''
	print(jv[0])
	#np.save('./inp2.npy', inputCorrelations[xyz_triu])
	#np.save('./jv1.npy', jv)
	#np.save('./sph1.npy', sph)
	plt.plot(jv)
	plt.figure(2)
	plt.plot(sph)
	plt.figure(3)
	plt.plot([val.imag for val in sph], c = 'r')
	#plt.twinx().plot(theta.T, c = 'g')
	plt.show()
	print(sph[0])
	plt.plot(ores)
	plt.show()
	'''
	'''
	results = np.zeros(almMap.getsize(lMax), dtype = complex)
	for l in np.arange(lMax): #increase lmax by 1 to account for starting from 0
		#jvVals = np.reshape(sphBj(l, kr.flatten(), autos=True), kr.shape) #compute Bessel function radius values
		jvVals = np.repeat(SWHT.swht.sphBj(l, kRVec.flatten(), autos=False)[np.newaxis], phi.shape[0], axis = 0) #compute Bessel function radius values
		lRev = (4 * np.pi * (- 1j)** -l) # (1,)
		for m in np.arange(l, l+1):
			print(l, m)
			#Compute visibility spherical harmonic coefficients according to SWHT, i.e. multiply visibility by spherical wave harmonics for each L&M and sum over all baselines.
			#Note that each non-zero baseline is effectively summed twice in the precedi(In the MNRAS letter image the NZ baselines were only weighted once, i.e. their conjugate baselines were not summed.)
			#spharmlm = np.repeat(np.conj(Ylm.Ylm(l, m, phi[:, sbIdx:sbIdx+1], theta[:, sbIdx:sbIdx+1])), nsbs, axis=1) #spherical harmonics only needed to be computed once for all baselines, independent of observing frequency
			spharmlm = np.conj( SWHT.Ylm.Ylm( l, m, phi, theta))
			resultsIndex = almMap.getidx(lMax, l, l + m)
			print(inputCorrelations[xyz_triu].shape, jvVals.shape, spharmlm.shape)
			results[resultsIndex] = np.mean(preFac * np.nansum(inputCorrelations[xyz_triu].T * jvVals * spharmlm, axis = 1) / lRev)
	'''
	#saveCache(True)

	print('Map Start')
	#results = SWHT.util.array2almVec(ores)
	#hpMap = healpy.ma(healpy.alm2map(results, 64))
	hpMap = np.log2(np.abs(healpy.alm2map(results, 64)))

	galRotator = healpy.rotator.Rotator(coord = ['C','G'])
	hpMap = galRotator.rotate_map(hpMap)
	hpMapMasked = healpy.ma(hpMap)
	hpMapMasked.mask = np.logical_not(hpMapMasked)
	print('Map End')

	#np.save('./hpMap1.npy', hpMap)
	#np.save('./res1.npy', results)
	#plt.plot(hpMap)
	#plt.show()

	mask = np.zeros(healpy.nside2npix(64), dtype=np.bool)
	print(mask.shape)
	pixel_theta, pixel_phi = healpy.pix2ang(64, np.arange(healpy.nside2npix(64)), lonlat = True)
	

	# Fuck everyhting about this coordinate system...
	invLoc = pixel_theta > 180.
	# We cannot simply do reverse orders as pixels are allocated on a per-area basis, not a numbers basis
	
	pixel_theta[invLoc] =  (360. - pixel_theta[invLoc]) * -1.
	pixel_theta[~invLoc] = pixel_theta[~invLoc]

	#healpy.mollview(pixel_theta)
	#healpy.mollview(pixel_phi)
	#plt.show()
	#pixel_theta = 360. - pixel_theta

	galCoordLon = np.array([skyCoord.l.deg for skyCoord in zenithArr])
	galCoordLat = np.array([skyCoord.b.deg for skyCoord in zenithArr])

	healpy.mollview(pixel_theta)
	healpy.projplot(galCoordLon, galCoordLat, lonlat = True)
	healpy.visufunc.graticule()

	healpy.mollview(pixel_phi)
	healpy.projplot(galCoordLon, galCoordLat, lonlat = True)
	healpy.visufunc.graticule()


	print(zip(list(galCoordLon), list(galCoordLat)))
	plt.figure()
	plt.plot(galCoordLon)
	plt.twinx().plot(galCoordLat, c = 'r')
	plt.show()

	galCoordSampled = np.deg2rad(np.vstack([galCoordLon, galCoordLat])[:, np.newaxis])
	pixelCoord = np.deg2rad(np.vstack([pixel_theta, pixel_phi])[..., np.newaxis])

	a_1 = np.square(np.sin(0.5 * (galCoordSampled[0] - pixelCoord[0]))) + np.cos(galCoordSampled[0]) * np.cos(pixelCoord[0]) * np.square(np.sin(0.5 * (galCoordSampled[1] - pixelCoord[1])))
	deltaLoc = 2. * np.arctan2(np.sqrt(a_1), np.sqrt(1 - a_1))

	healpy.mollview(np.min(deltaLoc, axis = 1))
	healpy.projplot(galCoordLon, galCoordLat, lonlat = True)
	plt.show()
	#deltaLoc = np.arccos(np.sin(galCoordSampled[0]) * np.sin(pixelCoord[0]) + np.cos(galCoordSampled[0]) * np.cos(pixelCoord[0]) * np.cos(galCoordSampled[1] - pixelCoord[1]))
	
	mask = np.logical_not(np.any(np.rad2deg(deltaLoc) < 90, axis = 1)) # Assume 70 degrees field of view in each direction. phi modified to range from0  -> 360 as well for equal weighting.
	#mask = np.logical_not(pixel_phi > -30.)
	hpMapMasked.mask = mask

	#print(mask.shape, hpMap.shape, deltaLoc.shape, galCoordSampled.shape, pixelCoord.shape)
	#print(hpMap.shape, hpMap.shape)
	#hpMap.mask = mask

	#print(hpMap.shape, hpMap.shape, hpMap.real.shape)


	print('Map Reals')

	'''
	for i in range(galCoordSampled.shape[2]):
		healpy.mollview(deltaLoc[:, i])
		print(galCoordLon.shape, galCoordLon[i])
		healpy.projscatter([galCoordLon[i]], [galCoordLat[i]], lonlat = True, c = 'g')
		plt.savefig('./{0}.png'.format(i))
 		plt.clf()

 	plt.close('all')
	'''
	'''
	skyCoord = cachedSkyCoords('Cas A')
	print(skyCoord, skyCoord.galactic)
	for bodyArr in plotLists:
		print(bodyArr)
		healpy.projplot(np.array(bodyArr).T, lonlat = True)
	'''
	healpy.mollview(hpMap.real)
	healpy.visufunc.graticule()
	plt.show()
	healpy.mollview(hpMapMasked.mask)
	healpy.projplot(galCoordLon, galCoordLat, lonlat = True, coord = 'G')
	healpy.visufunc.graticule()
	plt.show()
	healpy.mollview(hpMapMasked)
	healpy.visufunc.graticule()

	skyCoord = astropy.coordinates.SkyCoord.from_name('Cas A')
	print(skyCoord.galactic)
	healpy.projscatter([skyCoord.galactic.l.deg], [skyCoord.galactic.b.deg], lonlat = True, c = 'r')

	skyCoord = astropy.coordinates.SkyCoord.from_name('Cyg A')
	print(skyCoord.galactic)
	healpy.projscatter([skyCoord.galactic.l.deg], [skyCoord.galactic.b.deg], lonlat = True, c = 'g')
	
	skyCoord = astropy.coordinates.SkyCoord.from_name('Vir A')
	print(skyCoord.galactic)
	healpy.projscatter([skyCoord.galactic.l.deg], [skyCoord.galactic.b.deg], lonlat = True, c = 'b')

	skyCoord = astropy.coordinates.SkyCoord.from_name('Polaris')
	print(skyCoord.galactic)
	healpy.projscatter([skyCoord.galactic.l.deg], [skyCoord.galactic.b.deg], lonlat = True, c = 'm')

	plt.show()
	#healpy.mollview(hpMap.real,deg=True, rot = [90, 0], coord = 'CG')
	#healpy.mollview(hpMap.real, deg=True, rot = [0, 90], coord = 'CG')

	exit;

'''
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _cython_loop(double[:] kRVec, double[:] kZeroBool, double[:, :] phi, double[:, :] theta, int lMax, complex[:] results, float preFac, int[:] index) nogil:

	cdef int l, m
	cdef complex lRev
	cdef complex[:, :, :] j_l, y_lm_star

	for l in prange(lMax):
		j_l = csc.spherical_jn(l, kRVec)
		j_l[kZeroBool] = 0.
		lRev = (4 * np.pi * (-1j)**l)
		for m, in range(i + 1):
			y_lm_star = (csc.sph_harm(m, l, phi, theta)).conjugate()
			results[index[l,m]] = preFac * (inputDropped * j_l * y_lm_star.transpose()).sum() / lRev
'''

def cartToSpherical(uvw):
	print(uvw.shape)

	# r = sqrt(sq(x) + sq(y) + sq(z))
	# theta = acos(z / r)
	# phi = atan(y / x)
	
	#r = np.sqrt(np.sum(np.square(uvw), axis = 2))
	r = np.sqrt(np.sum(np.square(uvw[0]), axis = (1)))
	r = np.repeat(r[np.newaxis], uvw.shape[0], axis = 0)

	theta = np.arccos(uvw[..., 2] / r)
	phi = np.arctan2(uvw[..., 1], uvw[..., 0]) + np.pi # Scipy requires 0 -> 2pi

	zerothBaselines = np.where(r == 0.)
	theta[zerothBaselines] = np.pi / 2.
	phi[zerothBaselines] = np.pi
	return r, theta, phi

def plotAllSkyImage(allSkyImage, plotOptions, labelOptions, pixels, stationLocation, lVec, mVec):
	logPlot, skyObjColor, gridThickness, backgroundColor, foregroundColor, saveImage, radialLabelAngle, colorBar, obsSite, outputFolder = plotOptions
	dateTime, rcuMode, subband, frequency, polarity, figNum = labelOptions

	telescopeLoc = astropy.coordinates.EarthLocation( lat = stationLocation[0] * u.deg, lon = stationLocation[1] * u.deg, height = stationLocation[2] * u.m )
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

	print(figNum)
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

	global vmaxCache
	global vminCache
	
	lCoord, mCoord = np.meshgrid(lVec.squeeze(), mVec.squeeze())

	print(lCoord, lCoord.shape, mCoord.shape, allSkyImage.shape)
	if logPlot:
		allSkyImageLog = np.log10(allSkyImage)
		#vmaxVar = np.nanpercentile(allSkyImageLog, 99)
		#vminVar = np.nanpercentile(allSkyImageLog, 33)

		#vminCache.append(vminVar)
		#vmaxCache.append(vmaxVar)

		#vminVar = np.mean(vminCache)
		#vmaxVar = np.mean(vmaxCache)

		pltIm = axImage.imshow(allSkyImageLog, cmap='jet', label = 'ax_image', interpolation = 'nearest')
	else:
		vmaxVar = np.nanpercentile(allSkyImage, 99)
		vminVar = np.nanpercentile(allSkyImage, 5)

		#vminCache.append(vminVar)
		#vmaxCache.append(vmaxVar)

		#vminVar = np.mean(vminCache)
		#vmaxVar = np.mean(vmaxCache)
		pltIm = axImage.imshow(allSkyImage, cmap='jet', label = 'ax_image', interpolation = 'nearest')

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
		plotFilename = plotFilename.replace(' ', '_').replace(':', '')
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
