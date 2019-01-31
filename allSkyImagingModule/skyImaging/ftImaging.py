"""Process correlations through Fourier Transform-based method with these functions.
"""
import scipy.constants
import scipy.ndimage

import multiprocessing as mp
import numpy as np
from skyPlotter import ftPlot

from workers import ftWorkers

def ftProcessCorrelations(procMethod, rfiFlag, corrDict, processDict, options, posArr, frequency, stationRotation, stationLocation, labelOptions, metaDataArr):
	"""Head function that sets up values needed by both the FFT and DFT.
	
	Args:
	    procMethod (procMethod): Processing method (DFT, FFt, FFT/gauss/0.5 for FFT/method/kernalextrafeature)
	    rfiFlag (bool): Bool to enable RFI mode (interactive plotting)
	    corrDict (dict): Dictionary of elements in the format 'CORRELATION': (nAnts x nAnts x nSamples array of correlations)
	    processDict (dict): Dictionary of correlation times ('XX', 'XY'...) and values of whether or not to process them
	    options (dict): Standard options dictionary
	    posArr (np.ndarray): Antenna locations in station coordinates
	    frequency (float): Frequency of observation (in Hz)
	    stationRotation (float): Rotation of the station with respect to North
	    stationLocation (list): [lon, lat, alt] of station
	    labelOptions (list): List of useful values for plotting
	    metaDataArr (list): List of str(dict) from observation group attributes
	
	Returns:
	    dict: Output sky images
	"""

	# No w-term used for FT operations for speedup / due to negligible gains
	posX, posY, __ = posArr

	lVec, mVec = options['imagingOptions']['pixelCount']
	fov = options['imagingOptions']['fieldOfView']
	rfiFlag = options['rfiMode']

	if rfiFlag:
		print('RFI Mode enabled, we will be forcing the FT method of choice (or FFT) and providing a GUI.')

	# Get the observing angle + generate the l, m values to FT
	vecRange = np.sin([-1. * fov / 2., fov / 2.])

	lVec = np.arange(vecRange[0], vecRange[1] + 0.01, (2. / lVec) or 0.008 ) # Default to 256 pixels
	mVec = np.arange(vecRange[0], vecRange[1] + 0.01, (2. / mVec) or 0.008 )

	# If needed, mask outside the fov
	mask = np.zeros([lVec.size, mVec.size], dtype = bool)
	if options['imagingOptions']['maskOutput']:
		dist = np.sqrt(np.square(np.meshgrid(lVec)) + np.square(np.meshgrid(mVec)).T)
		mask[dist >= vecRange[1] * 1.0005] = True # Offset to account for pixelisation


	if 'dft' in procMethod:
		imagingFunc = dftImage

	elif 'fft' in procMethod:
		imagingFunc = fftImage

	else:
		print('Unknown Processing Method {0}, forcing DFT. (Debug: rfiFlag={1})'.format(procMethod, rfiFlag))
		options['imagingOptions']['method'] = 'dft'
		imagingFunc = dftImage

	# Process the correlations through the provided method
	procDict = {}
	for key, value in processDict.items():
		if value:
			allSkyImage = processAllSkyIm(imagingFunc, corrDict[key], posX, posY, frequency, lVec, mVec, stationRotation, mask, options)
			print("{0} Polarisation Processed.".format(key))
			procDict[key] = allSkyImage

	# Process the requirested all sky output types
	for method in options['imagingOptions']['correlationTypes']:
		if len(method) == 2:
			allSkySave = procDict[method]
			allSkyImage = procDict[method].real
		elif 'I' is method:
			allSkySave = (procDict['XX'] + procDict['YY'])
			allSkyImage = allSkySave.real
		elif 'Q' is method:
			allSkySave = (procDict['XX'] - procDict['YY'])
			allSkyImage = allSkySave.real
		elif 'U' is method:
			allSkySave = (procDict['XY'] + procDict['YX'])
			allSkyImage = allSkySave.real
		elif 'V' is method:
			allSkySave = (procDict['YX'] - procDict['XY'])
			allSkyImage = allSkySave.imag
		else:
			print('How did we get this far with a broken method \'{0}\'? Processing the first value in procDict so it hasn\'t gone to waste...'.format(method))
			allSkyImage = procDict.values()[0]
			method = procDict.keys()[0]

		if options['imagingOptions']['ftSubtractBackground']:
			allSkyImage -= np.nanmin(allSkyImage) * 0.999999

		procDict[method] = allSkySave

		if options['plottingOptions']['plotImages']:
			ftPlot(allSkyImage, options, labelOptions + [method, 0], stationLocation, lVec, mVec, metaDataArr)

	return procDict

def processAllSkyIm(imagingFunc, inputCorrelations, posX, posY, frequency, lVec, mVec, stationRotation, mask, options):
	"""Generate an all sky image and rotate/mask as needed
	
	Args:
	    inputCorrelations (TYPE): Description
	    posX (TYPE): Description
	    posY (TYPE): Description
	    frequency (TYPE): Description
	    lVec (TYPE): Description
	    mVec (TYPE): Description
	    stationRotation (TYPE): Description
	    mask (TYPE): Description
	    plotOptions (TYPE): Description
	    labelOptions (TYPE): Description
	
	Returns:
	    np.ndarray: All Sky Image
	"""
	allSkyIm = imagingFunc(inputCorrelations, posX, posY, frequency, lVec, mVec, options)
	allSkyIm = np.rot90(allSkyIm, 3)
	allSkyIm = np.flipud(allSkyIm)
	allSkyIm.real = scipy.ndimage.rotate(allSkyIm.real, stationRotation, mode = 'constant', cval = 100., reshape = False)
	allSkyIm.imag = scipy.ndimage.rotate(allSkyIm.imag, stationRotation, mode = 'constant', cval = 100., reshape = False)
	
	allSkyIm[mask] = np.nan
	allSkyIm[allSkyIm == 100.] = np.nan

	return allSkyIm

def dftImage(correlationMatrix, posX, posY, obsFreq, lVec, mVec, options):
	"""DFT methodology for generating all sky image
	
	Args:
	    correlationMatrix (np.ndarray): Input correlations (nAnts x nAnts x nSamples array of correlations)
	    posX (np.array): Antenna location x components
	    posY (np.array): Antenna location y components
	    obsFreq (float): Frequency of observation (Hz)
	    lVec (np.array): l-Values to sample (cosine of theta-x)
	    mVec (np.array): l-Values to sample (cosine of theta-y)
	
	Returns:
	    np.ndarray: Output sky image/
	"""
	# corrMat: (96,96,nChan)
	frameCount = correlationMatrix.shape[2]

	# Ensure right shape for broadcasting
	lVec, mVec, posX, posY = __ftVectorChecks(lVec, mVec, posX, posY)

	# Initialise an array to fill with out results
	skyView = np.zeros([lVec.size, mVec.size, frameCount], dtype = np.complex128) # (l,m,nChan)
	
	wavelength = scipy.constants.c / obsFreq
	k = (2 * np.pi) / wavelength

	wX = np.exp(-1j * k * posX * lVec) # (l, 96)
	wY = np.exp(-1j * k * posY * mVec) # (m, 96)
	weight = np.multiply(wX[:, np.newaxis, :], wY[:, :, np.newaxis]).transpose((1,2,0))[..., np.newaxis] # (l,m,96,1)
	conjWeight = np.conj(weight).transpose((0,1,3,2)) # (l,m,1,96)

	# Should be able to fully vectorise this over all channels, should give a nice speedup if we pass frames as alternative channels... for another day.
	if options['multiprocessing']:
		processCount = int(mp.cpu_count()-1) # Leave a thread for OS / disk I/O, etc.
		mpPool = mp.Pool(processes = processCount)
		
		# Attempting to pass more than ~350 correlations per thread in a 3 thread/16GB system is a OOM situation -- limit each set to 200
		altCalc = int(np.ceil(frameCount / 200.))

		# In the case that we do use the alternative limit, spin up a multiple of n_proc threads to ensure a minimal amount of cores are idle.
		if (altCalc > processCount):
			altCalc += altCalc % processCount
			processCount = altCalc - 1

		chunkCount = max([altCalc, processCount])
		fragments = np.array_split(np.arange(frameCount), chunkCount) # int 200: err on the side of caution and use another worker thread if needed.
		callBacks = [mpPool.apply_async(ftWorkers.dftWorker, args = ([idx, correlationMatrix[..., fragmentChunk], obsFreq, weight, conjWeight, skyView[..., fragmentChunk]])) for idx, fragmentChunk in enumerate(fragments)]

		mpPool.close()
		mpPool.join()

		for asyncResult in callBacks:
			idx, skyViewChunk = asyncResult.get()
			skyView[..., fragments[idx]] = skyViewChunk

	else:
		skyView = ftWorkers.dftWorker(0, correlationMatrix, obsFreq, weight, conjWeight, skyView)[1]

	# Swap x/y axis (why when we just do an ndflip later?)
	return skyView.transpose((1,0,2))

def fftImage(correlationMatrix, posX, posY, obsFreq, lVec, mVec, options):
	"""Currently nonfunctional; will return later if I can.
	
	Args:
	    correlationMatrix (TYPE): Description
	    posX (TYPE): Description
	    posY (TYPE): Description
	    obsFreq (TYPE): Description
	    lVec (TYPE): Description
	    mVec (TYPE): Description
	    kernel (str, optional): Description
	    generalSettings (None, optional): Description
	
	Returns:
	    TYPE: Description
	"""

	raise RuntimeError("Non-functional.")
	kernel = options['imagingOptions']['method']

	# corrMat: (96,96,nChan)
	frameCount = correlationMatrix.shape[2]


	# Ensure right shape for broadcasting
	lVec, mVec, posX, posY = __ftVectorChecks(lVec, mVec, posX, posY)

	wavelength = scipy.constants.c / obsFreq
	k = (2 * np.pi) / wavelength
	xy = np.dstack([posX, posY]) # (96, 1, 2)
	uv = (xy - xy.transpose((1,0,2))) * k # (96, 96, 2)
	uv = uv.transpose((2,0,1)).reshape(2, -1) # (2, 96 * 96)

	uv += np.abs(np.min(uv))
	uvMax = np.max(np.abs(uv))
	uv /= uvMax

	skyViewPixels = int(np.max([lVec.size, mVec.size]))
	pixelBuffer = int(0.1 * skyViewPixels)
	pixelBuffer += pixelBuffer % 2 # Ensure it's an even number
	print(pixelBuffer)
	uv *= skyViewPixels
	uv += pixelBuffer / 2 # 5% buffer

	skyViewPixels += pixelBuffer


	uv = uv.reshape(2, -1)
	correlationMatrix = correlationMatrix.reshape(1, -1, frameCount) # (1, 96 * 96, nChan)
	skyView = np.zeros([skyViewPixels, skyViewPixels, frameCount], dtype = np.complex128)


	if 'gaus' in kernel:
		gaussCoeff = float(kernel.split('/')[-1])

		pixelBoundary = int((skyView.shape[0] / 1.05) * 0.05)
		sampleArea = min(max(int(gaussCoeff * 10), 1), pixelBoundary) # Samples an area effectively ~5x the length/width of the Gaussian width

		print(sampleArea, gaussCoeff, pixelBoundary)
		guassMatrix = np.square(np.array([[gaussCoeff, 0.], [0., gaussCoeff]]))
		convKernel = lambda sampledPoints: (1. / np.sqrt(8. * np.pi) / gaussCoeff ) * np.sum(np.exp(-1. * np.dot(guassMatrix, np.square(sampledPoints))), axis = (0))

	else:
		convKernel = lambda sampledPoints: np.ones([sampledPoints.shape[0], sampledPoints.shape[2]])

	if 'rect' in kernel:
		sampleArea = 1

	sampleCache = np.mgrid[-sampleArea + 1:sampleArea, -sampleArea + 1:sampleArea].reshape(2, -1).astype(float).T[..., np.newaxis] # (nSamples, 2, 1)

	# Find the closest points in the sampled array
	offsets = (uv % 1.)[np.newaxis] # (1, 2, 96 * 96)

	offsets[offsets > 0.5] += 0.5

	sampleOffsets = sampleCache + offsets # (nSamples, 2, 96 * 96)
	uvIdx = (sampleOffsets + uv).astype(int)


	convolvedWeight = convKernel(sampleOffsets) # (nSamples, 96 * 96)


	skyView = ftWorkers.fftWorker(0, correlationMatrix, uvIdx, convolvedWeight, skyView)[1]

	skyView = skyView.transpose((1,0,2))[pixelBuffer / 2: -pixelBuffer / 2, pixelBuffer / 2 : -pixelBuffer / 2, ...]
	return skyView




def __ftVectorChecks(lVec, mVec, posX, posY):
	"""Summary
	
	Args:
	    lVec (TYPE): Description
	    mVec (TYPE): Description
	    posX (TYPE): Description
	    posY (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	if len(lVec.shape) != 2:
		lVec = lVec[np.newaxis, :] # (1, l)

	if len(mVec.shape) != 2:
		mVec = mVec[np.newaxis, :] # (1, m)

	if len(posX.shape) != 2:
		posX = posX[:, np.newaxis] # (96, 1)

	if len(posY.shape) != 2:
		posY = posY[:, np.newaxis] # (96, 1)

	return lVec, mVec, posX, posY
