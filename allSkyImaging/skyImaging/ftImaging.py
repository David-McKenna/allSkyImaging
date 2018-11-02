"""Summary
"""
import multiprocessing as mp
import numpy as np

import ftWorkers
import skyPlotter

def ftprocessCorrPlot():
	"""Summary
	"""
	labelOptions = labelOptions + [polChar, 0]
	allSkyIm, __, __ = __processAllSkyIm(inputCorrelationsArr, posX, posY, frequency, lVec, mVec, stationRotation, mask, plotOptions, labelOptions)
	print("{0} Polarisation Processed, begining plotting".format(polChar))

	fileLoc = []
	for i in range(allSkyIm.shape[-1]):
		labelOptions[0] = dateArr[i] + 1
		labelOptions[-1] = figNum
		figNum += 1
		fileLoc.append(skyPlotter.plotAllSkyImage(allSkyIm[..., i], plotOptions, labelOptions, pixels, stationLocation, lVec, mVec))

	if plotOptions[5] and allSkyIm.shape[-1] > 20:
		filePrefix = fileLoc[0].split(' ')[0]
		fileSuffix = '_'.join(fileLoc[0].split('/')[-1].split('_')[1:])[:-4]
		print("Exporting frames to video at " + "./{0}{1}.mpg".format(filePrefix, fileSuffix))
		subprocess.call([ffmpegLoc, '-y',  '-r',  '20', '-pattern_type', 'glob',  '-i',  '{0}*{1}.png'.format(filePrefix, fileSuffix), '-vb', '50M', "{0}{1}.mpg".format(filePrefix, fileSuffix)])

	

def __processAllSkyIm(inputCorrelations, posX, posY, frequency, lVec, mVec, stationRotation, mask, plotOptions, labelOptions):
	"""Summary
	
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
	    TYPE: Description
	"""
	allSkyIm = corrToSkyImage(inputCorrelations, posX, posY, frequency, lVec, mVec)
	allSkyIm = np.rot90(allSkyIm, 3)
	allSkyIm = np.flipud(allSkyIm)
	allSkyIm = scipy.ndimage.rotate(allSkyIm, stationRotation, mode = 'constant', cval = 100., reshape = False)
	
	allSkyIm[mask] = np.nan
	allSkyIm[allSkyIm == 100.] = np.nan

	return allSkyIm, plotOptions, labelOptions

def dftImage(correlationMatrix, posX, posY, obsFreq, lVec, mVec):
	"""Summary
	
	Args:
	    correlationMatrix (TYPE): Description
	    posX (TYPE): Description
	    posY (TYPE): Description
	    obsFreq (TYPE): Description
	    lVec (TYPE): Description
	    mVec (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	# corrMat: (96,96,nChan)
	frameCount = correlationMatrix.shape[2]

	# Ensure right shape for broadcasting
	lVec, mVec, posX, posY = __ftVectorChecks(lVec, mVec, posX, posY)

	# Initialise an array to fill with out results
	skyView = np.zeros([lVec.size, mVec.size, frameCount]) # (l,m,nChan)
	
	wavelength = scipy.constants.c / obsFreq
	k = (2 * np.pi) / wavelength

	wx = np.exp(-1j * k * posX * lVec) # (l, 96)
	wy = np.exp(-1j * k * posY * mVec) # (m, 96)
	weight = np.multiply(wx[:, np.newaxis, :], wy[:, :, np.newaxis]).transpose((1,2,0))[..., np.newaxis] # (l,m,96,1)
	conjWeight = np.conj(weight).transpose((0,1,3,2)) # (l,m,1,96)

	# Should be able to fully vectorise this over all channels, should give a nice speedup if we pass frames as alternative channels... for another day.
	if multiprocessing:
		fragments = np.array_split(np.arange(frameCount))
		processCount = int(mp.cpu_count()-1)
		mpPool = mp.Pool(processes = processCount)
		callBacks = [mpPool.apply_async(ftWorkers.dftWorker, args = ([idx, correlationMatrix[..., fragmentChunk], obsFreq, conjWeight, skyView[..., fragmentChunk]])) for idx, fragmentChunk in enumerate(fragments)]

		mpPool.close()
		mpPool.join()

		for asyncResult in callBacks:
			idx, skyViewChunk = asyncResult.get()
			skyView[..., fragments[idx]] = skyViewChunk

	else:
		skyView = ftWorkers.dftWorker(0, correlationMatrix, obsFreq, conjWeight, skyView)[1]

	return skyView.transpose((1,0,2))

def fftImage(correlationMatrix, posX, posY, obsFreq, lVec, mVec, kernel = 'gaus/0.5', generalSettings = None):
	"""Summary
	
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
	generalSettings = generalSettings or None

	multiprocessing = False
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

	print('uvminmax', np.min(uv), np.max(uv))

	skyViewPixels = int(np.max([lVec.size, mVec.size]))
	pixelBuffer = int(0.1 * skyViewPixels)
	pixelBuffer += pixelBuffer % 2 # Ensure it's an even number
	print(pixelBuffer)
	uv *= skyViewPixels
	uv += pixelBuffer / 2 # 5% buffer

	print('uvminmax', np.min(uv), np.max(uv))

	skyViewPixels += pixelBuffer

	# Ensure we have an odd-sized grid for efficient fft operations
	if skyViewPixels % 2:
		skyViewPixels += 1


	uv = uv.reshape(2, -1)
	correlationMatrix = correlationMatrix.reshape(1, -1, frameCount) # (1, 96 * 96, nChan)
	skyView = np.zeros([skyViewPixels, skyViewPixels, frameCount])
	if multiprocessing:
		fragments = np.array_split(np.arange(frameCount))
		processCount = int(mp.cpu_count()-1)
		mpPool = mp.Pool(processes = processCount)
		callBacks = [mpPool.apply_async(ftWorkers.fftWorker, args = ([idx, correlationMatrix[..., fragmentChunk], uv, kernel, skyView[..., fragmentChunk]])) for idx, fragmentChunk in enumerate(fragments)]

		mpPool.close()
		mpPool.join()

		for asyncResult in callBacks:
			idx, skyViewChunk = asyncResult.get()
			skyView[..., fragments[idx]] = skyViewChunk

	else:
		skyView = fftWorker(0, correlationMatrix, uv, kernel, skyView)[1]

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