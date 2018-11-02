
from ..dataTools import usefulFunctions

def generatePlots(inputCorrelations, antPosArr, processingOptions, plotOptions, dateArr, rcuMode, subband, stationRotation, stationLocation = None, calibrationArr = None):
	"""Summary
	
	Args:
	    inputCorrelations (TYPE): Description
	    antPos (TYPE): Description
	    plotOptions (TYPE): Description
	    dateArr (TYPE): Description
	    rcuMode (TYPE): Description
	    subband (TYPE): Description
	    multiprocessing (bool, optional): Description
	    stationRotation (TYPE, optional): Description
	    plotX (bool, optional): Description
	    plotY (bool, optional): Description
	    mask (bool, optional): Description
	    lVec (None, optional): Description
	    mVec (None, optional): Description
	    calibrationX (None, optional): Description
	    calibrationY (None, optional): Description
	    baselineLimits (None, optional): Description
	    stationLocation (None, optional): Description
	"""
	crossCorr, posList, rawAnts = __initialiseImaging(inputCorrelations, antPosArr, calibrationArr, subband, options['baselineLimits'])
	xxCorr, yyCorr, xyCorr, yxCorr = crossCorr
	posX, posY, posZ = posList

	frequency = usefulFunctions.calcFreq(rcuMode, subband) * 1e6

	labelOptions = [dateArr, rcuMode, subband, frequency]

	lVec, mVec = processingOptions['pixelCount']
	fov = processingOptions['fieldOfView']

	vecRange = np.sin([-1. * fov / 2., fov / 2.])

	lVec = np.arange(vecRange[0], vecRange[1] + 0.01, (2. / lVec) or 0.008 ) # Default to 250 pixels for a 180degree field of view
	mVec = np.arange(vecRange[0], vecRange[1] + 0.01, (2. / mVec) or 0.008 )

	pixels = np.max([lVec.size, mVec.size])

	mask = np.zeros([lVec.size, mVec.size]).astype(bool)
	if mask:
		dist = np.sqrt(np.square(np.meshgrid(lVec)) + np.square(np.meshgrid(mVec)).T)
		mask[dist >= vecRange[1]] = True

	figNum = 0

def __initialiseImaging(inputCorrelations, antPosArr, calibrationArr, subband, baselineLimits):

	antPos, rawAnts = antPosArr
	rawAnts = rawAnts[:, 0, :][:, np.newaxis]
	posX = antPos[..., 0]
	posY = antPos[..., 1]
	posZ = antPos[..., 2]

	inputCorrelationsX = inputCorrelations[..., 0]
	inputCorrelationsY = inputCorrelations[..., 1]

	if calibrationArr is not None:
		print('Calibrating data for subband {0}'.format(subband))

		calSubbandX = calibrationArr[:, subband, 0]
		calSubbandY = calibrationArr[:, subband, 1]

		calMatrixXArr = np.outer(calSubbandX, np.conj(calSubbandX).T)[..., np.newaxis]
		inputCorrelationsX = np.multiply(np.conj(calMatrixXArr), inputCorrelationsX)

		calMatrixYArr = np.outer(calSubbandY, np.conj(calSubbandY).T)[..., np.newaxis]
		inputCorrelationsY = np.multiply(np.conj(calMatrixYArr), inputCorrelationsY)

	if baselineLimits:
		print('Test baseline limits.')
		baselines = posX - posY.T
		if baselineLimits[0]:
			baselineLimits[0] = np.min(baselines) * 1.0001
			print('By not using autocorrelations the sky is a lot quieter: we are disabling log plotting as a result.')
			plotOptions[0] = False
		elif not baselineLimits[1]:
			baselineLimits[1] = np.max(baselines) * 2.
		minBaseline, maxBaseline = baselineLimits

		baselines = np.logical_or(np.abs(baselines) > maxBaseline, np.abs(baselines) < minBaseline)

		inputCorrelationsX[baselines] = 0.
		inputCorrelationsY[baselines] = 0.
	# Reference
	stackedArr = np.zeros(np.array(inputCorrelationsX.shape[:2]) * 2, dtype = 'complex')
	stackedArr[::2, ::2] = inputCorrelationsX[..., 0]
	stackedArr[1::2, 1::2] = inputCorrelationsY[..., 0]

	xxCorr = inputCorrelationsX
	yyCorr = inputCorrelationsY
	xyCorr = stackedArr[1::2, ::2]
	yxCorr = stackedArr[::2, 1::2]

	return [xxCorr, yyCorr, xyCorr, yxCorr], [posX, posY, posZ], rawAnts

def __processCorrPlot():
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

@cache
def cachedSkyCoords(name):
	return astropy.coordinates.SkyCoord.from_name(name)