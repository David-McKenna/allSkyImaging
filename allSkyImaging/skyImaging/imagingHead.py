"""Summary
"""
import numpy as np
import scipy.ndimage

from . import ftImaging
from . import swhtImaging

from dataTools.usefulFunctions import calcFreq

def generatePlots(inputCorrelations, antPosArr, options, dateArr, rcuMode, subband, stationRotation, stationLocation = None, calibrationArr = None):
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
	processingOptions = options['imagingOptions']

	crossCorr, posArr, rawAnts = __initialiseImaging(inputCorrelations, antPosArr, calibrationArr, subband, options['imagingOptions']['baselineLimits'])
	xxCorr, yyCorr, xyCorr, yxCorr = crossCorr

	corrDict = {
			'XX': xxCorr,
			'YY': yyCorr,
			'XY': xyCorr,
			'YX': yxCorr
	}

	frequency = calcFreq(rcuMode, subband) * 1e6

	labelOptions = [dateArr, rcuMode, subband, frequency]

	processDict = {'XX': False, 'YY': False, 'XY': False, 'YX': False}
	for value in processingOptions['correlationTypes']:
		if len(value) == 2:
			processDict[value] = True

		elif value in ['I', 'Q']:
			processDict['XX'] = True
			processDict['YY'] = True

		elif value in ['U', 'V']:
			processDict['XY'] = True
			processDict['YX'] = True


	rfiFlag = options['rfiMode']
	if rfiFlag:
		print('RFI Mode enabled, we will be forcing the FT method of choice (or FFT) and providing a GUI.')

	procMethod = processingOptions['method']
	if 'ft' in procMethod or rfiFlag:
		imgData = ftImaging.ftProcessCorrelations(procMethod, rfiFlag, corrDict, processDict, options, posArr, frequency, stationRotation, stationLocation, labelOptions, dateArr)
	elif 'swht' in procMethod:
		imgData = swhtImaging.swhtProcessCorrelations(corrDict, options, processDict, rawAnts, stationLocation, dateArr, frequency, labelOptions)
	else:
		raise RuntimeError('Unknown processing method "{0}".'.format(procMethod))

	return imgData

def __initialiseImaging(inputCorrelations, antPosArr, calibrationArr, subband, baselineLimits):
	"""Summary
	
	Args:
	    inputCorrelations (TYPE): Description
	    antPosArr (TYPE): Description
	    calibrationArr (TYPE): Description
	    subband (TYPE): Description
	    baselineLimits (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
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

	if baselineLimits[0] or baselineLimits[1]:
		baselines = posX - posY.T
		
		if baselineLimits[0] == 'noAuto':
			baselineLimits[0] = np.min(baselines[baselines > 0.]) * 1.0001

		minBaseline, maxBaseline = baselineLimits

		baselines = np.logical_or(np.abs(baselines) > maxBaseline, np.abs(baselines) < minBaseline)

		inputCorrelationsX[baselines] = 0.
		inputCorrelationsY[baselines] = 0.

	# Reference
	stackedArr = np.zeros(list(np.array(inputCorrelationsX.shape[:2]) * 2) +  [inputCorrelationsX.shape[2]], dtype = 'complex')
	stackedArr[::2, ::2] = inputCorrelationsX[...]
	stackedArr[1::2, 1::2] = inputCorrelationsY[...]

	xxCorr = inputCorrelationsX
	yyCorr = inputCorrelationsY
	xyCorr = stackedArr[1::2, ::2]
	yxCorr = stackedArr[::2, 1::2]

	return [xxCorr, yyCorr, xyCorr, yxCorr], [posX, posY, posZ], rawAnts
