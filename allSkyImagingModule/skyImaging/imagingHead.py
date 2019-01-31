"""Main Imaging Functions
"""
import numpy as np

from . import ftImaging
from . import swhtImaging

from dataTools.usefulFunctions import calcFreq

def generatePlots(inputCorrelations, antPosArr, options, metaDataArr, rcuMode, subband, stationRotation, stationLocation = None, calibrationArr = None):
	"""Head function for generating output data
	
	Args:
	    inputCorrelations (np.ndarray): Input correlations (nAnts x nAnts x nSamples array of correlations)
	    antPosArr (list): Raw/StationCoord antenna location
	    options (dict): Standard options dictionary
	    metaDataArr (list): List of str(dict) from observation group attributes
	    rcuMode (int): RCU Mode of observation
	    subband (int): Subband observed
	    stationRotation (float, optional): Station rotation, degrees east of north
	    stationLocation (list, optional): [lon, lat, alt] of station
	    calibrationArr (np.ndarray, optional): (nAnts, 512) calibration array
	
	Returns:
	    np.ndarray: Output sky data
	
	Raises:
	    RuntimeError: If an unknown processing method is provided
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

	labelOptions = [metaDataArr, rcuMode, subband, frequency]

	# Check what output types are requested: Only process correlations if needed
	# This is done by generating a dictionary of 'corrType': boolProcOrNot pairs.
	processDict = {'XX': False, 'YY': False, 'XY': False, 'YX': False}

	if isinstance(processingOptions['correlationTypes'], str):
		processingOptions['correlationTypes'] = [processingOptions['correlationTypes']]
	for value in processingOptions['correlationTypes']:
		if len(value) == 2:
			processDict[value] = True

		elif value in ['I', 'Q']:
			processDict['XX'] = True
			processDict['YY'] = True

		elif value in ['U', 'V']:
			processDict['XY'] = True
			processDict['YX'] = True

	# RFI mode is an interactive plot allowing you to see the sky location / intensity by clicking on the plot
	rfiFlag = options['rfiMode']
	if rfiFlag:
		print('RFI Mode enabled, we will be forcing the FT method of choice (or FFT) and providing a GUI.')

	# Pass our variables onto the processing type class.
	procMethod = processingOptions['method']
	if 'ft' in procMethod or rfiFlag:
		imgData = ftImaging.ftProcessCorrelations(procMethod, rfiFlag, corrDict, processDict, options, posArr, frequency, stationRotation, stationLocation, labelOptions, metaDataArr)
	elif 'swht' in procMethod:
		imgData = swhtImaging.swhtProcessCorrelations(corrDict, options, processDict, rawAnts, stationLocation, metaDataArr, frequency, labelOptions)
	else:
		raise RuntimeError('Unknown processing method "{0}".'.format(procMethod))

	return imgData

def __initialiseImaging(inputCorrelations, antPosArr, calibrationArr, subband, baselineLimits):
	"""Extract the starting data from input variables
	
	Args:
	    inputCorrelations (np.ndarray): Input correlations (nAnts x nAnts x nSamples array of correlations)
	    antPosArr (np.ndarray): Raw/StationCoord atenna positions
	    calibrationArr (np.ndarray, optional): (nAnts, 512) calibration array
	    subband (int): Observing subband
	    baselineLimits (list): Limit processed correlations based on provided baseline limits
	
	Returns:
	    list: Processed results
	"""
	antPos, rawAnts = antPosArr
	rawAnts = rawAnts[:, 0, :][:, np.newaxis]
	posX = antPos[..., 0]
	posY = antPos[..., 1]
	posZ = antPos[..., 2]

	baselines = np.sqrt(np.square(posX - posX.T) + np.square(posY - posY.T))
	print(np.max(baselines), np.min(baselines[baselines > 0.]))
	np.save('./bl.npy', baselines)
	
	if baselineLimits[0] or baselineLimits[1]:
		baselines = np.sqrt(np.square(posX - posX.T) + np.square(posY - posY.T))
		
		if baselineLimits[0] == 'noAuto':
			baselineLimits[0] = np.min(baselines[baselines > 0.]) * 1.0001

		minBaseline, maxBaseline = baselineLimits

		baselines = np.logical_or(np.abs(baselines) > maxBaseline, np.abs(baselines) < minBaseline)

		inputCorrelations[::2, ::2, ...][baselines] = 0.
		inputCorrelations[1::2, 1::2, ...][baselines] = 0.


	if calibrationArr is not None:
		print('Calibrating data for subband {0}'.format(subband))

		calSubbandX = calibrationArr[:, subband, 0]
		calSubbandY = calibrationArr[:, subband, 1]

		calMatrixXArr = np.outer(calSubbandX, np.conj(calSubbandX).T)[..., np.newaxis]
		inputCorrelations[::2, ::2, ...] = np.multiply(np.conj(calMatrixXArr), inputCorrelations[::2, ::2, ...])

		calMatrixYArr = np.outer(calSubbandY, np.conj(calSubbandY).T)[..., np.newaxis]
		inputCorrelations[1::2, 1::2, ...] = np.multiply(np.conj(calMatrixYArr), inputCorrelations[1::2, 1::2, ...])


	xxCorr = inputCorrelations[::2, ::2, ...]
	yyCorr = inputCorrelations[1::2, 1::2, ...]
	xyCorr = inputCorrelations[::2, 1::2, ...]
	yxCorr = inputCorrelations[1::2, ::2, ...]

	return [xxCorr, yyCorr, xyCorr, yxCorr], [posX, posY, posZ], rawAnts
