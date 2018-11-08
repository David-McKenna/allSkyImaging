"""Summary
"""
import h5py
import numpy as np
import multiprocessing as mp
import os
import time

from dataTools import genericImportTools as importTools
from dataTools import antennaHandler
from dataTools import xstImporter
from dataTools import accImporter
from dataTools import defaultDict
from skyImaging import imagingHead

reload(importTools)
reload(antennaHandler)
reload(xstImporter)
reload(defaultDict)
reload(imagingHead)

def image(fileLocation, obsType = 'XST', options = None):
	"""Summary
	
	Args:
		fileLocation (TYPE): Description
		obsType (str, optional): Description
		options (dict, optional): Description
	
	Returns:
		TYPE: Description
	
	Raises:
		RuntimeError: Description
	"""
	options, rcuMode, __, rotation, activationPattern, calibrationLocation, stationPackage = basicInitialisation(options)
	posXPol, stationLocation, __, antLocs = stationPackage

	posX = posXPol[:, 0, np.newaxis]
	posY = posXPol[:, 1, np.newaxis]
	posZ = posXPol[:, 2, np.newaxis]
	antPos = np.dstack([posX, posY, posZ])

	if obsType.lower() == 'xst':
		outputFile, groupPrefix = xstImporter.importXST(fileLocation, outputFile = options['fileLocations']['outputH5Location'], groupNamePrefix = options['h5GroupName'], rcuMode = rcuMode, calibrationFile = calibrationLocation, activationPattern = activationPattern)
	elif obsType.lower() == 'acc':
		'''
		# In the case we have been provided a folder of ACC files, run recursively to process them. 
		if os.path.isdir(fileName):
			fileList, __, __ = processInputLocation(fileName, 'acc')
			for fileName in fileList:
				main(fileName, 'ACC', breakThings, rcuMode, subbandArr, deltasLoc, fieldLoc, plotOptions, activation, calLoc, outputH5Loc, baselineLimits)

			return

		outputFile, groupPrefix, rcuMode = importXST.importACC(fileLocation, rcuMode, calibrationFile = calLoc, outputFile = outputH5Loc)
		'''
	else:
		raise RuntimeError('Unknown observation file type.')

	if not options['plottingOptions']['outputFolder']:
		options['plottingOptions']['outputFolder']= '/'.join(fileLocation.split('/')[:-1]) + '/allSkyOutput/'


	with h5py.File(outputFile, 'r+') as corrRef:
		if obsType != 'acc':
			subbandArr = [subbandInt[0][2:] for subbandInt in corrRef[groupPrefix].items() if 'sb' in subbandInt[0]]
			#totalFrames = np.sum([corrRef[groupPrefix][subbandName]['correlationArray'].shape[-1] for subbandName in corrRef[groupPrefix].items() if 'sb' in subbandName[0]])
			
			# Set a (two originally) threshold(s) for higher level array splitting:
			# 	#1: If we have more than 5 subbands to image (Unoptomised for parsing single subbands at a time)
			# 	#2: If we have more than 127 frames to image (to prevent OutOfMemory conditions) # Removed for now
			splitThreshold = (len(subbandArr) > 5) #  or (totalFrames > 127)
		else:
			'''
			subbandArr = np.arange(512)
			'''

		if 'calibrationArray' in corrRef[groupPrefix]:
			print('Extracting Calibrations')
			corrGroup = corrRef[groupPrefix]['calibrationArray']
			calibrationArr = corrGroup[...]
		else:
			calibrationArr = None


		if options['multiprocessing'] and splitThreshold:
			options['multiprocessing'] = False

			allSkyDatasetNames = np.zeros([len(subbandArr)], dtype = str)
			allSkyDataArr = np.zeros([len(subbandArr)], dtype = object)

			processCount = int(mp.cpu_count() * options['multiprocessingCoreFrac'])
			mpPool = mp.Pool(processes = processCount)
			callBacks = [mpPool.apply_async(processData, args = ([idx, corrRef, calibrationArr, antPos, antLocs, options, rcuMode, stationLocation, rotation, groupPrefix, [subband]])) for idx, subband in enumerate(subbandArr)]
	
			mpPool.close()
			mpPool.join()
	
			for asyncResult in callBacks:
				idx, allSkyDatasetName, allSkyDataArr = asyncResult.get()
				allSkyDatasetNames[idx] = allSkyDatasetName[0]
				allSkyDataArr[idx] = allSkyDataArr[0]
			
			allSkyDatasetNames = list(allSkyDatasetNames)
			allSkyDataArr = list(allSkyDataArr)

		else:
			__, allSkyDatasetNames, allSkyDataArr = processData(1, corrRef, calibrationArr, antPos, antLocs, options, rcuMode, stationLocation, rotation, groupPrefix, subbandArr)

		for allSkyDatasetName, allSkyData in zip(allSkyDatasetNames, allSkyDataArr):
			for key, allSkyDataArr in allSkyData.items():
				# h5py doesn't support complex values?!
				if isinstance(allSkyDataArr, np.ma.MaskedArray):
					outputDataGroup = corrRef.require_dataset('{0}-{1}-{2}'.format(allSkyDatasetName, key, 'mask') , allSkyDataArr.shape, compression = 'lzf', dtype = np.float64)
					outputDataGroup[...] = allSkyDataArr.mask
					outputDataGroup.attrs.create('options', str(options))
					outputDataGroup.attrs.create('fillValue', 100.)

					allSkyDataArr = allSkyDataArr.filled(100.)


				outputDataGroup = corrRef.require_dataset('{0}-{1}-{2}'.format(allSkyDatasetName, key, 'real') , allSkyDataArr.shape, compression = 'lzf', dtype = np.float64)
				outputDataGroup[...] = allSkyDataArr.real
				outputDataGroup.attrs.create('options', str(options))

				if isinstance(allSkyDataArr, complex):
					outputDataGroup = corrRef.require_dataset('{0}-{1}-{2}'.format(allSkyDatasetName, key, 'imag') , allSkyDataArr.shape, compression = 'lzf', dtype = np.float64)
					outputDataGroup[...] = allSkyDataArr.imag
					outputDataGroup.attrs.create('options', str(options))

	return outputFile, allSkyDatasetNames

def monitor(fileLocation, obsType = 'XST', options = None, processedOutput = '/processed/', fileThreshold = 1, checkEvery = 10.):
	"""Summary
	
	Returns:
	    TYPE: Description
	"""
	totalNewFiles = 0
	while True:
		currentFiles = os.listdir(fileLocation)
		obsFiles = [fileLocation + fileName for fileName in currentFiles if obsType.lower() in fileName]

		if not os.path.isdir(fileLocation + processedOutput):
			os.mkdir(fileLocation + processedOutput)

		if len(obsFiles) > fileThreshold:
			print("Succifient new files found! Begining to process {0} files.".format(len(obsFiles)))
			image(fileLocation, obsType, options)

			totalNewFiles += len(obsFiles)
			print("Processing complete, moving files to output folder at {0}".format(fileLocation + processedOutput))
			for fileName in obsFiles:
				os.rename(fileName, fileName.replace(fileLocation, fileLocation + processedOutput))
		else:
			print('Insufficient new files found ({0}), sleeping for {1} seconds...'.format(len(obsFiles), checkEvery))
			time.sleep(checkEvery)

		print("We have currently processed {0} new files.".format(totalNewFiles))


def processData(idx, corrRef, calibrationArr, antPos, antLocs, options, rcuMode, stationLocation, rotation, groupPrefix, subbandArr):
	"""Summary
	
	Args:
		idx (TYPE): Description
		corrRef (TYPE): Description
		calibrationArr (TYPE): Description
		antPos (TYPE): Description
		antLocs (TYPE): Description
		options (TYPE): Description
		rcuMode (TYPE): Description
		stationLocation (TYPE): Description
		rotation (TYPE): Description
	
	Returns:
		TYPE: Description
	"""
	allSkyDatasetNames = []
	allSkyDataArr = []
	for subbandVal in subbandArr:
		if 'accObs' not in corrRef[groupPrefix]:
			corrArr = corrRef[groupPrefix]['sb' + str(subbandVal)]['correlationArray']
		else:
			corrArr = corrRef[groupPrefix]['accObs/correlationArray']

		dateArr = np.vstack(corrArr.attrs.values()).astype(str)[:, -1]

		allSkyData = imagingHead.generatePlots(corrArr, [antPos, antLocs], options, dateArr, rcuMode, int(subbandVal), stationRotation = rotation, stationLocation = stationLocation, calibrationArr = calibrationArr)

		# h5py only offers multithreaded writes through the use of MPI4py
		allSkyDatasetName = '{0}sb{1}/imageData/'.format(groupPrefix, subbandVal)
		allSkyDatasetNames.append(allSkyDatasetName)
		allSkyDataArr.append(allSkyData)


	return idx, allSkyDatasetNames, allSkyDataArr

def basicInitialisation(options):
	"""Summary
	
	Args:
	    options (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	if not options:
		options = defaultDict.default()
	elif 'fullStructure' not in options:
		options = defaultDict.patchDefault(options)

	rcuMode = options['rcuMode']
	stationName = options['stationID']
	activationPattern = options['activationPattern']
	options = defaultDict.updateLocation(options, stationName, rcuMode)
	importTools.checkRequiredFiles(options['fileLocations'], stationName, rcuMode)

	rotation = importTools.processTxt(options['fileLocations']['lbaRotation'], stationName)

	if options['imagingOptions']['calibrateData']:
		calibrationLocation = options['fileLocations']['calibrationLocation']
	else:
		calibrationLocation = None

	if rcuMode > 4:
		activationPattern = options['activationPattern']
	else:
		activationPattern = None

	stationPackage = antennaHandler.configureStation(options['fileLocations'], activationPattern, rcuMode, True)

	return options, rcuMode, stationName, rotation, activationPattern, calibrationLocation, stationPackage
