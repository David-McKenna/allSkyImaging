"""Module head function: call the functions here to run anything you need the moudle to do.
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
	"""Image a given file/folder. Options can be given a patch dictionary.

	Patch dictionaries are dictionaries where you want to keep most of the default values, with some modifications. 
	See defaultDict.py for details. For a simple rundown: change options by referencing them, or subdictioanry values
		by using the 'subdict:option' syntax as the dictionary key.
	Eg:
	{
		'rcuMode': 5,
		'imagingOptions:method': 'swht'
	}
	
	Args:
		fileLocation (str): File/Folder location
		obsType (str, optional): Observation Type
		options (dict, optional): Stadard options dictionary or patch value dictionary
	
	Returns:
		str, str: Output file name, output h5 group name
	
	Raises:
		RuntimeError: Passed an unknwon type of observation.
	"""

	# Get station information
	options, rcuMode, __, rotation, activationPattern, calibrationLocation, stationPackage = basicInitialisation(options)
	posXPol, stationLocation, __, antLocs = stationPackage

	posX = posXPol[:, 0, np.newaxis]
	posY = posXPol[:, 1, np.newaxis]
	posZ = posXPol[:, 2, np.newaxis]
	antPos = np.dstack([posX, posY, posZ])

	# Process the input correlations
	if obsType.lower() == 'xst':
		outputFile, groupPrefix = xstImporter.importXST(fileLocation, outputFile = options['fileLocations']['outputH5Location'], groupNamePrefix = options['h5GroupName'], rcuMode = rcuMode, calibrationFile = calibrationLocation, activationPattern = activationPattern)
	elif obsType.lower() == 'acc':
		raise RuntimeError('Reimplementation needed')
	else:
		raise RuntimeError('Unknown observation file type.')

	if not options['plottingOptions']['outputFolder']:
		options['plottingOptions']['outputFolder']= '/'.join(fileLocation.split('/')[:-1]) + '/allSkyOutput/'


	# Read data from the processed correlations from the h5 file and start imaging.
	with h5py.File(outputFile, 'r+') as corrRef:
		if obsType != 'acc':
			subbandArr = [subbandInt[0][2:] for subbandInt in corrRef[groupPrefix].items() if 'sb' in subbandInt[0]]
			#totalFrames = np.sum([corrRef[groupPrefix][subbandName]['correlationArray'].shape[-1] for subbandName in corrRef[groupPrefix].items() if 'sb' in subbandName[0]])
			
			# Set a (two originally) threshold(s) for higher level array splitting:
			# 	#1: If we have more than 5 subbands to image (Unoptomised for parsing single subbands at a time)
			# 	#2: If we have more than 127 frames to image (to prevent OutOfMemory conditions) # Removed for now
			splitThreshold = (len(subbandArr) > 5) #  or (totalFrames > 127)
		else:
			raise RuntimeError('Reimplementation needed')


		# If we had a correlation array and want to use it, apply it.
		if ('calibrationArray' in corrRef[groupPrefix]) and options['imagingOptions']['calibrateData']:
			print('Extracting Calibrations')
			corrGroup = corrRef[groupPrefix]['calibrationArray']
			calibrationArr = corrGroup[...]
		else:
			calibrationArr = None

		# If we have more than 5 subbands in the given processing folder, we probably are dealing with a subband scan. It'll be faster to parallelise
		#	on a subband level than per-observation level as a result.
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

		method = options['imagingOptions']['method']
		updatedDatasetNames = []
		# Store the output data in the same h5 file.
		for allSkyDatasetName, allSkyData in zip(allSkyDatasetNames, allSkyDataArr):
			for key, allSkyDataArr in allSkyData.items():
				# Store masks as separate arrays
				allSkyName = [allSkyDatasetName + method, key, 'mask']
				if isinstance(allSkyDataArr, np.ma.MaskedArray):
					writeDataset(corrRef, allSkyDataArr.mask, allSkyName, [('options', str(options)), ('fillValue', 100.)], np.float64)
					allSkyDataArr = allSkyDataArr.filled(100.)
					updatedDatasetNames.append(allSkyName)

				allSkyName[-1] = 'real'
				writeDataset(corrRef, allSkyDataArr.real.astype(np.float64), allSkyName, [('options', str(options))], np.float64)
				updatedDatasetNames.append(allSkyName)
				
				# h5py will fail to write if we have a complex array but zero values in the complex component, so...
				if isinstance(allSkyDataArr.flatten()[0], np.complex128):
					allSkyName[-1] = 'imag'
					writeDataset(corrRef, allSkyDataArr.imag.astype(np.float64),allSkyName, [('options', str(options))], np.float64)
					updatedDatasetNames.append(allSkyName)

	return outputFile, updatedDatasetNames

def monitor(fileLocation, obsType = 'XST', options = None, processedOutput = '/processed/', fileThreshold = 1, checkEvery = 10.):
	"""Monitor a folder for new observation and process them when thresholds are met
	
	Args:
	    fileLocation (str): Folder location to monitor
	    obsType (str, optional): Type of observation
	    options (dict, optional): Standard options / patch values dictionary
	    processedOutput (str, optional): Subfolder to store processed files in
	    fileThreshold (int, optional): Number of new files needed before we start a batch processing
	    checkEvery (float, optional): How often to check for new files

	"""
	totalNewFiles = 0
	while True:
		currentFiles = os.listdir(fileLocation)
		obsFiles = [fileLocation + fileName for fileName in currentFiles if obsType.lower() in fileName]

		if not os.path.isdir(fileLocation + processedOutput):
			os.makedirs(fileLocation + processedOutput)

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


def writeDataset(groupRef, dataset, nameList, attributes, dtypeVar):
	"""Handler for output dataset writes
	
	Args:
	    groupRef (h5Group): Group Reference to add dataset to
	    dataset (np.ndarray): Numpy array
	    nameList (list): Strings to format to give group name
	    attributes (list(tuple)): List of tuples of attributes to add to the dataset
	    dtypeVar (): Data type to save the dataset as (np.float64, np.complex128,...)
	"""
	print(nameList, type(dataset), type(dataset.flatten()[0]))
	outputDataGroup = groupRef.require_dataset('-'.join(nameList), dataset.shape, compression = 'lzf', dtype = dtypeVar)
	outputDataGroup[...] = dataset

	for key, val in attributes:
		outputDataGroup.attrs.create(key, val)




def processData(idx, corrRef, calibrationArr, antPos, antLocs, options, rcuMode, stationLocation, rotation, groupPrefix, subbandArr):
	"""Handle the image generation to give sane outputs for multiprocessing
	
	Args:
	    idx (int): Multiprocessing ID
	    corrRef (h5Group): Correlation h5 group reference
	    calibrationArr (np.array): Imported calibration data array
	    antPos (np.array): Raw antenna locations (GCRS)
	    antLocs (np.array): StationCoord Antenna location
	    options (dict): Options Dictionary
	    rcuMode (int): RCU mode for observation
	    stationLocation (list): [lon, lat, alt] for station
	    rotation (float): Station Rotation (east of north)
	    groupPrefix (list): Group name prefix
	    subbandArr (list): List of subbands to process
	
	Returns:
	    int, list, list: Multiprocessing id, list of dataset name, list of datasets
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
	"""Get station information from the passed dictionary.
	
	Args:
	    options (dict): Default/patched dictionary
	
	Returns:
	    list: Station / imaging useful parameters
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
