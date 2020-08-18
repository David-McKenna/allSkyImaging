"""Functions related to importing XST data from binary .dat files.
"""
import numpy as np
import h5py
import ast
import datetime

from .genericImportTools import processInputLocation, processRCUMode, includeCalibration

def importXST(fileName, outputFile, groupNamePrefix, rcuMode = None, calibrationFile = None, activationPattern = None):
	"""Import a/a folder of XST observations.
	
	Args:
	    fileName (str): Input file/folder name
	    outputFile (str): Output h5 name
	    groupNamePrefix (str): Output h5 group name
	    rcuMode (int, optional): RCU Observing Mode
	    calibrationFile (str, optional): Calibration file location
	    activationPattern (str, optional): HBA Activation pattern name
	
	Returns:
	    str, str: Output file location, output group name
	"""

	# Gather the XST files in a given location
	fileList, fileName, folderPath = processInputLocation(fileName, dataType = 'XST')
	fileList.sort(key = lambda f: int(''.join(l for l in f.split('sb')[-1] if l.isdigit()))) # Reorder by subband afterwards. Shouldn't be needed anymore, but it's nice to keep for peace of mind.

	# Check if we have logs provided to extract metadata, otherwise set some sane defaults.
	try:
		# Perform a quick test of RCU mode to see if we have a log file. Otherwise this will raise an IOError we can catch.
		testRef = open(fileList[0] + '.log', 'r')
		line = testRef.readline().strip('\n')
		rcuModeRead = int(line.split(' ')[-1])
		testRef.close()
		assert(rcuMode == rcuModeRead)

		logFiles = True

	except IOError:
		logFiles = False

		# Check if we were provided an rcuMode or try determine it from the input data folder or calibration file
		rcuMode = rcuMode or processRCUMode(folderPath, calibrationFile)

		# mode: [mode, subband, integration time, integrations]
		metadata = {'1': [1, 100, 5, 1], '2': [2, 100, 5, 1], '3': [3, 100, 5, 1], '4': [4, 100, 5, 1], '5': [5, 200, 10, 1], '6': [6, 200, 10, 1], '7': [7, 200, 10, 1]}
		metadata = metadata[str(rcuMode)]

		# Could now be a single print statement; used to take more asusmptions based on rcuMode that are now determined later on
		warnMessage = 'Unable to open log files, we will make assumptions for the observation\'s metadata (setting integration time and observation duration to {0}s)'.format(metadata[2])
		
		print(warnMessage)

	# Count the number of incomplete files
	nullFiles = 0

	# Create an output file containing the observations in a compressed h5 dataset
	with h5py.File(outputFile, 'a') as outputRef:

		# Initialise the file, observation group
		groupRef = outputRef.require_group(groupNamePrefix)

		dataArr = []
		for fileNameVar in fileList:
			with open(fileNameVar, 'rb') as dataRef:
				# Read the dataset fromt he binary file
				datasetComplex = np.fromfile(dataRef, dtype = np.complex128)
				reshapeSize = datasetComplex.size / (192 ** 2)

				# Check if the file is incomplete (missing bytes or corrupted, all files should be multiples of a 192, 192 array)
				if datasetComplex.size % (192 ** 2) > 0:
					print('INCOMPLETE FILE SKIPPED: {0}, SIZE ONLY {1}, MISSING {2} DATAPOINTS'.format(fileNameVar, datasetComplex.size, (192 ** 2) - datasetComplex.size % (192 ** 2)))
					nullFiles += 1
					continue
				reshapeSize = int(reshapeSize)
				datasetComplex = datasetComplex.reshape(192, 192, reshapeSize, order = 'F')

				if 'dropData' in outputFile:
					datasetComplex = datasetComplex[..., -1]
					reshapeSize = 1
				
				# Extract basic information from the filename
				fileNameMod = fileNameVar.split('/')[-1]
				fileNameExtract = fileNameMod.split('_')
				dateTimeObj = datetime.datetime.strptime(''.join(fileNameExtract[0:2]), '%Y%m%d%H%M%S')
				dateTime = str(datetime.datetime.strptime(''.join(fileNameExtract[0:2]), '%Y%m%d%H%M%S'))
				subbandStr = fileNameExtract[2]

				print('Processing {0} observations from {1} at subband {2}'.format(reshapeSize, dateTime, subbandStr))


				# If we have the log files, parse some data from them
				if logFiles:
					logName = fileNameVar + '.log'
					with open(logName, 'r') as logRef:
						lines = [line.strip('\n').split(' ')[-1] for line in logRef]

						# AST broke for limited number of cases; revert to manual filtering.
						logData = [[]] * 6

						logData[0] = int(lines[0])
						logData[1] = int(lines[1])
						logData[2] = int(lines[2].strip('s'))
						logData[3] = int(lines[3].strip('s')) / logData[2]
						logData[4] = dateTimeObj
						logData[5] = str(datetime.datetime.strptime(''.join(lines[4]), '%Y/%m/%d@%H:%M:%S'))
						# If an assert has been raised here, you have an observation taken while the activation pattern was not fully applied.
						assert(logData[0] == rcuMode)
						assert(logData[1] == int(subbandStr[2:]))
						if not (rcuMode < 5) and False in ['253' not in line for line in lines[6:]]:
							print('WARNGING: We have detected that the activation pattern was not properly applied. \nAs a result, we are skipping the observation at time {0}, subband {1}'.format(logData[5], subbandStr))
							continue

				# Otherwise, update the metadata with information we can gleam from the filename
				else:
					logData = metadata + [dateTimeObj, None]
					logData[1] = int(subbandStr[2:])
					logData[3] = reshapeSize
					logData[4] = dateTimeObj
				dataArr.append(np.array([subbandStr, dateTime, datasetComplex, reshapeSize, logData], dtype = object))

		if groupNamePrefix == 'allSkyObservation':
			groupNamePrefix = 'allSky-mode{0}-{1}'.format(rcuMode, str(dataArr[0][1]))

			if calibrationFile:
				groupNamePrefix += '-calibrated-'

		# Initialise the file, observation group
		groupRef = outputRef.require_group(groupNamePrefix)

		# Group our data by subband
		dataArr= np.vstack(dataArr)
		subbandArr = np.unique(dataArr[:, 0])
		#dateTimeArr = [dataArr[:, 1][dataArr[:, 0] == subbandVal] for subbandVal in subbandArr]
		datasetComplexArr = [dataArr[:, 2][dataArr[:, 0] == subbandVal] for subbandVal in subbandArr]
		logDataArr = [list(dataArr[:, 4][dataArr[:, 0] == subbandVal]) for subbandVal in subbandArr]
				
		for idx, subband in enumerate(subbandArr):
			datasetComplex = np.dstack(datasetComplexArr[idx])

			# h5py doesn't want to place nicely, fill the array after creating it.
			corrDataset = groupRef.require_dataset("{0}/correlationArray".format(subband), datasetComplex.shape, dtype = np.complex128, compression = "lzf")
			corrDataset[...] = datasetComplex

			# For each subband, iterate over every saved frame and fill in the group attributes.
			timeStep = 0
			print("Writing metadata to file. Note: this will be slow if you have previously processed this dataset, move or remove the old H5 file as needed.")
			for logData in logDataArr[idx]:
				print('Imported frame at time {0} in subband {1}'.format(logData[-2], subband))
				mode, subband, intTime, intCount, dateTimeObj, endTime = logData

				timeDelta = datetime.timedelta(seconds = intTime / 2.)
				intTime = datetime.timedelta(seconds = intTime)

				# Each integration gets it's own frame, so it needs a unique attribute to corretly timestamp it.
				# We log the time of the integration as the central time to better account for the sky during long integrations
				for i in range(int(intCount)):
					corrDataset.attrs.create('{0:04d}'.format(timeStep), str({'mode': mode, 'subband': subband, 'integrationTime': intTime, 'activationPattern': str(activationPattern), 'integrationMidpoint': str(dateTimeObj + timeDelta * (i + 1)), 'integrationEnd': endTime})) # Can be decoded with ast
					timeStep += 1

		# If provided a calibration file, include it in the dataset. This call is skipped if a calibration is already included
		#	as we do not expect the calibration have majour changes over time.
		if (calibrationFile is not None) and "calibrationArray" not in groupRef:
			includeCalibration(calibrationFile, groupRef)

	return outputFile, groupNamePrefix
