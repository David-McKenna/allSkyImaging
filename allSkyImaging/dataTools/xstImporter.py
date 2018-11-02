"""Functions related to importing XST data from binary .dat files.
"""
import numpy as np
import h5py
import ast
import datetime

def importXST(fileName, outputFile, groupNamePrefix, rcuMode = None, calibrationFile = None, activationPattern = None):
	"""Summary
	
	Args:
	    fileName (TYPE): Description
	    outputFile (None, optional): Description
	    groupNamePrefix (None, optional): Description
	    rcuMode (None, optional): Description
	    calibrationFile (None, optional): Description
	
	Returns:
	    TYPE: Description	
	"""

	# Gather the XST files in a given location
	fileList, fileName, folderPath = processInputLocation(fileName, dataType = 'XST')
	fileList.sort(key = lambda f: int(filter(str.isdigit, f.split('sb')[-1]))) # Reorder by subband afterwards. Shouldn't be needed anymore, but it's nice to keep for peace of mind.

	# Check if we have logs provided to extract metadata, otherwise set some sane defaults.
	try:
		# Perform a quick test of RCU mode to see if we have a log file. Otherwise this will raise an IOError we can catch.
		testRef = open(fileList[0] + '.log', 'r')
		rcuModeRead = int(testRef[0].split(' ')[-1])
		testRef.close()

		if not rcuMode:
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

	if groupNamePrefix == 'allSkyObservation':
		groupNamePrefix = 'allSky-mode{0}'.format(rcuMode)

		if calibrationFile:
			groupNamePrefix += '-calibrated'

	# Create an output file containing the observations in a compressed h5 dataset
	with h5py.File(outputFile, 'a') as outputRef:

		# Initialise the file, observation group
		datasetComplexHead = {'processTime': str(datetime.datetime.utcnow())}
		groupRef = outputRef.require_group(groupNamePrefix)

		dataArr = []
		for fileName in fileList:
			with open(fileName, 'rb') as dataRef:
				# Read the dataset fromt he binary file
				datasetComplex = np.fromfile(dataRef, dtype = np.complex128)
				reshapeSize = datasetComplex.size / (192 ** 2)

				# Check if the file is incomplete (missing bytes or corrupted, all files should be multiples of a 192, 192 array)
				if not datasetComplex.size % (192 ** 2):
					print('INCOMPLETE FILE SKIPPED: {0}'.format(fileName))
					nullFiles += 1
					continue

				datasetComplex = datasetComplex.reshape(192, 192, reshapeSize, order = 'F')

				# Extract basic information from the filename
				fileNameMod = fileName.split('/')[-1]
				fileNameExtract = fileNameMod.split('_')
				dateTime = str(datetime.datetime.strptime(''.join(fileNameExtract[0:2]), '%Y%m%d%H%M%S'))
				subbandStr = fileNameExtract[2]


				# If we have the log files, parse some data from them
				if logFiles:
					logName = fileName + '.log'
					with open(logName, 'r') as logRef:
						logData = [ast.literal_eval(line.split(' ')[-1].strip('\n')) for line in logRef[:6]]

						assert(logData[0] == rcuMode)
						assert(logData[1] == int(subband[2:]))
						logData[2] = int(logData[2].strip('s'))
						logData[3] = int(logData[3].strip('s')) / logData[2]

						# If an assert has been raised here, you have an observation taken while the activation pattern was not fully applied.
						assert((rcuMode < 5) or not False in ['253' not in line for line in logRef[6:]])

				# Otherwise, update the metadata with information we can gleam from the filename
				else:
					logData = metadata + [None, None]
					logData[1] = int(subbandStr[2:])
					logData[3] = reshapeSize

				dataArr.append(np.array([subbandStr, dateTime, datasetComplex, reshapeSize, logData], dtype = object))


		# Group our data by subband
		dataArr= np.vstack(dataArr)
		subbandArr = np.unique(dataArr[:, 0])
		dateTimeArr = [dataArr[:, 1][dataArr[:, 0] == subbandVal] for subbandVal in subbandArr]
		datasetComplexArr = [dataArr[:, 2][dataArr[:, 0] == subbandVal] for subbandVal in subbandArr]
		reshapeSizeArr = [list(dataArr[:, 3][dataArr[:, 0] == subbandVal]) for subbandVal in subbandArr]
		logDataArr = [list(dataArr[:, 4][dataArr[:, 0] == subbandVal]) for subbandVal in subbandArr]
				
		for idx, subband in enumerate(subbandArr):
			# Extract the XX and YY correlations and stack them before saving them to disk
			datasetComplex = np.dstack(datasetComplexArr[idx])
			datasetComplexX = datasetComplex[::2, ::2, ...]
			datasetComplexY = datasetComplex[1::2, 1::2, ...]
			datasetComplex = np.stack([datasetComplexX, datasetComplexY], axis = -1)

			# h5py doesn't want to place nicely, fill the array after creating it.
			corrDataset = groupRef.require_dataset("{0}/correlationArray".format(subband), datasetComplex.shape, dtype = np.complex128, compression = "lzf")
			corrDataset[...] = datasetComplex

			# Save the activation pattern used for the observation as an attribute
			corrDataset.attrs.create('activationPattern', str(activationPattern))

			# For each subband, iterate over every saved frame and fill in the group attributes.
			timeStep = 0
			for dtIdx, logData in enumerate(logDataArr[idx]):
					currDelta = reshapeSizeArr[idx][dtIdx]
					mode, subband, intTime, intCount = logData[idx][dtIdx]

					timeDelta = datetime.timedelta(seconds = intTime / 2.)
					intTime = datetime.timedelta(seconds = intTime)

					# Each integration gets it's own frame, so it needs a unique attribute to corretly timestamp it.
					# We log the time of the integration as the central time to better account for the sky during long integrations
					for i in range(intCount):
						corrDataset.attrs.create(str(timeStep), [mode, subband, intTime, dateTimeArr[idx][dtIdx] + timeDelta * (i + 1)])
						timeStep += 1

		# If provided a calibration file, include it in the dataset. This call is skipped if a calibration is already included
		#	as we do not expect the calibration have majour changes over time.
		if (calibrationFile is not None) and "calibrationArray" not in groupRef:
			includeCalibration(calibrationFile, groupRef)

	return outputFile, groupNamePrefix