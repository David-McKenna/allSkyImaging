"""Summary
"""

def importACC(fileName, outputFile, groupNamePrefix, rcuMode = None, calibrationFile = None, activationPattern = None): 
	"""Summary
	
	Args:
	    fileName (TYPE): Description
	    rcuMode (None, optional): Description
	    calibrationFile (None, optional): Description
	    outputFile (None, optional): Description
	    groupNamePrefix (None, optional): Description
	    integrationTime (None, optional): Description
	
	Returns:
	    TYPE: Description
	"""
	fileList, fileName, folderPath = processInputLocation(fileName, dataType = 'ACC')
	rcuMode = rcuMode or processRCUMode(rcuMode, folderPath)

	
			
	outputFile, groupNamePrefix = h5PrepNames(outputFile, groupNamePrefix)

	with h5py.File(outputFile, 'a') as outputRef:
		datasetComplexHead = {'processTime': str(datetime.datetime.utcnow())}
		groupRef = outputRef.require_group(groupNamePrefix)

		dataArr = []
		dateArr = []
		for fileName in fileList:
			with open(fileName, 'rb') as dataRef:
				datasetComplex = np.fromfile(dataRef, dtype = np.complex128)
				datasetComplex = datasetComplex.reshape(192, 192, 512)

				fileNameMod = fileName.split('/')[-1]
				fileNameExtract = fileNameMod.split('_')
				dateTime = datetime.datetime.strptime(''.join(fileNameExtract[0:2]), '%Y%m%d%H%M%S')


				dataArr.append(datasetComplex)
				dateArr.append(dateTime)
				
		for idx, observation in enumerate(dataArr):

			observationX = observation[::2, ::2, ...]

			observationY = observation[1::2, 1::2, ...]
			observation = np.stack([observationX, observationY], axis = -1)

			corrDataset = groupRef.require_dataset("correlationArray", observation.shape, dtype = np.complex128, compression = "lzf")
			corrDataset[...] = observation

			dateTimeVar = np.datetime64(dateArr[idx])
			allTimes = dateTimeVar + np.arange(-511, 1) * np.timedelta64(1, 's')

			corrDataset.attrs.create('rcuMode', rcuMode)
			corrDataset.attrs.create('subband', np.arange(512))
			corrDataset.attrs.create('intTime', datetime.timedelta(seconds = 1.))
			corrDataset.attrs.create('sampleTimes', allTimes)
			corrDataset.attrs.create('activationPattern', activationPattern)


			if calibrationFile is not None:
				includeCalibration(calibrationFile, groupRef)

	return outputFile, groupNamePrefix, rcuMode 
