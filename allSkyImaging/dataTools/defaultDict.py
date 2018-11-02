files = [[defaultField.replace('IE613', stationName), 'https://raw.githubusercontent.com/griffinfoster/SWHT/hbaDeltas/SWHT/data/LOFAR/StaticMetaData/{0}-AntennaField.conf'.format(stationName)], [defaultDeltas.replace('IE613', stationName), 'https://raw.githubusercontent.com/griffinfoster/SWHT/hbaDeltas/SWHT/data/LOFAR/StaticMetaData/iHBADeltas/{0}-iHBADeltas.conf'.format(stationName)], [lbaRotation, 'https://raw.githubusercontent.com/cosmicpudding/lofarimaging/master/stationrotations.txt']]


def default():

	fileLocations = {
		'antennaField': './{0}-AntennaField.conf',
		'hbaDeltas': './{0}-iHBADeltas.conf',
		'lbaRotation': './stationrotations.txt',
		'outputH5Location': './{0}.h5',
		'calibrationLocation': './{0}-mode{1}dat', 

		'remoteURL': {
			'antennaField': 'https://raw.githubusercontent.com/griffinfoster/SWHT/hbaDeltas/SWHT/data/LOFAR/StaticMetaData/',
			'hbaDeltas': 'https://raw.githubusercontent.com/griffinfoster/SWHT/hbaDeltas/SWHT/data/LOFAR/StaticMetaData/',
			'lbaRotation': 'https://raw.githubusercontent.com/cosmicpudding/lofarimaging/master/',
			'calibrationLocation': 'https://svn.astron.nl/Station/trunk/CalTables/{0}/'
		},
	}

	imagingOptions = {
		'calibrateData': True,
		'baselineLimits': [None, None],
		'pixelCount': [255, 255],
		'fieldOfView': np.pi, # 180 degrees

	}

	plottingOptions = {
		'outputFolder': None,

	}

	defaultDictionary= {
		'fullStructure': True,
		'stationID': 'IE613',
		'rcuMode': 3,
		'activationPattern': 'Generic_Int_201512',
		'h5GroupName': 'allSkyObservation',

		'multiprocessing': False,
		'multiprocessingCoreFrac': 0.66,

		'fileLocations': fileLocations,
		'imagingOptions': imagingOptions,
		'plottingOptions': plottingOptions
	}

def patchDefault(newOptions):
	
	defaultDictionary = default()
	for option, newValue in newOptions.items():
		if ':' in option:
			subDict, subOpt = option.split(':')
			defaultDictionary[subDict][subOpt] = newValue
		else:
			defaultDictionary[option] = newValue

	return defaultDictionary

def updateLocation(currDict, stationID, rcuMode):

	fileLocations = currDict['fileLocations']

	fileLocations['antennaField'] = fileLocations['antennaField'].format(stationID)
	fileLocations['hbaDeltas'] = fileLocations['hbaDeltas'].format(stationID)
	fileLocations['calibrationLocation'] = fileLocations['calibrationLocation'].format(stationID, rcuMode)
	
	if '{' in fileLocations['outputH5Location']:
		fileLocations['outputH5Location'] = fileLocation['outputH5Location'].format(stationID)

	currDict['fileLocations'] = fileLocations

	return currDict