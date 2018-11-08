"""Summary
"""
import numpy as np

def default():
	"""Generate the default dictionary.
	"""
	fileLocations = {
		'antennaField': './{0}-AntennaField.conf',
		'hbaDeltas': './{0}-iHBADeltas.conf',
		'lbaRotation': './stationrotations.txt',
		'outputH5Location': './{0}-mode{1}.h5',
		'calibrationLocation': './{0}-mode{1}.dat', 

		'remoteURL': {
			'antennaField': 'https://raw.githubusercontent.com/griffinfoster/SWHT/master/SWHT/data/LOFAR/StaticMetaData/',
			'hbaDeltas': 'https://raw.githubusercontent.com/griffinfoster/SWHT/master/SWHT/data/LOFAR/StaticMetaData/',
			'lbaRotation': 'https://raw.githubusercontent.com/cosmicpudding/lofarimaging/master/',
			'calibrationLocation': 'https://svn.astron.nl/Station/trunk/CalTables/{0}/'
		}
	}

	imagingOptions = {
		'method': 'dft',
		'calibrateData': True,
		'baselineLimits': [None, None], # Defaults to no culling. To remove autocorrelations, set element 0 to 'noAuto'
		'pixelCount': [255, 255],
		'fieldOfView': np.pi, # 180 degrees
		'maskOutput': True,

		'correlationTypes': ['I', 'YY'],

		'subtractBackground': True,

		'lMax': 32

	}

	plottingOptions = {
		'outputFolder': None,
		'plotImages': True,
		'displayImages': False,
		'generateVideo': True,
		'videoFramerate': 8,

		'figureShape': [18, 14],

		'interstellarSources': ['Polaris', 'Cas A', 'Cyg A', 'Sgr A', 'Tau A', 'Vir A', 'Cen A', 'Vela'],
		'solarSystemSources': ['Sun', 'Jupiter', 'Moon', 'Uranus', 'Neptune'],

		'colorBarLimits': 16, #n time steps, maxmin (defaults to 16, 2 seconds of video playback)
		'maxPercentile': 99,
		'minPercentile': 5,

		'logPlot': True,
		'skyObjColor': 'black',
		'gridThickness': 0.5,
		'backgroundColor': 'black',
		'foregroundColor': 'white',
		'radialLabelAngle': 0,
		'fontSizeFactor': None,
		'colorBar': True,
		'graticule': True

	}

	defaultDictionary = {
		'fullStructure': True,
		'stationID': 'IE613',
		'rcuMode': 3,
		'activationPattern': 'Generic_Int_201512',
		'h5GroupName': 'allSkyObservation',

		'multiprocessing': False,
		'multiprocessingCoreFrac': 0.66,

		'rfiMode': False,

		'fileLocations': fileLocations,
		'imagingOptions': imagingOptions,
		'plottingOptions': plottingOptions
	}

	return defaultDictionary

def patchDefault(newOptions):
	"""Patch the default dictionary with the input values.

	Input values should be in a dictionary in the format 'option': newValue, or for a sub-dictionary 'dictName:subOption': newValue.
	
	Args:
	    newOptions (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	defaultDictionary = default()
	for option, newValue in newOptions.items():
		if ':' in option:
			subDict, subOpt = option.split(':')
			defaultDictionary[subDict][subOpt] = newValue
		else:
			defaultDictionary[option] = newValue

	return defaultDictionary

def updateLocation(currDict, stationID, rcuMode):
	"""Summary
	
	Args:
	    currDict (TYPE): Description
	    stationID (TYPE): Description
	    rcuMode (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	fileLocations = currDict['fileLocations']

	fileLocations['antennaField'] = fileLocations['antennaField'].format(stationID)
	fileLocations['hbaDeltas'] = fileLocations['hbaDeltas'].format(stationID)
	fileLocations['calibrationLocation'] = fileLocations['calibrationLocation'].format(stationID, rcuMode)
	
	if '{' in fileLocations['outputH5Location']:
		fileLocations['outputH5Location'] = fileLocations['outputH5Location'].format(stationID, rcuMode)

	currDict['fileLocations'] = fileLocations

	return currDict
