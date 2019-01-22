"""Options dictionary generator / helper.
"""
import numpy as np

def default():
	"""Generate the default dictionary.
	"""
	fileLocations = {
		'antennaField': './config/{0}-AntennaField.conf',
		'hbaDeltas': './config/{0}-iHBADeltas.conf',
		'lbaRotation': './config/stationrotations.txt',
		'outputH5Location': './{0}-mode{1}.h5',
		'calibrationLocation': './config/{0}-mode{1}.dat', 

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

		'ftSubtractBackground': True,

		'swhtlMax': 32

	}

	plottingOptions = {
		'outputFolder': None,
		'plotImages': True,
		'displayImages': False,
		'generateVideo': True,
		'videoFramerate': 8,
		'ffmpegLoc': 'ffmpeg',

		'figureShape': [18, 14],

		'plotSkyObjects': True,
		'interstellarSources': ['Polaris', 'Cas A', 'Cyg A', 'Sgr A', 'Tau A', 'Vir A', 'Cen A', 'Vela'],
		'solarSystemSources': ['Sun', 'Jupiter', 'Moon', 'Uranus', 'Neptune'],

		'colorBarMemory': 1, #n time steps, maxmin (defaults to 16, 2 seconds of video playback)
		'maxPercentile': 99,
		'minPercentile': 33,

		'logPlot': True,
		'skyObjColor': 'black',
		'gridThickness': 0.5,
		'backgroundColor': 'black',
		'foregroundColor': 'white',
		'radialLabelAngle': -20,
		'fontSizeFactor': None,
		'colorBar': True,
		'graticule': True,

		'swhtZenithPointings': True,

	}

	defaultDictionary = {
		'stationID': 'IE613',
		'rcuMode': 3,
		'activationPattern': 'Generic_Int_201512', #Effelsberg_elements_20091110, Generic_International_Station_20091110, Generic_Core_201512, Generic_Remote_201512
		'h5GroupName': 'allSkyObservation',

		'multiprocessing': True,

		'rfiMode': False,

		'fileLocations': fileLocations,
		'imagingOptions': imagingOptions,
		'plottingOptions': plottingOptions
	}

	return defaultDictionary

def patchDefault(newOptions):
	"""Patch the default dictionary with the input values.

	Input values should be in a dictionary in the format 
	
		{'option': newValue}

	or for a sub-dictionary 

		{'dictName:subOption': newValue}

	This function is autmatically called if you pass a dictionary to the head skyImager module.
	
	Args:
	    newOptions (dict): Dictionary of new options in the syntax described above.
	
	Returns:
	    dict: A patched version of the default dictionary.
	"""

	defaultDictionary = default()
	for option, newValue in newOptions.items():
		if ':' in option:
			subDict, subOpt = option.split(':')
			if (subDict in defaultDictionary) and (subOpt in defaultDictionary[subDict]):
				defaultDictionary[subDict][subOpt] = newValue
			else:
				raise LookupError("The option '{0}' was not found".format(option))
		else:
			if option in defaultDictionary:
				defaultDictionary[option] = newValue
			else:
				raise LookupError("The option '{0}' was not found.".format(option))
	return defaultDictionary

def updateLocation(currDict, stationID, rcuMode):
	"""Update the stationary values to reflect the station being used.
	
	Args:
	    currDict (dict): Current default/patched by patchDefault dictionary
	    stationID (string): Station ID
	    rcuMode (int): RCU mode used for the observation.
	
	Returns:
	    dict: Patched version of the input dictionary.
	"""
	fileLocations = currDict['fileLocations']

	fileLocations['antennaField'] = fileLocations['antennaField'].format(stationID)
	fileLocations['hbaDeltas'] = fileLocations['hbaDeltas'].format(stationID)
	fileLocations['calibrationLocation'] = fileLocations['calibrationLocation'].format(stationID, rcuMode)
	
	if '{' in fileLocations['outputH5Location']:
		fileLocations['outputH5Location'] = fileLocations['outputH5Location'].format(stationID, rcuMode)

	currDict['fileLocations'] = fileLocations

	return currDict
