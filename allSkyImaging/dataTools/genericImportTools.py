"""Generic functions used for importing multiple data types.
"""
import numpy as np
import os
import sys
import select
import ast
import urllib
import datetime


def processInputLocation(fileName, dataType):
	"""Get all files continaing a string in a given directory, or return the folder location if the input is a file.
	
	Args:
	    fileName (string): Folder / file location
	    dataType (string): String to check for (ACC, XST...)
	
	Returns:
	    list(fileList, fileName, folderPath): List of file(s), a sample filename and the folder location.,
	
	Raises:
	    IOError: IOError raised if the location cannot be found.
	"""
	if os.path.isdir(fileName):
		fileList = [os.path.join(fileName, fileVar) for fileVar in os.listdir(fileName) if (fileVar.endswith('.dat') and (dataType.lower() in fileVar))]
		fileList.sort(key = lambda f: int(filter(str.isdigit, f)))
		if not len(fileList):
			raise IOError('No {0} files found in provided directory.'.format(dataType.upper()))
		folderPath = os.path.dirname(os.path.abspath(fileName))
		fileName = fileList[0]

	else:
		fileList = [fileName]
		folderPath = os.path.dirname(os.path.abspath(fileName))

	return fileList, fileName, folderPath


def processRCUMode(folderPath, calFile = None, modeApprox = 3):
	"""Given the calibration file location and folder location for the data, attempt to determine the rcuMode used for the observation.
	
	Args:
	    folderPath (string): Data folder location.
	    calFile (string, optional): Calibration file location.
	    modeApprox (int, optional): Fallback mode
	
	Returns:
	    rcuMode (int): Predicted rcuMode

	"""
	try:
		if not calFile:
			print('Attempting to determine mode from path name')
			modeApprox = folderPath.split('mode')[1][:1]
			modeApprox = int(modeApprox)

		else:
			print('Determining rcuMode from provided Calibration File: {0}'.format(calFile.split('/')[-1]))
			modeDict = {'110': 5, '170': 6, '210': 7, '10': 3, '30': 4}
			modeStr = calFile.split('-')[-1].split('_')[0]
			rcuMode = modeDict[modeStr]

			if rcuMode < 5:
				if 'INNER' in calFile:
					rcuMode -= 2

	except (IndexError, ValueError):
		print("We have been unable to determine the rcuMode form the input variables. Type 'exit' in the next 10 seconds to kill the script, or wait for the script to continue with the assumption of mode {0}.".format(modeApprox))
		inputVar, __, __ = select.select([sys.stdin], [], [], 10)

		if inputVar != 'exit':
			rcuMode = modeApprox
		else:
			exit()

	return rcuMode

def h5PrepNames(outputFile, fileName, groupNamePrefix):
	"""Generate names for the h5 groups / output file.
	
	Args:
	    outputFile (string): Output filename (or None) to cleanup
	    fileName (string): Sample filename
	    groupNamePrefix (string): Sample groupname (or None)
	
	Returns:
	    TYPE: Description
	"""
	if outputFile is None:
		outputFile = fileName.split('.dat')[0] + '.h5'
	elif not '/' in outputFile:
		outputFile = '/'.join(fileName.split('/')[:-1]) + '/' + outputFile

	if groupNamePrefix is None:
		groupNamePrefix = '-'.join(fileName.split('/')[-1].split('_')[:2]) + '/'

	return outputFile, groupNamePrefix

def includeCalibration(calibrationFile, groupRef):
	"""Parse a calibration file and store it in the given h5 group.
	
	Args:
	    calibrationFile (stirng): Calibration file location
	    groupRef (h5Group): h5 Group Reference
	"""
	with open(calibrationFile, 'rb') as calRef:
		headerArray = []
		currOffset = calRef.tell()

		#for line in calRef: # 8192 binary chunks per read, skips over some of our data so we can't break when we find the end of the header.
		while True:
			line = calRef.readline()
	
			# Read until the header
			if 'HeaderStart' in line:
				continue

			# Break at the end of the header
			if 'HeaderStop' in line:
				break

			# Store the header values
			headerArray.append(line)
		# After HeaderStop, we can import the rest of the file as our input data
		calVar = np.fromfile(calRef, dtype = np.complex128)
		rcuCount = calVar.size / 512
		calData = calVar.reshape(rcuCount, 512, order = 'F')

		# X data: even values, y data: odd values.
		calXData = calData[::2]
		calYData = calData[1::2]

		calData = np.dstack([calXData, calYData])

		calDataset = groupRef.require_dataset("calibrationArray", calData.shape, dtype = np.complex128, compression = "lzf")
		calDataset[...] = calData
		headerArray = [headerStr.strip('CalTableHeader').split(' = ') for headerStr in headerArray]

		# Convert the header to a dictionary
		keyValDict = {}
		for key, val in headerArray:
			key = ''.join(key.split('.'))
			val = val.strip('\n')

			keyValDict[key] = val

		# Convert some of the header elements to their natural forms and store them as attributes.
		keyValDict = patchKeyValDict(keyValDict)
		for key, val in keyValDict.items():
			calDataset.attrs.create(key, val)

def patchKeyValDict(keyValDict):
	"""Get header elements and store them in their natural forms.
	
	Args:
	    keyValDict (dict): Dictionary of values generated from a calfile header
	
	Returns:
	    dict: Dictionary with modified values
	"""
	keyValDict['ObservationMode'] = int(keyValDict['ObservationMode'])
	keyValDict['ObservationDate'] = str(datetime.datetime.strptime(keyValDict['ObservationDate'], '%Y%m%d%H%M'))
	keyValDict['CalibrationDate'] = str(datetime.datetime.strptime(keyValDict['CalibrationDate'], '%Y%m%d'))
	keyValDict['CalibrationVersion'] = int(keyValDict['CalibrationVersion'])
	keyValDict['CalibrationPPSDelay'] = ast.literal_eval(keyValDict['CalibrationPPSDelay'].replace(' ', ',')[:-2] + ']')

	return keyValDict 

# Double check rcu mode 2 / 4.... for calibration files
# Store hardcoded naming schemes for the different RCU mode calibration files.
rcuMode2Str = [None, 							# 0
			'CalTable-{0}-LBA_OUTER-10_90.dat', # 1
			'CalTable-{0}-LBA_OUTER-30_90.dat', # 2
			'CalTable-{0}-LBA_INNER-10_90.dat', # 3
			'CalTable-{0}-LBA_INNER-30_90.dat', # 4
			'CalTable-{0}-HBA-110_190.dat', 	# 5
			'CalTable-{0}-HBA-170_230.dat',		# 6
			'CalTable-{0}-HBA-210_250.dat']		# 7

def checkRequiredFiles(fileLocations, stationName, rcuMode):
	"""Ensure we have the files needed to continue, or download them.
	
	Args:
	    fileLocations (dict): Dictionary of expected file locations
	    stationName (string): Station ID
	    rcuMode (int): Observing RCU Mode
	"""
	remoteURLDict = fileLocations['remoteURL']

	for key, value in fileLocations.items():
		if key not in ['outputH5Location', 'remoteURL', 'calibrationLocation']:
			# If the local file doesn't exist, download a copy.
			folderLoc = '/'.join(value.split('/')[:-1])
			fileName = value.split('/')[-1]
			if not os.path.exists(folderLoc):
				os.makedirs(folderLoc)
			if not os.path.exists(value):
				print('Downloading {0} from {1}'.format(fileName, remoteURLDict[key] + fileName))
				urllib.urlretrieve(remoteURLDict[key] + fileName, value)

		elif key == 'calibrationLocation': # Grabbing the calibration file needs a few changes to the call
			if os.path.exists(value):
				continue
			else:
				fileName = rcuMode2Str[rcuMode].format(stationName[2:])
				remoteName = remoteURLDict[key].format(stationName) + fileName

				urllib.urlretrieve(remoteName, value)

def processTxt(fileLoc, stationName):
	"""Used to process the stationrotation.txt file
	
	Args:
	    fileLoc (string): File location
	    stationName (string): Station ID
	
	Returns:
	    float: Station rotation (degrees east of north)
	"""
	with open(fileLoc, 'r') as fileRef:
		stationLine = [line for line in fileRef if stationName.upper() in line][0]
		stationRotation = -1. * float(stationLine.split(' ')[1][:-2])

	return stationRotation

def parseBlitzFile(linesArray, keyword, refLoc = False):
	"""Parse an element of a blitz config given the name of the array
	
	Args:
	    linesArray (list): List of lines in a given config file
	    keyword (string): Keyword to trigger extracting the next array
	    refLoc (bool, optional): if true, add an offset to account for an extra line prior to the array being read
	
	Returns:
	    np.ndarray: Blitz config array converted to numpy array
	"""

	# Strip all the newline characters
	linesArray = [line.strip('\n') for line in linesArray]

	# Find the provided trigger word
	triggerLine = [idx for idx, line in enumerate(linesArray) if line.startswith(keyword)][0]
	triggerLine += refLoc + 1 # refLoc = add a buffer line to account for a location refernece before processing.

	# If there is a ] on the trigger line, we are specifying a data shape. Parse it so we know the correct output shape.
	if ']' in linesArray[triggerLine]:
		line = linesArray[triggerLine]
		startIdx = line.find('[') + 1 # Don't include it as a character
		endIdx = line.find(']')

		dataElements = __processLine(line[startIdx:endIdx])
		return np.array(dataElements)

	# Find the close for the given data chunk
	endLine = [idx + triggerLine for idx, line in enumerate(linesArray[triggerLine:]) if ']' in line][0]
	iterateLines = range(triggerLine + 1, endLine)

	# Determine the number of dimensions in the output shape
	arrayShapeParts = linesArray[triggerLine].count('x') + 1

	# Get the line index of the start line
	startArrayLoc = linesArray[triggerLine].index('[') 

	# Get the elements defining the array shape
	splitTuples =  linesArray[triggerLine][:startArrayLoc].split('x')

	# Filter out spaces / null characters
	arrayShape = filter(None, splitTuples)[:arrayShapeParts]

	# Use ast to parse the values to their true type (int, tuple, string...), but it should always return an int for us here.
	arrayShape = [ast.literal_eval(strEle.strip(' ')) for strEle in arrayShape]
	arrayShape = [tuplePair[1] - tuplePair[0] + 1 for tuplePair in arrayShape]

	# Parse in the data, store in an array and then reshape it to the correct shape.
	# Process line splits up the elements and filters off empty strings, neede as there is an arbitrary number
	#	of spaces between data elements as tab characters are not used to separate them.
	arrayData = []
	for lineIdx in iterateLines:
		procLine = __processLine(linesArray[lineIdx])
		arrayData.append(procLine)

	arrayData = np.vstack(arrayData).reshape(arrayShape)

	return arrayData

def __processLine(line):
	"""Convert a blitz line to be purely float elements
	
	Args:
	    line (string): Blitz config sub-array line
	
	Returns:
	    list(float): Processed line
	"""
	line = filter(None, line.split(' '))
	return [float(element) for element in line] 