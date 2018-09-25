

import ast

def CalcFreq(rcuMode, subband):
	if rcuMode == 5:
		freqoff = 100e6;
		basefreq = 200.0;
	elif rcuMode == 6:
		freqoff = 160e6;
		basefreq = 160.0;
	elif rcuMode == 7:
		freqoff = 200e6;
		basefreq = 200.0;
	else:
		freqoff = 0;
		basefreq = 200.0;
	frequency = ((basefreq / 1024 ) * subband + (freqoff/1e6))
	return frequency, freqoff


def parseiHBAField(afilename, hbaDeltasFile,active_elems, rcuMode, EU):
	# [posxpol, posypol, lon, lat, refpos, rotmat] =
	#     parseiHBAField(filename, hbaDeltasFile, selection, rcuMode, EU)
	#
	# Parser for AntennaFields.conf and iHBADeltas.conf files.
	#
	# arguments
	# filename : filename (including path if necessary) of the AntennaFields.conf
	#            file
	# hbaDeltasFile: filename (including path if necessary) of the iHBADeltas.conf
	#            file
	# selection: element number per tile [0..15]
	# rcuMode  : RCU mode number (determines array configuration)
	# EU       : true for European stations
	#
	# return values
	# posxpol : nElem x 3 matrix with (x, y, z)-positions of x-dipoles (in m)
	#           w.r.t. the reference position in ITRF coordinates
	# posypol : nElem x 3 matrix with (x, y, z)-positions of y-dipoles (in m)
	#           w.r.t. the reference position in ITRF coordinates
	# lon     : geographic longitude of the station (in degrees)
	# lat     : geographic latitude of the station (in degrees)
	# refpos  : nElem x 3 vector with (x, y, z)-coordinate of ITRF reference
	#           position of the specified antenna array (in m)
	# rotmat  : 3 x 3 rotation matrix to convert ITRF to local coordinates
	#
	# SJW, January 2011
	# modified by SJW, January 2012
	# including iHBADeltas.conf by MJN, November 2015
	# update new format Antenna Field size by MJN, June 2017

	# added by Pearse Murphy, August 2018 to account for sparse HBA antennas 
	elemsDict={'Effelsberg_elements_20091110' : [1,4,13,15,11,9,14,1,15,0,8,2,11,3,14,0,2,4,3,0,0,2,12,12,12,12,15,11,14,15,7,5,1,0,3,10,1,11,0,12,12,1,6,7,0,10,9,6,15,14,11,7,2,0,7,12,15,8,13,3,7,6,3,15,11,1,4,11,8,1,8,15,4,0,5,6,12,0,12,15,3,7,14,8,3,12,12,2,9,8,14,2,5,6,12,0], 
	'Generic_International_Station_20091110' : [15,0,15,3,9,15,14,2,0,3,4,14,10,8,5,15,12,0,2,11,3,12,12,1,5,4,4,8,6,3,0,5,3,11,3,2,8,15,13,8,3,2,9,1,14,8,8,0,12,13,0,11,15,3,12,3,13,3,10,5,0,10,1,6,4,10,3,15,3,14,0,12,0,7,0,12,7,3,13,0,7,3,15,4,14,4,3,8,4,9,12,0,14,9,3,11],
	'Generic_Int_201512' : [0,5,3,1,8,3,12,15,10,13,11,5,12,12,5,2,10,8,0,3,5,1,4,0,11,6,2,4,9,14,15,3,7,5,13,15,5,6,5,12,15,7,1,1,14,9,4,9,3,9,3,13,7,14,7,14,2,8,8,0,1,4,2,2,12,15,5,7,6,10,12,3,3,12,7,4,6,0,5,9,1,10,10,11,5,11,7,9,7,6,4,4,15,4,1,15],
	'Generic_Core_201512' : [0,10,4,3,14,0,5,5,3,13,10,3,12,2,7,15,6,14,7,5,7,9,0,15,0,10,4,3,14,0,5,5,3,13,10,3,12,2,7,15,6,14,7,5,7,9,0,15],
	'Generic_Remote_201512' : [0,13,12,4,11,11,7,8,2,7,11,2,10,2,6,3,8,3,1,7,1,15,13,1,11,1,12,7,10,15,8,2,12,13,9,13,4,5,5,12,5,5,9,11,15,12,2,15]}
	
	if activeElems:
		activeElements = elemsDict[activeElems]

		with open(hbaDeltasFile, 'r') as hbaDeltasRef:
			deltaLines = [line for line in hbaDeltasRef]

		arrayDeltas = parseBlitzFile(deltaLines, 'HBADeltas', True)
	
	if rcuMode in [1, 2, 3, 4]:
		arrayName = 'LBA'
	else:
		arrayName = 'HBA'

	with open(afilename, 'r') as arrayRef:
		arrayLines = [line for line in arrayRef]
		
	arrayLoc = parseBlitzFile(arrayLines, arrayName, False)
	antLocs = parseBlitzFile(arrayLines, arrayName, True)

	posX = antLocs[:, 0, :]
	poxY = antLocs[:, 1, :]

	if activeElems:
		poxX += arrayDeltas[activeElements]
		poxY += arrayDeltas[activeElements]
	
	# select the right set of antennas
	if not EU:
		nElem = posX.shape[0]
		if rcuMode in [1, 2]:
			sel = range(nElem/2,nElem) #sel=25:48
		elif rcuMode in [3, 4]:
			sel = range(0,nElem/2)         #sel=1:24
		else:
			sel = range(0,nElem)
		
		posX = posX[sel, :]
		posY = posY[sel, :]
	
#    # obtain rotation matrix
	if 'HBA0\n' in arrayLines:
		####Matlab code. Needs porting for completeness
		#        while ~((strcmp(token, 'HBA1') && strcmp(prevtoken, 'ROTATION_MATRIX')) || isempty(token))
		#            prevtoken = token;
		#            token = fscanf(fid, '#s', 1);
		#        end
		#        fscanf(fid, '#s', 4);
		#        rotmat2 = zeros(3, 3);
		#        rotmat2(:) = fscanf(fid, '#f', 9);
		#        posxpol(1:nElem/2, :) = posxpol(1:nElem/2, :) / rotmat;
		#        posxpol(nElem/2+1:nElem, :) = posxpol(nElem/2+1:nElem, :) / rotmat2;
		#        posypol(1:nElem/2, :) = posypol(1:nElem/2, :) / rotmat;
		#        posypol(nElem/2+1:nElem, :) = posypol(nElem/2+1:nElem, :) / rotmat2;
		####Matlab code end

		# David McKenna: Port attempted based on python implementation of else statement by Joe McCauley. Untested as I don't have matlab to test the actual intended output and haven't pointed as a core/remote station yet.
		rotationMatrixOne = parseBlitzFile(arrayLines, 'ROTATION_MATRIX ' + arrayName + '0', False).T
		rotationMatrixTwo = parseBlitzFile(arrayLines, 'ROTATION_MATRIX ' + arrayName + '0', False).T

		posX = np.matrix(posX)
		posY = np.matrix(posY)

		posX[:, :nElem/2+1] = posX[:, :nElem/2+1] * np.linalg.pinv(rotationMatrixOne)
		posY[:, :nElem/2+1] = posY[:, :nElem/2+1] * np.linalg.pinv(rotationMatrixOne)

		posX[:, nElem/2+1:] = posX[:, nElem/2+1:] * np.linalg.pinv(rotationMatrixTwo)
		posY[:, nElem/2+1:] = posY[:, nElem/2+1:] * np.linalg.pinv(rotationMatrixTwo)

		roationMatrix = np.dstack([rotationMatrixOne, rotationMatrixTwo])

	else:
		rotationMatrix = parseBlitzFile(arrayLines, 'ROTATION_MATRIX ' + arrayName, False).T #has to be transposed to get into correct format for further operations below
		rotationMatrix = np.matrix(rotationMatrix)

		posX = np.matrix(posX)
		posX = posX * np.linalg.pinv( rotationMatrix ) #right matrix division
		posY = posY * np.linalg.pinv( rotationMatrix )

	
	# obtain longitude and latitude of the station
	# David McKenna: I don't want to touch this wizardry.
	wgs84_f = 1 / 298.257223563
	wgs84_a = 6378137
	wgs84_e2 = wgs84_f * (2 - wgs84_f)
	lon = np.arctan2(arrayLoc[1], arrayLoc[0])
	lon=lon* 180 / np.pi
	r = np.sqrt(arrayLoc[0]**2 + arrayLoc[1]**2)
	prev_lat = 100000
	lat = np.arctan2(arrayLoc[2], r)
	while (abs(lat - prev_lat) >= 1e-12):
		prev_lat = lat
		normalized_earth_radius = 1 / np.sqrt((1-wgs84_f)**2 * np.sin(lat)**2 + np.cos(lat)**2)
		lat = np.arctan2(wgs84_e2 * wgs84_a * normalized_earth_radius * np.sin(lat) + arrayLoc[2], r)
	lat = lat * 180 /np.pi

	return posX, posY, lon, lat, arrayLoc, rotationMatrix

def parseBlitzFile(linesArray, keyword, refLoc = False):
	linesArray = [line.strip('\n') for line in linesArray]

	triggerLine = [idx for idx, line in enumerate(linesArray) if line == keyword][0]
	triggerLine += refLoc + 1 # refLoc = add a buffer line to account for a location refernece before processing.


	if ']' in linesArray[triggerLine]:
		line = linesArray[triggerLine]
		startIdx = line.find('[') + 1 # Don't include it as a character
		endIdx = line.find(']')

		dataElements = __processLine(line[startIdx:endIdx])
		return np.array(dataElements)

	endLine = [idx + triggerLine for idx, line in enumerate(linesArray[triggerLine:]) if ']' in line][0]
	print(triggerLine, endLine)
	iterateLines = range(triggerLine + 1, endLine)

	arrayShapeParts = linesArray[triggerLine].count('x') + 1

	startArrayLoc = linesArray[triggerLine].index('[') 
	splitTuples =  linesArray[triggerLine][:startArrayLoc].split('x')

	arrayShape = filter(None, splitTuples)[:arrayShapeParts]
	arrayShape = [ast.literal_eval(strEle.strip(' ')) for strEle in arrayShape]
	arrayShape = [tuplePair[1] - tuplePair[0] + 1 for tuplePair in arrayShape]

	arrayData = []
	for lineIdx in iterateLines:
		procLine = __processLine(linesArray[lineIdx])
		arrayData.append(procLine)

	arrayData = np.vstack(arrayData).reshape(arrayShape)

	return arrayData

def __processLine(line):
	line = filter(None, line.split(' '))
	return [float(element) for element in line]