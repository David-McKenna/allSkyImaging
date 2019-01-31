"""Generate all of the information required to describe a given LOFAR station.

Originally ported to Python by Joe McCauley, slightly modified for this module. The source contained the following header:

@author: Joe McCauley (joe.mccauley@tcd.ie)
Written for Python 2.7
Based on a translated matlab script originally from ASTRON for processing 
xst data from an international LOFAR station.

Attributes:
    elemsDict (TYPE): Description
"""

import numpy as np

from .genericImportTools import parseBlitzFile

elemsDict={'Effelsberg_elements_20091110' : [1,4,13,15,11,9,14,1,15,0,8,2,11,3,14,0,2,4,3,0,0,2,12,12,12,12,15,11,14,15,7,5,1,0,3,10,1,11,0,12,12,1,6,7,0,10,9,6,15,14,11,7,2,0,7,12,15,8,13,3,7,6,3,15,11,1,4,11,8,1,8,15,4,0,5,6,12,0,12,15,3,7,14,8,3,12,12,2,9,8,14,2,5,6,12,0], 
	'Generic_International_Station_20091110' : [15,0,15,3,9,15,14,2,0,3,4,14,10,8,5,15,12,0,2,11,3,12,12,1,5,4,4,8,6,3,0,5,3,11,3,2,8,15,13,8,3,2,9,1,14,8,8,0,12,13,0,11,15,3,12,3,13,3,10,5,0,10,1,6,4,10,3,15,3,14,0,12,0,7,0,12,7,3,13,0,7,3,15,4,14,4,3,8,4,9,12,0,14,9,3,11],
	'Generic_Int_201512' : [0,5,3,1,8,3,12,15,10,13,11,5,12,12,5,2,10,8,0,3,5,1,4,0,11,6,2,4,9,14,15,3,7,5,13,15,5,6,5,12,15,7,1,1,14,9,4,9,3,9,3,13,7,14,7,14,2,8,8,0,1,4,2,2,12,15,5,7,6,10,12,3,3,12,7,4,6,0,5,9,1,10,10,11,5,11,7,9,7,6,4,4,15,4,1,15],
	'Generic_Core_201512' : [0,10,4,3,14,0,5,5,3,13,10,3,12,2,7,15,6,14,7,5,7,9,0,15,0,10,4,3,14,0,5,5,3,13,10,3,12,2,7,15,6,14,7,5,7,9,0,15],
	'Generic_Remote_201512' : [0,13,12,4,11,11,7,8,2,7,11,2,10,2,6,3,8,3,1,7,1,15,13,1,11,1,12,7,10,15,8,2,12,13,9,13,4,5,5,12,5,5,9,11,15,12,2,15]}

def configureStation(fileLocations, activeElems, rcuMode, EU):
	"""Previously known as parseiHBAField, rewritten to be a lot cleaner and easier to maintain. Renamed to more clearly describe it's purpose.

	Initially ported from ASTRON's MATLab code to Python by Joe McCauley.
	
	Args:
	    fileLocations (list): Location of the *-AntennaField.conf and *-iHBADeltas.conf files
	    activeElems (string): Name of the activation pattern ()
	    rcuMode (int): RCU Mode during observation
	    EU (bool): Are we using an international station? Set this to True.
	
	
	Returns:
	    posXPol, [lon, lat, height], arrayLoc, antLocs: X-Polarisation antenna locations (matches Y), station Location (lon, lat, alt (m)), array location in GCRS, Raw GCRC and Station-Coordinate antenna locations
	"""
	afilename, hbaDeltasFile = fileLocations['antennaField'], fileLocations['hbaDeltas']

	# If we specify an activation pattern, grab it and cache the offsets from tile centers.
	if activeElems:
		activeElements = elemsDict[activeElems]

		with open(hbaDeltasFile, 'r') as hbaDeltasRef:
			deltaLines = [line for line in hbaDeltasRef]

		arrayDeltas = parseBlitzFile(deltaLines, 'HBADeltas', False)

	
	# Classify antenna type by rcuMode
	if rcuMode in [1, 2, 3, 4]:
		arrayName = 'LBA'
	else:
		arrayName = 'HBA'

	# Load in a copy of the array config and parse
	with open(afilename, 'r') as arrayRef:
		arrayLines = [line for line in arrayRef]
		
	arrayLoc = parseBlitzFile(arrayLines, arrayName, False)
	antLocs = parseBlitzFile(arrayLines, arrayName, True)

	# Select the X/Y antenna locations (should be the same) then apply an activation pattern if we defined one
	posX = antLocs[:, 0, :]
	posY = antLocs[:, 1, :]

	if activeElems:
		posX += arrayDeltas[activeElements]
		posY += arrayDeltas[activeElements]

	# For incomplete antenna sets (Core / remote?) only include the relevant antenna due to their reduced RSP count.
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
	
    # Parse and apply the rotation matrix from the blitz file to convert to station coordinates from real coordinates.
    # If 'HBA0' is in the array config, we have a split array defined and thus need to rotate the elements individually.
    # This should only be triggered when EU is false.

    # Use matrices for easier linear algebra operations
	posX = np.matrix(posX)
	posY = np.matrix(posY)
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

		# David McKenna: Port attempted based on python implementation of else statement by Joe McCauley (don't have access to orignal source outside of the block left above). 
		# Untested as I don't have matlab to test the actual intended output and haven't pointed from a core/remote station yet.
		rotationMatrixOne = parseBlitzFile(arrayLines, 'ROTATION_MATRIX ' + arrayName + '0', False).T
		rotationMatrixTwo = parseBlitzFile(arrayLines, 'ROTATION_MATRIX ' + arrayName + '0', False).T

		posX[:, :nElem/2+1] = posX[:, :nElem/2+1] * np.linalg.pinv(rotationMatrixOne)
		posY[:, :nElem/2+1] = posY[:, :nElem/2+1] * np.linalg.pinv(rotationMatrixOne)

		posX[:, nElem/2+1:] = posX[:, nElem/2+1:] * np.linalg.pinv(rotationMatrixTwo)
		posY[:, nElem/2+1:] = posY[:, nElem/2+1:] * np.linalg.pinv(rotationMatrixTwo)

		rotationMatrix = np.dstack([rotationMatrixOne, rotationMatrixTwo])

	else:
		rotationMatrix = parseBlitzFile(arrayLines, 'ROTATION_MATRIX ' + arrayName, False).T #has to be transposed to get into correct format for further operations below
		rotationMatrix = np.matrix(rotationMatrix)

		posX = posX * np.linalg.pinv(rotationMatrix)
		posY = posY * np.linalg.pinv(rotationMatrix)

	posX = np.array(posX)
	posY = np.array(posY)

	# David McKenna: I don't want to touch this wizardry. But I recognise the code from Michiel Brentjens' lofar-antenna-positions.
	# Only changes to below this this block are supplimeneting it with Michiel's height code and variable renaming to make pylint happy.
	wgs84F = 1 / 298.257223563
	wgs84A = 6378137
	wgs84E2 = wgs84F * (2 - wgs84F)
	lon = np.arctan2(arrayLoc[1], arrayLoc[0])
	lon=lon* 180 / np.pi
	wgsRad = np.sqrt(arrayLoc[0]**2 + arrayLoc[1]**2)
	prevLat = 100000
	lat = np.arctan2(arrayLoc[2], wgsRad)
	while (abs(lat - prevLat) >= 1e-12):
		prevLat = lat
		normalisedEarthRadius = 1 / np.sqrt((1-wgs84F)**2 * np.sin(lat)**2 + np.cos(lat)**2)
		lat = np.arctan2(wgs84E2 * wgs84A * normalisedEarthRadius * np.sin(lat) + arrayLoc[2], wgsRad)

	lat = lat * 180 /np.pi
	height = wgsRad * np.cos(lat * np.pi / 180.) + arrayLoc[2] * np.sin(lat * np.pi / 180.) - wgs84A * np.sqrt(1. - wgs84E2 * np.sin(lat * np.pi / 180.) ** 2)

	print(lat, lon, height, arrayLoc)
	return posX, [lon, lat, height], arrayLoc, antLocs
