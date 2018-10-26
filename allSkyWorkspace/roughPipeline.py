
import h5py
import numpy as np
#import allSkyImager_wproj as allSkyImager
import allSkyImager as allSkyImager
import importXST
import os

global defaultDeltas
global defaultField
defaultDeltas = './IE613-iHBADeltas.conf'
defaultField = './IE613-AntennaField.conf'
lbaRotation = './lbaRotations.txt'
plotOptions = [True, 'black', .5, 'black', 'white', True, 0, True, 'IE613', 'None']


reload(importXST)
reload(allSkyImager)
# Extract data from blitz so we don't have to keep referencing them? Store them in the h5 on initial processing?
def main(fileLocation, obsType = 'XST', breakThings = False, rcuMode = None, subbandArr = None, deltasLoc = defaultDeltas, fieldLoc = defaultField, plotOptions = plotOptions, activation = None, calLoc = None, outputH5Loc = None, baselineLimits = None):
	fieldLoc, deltasLoc, lbaRotLoc = checkRequiredFiles(plotOptions[-2])


	rotation = processTxt(lbaRotLoc, plotOptions[-2])
	print(rotation)

	#if plotOptions[-2] != 'Birr':
	#	print('Attempting to change station to {0}'.format(plotOptions[-2]))
	#	deltasLoc = plotOptions[-2].join(defaultDeltas.split('IE613'))
	#	fieldLoc = plotOptions[-2].join(defaultField.split('IE613'))


	if obsType.lower() == 'xst':
		outputFile, groupPrefix, rcuMode = importXST.importXST(fileLocation, rcuMode, calibrationFile = calLoc, outputFile = outputH5Loc)
	elif obsType.lower() == 'acc':
		
		# In the case we have been provided a folder of ACC files, run recursively to process them. 
		if os.path.isdir(fileName):
			fileList, __, __ = processInputLocation(fileName, 'acc')
			for fileName in fileList:
				main(fileName, 'ACC', breakThings, rcuMode, subbandArr, deltasLoc, fieldLoc, plotOptions, activation, calLoc, outputH5Loc, baselineLimits)

			return

		outputFile, groupPrefix, rcuMode = importXST.importACC(fileLocation, rcuMode, calibrationFile = calLoc, outputFile = outputH5Loc)
	else:
		raise RuntimeError('Unknown observation file type.')
	print(deltasLoc)
	posXPol, posYPol, lon, lat, height, arrayLoc, antLocs = allSkyImager.parseiHBAField(fieldLoc, deltasLoc, activation, rcuMode, True)

	print(lon, lat, height)
	
	stationLocation = [lat, lon, height]
	np.save('antLoc.npy', posXPol)

	posX = posXPol[:, 0, np.newaxis]
	posY = posYPol[:, 1, np.newaxis]
	posZ = posXPol[:, 2, np.newaxis]
	antPos = np.dstack([posX, posY, posZ])

	if not plotOptions[-1]:
		plotOptions[-1] = '/'.join(outputFile.split('/')[:-1]) + '/'

	with h5py.File(outputFile, 'r+') as corrRef:
		if obsType != 'acc':
			if subbandArr == None:
				subbandArr = [subbandInt[0][2:] for subbandInt in corrRef[groupPrefix].items() if 'sb' in subbandInt[0]]
			elif instanceof(subbandArr, int):
				subbandArr = [subband]
		else:
			subbandArr = np.arange(512)

		if 'calibrationArray' in corrRef[groupPrefix]:
			print('Extracting Calibrations')
			corrGroup = corrRef['{0}calibrationArray'.format(groupPrefix)]
			print(corrGroup.shape)
			calibrationX, calibrationY = corrGroup[..., 0], corrGroup[..., 1]
		else:
			calibrationX, calibrationY = None, None

		for subbandVal in subbandArr:
			if obsType != 'acc':
				corrArr = corrRef['{0}sb{1}/correlationArray'.format(groupPrefix, subbandVal)]
			else:
				corrArr = corrRef['{0}/correlationArray'.format(groupPrefix)]
			datesArr = np.vstack(corrArr.attrs.values()).astype(str)[:, -1]

			allSkyData = allSkyImager.generatePlots(corrArr, [antPos, antLocs], plotOptions, datesArr, rcuMode, int(subbandVal), calibrationX = calibrationX, calibrationY = calibrationY, baselineLimits = baselineLimits, stationLocation = stationLocation, stationRotation = rotation)
			#print(allSkyData)
			# Currently assuming we will always generate plots for both X and Y polarisations
			allSkySubband = np.stack([allSkyData['X'], allSkyData['Y']], axis = -1)
			imageArr = corrRef.require_dataset('{0}sb{1}/imageData'.format(groupPrefix, subbandVal), allSkySubband.shape, compression = 'lzf', dtype = np.float64)
			imageArr[...] = allSkySubband

def checkRequiredFiles(stationName):
	files = [[defaultField.replace('IE613', stationName), 'https://raw.githubusercontent.com/David-McKenna/SWHT/hbaDeltas/SWHT/data/LOFAR/StaticMetaData/{0}-AntennaField.conf'.format(stationName)], [defaultDeltas.replace('IE613', stationName), 'https://raw.githubusercontent.com/David-McKenna/SWHT/hbaDeltas/SWHT/data/LOFAR/StaticMetaData/iHBADeltas/{0}-iHBADeltas.conf'.format(stationName)], [lbaRotation, 'https://raw.githubusercontent.com/cosmicpudding/lofarimaging/master/stationrotations.txt']]
	currFiles = []
	for filePath, defaultLocation in files:
		if not os.path.exists(filePath):
			newLoc = './' + defaultLocation.split('/')[-1]
			if not os.path.exists(newLoc):
				import urllib
				urllib.urlretrieve(defaultLocation, newLoc)
			currFiles.append(newLoc)
		else:
			currFiles.append(filePath)

	return currFiles

def processTxt(fileLoc, stationName):
	with open(fileLoc, 'r') as fileRef:
		stationLine = [line for line in fileRef if stationName.upper() in line][0]
		stationRotation = -1. * float(stationLine.split(' ')[1][:-2])

	return stationRotation