
import h5py
import numpy as np
import os
import multiprocessing as mp

from .dataTools import genericImportTools as importTools
from .dataTools import antennaHandler
from .dataTools import xstImporter
from .dataTools import accImporter
from .dataTools import defaultDict
from .skyImaging import imagingHead


def image(fileLocation, obsType = 'XST', rcuMode = None, deltasLoc = defaultDeltas, fieldLoc = defaultField, plotOptions = plotOptions, activation = None, calLoc = None, outputH5Loc = None, baselineLimits = None):
	

	options, rcuMode, stationName, rotation, activationPattern, calibrationLocation, stationPackage = basicInitialisation(options)
	posXPol, posYPol, stationLocation, __, antLocs = stationPackage

	posX = posXPol[:, 0, np.newaxis]
	posY = posXPol[:, 1, np.newaxis]
	posZ = posXPol[:, 2, np.newaxis]
	antPos = np.dstack([posX, posY, posZ])

	if obsType.lower() == 'xst':
		outputFile, groupPrefix = xstImporter.importXST(fileLocation, outputFile = options['fileLocations']['outputH5Location'], groupNamePrefix = options['h5GroupName'], rcuMode = rcuMode, calibrationFile = calibrationLocation, activationPattern = activationPattern)
	elif obsType.lower() == 'acc':
		'''
		# In the case we have been provided a folder of ACC files, run recursively to process them. 
		if os.path.isdir(fileName):
			fileList, __, __ = processInputLocation(fileName, 'acc')
			for fileName in fileList:
				main(fileName, 'ACC', breakThings, rcuMode, subbandArr, deltasLoc, fieldLoc, plotOptions, activation, calLoc, outputH5Loc, baselineLimits)

			return

		outputFile, groupPrefix, rcuMode = importXST.importACC(fileLocation, rcuMode, calibrationFile = calLoc, outputFile = outputH5Loc)
		'''
	else:
		raise RuntimeError('Unknown observation file type.')

	if not options['plottingOptions']['outputFolder']:
		plotOptions[-1] = '/'.join(outputFile.split('/')[:-1]) + '/allSkyOutput/'


	with h5py.File(outputFile, 'r+') as corrRef:
		if obsType != 'acc':
			subbandArr = [subbandInt[0][2:] for subbandInt in corrRef[groupPrefix].items() if 'sb' in subbandInt[0]]
		else:
			'''
			subbandArr = np.arange(512)
			'''

		if 'calibrationArray' in corrRef[groupPrefix]:
			print('Extracting Calibrations')
			corrGroup = corrRef[groupPrefix]['calibrationArray']
			calibrationArr = corrGroup[...]
		else:
			calibrationArr = None

		if options['multiprocessing'] and len(subbandArr) > 3:
			options['multiprocessing'] = False

			allSkyDatasetNames = np.zeros([len(subbandArr)], dtype = str)
			allSkyDataArr = np.zeros([len(subbandArr)], dtype = object)

			processCount = int(mp.cpu_count() * options['multiprocessingCoreFrac'])
			mpPool = mp.Pool(processes = processCount)
			callBacks = [mpPool.apply_async(processData, args = ([idx, corrRef, calibrationArr, antPos, antLocs, options, rcuMode, stationLocation, rotation])) for idx, fragmentChunk in enumerate(subbandArr)]
	
			mpPool.close()
			mpPool.join()
	
			for asyncResult in callBacks:
				idx, allSkyDatasetName, allSkyDataArr = asyncResult.get()
				allSkyDatasetNames[idx] = allSkyDatasetName[0]
				allSkyDataArr[idx] = allSkyDataArr[0]
			
			allSkyDatasetNames = list(allSkyDatasetNames)
			allSkyDataArr = list(allSkyDataArr)

		else:
			__, allSkyDatasetNames, allSkyDataArr = processData(1, corrRef, calibrationArr, antPos, antLocs, options, rcuMode, stationLocation, rotation)

		for allSkyDatasetName, allSkyData in zip(allSkyDatasetNames, allSkyDataArr):
			outputDataGroup = corrRef.require_dataset(allSkyDatasetName, allSkyData.shape, compression = 'lzf', dtype = np.float64)
			outputDataGroup[...] = datasetData
			outputDataGroup.attrs.create('options', options)

	return outputFile, allSkyDatasetNames

def monitor():


	return

def processData(idx, corrRef, calibrationArr, antPos, antLocs, options, rcuMode, stationLocation, rotation):
	allSkyDatasetNames = []
	allSkyDataArr = []
	for subbandVal in subbandArr:
		if obsType != 'acc':
			corrArr = corrRef[groupPrefix]['sb' + str(subbandVal)]['correlationArray']
		else:
			corrArr = corrRef[groupPrefix]['correlationArray']

		dateArr = np.vstack(corrArr.attrs.values()).astype(str)[:, -1]

		allSkyData = imagingHead.generatePlots(corrArr, [antPos, antLocs], processingOptions, plottingOptions, dateArr, rcuMode, int(subbandVal), stationRotation = rotation, stationLocation = stationLocation, calibrationArr = calibrationArr)

		allSkyDatasetName = '{0}sb{1}/imageData'.format(groupPrefix, subbandVal)
		
		allSkyDatasetNames.append(allSkyDatasetName)
		allSkyDataArr.append(allSkyData)


	return idx, allSkyDatasetNames, imageArr
def basicInitialisation(options):

	if not options:
		options = defaultDict.default()
	elif 'fullStructure' not in options:
		options = defaultDict.patchDefault(options)

	rcuMode = options['rcuMode']
	stationName = options['stationID']

	options = updateLocation(options, stationName, rcuMode)
	importTools.checkRequiredFiles(options['fileLocations'], stationName, rcuMode)

	rotation = importTools.processTxt(lbaRotLoc, stationName)

	if options['imagingOptions']['calibrateData']:
		calibrationLocation = options['fileLocations']['calibrationLocation']
	else:
		calibrationLocation = None

	if rcuMode > 4:
		activation = options['activationPattern']
	else:
		activation = None

	stationPackage = antennaHandler.configureStation(options['fileLocations'], activationPattern, rcuMode, True)

	return options, rcuMode, stationName, rotation, calibrationLocation, stationPackage