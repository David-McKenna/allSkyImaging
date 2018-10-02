import h5py
import numpy as np
import os
import ast
import csv
import datetime

def importXST(fileName, rcuMode = None, calibrationFile = None, outputFile = None, groupNamePrefix = None, integrationTime = None): # TODO: optional subband split
	if os.path.isdir(fileName):
		fileList = [os.path.join(fileName, fileVar) for fileVar in os.listdir(fileName) if ('xst.dat' in fileVar) and not ('.log' in fileVar)]
		fileList.sort(key = lambda f: int(filter(str.isdigit, f)))
		fileList.sort(key = lambda f: int(filter(str.isdigit, f.split('sb')[-1]))) # Reorder by subband afterwards. Shouldn't be needed anymore, but it's nice to keep for peace of mind.
		if not len(fileList):
			raise IOError('No files found in provided directory.')
		folderPath = fileName
		fileName = fileList[0]
		#print(fileList)

	else:
		fileList = [fileName]
		folderPath = os.path.dirname(os.path.abspath(fileName))


	try:
		testRef = open(fileList[0] + '.log', 'rb')
		logFiles = False # Debug ast later...
		raise IOError
	except IOError:
		print('Unable to open log files, we will make assumptions for the observation\'s metadata.')
		logFiles = False

		if not rcuMode:
			try:
				print('Attempting to guess mode from path name')
				modeApprox = fileList[0].split('mode')[1][:1]
				modeApprox = int(modeApprox)
				print(modeApprox)
			except (IndexError, ValueError):
				print('Assuming Mode 3 observation.')
				modeApprox = 3
			rcuMode = modeApprox
		else:
			print('Using provided rcu mode, {0}'.format(rcuMode))
			modeApprox = rcuMode

		# mode: [mode, subband, integration time]
		# Assuming one frame per observation for non-logged files
		metadata = {'1': [1, 100, 5], '2': [2, 100, 5], '3': [3, 100, 5], '4': [4, 100, 5], '5': [5, 200, 10], '6': [6, 200, 10], '7': [7, 200, 10]}
		metadata = metadata[str(modeApprox)]

	if 'HBA_elements.log' in os.listdir(folderPath):
		with open(os.path.join(folderPath, 'HBA_elements.log')) as actRef:
			fileLines = [line for line in actRef]
			#assert(np.array(['Traceback' not in line for line in fileLines]).all())
			
			## TODO
			



	if outputFile is None:
		outputFile = fileName.split('.dat')[0] + '.h5'
	elif not '/' in outputFile:
		outputFile = '/'.join(fileName.split('/')[:-1]) + '/' + outputFile

	if groupNamePrefix is None:
		groupNamePrefix = '-'.join(fileName.split('/')[-1].split('_')[:2]) + '/'

	with h5py.File(outputFile, 'a') as outputRef:
		datasetComplexHead = {'processTime': str(datetime.datetime.utcnow())}
		groupRef = outputRef.require_group(groupNamePrefix)

		dataArr = []
		for fileName in fileList:
			with open(fileName, 'rb') as dataRef:
				datasetComplex = np.fromfile(dataRef, dtype = np.complex128)
				reshapeSize = datasetComplex.size / (192 ** 2)
				if datasetComplex.size < 1000:
					print('IMCOMPLETE FILE SKIPPED: {0}'.format(fileName))
					continue
				datasetComplex = datasetComplex.reshape(192, 192, reshapeSize, 
order = 'F')

				fileNameMod = fileName.split('/')[-1]
				fileNameExtract = fileNameMod.split('_')
				dateTime = str(datetime.datetime.strptime(''.join(fileNameExtract[0:2]), '%Y%m%d%H%M%S'))
				subbandStr = fileNameExtract[2]
				if False:
				#if logFiles:
					logName = fileName + '.log'
					with open(logName, 'rb') as logRef:
						logData = [] # AST debug...
						print(logRef)
						for line in logRef:
							print(line)
							#logData.append(ast.literal_eval(line.split(' ')[-1]))
						print(logName)
						#logData = [ast.literal_eval(line.split(' ')[-1].strip('\n')) for line in logRef]
						#assert(logData[0] == rcuMode)
						#assert(logData[1] == int(subband[2:]))
						#logData[2] = datetime.timedelta(seconds = logData[2])
						#logData[3] = datetime.datetime.strptime(logData[3], '%Y/%m/%d@%H:%M:%S')
						#logData[4] = datetime.datetime.strptime(logData[4], '%Y/%m/%d@%H:%M:%S')

				else:
					logData = metadata
					logData[1] = int(subbandStr[2:])

				dataArr.append(np.array([subbandStr, dateTime, datasetComplex, reshapeSize, logData], dtype = object))


		dataArr= np.vstack(dataArr)
		subbandArr = np.unique(dataArr[:, 0])
		dateTimeArr = [dataArr[:, 1][dataArr[:, 0] == subbandVal] for subbandVal in subbandArr]
		datasetComplexArr = [dataArr[:, 2][dataArr[:, 0] == subbandVal] for subbandVal in subbandArr]
		reshapeSizeArr = [list(dataArr[:, 3][dataArr[:, 0] == subbandVal]) for subbandVal in subbandArr]
		logDataArr = [list(dataArr[:, 4][dataArr[:, 0] == subbandVal]) for subbandVal in subbandArr]
				
		for idx, subband in enumerate(subbandArr):
			datasetComplex = np.dstack(datasetComplexArr[idx])

			datasetComplexX = datasetComplex[::2, ::2, ...]

			datasetComplexY = datasetComplex[1::2, 1::2, ...]
			datasetComplex = np.stack([datasetComplexX, datasetComplexY], axis = -1)
			corrDataset = groupRef.require_dataset("{0}/correlationArray".format(subband), datasetComplex.shape, dtype = np.complex128, compression = "lzf")
			
			corrDataset[...] = datasetComplex

			offset = 0
			if logFiles:
				for dtIdx, logData in enumerate(logDataArr[idx]):
					currDelta = reshapeSizeArr[idx][dtIdx]
					mode, subband, intTime, startTime, endTime = logData[idx][dtIdx]

					if currDelta > 1:
						timeArr = [startTime + nTimes * intTime for nTimes in range(currDelta)]
					else:
						timeArr = [startTime]

					for i in range(currDelta):
						corrDataset.attrs.create(str(offset + i), [mode, subband, intTime, timeArr[i]])

					offset += currDelta

			else:
				# Todo: time interpolation for multiple frames in a single file witohut metadata
				for dtIdx, time in enumerate(dateTimeArr[idx]):
					mode, subband, intTime = logDataArr[idx][dtIdx]
					corrDataset.attrs.create(str(offset), [mode, subband, intTime, time])

					offset += reshapeSizeArr[idx][dtIdx]


		if calibrationFile is not None:
			with open(calibrationFile, 'rb') as calRef:
				headerArray = []
				currOffset = calRef.tell()

				#for line in calRef: # 8192 binary chunks per read, skips over some of our data so we can't break when we find the end of the header.

				while True:
					line = calRef.readline()

					if 'HeaderStart' in line:
						continue

					if 'HeaderStop' in line:
						break

					headerArray.append(line)

				calVar = np.fromfile(calRef, dtype = np.complex128)
				rcuCount = calVar.size / 512
				calData = calVar.reshape(rcuCount, 512, order = 'F')
				calXData = calData[::2]
				calYData = calData[1::2]

				calData = np.dstack([calXData, calYData])

			calDataset = groupRef.require_dataset("calibrationArray", calData.shape, dtype = np.complex128, compression = "lzf")
			calDataset[...] = calData
			headerArray = [headerStr.strip('CalTableHeader').split(' = ') for headerStr in headerArray]

			keyValDict = {}
			for key, val in headerArray:
				key = ''.join(key.split('.'))
				val = val.strip('\n')

				keyValDict[key] = val

			keyValDict = patchKeyValDict(keyValDict)
			for key, val in keyValDict.items():
				calDataset.attrs.create(key, val)

	return outputFile, groupNamePrefix, rcuMode

def patchKeyValDict(keyValDict):
	keyValDict['ObservationMode'] = int(keyValDict['ObservationMode'])
	keyValDict['ObservationDate'] = str(datetime.datetime.strptime(keyValDict['ObservationDate'], '%Y%m%d%H%M'))
	keyValDict['CalibrationDate'] = str(datetime.datetime.strptime(keyValDict['CalibrationDate'], '%Y%m%d'))
	keyValDict['CalibrationVersion'] = int(keyValDict['CalibrationVersion'])
	keyValDict['CalibrationPPSDelay'] = ast.literal_eval(keyValDict['CalibrationPPSDelay'].replace(' ', ',')[:-2] + ']')

	return keyValDict




#importXST('./20180922_152156_sb200_xst.dat', 5, '/cphys/ugrad/2015-16/JF/MCKENND2/allSkyDump/allSkyDump/Config_Cal/caltables/Cal201804/CalTable-613-HBA-110_190.dat')
