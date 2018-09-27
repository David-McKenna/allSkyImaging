import h5py
import numpy as np
import os
import ast
import csv
import datetime

def importXST(fileName, rcuMode, calibrationFile = None, outputFile = None, groupNamePrefix = None, integrationTime = None): # TODO: optional subband split
	if os.path.isdir(fileName):
		fileList = [os.path.join(fileName, fileVar) for fileVar in os.listdir(fileName) if 'xst.dat' in fileVar]
		fileList.sort(key = lambda f: int(filter(str.isdigit, f)))
		fileList.sort(key = lambda f: int(filter(str.isdigit, f.split('sb')[-1]))) # Reorder by subband afterwards.
		fileName = fileList[0]
		print(fileList)
	else:
		fileList = [fileName]

	if outputFile is None:
		outputFile = fileName.split('.dat')[0] + '.h5'

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
				datasetComplex = datasetComplex.reshape(192, 192, reshapeSize, order = 'F')

				fileName = fileName.split('/')[-1]
				fileNameExtract = fileName.split('_')
				dateTime = str(datetime.datetime.strptime(''.join(fileNameExtract[0:2]), '%Y%m%d%H%M%S'))
				subbandStr = fileNameExtract[2]

				dataArr.append(np.array([subbandStr, dateTime, datasetComplex, reshapeSize], dtype = object))

		dataArr= np.vstack(dataArr)
		subbandArr = np.unique(dataArr[:, 0])
		dateTimeArr = [dataArr[:, 1][dataArr[:, 0] == subbandVal] for subbandVal in subbandArr]
		datasetComplexArr = [dataArr[:, 2][dataArr[:, 0] == subbandVal] for subbandVal in subbandArr]
		reshapeSizeArr = [list(dataArr[:, 3][dataArr[:, 0] == subbandVal]) for subbandVal in subbandArr]
				
		for idx, subband in enumerate(subbandArr):
			datasetComplex = np.dstack(datasetComplexArr[idx])

			datasetComplexX = datasetComplex[::2, ::2]
			datasetComplexY = datasetComplex[1::2, 1::2]

			datasetComplex = np.stack([datasetComplexX, datasetComplexY], axis = -1)

			corrDataset = groupRef.require_dataset("{0}/correlationArray".format(subband), datasetComplex.shape, dtype = np.complex128, compression = "lzf")
			
			corrDataset[...] = datasetComplex

			offset = 0
			# Todo: time interpolation for multiple frames in a single file
			for dtIdx, time in enumerate(dateTimeArr[idx]):
				corrDataset.attrs.create(str(offset), time)
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
	return outputFile, groupNamePrefix

def patchKeyValDict(keyValDict):
	keyValDict['ObservationMode'] = int(keyValDict['ObservationMode'])
	keyValDict['ObservationDate'] = str(datetime.datetime.strptime(keyValDict['ObservationDate'], '%Y%m%d%H%M'))
	keyValDict['CalibrationDate'] = str(datetime.datetime.strptime(keyValDict['CalibrationDate'], '%Y%m%d'))
	keyValDict['CalibrationVersion'] = int(keyValDict['CalibrationVersion'])
	keyValDict['CalibrationPPSDelay'] = ast.literal_eval(keyValDict['CalibrationPPSDelay'].replace(' ', ',')[:-2] + ']')

	return keyValDict




#importXST('./20180922_152156_sb200_xst.dat', 5, '/cphys/ugrad/2015-16/JF/MCKENND2/allSkyDump/allSkyDump/Config_Cal/caltables/Cal201804/CalTable-613-HBA-110_190.dat')