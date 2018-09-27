
import h5py
import numpy as np
import allSkyImager
import importXST

defaultDeltas = '/cphys/ugrad/2015-16/JF/MCKENND2/allSkyDump/allSkyDump/Config_Cal/iHBADeltas/IE613-iHBADeltas.conf'
defaultField = '/cphys/ugrad/2015-16/JF/MCKENND2/allSkyDump/allSkyDump/Config_Cal/AntennaFields/IE613-AntennaField.conf'
plotOptions = [True, 'black', .5, 'black', 'white', True, 0, True, 'Birr']


# Extract data from blitz so we don't have to keep referencing them? Store them in the h5 on initial processing?
def main(fileLocation, rcuMode, subbandArr = None, deltasLoc = defaultDeltas, fieldLoc = defaultField, plotOptions = plotOptions, activation = None):
	outputFile, groupPrefix = importXST.importXST(fileLocation, rcuMode)

	posXPol, posYPol, __, __, __, __ = allSkyImager.parseiHBAField(fieldLoc, deltasLoc, activation, rcuMode, True)
	posX = posXPol[:, 0, np.newaxis]
	posY = posYPol[:, 1, np.newaxis]

	antPos = np.dstack([posX, posY])

	with h5py.File(outputFile, 'r') as corrRef:
		if subbandArr == None:
			subbandArr = [subbandInt[0][2:] for subbandInt in corrRef[groupPrefix].items()]
		elif instanceof(subbandArr, int):
			subbandArr = [subband]

		for subbandVal in subbandArr:
			corrArr = corrRef['{0}sb{1}/correlationArray'.format(groupPrefix, subbandVal)]
			datesArr = corrArr.attrs.values()
			allSkyImager.generatePlots(corrArr, antPos, plotOptions, datesArr, rcuMode, subbandVal)

