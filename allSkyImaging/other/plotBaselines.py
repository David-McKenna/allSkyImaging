"""Summary

Attributes:
    custom_debug (list): Debug list of activated station, increasing in order from 0 to 15 to ensure correct selections
    Effelsberg_elements_20091110 (list): Effelsberg HBA tile activations
    Generic_International_Station_20091110 (list): Generic HBA tile activations
"""
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import copy
import csv
import os

from astropy.coordinates import SkyCoord
randWalk = [ 1, 14, 15, 15, 12, 11, 10,  9,  2, 15,  1,  4,  8,  4,  7,  0,  7,
        0, 14,  0, 11,  3,  2,  9,  5, 11,  3,  3, 11, 12,  3, 10, 12,  7,
       15, 13, 10,  9,  8,  2,  2, 13, 13,  0,  5, 14,  0,  9,  2,  0,  5,
        8,  6,  9,  8, 13, 13, 15,  2, 10, 13,  9, 13, 11,  7,  9, 11, 14,
       13, 13, 12,  0,  0,  2,  8,  7,  6,  2, 10,  0,  4,  9,  2,  4,  3,
        4, 11,  5,  6,  1, 15, 15, 14,  6, 14,  2]


Effelsberg_elements_20091110 = [1,4,13,15,11,9,14,1,15,0,8,2,11,3,14,0,2,4,3,0,0,2,12,12,12,12,15,11,14,15,7,5,1,0,3,10,1,11,0,12,12,1,6,7,0,10,9,6,15,14,11,7,2,0,7,12,15,8,13,3,7,6,3,15,11,1,4,11,8,1,8,15,4,0,5,6,12,0,12,15,3,7,14,8,3,12,12,2,9,8,14,2,5,6,12,0]
Generic_International_Station_20091110 = [15,0,15,3,9,15,14,2,0,3,4,14,10,8,5,15,12,0,2,11,3,12,12,1,5,4,4,8,6,3,0,5,3,11,3,2,8,15,13,8,3,2,9,1,14,8,8,0,12,13,0,11,15,3,12,3,13,3,10,5,0,10,1,6,4,10,3,15,3,14,0,12,0,7,0,12,7,3,13,0,7,3,15,4,14,4,3,8,4,9,12,0,14,9,3,11]
custom_debug = (range(16) * 8)[:len(Effelsberg_elements_20091110)]

def getAntOrder():
	"""Get a list of arrays of increasing order, with lengths corresponding to a row of HBAs
	
	Returns:
	    list: list of 97 elements, 96 tiles and the dummy tile denoted as None
	"""
	lengths =  np.array([5,7,9,11,11,10,11,11,9,7,5])
	initVar = 0
	antennaeSets = []

	for lenRow in lengths:
		antennaeSets.append(np.arange(initVar, initVar + lenRow))
		initVar += lenRow

	antennaeSets[5] = np.concatenate([antennaeSets[5][:5], [None], antennaeSets[5][5:]])
	return antennaeSets

def plotTiles(lengths, titleVar = "HBA Activations", separationLines = False):
	"""Quick HBA layout plotter, noting activations from input dataset
	
	Args:
	    lengths (TYPE): List of rows of activated HBA tile numbers.
	    titleVar (string): Plot title
	    separationLines (bool, optional): Toggle lines between each tile on the graph
	
	Returns:
	    antLoc (np.array): Locations of activated tiles in the arbitrary coordinate system centered about 0
	"""
	antLoc = []

	yOffset = halfEven(len(lengths)) * 5.15
	plt.figure(figsize = (12,12))
	plt.xlim([-36, 36])
	plt.ylim([-36, 36])

	internalOffset = (np.linspace(-0.5, 0.25, 4) + 0.125) * 5
	
	internalOffset = np.array(np.meshgrid(internalOffset, internalOffset)).reshape(1, 2, -1)[0]
	for idx, row in enumerate(lengths):
		xOffset = -5.15 * halfEven(len(row))

		if separationLines:
			plt.axhline(yOffset + 2.5, lw = 5., c = '#2a5289', alpha = 0.3)
		for element in row:
			if idx == 3 and separationLines:
				plt.axvline(xOffset - 2.5, lw = 5., c = '#2a5289', alpha = 0.3)

			if element is not None:

				colEle = element % 4
				rowEle = ((element - colEle) % 16)

				xOffsetEle = xOffset + internalOffset[0][colEle]
				yOffsetEle = yOffset - internalOffset[1][rowEle]

				plt.scatter(xOffset + internalOffset[0], yOffset + internalOffset[1], s = 95, marker = 's', edgecolors = 'k', c = '#327e91', alpha = 0.3)
				plt.scatter([xOffsetEle], [yOffsetEle], s = 100, marker = 's', c = '#00f756', edgecolors = 'k')

				antLoc.append([xOffsetEle, yOffsetEle])

			xOffset += 5.15

		yOffset -= 5.15

	plt.title(titleVar)
	plt.savefig('./{0}.png'.format(titleVar))

	#plt.savefig("{0}.png".format(titleVar.strip(" ")))
	return np.array(antLoc)



def getBaselines(antLoc):
	"""Generate an array of baselines for the given input locations
	
	Args:
	    antLoc (list): List of antenna locations
	
	Returns:
	    baselines (list): List of baselines
	"""
	baselines = []
	for idx, ant in enumerate(antLoc[:-1]):
		for ant2 in antLoc[idx + 1:]:
			diff = ant - ant2
			baselines.append(diff)
			baselines.append(-1. * diff)

	return np.array(baselines)

def getUVWPlane(baselines, astropyCoord):
	"""Convert baselines to a UVW sampling for a given point on the sky
	
	Args:
	    baselines (np.array): List of baselines between antenna
	    astropyCoord (Skycoord): Astropy coordinate on the sky.
	
	Returns:
	    TYPE: Description
	"""
	ra = astropyCoord.icrs.ra.to(u.rad)
	dec = astropyCoord.icrs.dec.to(u.rad)

	transformMaxtrix = np.array([
		[np.sin(ra), np.cos(ra), 0],
		[-1. * np.sin(dec) * np.cos(ra), np.sin(dec) * np.sin(ra), np.cos(dec)],
		[np.cos(dec) * np.cos(ra), -1. * np.cos(dec) * np.sin(ra), np.sin(dec)]
		])

	if baselines.shape[1] != 3:
		transformMaxtrix = transformMaxtrix[:2, :2]

	uvw = np.dot(transformMaxtrix, baselines.T)

	return uvw.T

def populateSelections(lengths, selections):
	"""Given a structure of tiles and a list of activations per-tile, return an array in the 
		same structure as the tile layout but with the selected antenna in place of their ID number
	
	Args:
	    lengths (list): HBA tile structure to repopulate
	    selections (list): List of each antenna activated in each HBA tile (0-15)
	
	Returns:
	    newLengths (list): Repopulated list with the structure of lengths but contents of selections
	
	"""
	idx = 0
	newLengths = copy.deepcopy(lengths)
	for lenIdx, row in enumerate(lengths):
		for rowIdx, element in enumerate(row):
			if element is not None:
				newLengths[lenIdx][rowIdx] = selections[idx]
				idx += 1

	return newLengths

def halfEven(inputVar):
	"""Helper function: return the rounded, half of the input value
	
	Args:
	    inputVar (float-like): Input variable
	
	Returns:
	    float-like: Input /2 rounded down.
	"""
	return (inputVar - (inputVar % 2)) / 2


###


def getAntMap(station = 'IE613'):
	"""Generate a list of array elements for a given station in the ETRS coordinate system (results in meters)
	
	Args:
	    station (str, optional): Name of the station of interest
	
	Returns:
	    fullArr (list): All elements relating to the station
	    hbaLoc (np.array): All HBA coordinates in the given station (x,y,z)
	    lbaLoc (np.array): All LBA coordinates in the given station (x,y,z)
	"""
	fullArr = []
	hbaLoc = []
	lbaLoc = []

	if not os.path.exists("./etrs-antenna-positions.csv"):
		print("Unable to find antenna database! Downloading from Github.")
		import urllib
		urllib.urlretrieve("https://raw.githubusercontent.com/lofar-astron/lofar-antenna-positions/master/share/lofarantpos/etrs-antenna-positions.csv", "./etrs-antenna-positions.csv")
		print("Download Complete.")

	with open("./etrs-antenna-positions.csv", 'rb') as fileRef:
		elementReader = csv.reader(fileRef, delimiter = ',')

		for row in elementReader:
			if row[0] == station:
				fullArr.append(row)

				antType = row[1]
				if antType == 'HBA':
					hbaLoc.append(row[3:6])
				elif antType == 'LBA':
					lbaLoc.append(row[3:6])
				else:
					print('Unknown antenna type {0}'.format(antType))

	hbaLoc = np.array(hbaLoc, dtype = float)
	lbaLoc = np.array(lbaLoc, dtype = float)
	return fullArr, hbaLoc, lbaLoc

def plotTitleSave(dataX, dataY, title, scatterSize = None):
	"""Plot data, give it a title, save it to disk.
	
	Args:
	    dataX (list-like): X data
	    dataY (list-like): Y data
	    title (string): Plot title, name of saved file (without .png suffix)
	    scatterSize (float-like, optional): Size of scatter points
	"""
	plt.figure(figsize = (20, 20))
	plt.title(title)
	if scatterSize:
		plt.scatter(dataX, dataY, s = scatterSize)
	else:
		plt.scatter(dataX, dataY)
	plt.savefig(title + ".png")



if __name__ == '__main__':
	listOfAnt = getAntOrder()
	
	effelsPlotVar = populateSelections(listOfAnt, Effelsberg_elements_20091110)
	genericPlotVar = populateSelections(listOfAnt, Generic_International_Station_20091110)
	debugPlotVar = populateSelections(listOfAnt, custom_debug)
	debugWalkPlotVar = populateSelections(listOfAnt, randWalk)

	antLocEffels = plotTiles(effelsPlotVar, "HBA Activations using the Effelsberg Scheme", True)
	antLocGeneric = plotTiles(genericPlotVar, "HBA Activations using the Generic Scheme", True)
	antLocDebug = plotTiles(debugPlotVar, "HBA Activations for the Debug Scheme", True)
	antLocDebugWalk = plotTiles(debugWalkPlotVar, "HBA Activations for the Debug Random Walk Scheme", True)

	baselinesEffels = getBaselines(antLocEffels)
	baselinesGeneric = getBaselines(antLocGeneric)
	baselinesWalk = getBaselines(antLocDebugWalk)

	#__, __, lbaLoc = getAntMap('DE601')
	#effelsLBABaselines = getBaselines(lbaLoc)

	plotTitleSave(baselinesEffels[:, 0], baselinesEffels[:, 1], "Snapshot Baselines for Effelsberg Scheme", 5)
	plotTitleSave(baselinesGeneric[:, 0], baselinesGeneric[:, 1], "Snapshot Baselines for Generic Scheme", 5)
	plotTitleSave(baselinesWalk[:, 0], baselinesWalk[:, 1], "Snapshot Baselines for Random Walk Scheme", 5)
	#plotTitleSave(lbaLoc[:, 0], lbaLoc[:, 1], "Effelsberg LBA Station")
	#plotTitleSave(effelsLBABaselines[:, 0], effelsLBABaselines[:, 1], "Effelsberg LBA Baselines", 2.)

	plt.show()
