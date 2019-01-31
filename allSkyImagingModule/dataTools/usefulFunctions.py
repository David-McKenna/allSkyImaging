"""Handy functions used througohut the module.
"""

from functools import wraps
import numpy as np

def cache(function):
	"""Wrapper to cache values of deterministic functions often used with repeating inputs.
	
	Args:
	    function (func): Function to cache
	"""
	cachedLookup = {}
	@wraps(function)
	def trueFunc(*args):
		"""Wrap the actual function with a lookup before calling the original function if the call hasn't been made previously
		"""
		if args in cachedLookup:
			return cachedLookup[args]

		result = function(*args)
		cachedLookup[args] = result

		return result

	return trueFunc

def calcFreq(rcuMode, subband):
	"""Find the frequency of a LOFAR observation given the rcuMode and subband
	
	Args:
	    rcuMode (int): RCU Mode during the observation
	    subband (int): Subband observed
	
	Returns:
	    float: Frequency of observation (in MHz)
	"""

	baseFreq = 100.0

	# Ensure we have cased to ints...
	rcuMode = int(rcuMode)
	subband = int(subband)
	if rcuMode == 5:
		freqOff = 100

	elif rcuMode == 6:
		freqOff = 160
		baseFreq = 80.0

	elif rcuMode == 7:
		freqOff = 200

	else:
		freqOff = 0

	frequency = ((baseFreq / 512 ) * (subband) + (freqOff))

	return frequency

def lonToHealpyLon(pixelTheta, rads = False):
	"""Healpy uses a different coordinate system to Galactic/standard latlon. Convert a given array to be compatible.
	
	Args:
	    pixelTheta (np.ndarray): Array of longitutde values to convert
	    rads (bool, optional): Assume values are in radians if true, degrees if not
	
	Returns:
	    np.ndarray: Converted coordinates
	"""
	# Ensure we have a numpy array
	pixelTheta = np.array(pixelTheta)
	if not rads:
		initThreshold, upperLimit = 180., 360.
	else:
		initThreshold, upperLimit = np.pi, 2. * np.pi

	invLoc = pixelTheta > initThreshold
	# We cannot simply do reverse orders as pixels are allocated on a per-area basis, not a numbers basis
	pixelTheta[invLoc] =  (upperLimit - pixelTheta[invLoc]) * -1.
	pixelTheta[~invLoc] = pixelTheta[~invLoc]

	return pixelTheta

def lonToContinuousGalacticLon(pixelTheta, rads = False):
	"""Galactic longitude is not continuously plotted on healpy maps. Convert our coordinates so they are.

	This is used for matplotlib plotting and human-readable outputs (such as the legend)
	
	Args:
	    pixelTheta (np.ndarray): Array of longitutde values to convert
	    rads (bool, optional): Assume values are in radians if true, degrees if not
	
	Returns:
	    np.ndarray: Converted coordinates
	"""
	pixelTheta = np.array(pixelTheta)
	if not rads:
		initThreshold, upperLimit = 180., 360.
	else:
		initThreshold, upperLimit = np.pi, 2. * np.pi

	invLoc = pixelTheta > initThreshold
	# We cannot simply do reverse orders as pixels are allocated on a per-area basis, not a numbers basis
	pixelTheta[invLoc] =  (upperLimit - pixelTheta[invLoc])
	pixelTheta[~invLoc] *= -1.

	return pixelTheta
