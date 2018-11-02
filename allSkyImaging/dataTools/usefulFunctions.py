"""Summary
"""

def cache(function):
	"""Summary
	
	Args:
	    function (TYPE): Description
	"""
	cachedLookup = {}
	@wraps(function)
	def trueFunc(*args):
		if args in cachedLookup:
			return cachedLookup[args]

		result = function(*args)
		cachedLookup[args] = result

		return result

	return trueFunc

def calcFreq(rcuMode, subband):
	"""Summary
	
	Args:
	    rcuMode (TYPE): Description
	    subband (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	baseFreq = 100.0

	rcuMode = int(rcuMode)
	subband = int(subband)
	if rcuMode == 5:
		freqOff = 100e6

	elif rcuMode == 6:
		freqOff = 160e6
		baseFreq = 80.0

	elif rcuMode == 7:
		freqOff = 200e6

	else:
		freqOff = 0

	frequency = ((baseFreq / 512 ) * (subband) + (freqOff/1e6))

	return frequency