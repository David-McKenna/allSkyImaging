"""Mathematical implementation of the SWHT for full sky images.
"""

import scipy.special
import numpy as np

from healpy.sphtfunc import Alm as almMap

def swhtWorker(idx, lProc, lMax, kRVec, kZeroBool, inputDropped, phi, theta, preFac, results, uvw):
	"""SWHT implementatino for full sky images.
	
	Based on the methodology described by Carrozi 2015
	
	
	Args:
	    idx (int): Multiprocessing ID
	    lProc (int or list-like): l Values to process (list or upto l)
	    lMax (int): Maximum overall l value (for MP reference)
	    kRVec (np.ndarray): Product of k . r for reference
	    kZeroBool (np.ndarray): Locations of r = 0 to fix issues
	    inputDropped (np.ndarray): Input correlations, named as we tend to drop the lower triangle for SWHT observations
	    phi (np.ndarray): Element of UVW plane in psh. coords
	    theta (np.ndarray): Element of UVW plane in psh. coords
	    preFac (float): Constant for output
	    results (np.ndarray): Reference expected output map shape
	    uvw (np.ndarray): Reference for shapes, can be removed if a better reference is found.
	
	Returns:
	    idx, outputHealPyMap: Multiprocessing ID, resulting healpy map
	
	"""
	
	# Multiprocessing gives a list, single process gives an int.
	if isinstance(lProc, int):
		iterL = range(lProc)
	else:
		iterL = lProc

	for counter, l in enumerate(iterL):
		print("Processing l value {0} ({1}/{2}) on Chunk {3}.".format(l, counter + 1, len(iterL), idx + 1))
		j_l = scipy.special.spherical_jn(l, kRVec) # (rLen,)
		j_l[kZeroBool] = 0. # nan otherwise
		
		lRev = (4 * np.pi * (-1j)**l) # (1,)
		visBessel = inputDropped.flatten() * np.repeat(j_l, uvw.shape[0], axis = 0) 

		# Healpy dumps m < 0, so don't both calculating them.
		for m in range(l + 1):
			y_lm_star = np.conj(scipy.special.sph_harm(m, l, phi.T.flatten(), theta.T.flatten()))# (timeSamples * nants ** 2)
			resultsIndex = almMap.getidx(lMax - 1, l, m)
			results[resultsIndex] = preFac * np.sum(visBessel * y_lm_star) / lRev

	print("Chunk {0} work completed for l values {1}".format(idx + 1, iterL))

	return idx, results


# Cython implementation? Couldn't install the libraries on the astro computer so I couldn't test it.
# Probably doesn't work at all considering how bad my C knowledge is, plus the changes made to the above function since I coded it (no longer hard-types inputs)
'''
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _cython_loop(double[:] kRVec, double[:] kZeroBool, double[:, :] phi, double[:, :] theta, int lMax, complex[:] results, float preFac, int[:] index) nogil:

	cdef int l, m
	cdef complex lRev
	cdef complex[:, :, :] j_l, y_lm_star

	for l in prange(lMax):
		j_l = csc.spherical_jn(l, kRVec)
		j_l[kZeroBool] = 0.
		lRev = (4 * np.pi * (-1j)**l)
		for m, in range(i + 1):
			y_lm_star = (csc.sph_harm(m, l, phi, theta)).conjugate()
			results[index[l,m]] = preFac * (inputDropped * j_l * y_lm_star.transpose()).sum() / lRev
'''