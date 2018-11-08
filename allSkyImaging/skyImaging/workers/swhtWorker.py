"""Summary
"""

import scipy.special
import numpy as np

from healpy.sphtfunc import Alm as almMap

def swhtWorker(idx, lProc, lMax, kRVec, kZeroBool, inputDropped, phi, theta, preFac, results, uvw):
	"""Summary
	
	Args:
	    lMax (TYPE): Description
	    kRVec (TYPE): Description
	    kZeroBool (TYPE): Description
	    inputDropped (TYPE): Description
	    phi (TYPE): Description
	    theta (TYPE): Description
	    preFac (TYPE): Description
	    results (TYPE): Description
	    uvw (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
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
		for m in range(l + 1):
			y_lm_star = np.conj(scipy.special.sph_harm(m, l, phi.T.flatten(), theta.T.flatten()))# (timeSamples * nants ** 2)
			resultsIndex = almMap.getidx(lMax - 1, l, m)
			results[resultsIndex] = preFac * np.sum(visBessel * y_lm_star) / lRev

	print("Chunk {0} work completed for l values {1}".format(idx + 1, iterL))

	return idx, results


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