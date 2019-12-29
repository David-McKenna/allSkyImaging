"""Generate full sky images with the SWHT methodology.
"""

import astropy.coordinates
import astropy.units as u
import numpy as np
import scipy.constants
import multiprocessing as mp

from healpy.sphtfunc import Alm as almMap
import healpy.newvisufunc
import matplotlib.pyplot as plt

from .workers import swhtWorker
from .skyPlotter import swhtPlot
from . import usefulFunctions

def swhtProcessCorrelations(corrDict, options, processDict, xyz, stationLocation, metaDataArr, frequency, labelOptions):
	"""Handler for converting an observation to a full sky map through the spherical wave harmonic transform.
	
	Args:
	    corrDict (dict): Dictionary of elements in the format 'CORRELATION': (nAnts x nAnts x nSamples array of correlations)
	    options (dict): Standard options dictionary
	    processDict (dict): Dictionary of correlation times ('XX', 'XY'...) and values of whether or not to process them
	    xyz (np.ndarray): Antenna locations nAnts x 3
	    stationLocation (list): Station location in (lon (deg), lat (deg) and elevation (meters))
	    metaDataArr (list): List of dictionaries formatted as strings from the input dataset attributes
	    frequency (float): Observing frequency (in Hz)
	    labelOptions (list): List of useful values for plotting
	
	Returns:
	    dict: Dictionary of processed outputs ('XX', 'I',....)
	"""

	# Borrowed from Griffith's method as an optomisisation, cuts processed values by a factor of 2.
	# Should look into implementing in the FT method? After all, we can recreate the lost correlations
	#	by getting their complex conjugate.
	xyz = xyz - xyz.transpose(1,0,2)
	xyzTriu = np.triu_indices(xyz.shape[0])
	xyz = xyz[xyzTriu]

	# Extract the integration midpoints from the metadata array.
	dateArr = [dateTime.split("'integrationMidpoint': '")[1].split("', ")[0] for dateTime in metaDataArr]

	# We need to use two sets of locations: The true location for the zenith ponting and an arbitrary lon/lat for the uvw plane.
	trueTelescopeLoc = astropy.coordinates.EarthLocation( lat = stationLocation[1] * u.deg, lon = stationLocation[0] * u.deg, height = stationLocation[2] * u.m )
	telescopeLoc = astropy.coordinates.EarthLocation(lat = 0 * u.deg, lon = -90 * u.deg, height = 0 * u.m) # Not sure why, used by Griffin as a baseline for calculations and it works.
	
	# Initialise the UVW array
	uvw = np.zeros([0] + list(xyz.shape))

	zenithArr = []
	timeArr = []
	for timestamp in dateArr:
		# Get two sets of times: One for getting the zenith, the other for plotting the observation title times in local time
		obsTime = astropy.time.Time(timestamp, scale = 'utc', location = telescopeLoc)
		utcObsTime = astropy.time.Time(timestamp, scale = 'utc')


		altAz = astropy.coordinates.AltAz(az = 0. * u.deg, alt = 90. * u.deg, location = trueTelescopeLoc, obstime = utcObsTime)
		zenith = altAz.transform_to(astropy.coordinates.Galactic)
		sidAngle = float(obsTime.sidereal_time('mean').radian)

		zenithArr.append(zenith)
		timeArr.append(obsTime)

		# Rotate the baselines to generate a uvw-like plane (still needs to be divded by wavelength)
		sinSA = np.sin(sidAngle)
		cosSA = np.cos(sidAngle)
		rotationMatrix = np.array([[sinSA, cosSA, 0.],
						[-1. * cosSA, sinSA, 0.],
						[0., 0., 1.]])

		uvw = np.concatenate([uvw, np.dot(rotationMatrix, xyz.T).T[np.newaxis]], axis = 0)
		print('UVW Successfully generated for timestamp {0} with zenith at {1:.2f}, {2:.2f} (GAL)'.format(timestamp[:-5], zenith.l.deg, zenith.b.deg))

	# Convert the UVW coordinates to r, theta, phi (referred to as rtp) for spherical coordinates
	r, theta, phi = cartToSpherical(uvw)
	r = r[0] # All samples have the same r values

	# Generate the k-vector to normalise the r elements
	k_0 = 2. * np.pi * frequency / scipy.constants.c
	preFac = 2. * (k_0 ** 2) / np.pi
	kRVec = k_0 * r

	# Get the maximum l value from the options dictionary and prepare the result arrays
	lMax = options['imagingOptions']['swhtlMax']
	arraySize = almMap.getsize(lMax - 1)
	results = np.zeros(arraySize, dtype = complex)
	maskVar = np.ones(healpy.nside2npix(max(64, 2 * lMax)), dtype=np.bool)

	# Account for autocorrelations (remove excess np.nans from results as r = 0 causes issues to say the least)
	kZeroBool = kRVec == 0.

	# Generate the mask based on the pointings: pixels are masked if they are never within 90 degrees of a zenith pointing
	if options['imagingOptions']['maskOutput']:
		pixelTheta, pixelPhi = healpy.pix2ang(max(64, 2 * lMax), np.arange(healpy.nside2npix(max(64, 2 * lMax))), lonlat = True)
	
		# Fuck everyhting about this coordinate system...
		# For further context: standard latitutde and longitude, healpy and galactic coordinates all have different
		# 	ways of wrapping angles. I created functions to map between the different coordinate systems
		pixelTheta = usefulFunctions.lonToHealpyLon(pixelTheta, rads = False)
	
	
		galCoordLon = np.array([skyCoord.l.deg for skyCoord in zenithArr])
		galCoordLat = np.array([skyCoord.b.deg for skyCoord in zenithArr])

		# Convert to radians for easier angular maths
		galCoordSampled = np.deg2rad(np.vstack([galCoordLon, galCoordLat])[:, np.newaxis])
		pixelCoord = np.deg2rad(np.vstack([pixelTheta, pixelPhi])[..., np.newaxis])
		
		# Check where the angular difference is less than 70 (ie., 35 degrees in any direction., save as used for all sky maps)
		# Not sure about where the factor of 2 is coming from, but I can visually confirm that this works.
		deltaLoc = greatCircleAngularDiff(galCoordSampled, pixelCoord)
		maskVar = np.logical_not(np.any(np.rad2deg(deltaLoc) < 35, axis = 1))


	procDict = {}
	# Create a rotate to convert from arbitrary healpy logic to cartesian-galactic coordinates
	galRotator = healpy.rotator.Rotator(coord = ['C','G'])
	for key, value in processDict.items():
		# processDict is a dictionary of keys/values of correlations and whether or not they are needed for the output images
		# So if we were imaging the Stokes I, we would have {'XX': True, 'YY': True, 'XY': False,...}, making it easy to only processed
		# 	the required correlations.
		if value:
			# Drop the autocorrelations; we have to remove the r = 0 elements to get rid of inf/nans
			corrDict[key][np.eye(corrDict[key].shape[0], dtype = bool)] = 0.
			inputDropped = corrDict[key][xyzTriu]

			# Multithreaded, might not be more efficient if I implement Schaeffer 2015
			if options['multiprocessing']:
				processCount = int(mp.cpu_count()-1)
				mpPool = mp.Pool(processes = processCount)

				fragments = [np.arange(lMax)[i::processCount] for i in range(processCount)]
				callBacks = [mpPool.apply_async(swhtWorker.swhtWorker, args = ([idx, fragmentChunk, lMax, kRVec, kZeroBool, inputDropped, phi, theta, preFac, results, uvw])) for idx, fragmentChunk in enumerate(fragments)]

				mpPool.close()
				mpPool.join()

				for asyncResult in callBacks:
					idx, resultsVar = asyncResult.get()
					results[resultsVar != 0.] = resultsVar[resultsVar != 0.]

				allSkyImageAlm = results

			else:
				allSkyImageAlm = swhtWorker.swhtWorker(0, lMax, lMax, kRVec, kZeroBool, inputDropped, phi, theta, preFac, results, uvw)[1]

			# Convert the results to a healpy map
			procDict['{0}-alm'.format(key)] = allSkyImageAlm
			hpMap = healpy.alm2map(allSkyImageAlm, max(64, 2 * lMax)).astype(complex)

			# Clone + generate imaginary map as well (as healpy seems to hate imaginary values results)
			# Not sure if this is the right way to go about doing things, but it seems to function at least.
			# Can't find a Stokes-V map of the sky to compare.
			clone = allSkyImageAlm.copy()
			clone.real = clone.imag
			clone.imag = np.zeros_like(allSkyImageAlm)
			hpMap.imag = healpy.alm2map(clone, max(64, 2 * lMax))

			# Rotate to the right coordinate system
			hpMap = galRotator.rotate_map(hpMap)

			# Store the results]
			print("{0} Polarisation Processed.".format(key))
			procDict[key] = hpMap

	# Process each of the requested outputs
	for method in options['imagingOptions']['correlationTypes']:
		if len(method) == 2:
			allSkyImage = procDict[method].real
			allSkySave = procDict[method]
		elif 'I' is method:
			allSkySave = (procDict['XX'] + procDict['YY'])
			allSkyImage = allSkySave.real
		elif 'Q' is method:
			allSkySave = (procDict['XX'] - procDict['YY'])
			allSkyImage = allSkySave.real
		elif 'U' is method:
			allSkySave = (procDict['XY'] + procDict['YX'])
			allSkyImage = allSkySave.real
		elif 'V' is method:
			# Healpy doesn't like dealing with imaginary values, take the 'V' output with a mountain of salt.
			allSkySave  = (procDict['YX'] - procDict['XY'])
			allSkyImage = allSkySave.imag
		else:
			print('How did we get this far with a broken method \'{0}\'? Processing the first value in procDict so it hasn\'t gone to waste...'.format(method))
			allSkyImage = procDict.values()[0]
			method = procDict.keys()[0]

		# Mask off values if needed. 
		# For note, doing it before this point breaks stuff as the masked arrays fail to mask on simple maths (values on the other of 1e40 anyone?)
		allSkyImage = healpy.ma(allSkyImage)
		allSkySave = healpy.ma(allSkySave)

		allSkyImage.mask = maskVar
		allSkySave.mask = maskVar

		# Don't save the result out of we are doing a pure-correlation map
		procDict[method] = allSkySave

		if options['plottingOptions']['plotImages']:
			swhtPlot(allSkyImage, options, labelOptions + [method, 0], healpy.newvisufunc.mollview, [timeArr, trueTelescopeLoc], zenithArr, metaDataArr)

	# Save a copy of the mask encase it gets mangled again.
	procDict['mask'] = maskVar

	return procDict



def cartToSpherical(uvw):
	"""Convert cartesian coordinates to spherical coordinates.
	
	Args:
	    uvw (np.ndarray): u[0], v[1], w[2] elements.
	
	Returns:
	    list: rtp version of uvw
	"""

	# Assuming a constant frequency, r doesn't change between samples, it only gets rotated
	# (magnitude is constant, so just calculate for first sample and repeat)
	# r = (x^2 + y^2 + z^2)^0.5
	r = np.sqrt(np.sum(np.square(uvw[0]), axis = (1)))
	r = np.repeat(r[np.newaxis], uvw.shape[0], axis = 0)

	theta = np.arccos(uvw[..., 2] / r) # Error appears due to autocorrelations at r = 0
	phi = np.arctan2(uvw[..., 1], uvw[..., 0]) + np.pi # Scipy requires 0 -> 2pi for sph_harm

	# Set the autocorrelations to sane defaults instead of inf/nan
	zerothBaselines = np.where(r == 0.)
	theta[zerothBaselines] = np.pi / 2.
	phi[zerothBaselines] = np.pi

	return r, theta, phi

def greatCircleAngularDiff(angCoord1, angCoord2):
	"""Get the difference between two sets of angles, where we have the latitutde and longitude stacked in input elements [0] and [1].
	
	Implemented from maths on https://en.wikipedia.org/wiki/Great-circle_distance
	
	Args:
	    angCoord1 (np.ndarray-like): Coordinate set 1 (radians)
	    angCoord2 (np.ndarray-like): Coordinate set 2 (radians)
	
	Returns:
	    np.ndarray: Array of angle differences
	"""
	lon1, lon2 = angCoord1[0], angCoord2[0]
	lat1, lat2 = angCoord1[1], angCoord2[1]

	dlon = lon2 - lon1

	nom = np.square(np.cos(lat2) * np.sin(dlon)) + np.square(np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon))
	denom = np.sin(lat1) * np.sin(lat2) + np.cos(lat1) + np.cos(lat2) * np.cos(dlon)

	dsig = np.arctan2(nom.flatten(), denom.flatten()).reshape(nom.shape)

	return dsig
