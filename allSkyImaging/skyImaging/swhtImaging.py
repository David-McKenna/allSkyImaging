"""Summary
"""


def swhtSkyImage(inputCorrelations, xyz, stationLocation, timesArr, lMax, frequency):
	"""Summary
	
	Args:
	    inputCorrelations (TYPE): Description
	    xyz (TYPE): Description
	    stationLocation (TYPE): Description
	    timesArr (TYPE): Description
	    lMax (TYPE): Description
	    frequency (TYPE): Description
	"""
	xyz_raw = xyz.copy()
	xyz = xyz - xyz.transpose(1,0,2)
	xyz_triu = np.triu_indices(xyz.shape[0])

	inputCorrelations[np.eye(inputCorrelations.shape[0], dtype = bool)] = 0.
	inputDropped = inputCorrelations[xyz_triu]
	xyz = xyz[xyz_triu]

	print(stationLocation)
	trueTelescopeLoc = astropy.coordinates.EarthLocation( lat = stationLocation[0] * u.deg, lon = stationLocation[1] * u.deg, height = stationLocation[2] * u.m )
	telescopeLoc = astropy.coordinates.EarthLocation(lat = 0 * u.deg, lon = -90 * u.deg, height = 0 * u.m) # Not sure why, used by Griffin as a baseline for calculations and it works.
	uvw = np.zeros([0] + list(xyz.shape))
	print(uvw.shape)

	zenithArr = []
	altAzArr = []
	for timestamp in timesArr:
		print(timestamp)
		obsTime = astropy.time.Time(timestamp, scale = 'utc', location = telescopeLoc)
		utcObsTime = astropy.time.Time(timestamp, scale = 'utc')

		altAz = astropy.coordinates.AltAz(az = 0. * u.deg, alt = 90. * u.deg, location = trueTelescopeLoc, obstime = utcObsTime)
		zenith = altAz.transform_to(astropy.coordinates.Galactic)

		print(zenith.l.deg, zenith.b.deg)
		sidAngle = float(obsTime.sidereal_time('mean').radian)
		zenithArr.append(zenith)
		altAzArr.append(altAz)

		sinSA = np.sin(sidAngle)
		cosSA = np.cos(sidAngle)
		rotationMatrix = np.array([[sinSA, cosSA, 0.],
						[-1. * cosSA, sinSA, 0.],
						[0., 0., 1.]])

		#raw_input([sidAngle, obsTime, sinSA, cosSA])
		uvw = np.concatenate([uvw, np.dot(rotationMatrix, xyz.T).T[np.newaxis]], axis = 0)

	r, theta, phi = cartToSpherical(uvw)
	r = r[0] # All samples have the same r values

	k_0 = 2. * np.pi * frequency / scipy.constants.c

	preFac = 2. * (k_0 ** 2) / np.pi
	kRVec = k_0 * r

	arraySize = almMap.getsize(lMax - 1)
	results = np.empty(arraySize, dtype = complex)
	ores = np.zeros((lMax, 2* lMax + 1), dtype = complex)
	jv = []
	sph = []

	kZeroBool = kRVec == 0.

	for l in range(lMax):
		j_l = scipy.special.spherical_jn(l, kRVec) # (rLen,)
		j_l[kZeroBool] = 0. # nan otherwise
		
		lRev = (4 * np.pi * (-1j)**l) # (1,)
		print(l)
		jv.append(np.sum(j_l))
		visBessel = inputDropped.flatten() * np.repeat(j_l, uvw.shape[0], axis = 0) 
		for m in range(l + 1):
			y_lm_star = np.conj(scipy.special.sph_harm(m, l, phi.T.flatten(), theta.T.flatten()))# (timeSamples * nants ** 2)
			print(phi.shape, theta.shape)
			sph.append(np.sum(y_lm_star))
			#print(j_l, y_lm_star)

			resultsIndex = almMap.getidx(lMax - 1, l, m)

			#print(inputCorrelations.shape, inputCorrelations[xyz_triu].shape)
			#print(inputCorrelations[xyz_triu].shape, j_l[:, np.newaxis].shape, y_lm_star.T.shape)

			#print(l, m, resultsIndex, lMax)
			results[resultsIndex] = preFac * np.sum(visBessel * y_lm_star) / lRev



def cartToSpherical(uvw):
	print(uvw.shape)

	# r = sqrt(sq(x) + sq(y) + sq(z))
	# theta = acos(z / r)
	# phi = atan(y / x)
	
	#r = np.sqrt(np.sum(np.square(uvw), axis = 2))
	r = np.sqrt(np.sum(np.square(uvw[0]), axis = (1)))
	r = np.repeat(r[np.newaxis], uvw.shape[0], axis = 0)

	theta = np.arccos(uvw[..., 2] / r)
	phi = np.arctan2(uvw[..., 1], uvw[..., 0]) + np.pi # Scipy requires 0 -> 2pi

	zerothBaselines = np.where(r == 0.)
	theta[zerothBaselines] = np.pi / 2.
	phi[zerothBaselines] = np.pi
	return r, theta, phi