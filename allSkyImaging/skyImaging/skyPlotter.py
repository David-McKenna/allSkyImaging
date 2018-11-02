"""Summary
"""

def __processCorrPlot():
	if calibrationArr is not None:
		print('Calibrating {0} for subband {1} from shape {2}'.format(polChar, subband, calibrationX.shape))
		calSubband = calibrationArr[:, subband]
		calMatrixArr = np.outer(calSubband, np.conj(calSubband).T)[..., np.newaxis]
		inputCorrelations = np.multiply(np.conj(calMatrixArr), inputCorrelations)

	labelOptions = labelOptions + [polChar, 0]
	allSkyIm, __, __ = __processAllSkyIm(inputCorrelationsArr, posX, posY, frequency, lVec, mVec, stationRotation, mask, plotOptions, labelOptions)
	print("{0} Polarisation Processed, begining plotting".format(polChar))

	fileLoc = []
	for i in range(allSkyIm.shape[-1]):
		labelOptions[0] = dateArr[i] + 1
		labelOptions[-1] = figNum
		figNum += 1
		fileLoc.append(plotAllSkyImage(allSkyIm[..., i], plotOptions, labelOptions, pixels, stationLocation, lVec, mVec))

	if plotOptions[5] and allSkyIm.shape[-1] > 20:
		filePrefix = fileLoc[0].split(' ')[0]
		fileSuffix = '_'.join(fileLoc[0].split('/')[-1].split('_')[1:])[:-4]
		print("Exporting frames to video at " + "./{0}{1}.mpg".format(filePrefix, fileSuffix))
		subprocess.call([ffmpegLoc, '-y',  '-r',  '20', '-pattern_type', 'glob',  '-i',  '{0}*{1}.png'.format(filePrefix, fileSuffix), '-vb', '50M', "{0}{1}.mpg".format(filePrefix, fileSuffix)])


def plotAllSkyImage(allSkyImage, plotOptions, labelOptions, pixels, stationLocation, lVec, mVec):
	logPlot, skyObjColor, gridThickness, backgroundColor, foregroundColor, saveImage, radialLabelAngle, colorBar, obsSite, outputFolder = plotOptions
	dateTime, rcuMode, subband, frequency, polarity, figNum = labelOptions

	telescopeLoc = astropy.coordinates.EarthLocation( lat = stationLocation[0] * u.deg, lon = stationLocation[1] * u.deg, height = stationLocation[2] * u.m )
	latLon = telescopeLoc.to_geodetic()

	if len(dateTime) == 15:
		dateTime = str(datetime.datetime.strptime(dateTime, '%Y%m%d-%H%M%S'))

	obsTime = astropy.time.Time(dateTime)

	try:
		knownSources = ['Polaris', 'Cas A', 'Cyg A', 'Sgr A', 'Tau A', 'Vir A', 'Cen A', 'Vela']
		referenceObject = []
		referenceObject = [cachedSkyCoords(name) for name in knownSources]
	
		knownBodies = ['Sun', 'Jupiter', 'Moon', 'Uranus', 'Neptune']
		referenceObject.extend([astropy.coordinates.get_body(sourceName, obsTime, telescopeLoc) for sourceName in knownBodies])
		knownSources[0] = 'NCP (Polaris)'
		knownSources.extend(knownBodies)

	except astropy.coordinates.name_resolve.NameResolveError as timeoutException:
		print("Unable to resolve all sources (likely timed out on database lookup), skipping plotting some sources.")
		print(timeoutException)

	altAzRef = astropy.coordinates.AltAz(obstime = obsTime, location = telescopeLoc)
	altAzObjects = [skyObj.transform_to(altAzRef) for skyObj in referenceObject]

	gs = mpl.gridspec.GridSpec(1, 2, width_ratios = [18, 1])
	gs2 = mpl.gridspec.GridSpec(1, 2, width_ratios = [18, 1])

	print(figNum)
	fig = plt.figure(figNum, figsize = (18, 14))
	fig.patch.set_facecolor(backgroundColor)
	plt.suptitle( 'LOFAR mode {0}{1} all sky plot at {2}MHz (sb{3}) for {4}\n'.format(rcuMode, polarity, round(frequency / 1e6, 2), subband, obsSite), fontsize = 28, color=foregroundColor )#, va = 'top') 
	plt.rcParams["text.color"] = foregroundColor
	plt.rcParams["axes.labelcolor"] = foregroundColor
	plt.rcParams["xtick.color"] =  foregroundColor
	plt.rcParams["ytick.color"] = foregroundColor
	plt.rcParams['axes.edgecolor'] = foregroundColor
	plt.rcParams['axes.linewidth'] = gridThickness

	axImage = fig.add_subplot(gs[0], label = 'ax_image')
	axImage.axis('off')

	global vmaxCache
	global vminCache
	
	lCoord, mCoord = np.meshgrid(lVec.squeeze(), mVec.squeeze())

	print(lCoord, lCoord.shape, mCoord.shape, allSkyImage.shape)
	if logPlot:
		allSkyImageLog = np.log10(allSkyImage)
		#vmaxVar = np.nanpercentile(allSkyImageLog, 99)
		#vminVar = np.nanpercentile(allSkyImageLog, 33)

		#vminCache.append(vminVar)
		#vmaxCache.append(vmaxVar)

		#vminVar = np.mean(vminCache)
		#vmaxVar = np.mean(vmaxCache)

		pltIm = axImage.imshow(allSkyImageLog, cmap='jet', label = 'ax_image', interpolation = 'nearest')
	else:
		vmaxVar = np.nanpercentile(allSkyImage, 99)
		vminVar = np.nanpercentile(allSkyImage, 5)

		#vminCache.append(vminVar)
		#vmaxCache.append(vmaxVar)

		#vminVar = np.mean(vminCache)
		#vmaxVar = np.mean(vmaxCache)
		pltIm = axImage.imshow(allSkyImage, cmap='jet', label = 'ax_image', interpolation = 'nearest')

	axImage.axis('off')
	if colorBar:
		axColorBar = plt.subplot(gs[1])
		colorBarObj = plt.colorbar(pltIm, axColorBar)

		axColorBar.tick_params(which = 'minor', length = 2)
		axColorBar.tick_params(which = 'major', length = 4, width = 1)      
		axColorBar.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(10))
	   
	pltObj = fig.add_subplot(gs2[0], label = 'ax', polar = True)

	axImage.set_xlim((0, pixels))
	axImage.set_ylim((pixels, 0))
	pltObj.set_theta_zero_location("N")
	pltObj.set_theta_direction(1)

	plotStatus = [[True, __plotSkyObject(axImage, skyObj, pixels, skyObjColor, knownSources[idx])] if skyObj.alt.deg > 20. else [False, __plotSkyObject(axImage, skyObj, pixels, skyObjColor, knownSources[idx], offset = True)]  for idx, skyObj in enumerate(altAzObjects)]
	
	plt.rcParams["text.color"] = foregroundColor
	legend = axImage.legend( loc = 8, bbox_to_anchor = ( 0.5, -0.128 ), ncol = 4, framealpha = 0.0, fontsize = 14, title = str(obsTime)[:-4])
	legend.get_title().set_fontsize('22')

	for idx, skyText in enumerate(legend.get_texts()):
		if not plotStatus[idx][0]:
			plt.setp(skyText, color = 'red')
	radii = []

	for r in range(0, 90, 15): # r grid at 15 degree intervals
		radii.append(180 * np.cos(r * np.pi/180)) # plot the radii so as to display as an orthographic grid
		pltObj.set_rgrids(radii)
	if radialLabelAngle: # you would not want to put y ticks on 0 anyhow as it would be messy
		yLabel = [ '', '15' + u'\xb0', '30' + u'\xb0', '45' + u'\xb0', '60' + u'\xb0', '75' + u'\xb0' ]
		pltObj.set_yticklabels(yLabel, color = skyObjColor)
		pltObj.set_rlabel_position(radialLabelAngle)
	else:
		yLabel = []
		pltObj.set_yticklabels(yLabel)

	thetaticks = np.arange(0, 360, 45)
	pltObj.set_thetagrids(thetaticks, weight = 'bold', color = skyObjColor, fontsize = 18)
	pltObj.tick_params('x', pad = -35, rotation = 'auto')

	pltObj.grid(False, 'both', color = skyObjColor, linewidth = gridThickness)
	pltObj.patch.set(alpha = 0.0)
	plt.sca(axImage)

	if saveImage:
		if not os.path.exists(outputFolder):
			os.makedirs(outputFolder)

		plotFilename = "{6}{0}_{1}_sb{2}_mode{3}{4}_{5}MHz.png".format(dateTime, obsSite, subband, rcuMode, polarity, int(frequency/1e6), outputFolder)
		plotFilename = plotFilename.replace(' ', '_').replace(':', '')
		print("Saving output to {0}".format(plotFilename))
		
		fig.savefig(plotFilename, facecolor=fig.get_facecolor(), edgecolor='none')
		plt.close(figNum)
		return plotFilename
	else:
		fig.show()
		plt.close(figNum)
		return


def __plotSkyObject(axIm, skyObj, pixels, skyObjColor, sourceName, offset = False):
	if not offset:
		rho = np.sin(np.pi / 2. - skyObj.alt.rad )
		phi = skyObj.az.rad

		y, x = rho * np.cos(phi), rho * np.sin(phi)
		x, y = (pixels/2) - (pixels/2) * x, (pixels/2) - (pixels/2) * y
	else:
		x = 0. 
		y = 0.

	skyPlt = axIm.scatter(x, y, color = skyObjColor, marker = 'D', s = 50, label = u"{0} - Az={1}\xb0, El={2}\xb0".format(sourceName, round(skyObj.az.deg, 1), round(skyObj.alt.deg, 1)), alpha = 1)
	if not offset:
		textObj = axIm.annotate(sourceName, xy = (x,y), xytext = (x+2,y+2), color = skyObjColor, fontsize = 18)
		textObj.set_path_effects([mplPe.withStroke(linewidth = 5, foreground = 'w')])


def swhtPlotter():
	print('Map Start')
	#results = SWHT.util.array2almVec(ores)
	#hpMap = healpy.ma(healpy.alm2map(results, 64))
	hpMap = healpy.alm2map(results, 64)

	galRotator = healpy.rotator.Rotator(coord = ['C','G'])
	hpMap = galRotator.rotate_map(hpMap)
	hpMapMasked = healpy.ma(hpMap)
	hpMapMasked.mask = np.logical_not(hpMapMasked)
	print('Map End')

	#np.save('./hpMap1.npy', hpMap)
	#np.save('./res1.npy', results)
	#plt.plot(hpMap)
	#plt.show()

	mask = np.zeros(healpy.nside2npix(64), dtype=np.bool)
	print(mask.shape)
	pixel_theta, pixel_phi = healpy.pix2ang(64, np.arange(healpy.nside2npix(64)), lonlat = True)
	

	# Fuck everyhting about this coordinate system...
	invLoc = pixel_theta > 180.
	# We cannot simply do reverse orders as pixels are allocated on a per-area basis, not a numbers basis
	
	pixel_theta[invLoc] =  (360. - pixel_theta[invLoc]) * -1.
	pixel_theta[~invLoc] = pixel_theta[~invLoc]

	#healpy.mollview(pixel_theta)
	#healpy.mollview(pixel_phi)
	#plt.show()
	#pixel_theta = 360. - pixel_theta

	galCoordLon = np.array([skyCoord.l.deg for skyCoord in zenithArr])
	galCoordLat = np.array([skyCoord.b.deg for skyCoord in zenithArr])

	healpy.mollview(pixel_theta)
	healpy.projplot(galCoordLon, galCoordLat, lonlat = True)
	healpy.visufunc.graticule()

	healpy.mollview(pixel_phi)
	healpy.projplot(galCoordLon, galCoordLat, lonlat = True)
	healpy.visufunc.graticule()


	print(zip(list(galCoordLon), list(galCoordLat)))
	plt.figure()
	plt.plot(galCoordLon)
	plt.twinx().plot(galCoordLat, c = 'r')
	plt.show()

	galCoordSampled = np.deg2rad(np.vstack([galCoordLon, galCoordLat])[:, np.newaxis])
	pixelCoord = np.deg2rad(np.vstack([pixel_theta, pixel_phi])[..., np.newaxis])

	print(galCoordSampled.shape, pixelCoord.shape)
	#a_1 = np.square(np.sin(0.5 * (galCoordSampled[0] - pixelCoord[0]))) + np.cos(galCoordSampled[0]) * np.cos(pixelCoord[0]) * np.square(np.sin(0.5 * (galCoordSampled[1] - pixelCoord[1])))
	#a_2 = np.square(np.sin(0.5 * (pixelCoord[0] - galCoordSampled[0]))) + np.cos(galCoordSampled[0]) * np.cos(pixelCoord[0]) * np.square(np.sin(0.5 * (pixelCoord[1] - galCoordSampled[1])))
	
	#deltaLocAntiClockwise = (np.arctan2(np.sqrt(a_1), np.sqrt(1. - a_1)))[..., np.newaxis]
	#deltaLocClockwise = (np.arctan2(np.sqrt(a_2), np.sqrt(1. - a_2)))[..., np.newaxis]
	
	deltaLocAntiClockwise = greatCircleAngular(galCoordSampled, pixelCoord)[..., np.newaxis]
	deltaLocClockwise = greatCircleAngular(galCoordSampled, pixelCoord)[..., np.newaxis]
	print(deltaLocAntiClockwise.shape, deltaLocClockwise.shape)
	deltaLoc = np.stack([deltaLocAntiClockwise, deltaLocClockwise], axis = -1)
	print(deltaLoc.shape)
	deltaLoc = np.nanmin(deltaLoc, axis = -1)[..., 0]

	# When we get close to the equator every forumla for angle I've tried we lose precision, or introduces other errors.
	# 	to account for this issue, if we detect a frame with a sub 0.9 * pi maxima, we will add on the difference
	#	Another approach to this could be a scaling, but on my sample I got better results with an offset as the furthest point approached a delta function
	# 	making scaling have a minimal effect apart form the lat = 0 case
	maxima = np.max(deltaLoc, axis = 0)
	deltaLoc[:, maxima < 0.9 * np.pi] += np.pi - maxima[maxima < 0.9 * np.pi]

	print(deltaLoc.shape)

	healpy.mollview(np.min(deltaLoc, axis = 1))
	healpy.projplot(galCoordLon, galCoordLat, lonlat = True)
	healpy.visufunc.graticule()
	plt.show()
	#deltaLoc = np.arccos(np.sin(galCoordSampled[0]) * np.sin(pixelCoord[0]) + np.cos(galCoordSampled[0]) * np.cos(pixelCoord[0]) * np.cos(galCoordSampled[1] - pixelCoord[1]))
	
	mask = np.logical_not(np.any(np.rad2deg(deltaLoc) < 45, axis = 1)) # Assume 70 degrees field of view in each direction. phi modified to range from0  -> 360 as well for equal weighting.
	#mask = np.logical_not(pixel_phi > -30.)
	hpMapMasked.mask = mask

	#print(mask.shape, hpMap.shape, deltaLoc.shape, galCoordSampled.shape, pixelCoord.shape)
	#print(hpMap.shape, hpMap.shape)
	#hpMap.mask = mask

	#print(hpMap.shape, hpMap.shape, hpMap.real.shape)


	print('Map Reals')

	'''
	for i in range(galCoordSampled.shape[2]):
		healpy.mollview(deltaLoc[:, i])
		print(galCoordLon.shape, galCoordLon[i])
		healpy.projscatter([galCoordLon[i]], [galCoordLat[i]], lonlat = True, c = 'g')
		healpy.visufunc.graticule()
		plt.savefig('./{0}.png'.format(i))
 		plt.clf()

 	plt.close('all')
	'''
	'''
	skyCoord = cachedSkyCoords('Cas A')
	print(skyCoord, skyCoord.galactic)
	for bodyArr in plotLists:
		print(bodyArr)
		healpy.projplot(np.array(bodyArr).T, lonlat = True)
	'''
	healpy.mollview(hpMap.real)
	healpy.visufunc.graticule()
	plt.show()
	healpy.mollview(hpMapMasked.mask)
	healpy.projplot(galCoordLon, galCoordLat, lonlat = True, coord = 'G')
	healpy.visufunc.graticule()
	plt.show()
	healpy.mollview(hpMapMasked)
	healpy.projplot(galCoordLon, galCoordLat, lonlat = True)
	healpy.visufunc.graticule()

	skyCoord = astropy.coordinates.SkyCoord.from_name('Cas A')
	print(skyCoord.galactic)
	healpy.projscatter([skyCoord.galactic.l.deg], [skyCoord.galactic.b.deg], lonlat = True, c = 'r')

	skyCoord = astropy.coordinates.SkyCoord.from_name('Cyg A')
	print(skyCoord.galactic)
	healpy.projscatter([skyCoord.galactic.l.deg], [skyCoord.galactic.b.deg], lonlat = True, c = 'g')
	
	skyCoord = astropy.coordinates.SkyCoord.from_name('Vir A')
	print(skyCoord.galactic)
	healpy.projscatter([skyCoord.galactic.l.deg], [skyCoord.galactic.b.deg], lonlat = True, c = 'b')

	skyCoord = astropy.coordinates.SkyCoord.from_name('Polaris')
	print(skyCoord.galactic)
	healpy.projscatter([skyCoord.galactic.l.deg], [skyCoord.galactic.b.deg], lonlat = True, c = 'm')

	plt.show()
	#healpy.mollview(hpMap.real,deg=True, rot = [90, 0], coord = 'CG')
	#healpy.mollview(hpMap.real, deg=True, rot = [0, 90], coord = 'CG')



	R = np.cos(np.square(pixel_phi))
	X = R * np.cos(pixel_phi) * np.cos(pixel_theta)
	Y = R * np.sin(pixel_phi) * np.sin(pixel_theta)
	Z = R * np.cos(pixel_phi)[..., np.newaxis]

	X = X.reshape(32, 64)
	Y = Y.reshape(32, 64)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection = '3d')

	minVal = -1. * np.ma.min(hpMap.real)
	norm_values = (hpMapMasked.real + minVal) / (np.ma.max(hpMapMasked.real) + minVal)
	ax.plot_surface(X, Y, Z, facecolors = mpl.cm.jet(norm_values))

	plt.show()
	exit;

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


def greatCircleAngular(angCoord1, angCoord2):
	# Wikipedia 
	lon1, lon2 = angCoord1[0], angCoord2[0]
	lat1, lat2 = angCoord1[1], angCoord2[1]

	dlat = lat2 - lat1
	dlon = lon2 - lon1

	nom = np.square(np.cos(lat2) * np.sin(dlon)) + np.square(np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon))
	denom = np.sin(lat1) * np.sin(lat2) + np.cos(lat1) + np.cos(lat2) * np.cos(dlon)

	inputShape = nom.shape
	dsig = np.arctan2(nom.flatten(), denom.flatten()).reshape(inputShape)
	return dsig
