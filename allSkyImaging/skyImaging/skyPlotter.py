"""Handle matplotlib-based plotting of our sky images.
"""
import numpy as np
import astropy.units as u
import astropy.coordinates
import datetime
import collections
import subprocess
import os
import healpy.newvisufunc

import shutil

from dataTools.usefulFunctions import cache
from dataTools.usefulFunctions import lonToHealpyLon
from dataTools.usefulFunctions import lonToContinuousGalacticLon

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as mplPe

informationArr = {}

# Set some sane default lengths for the vmax/min queues so the global variables are initialised
vmaxCache = collections.deque([], 20)
vminCache = collections.deque([], 20)

def ftPlot(allSkyIm, options, labelOptions, stationLocation, lVec, mVec, dateArr):
	"""Summary
	
	Args:
	    allSkyIm (TYPE): Description
	    options (TYPE): Description
	    labelOptions (TYPE): Description
	    stationLocation (TYPE): Description
	    lVec (TYPE): Description
	    mVec (TYPE): Description
	    dateArr (TYPE): Description
	"""

	#Ensure we are on the right backend
	if options['rfiMode'] or options['plottingOptions']['displayImages']:
		plt.switch_backend('TkAgg')
	else:
		plt.switch_backend('Agg')

	global vmaxCache
	global vminCache

	if isinstance(options['plottingOptions']['colorBarMemory'], int):
		vmaxCache = collections.deque([], options['plottingOptions']['colorBarMemory'])
		vminCache = collections.deque([], options['plottingOptions']['colorBarMemory'])

		print('Colour bars will be averaged over {0} time steps ({1}s of video output if generated)'.format(options['plottingOptions']['colorBarMemory'], options['plottingOptions']['colorBarMemory'] / float(options['plottingOptions']['videoFramerate'])))
	else:
		if 'max' in options['plottingOptions']['colorBarMemory'].lower():
			vmaxCache = np.nanmax(allSkyIm)
			vminCache = np.nanmin(allSkyIm)

			print('Set the overall imaging limits (by max/min) to be between {0} and {1}'.format(vmaxCache, vminCache))
		elif 'percent' in options['plottingOptions']['colorBarMemory'].lower():
			vmaxCache = float(np.nanpercentile(allSkyIm, options['plottingOptions']['maxPercentile']))
			vminCache = float(np.nanpercentile(allSkyIm, options['plottingOptions']['minPercentile']))

			print('Set the overall imaging limits (by percentiles {2}/{3}) to be between {0} and {1}'.format(vmaxCache, vminCache, options['plottingOptions']['minPercentile'], options['plottingOptions']['maxPercentile']))

	fileLoc = []
	figNum = np.random.randint(0, 65000) # If multiprocessed we don't want our figure IDs to overlap.
	print('Begining Plotting: The first frame will take longer than expected as we need to lookup the sky objects.')
	for i in range(allSkyIm.shape[-1]):
		labelOptions[0] = dateArr[i].split("'integrationMidpoint': '")[1].split("', ")[0]
		labelOptions[-1] = figNum
		figNum += 1
		fileLoc.append(ftAllSkyImage(allSkyIm[..., i], options, labelOptions, stationLocation, lVec, mVec))

	if options['plottingOptions']['generateVideo'] and allSkyIm.shape[-1] > 10:
		dateTime, rcuMode, subband, frequency, polarity, figNum = labelOptions
		filePrefix = '/sb_{0}/{1}'.format(subband, dateTime[:2]) # Century of observation as a prefix
		fileSuffix = "mode{0}{1}_{2}_{3}MHz_sb{4}.png".format(rcuMode, polarity, options['imagingOptions']['method'].replace('/', '-'), int(frequency/1e6), subband)

		print("Exporting frames to video at " + "./{0}{1}.mpg".format(dateTime, fileSuffix))
		subprocess.call([options['plottingOptions']['ffmpegLoc'], '-y',  '-framerate',  str(options['plottingOptions']['videoFramerate']), '-pattern_type', 'glob',  '-i',  options['plottingOptions']['outputFolder'] + '{0}*{1}'.format(filePrefix, fileSuffix), '-r', str(max(options['plottingOptions']['videoFramerate'], 30)), options['plottingOptions']['outputFolder'] + "{0}{1}.mp4".format(dateTime, fileSuffix)])


def ftAllSkyImage(allSkyImage, options, labelOptions, stationLocation, lVec, mVec):
	"""Summary
	
	Args:
	    allSkyImage (TYPE): Description
	    options (TYPE): Description
	    labelOptions (TYPE): Description
	    stationLocation (TYPE): Description
	    lVec (TYPE): Description
	    mVec (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	plotOptions = options['plottingOptions']

	logPlot, skyObjColor, gridThickness, backgroundColor, foregroundColor, radialLabelAngle, colorBar, obsSite, outputFolder, maxPercentile, minPercentile, figureShape, fontSizeFactor = plotOptions['logPlot'], plotOptions['skyObjColor'], plotOptions['gridThickness'], plotOptions['backgroundColor'], plotOptions['foregroundColor'], plotOptions['radialLabelAngle'], plotOptions['colorBar'], options['stationID'], plotOptions['outputFolder'], plotOptions['maxPercentile'], plotOptions['minPercentile'], plotOptions['figureShape'], plotOptions['fontSizeFactor']
	dateTime, rcuMode, subband, frequency, polarity, figNum = labelOptions

	if fontSizeFactor is None:
		fontSizeFactor = figureShape[0] / 18.
	if options['rfiMode']:
		fontSizeFactor = min(fontSizeFactor, 0.66)

	telescopeLoc = astropy.coordinates.EarthLocation(lon = stationLocation[0] * u.deg, lat = stationLocation[1] * u.deg, height = stationLocation[2] * u.m)

	if len(dateTime) == 15:
		dateTime = str(datetime.datetime.strptime(dateTime, '%Y%m%d-%H%M%S'))

	obsTime = astropy.time.Time(dateTime)

	if plotOptions['plotSkyObjects']:
		knownSources, referenceObject = getSkyObjects(options, obsTime, telescopeLoc)

	altAzRef = astropy.coordinates.AltAz(obstime = obsTime, location = telescopeLoc)
	altAzObjects = [skyObj.transform_to(altAzRef) for skyObj in referenceObject]

	gs = mpl.gridspec.GridSpec(1, 2, width_ratios = [18, 1])
	gs2 = mpl.gridspec.GridSpec(1, 2, width_ratios = [18, 1])

	fig = plt.figure(figNum, figsize = (figureShape[0], figureShape[1]))
	fig.patch.set_facecolor(backgroundColor)
	plt.suptitle('LOFAR mode {0}{1} all sky plot at {2}MHz (sb{3}) for {4}\n'.format(rcuMode, polarity, round(frequency / 1e6, 2), subband, obsSite), fontsize = int(28 * fontSizeFactor), color=foregroundColor)#, va = 'top') 
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

	if logPlot:
		allSkyImageLog = np.log10(allSkyImage)

		if isinstance(vmaxCache, float):
			vmaxVar = np.log10(vmaxCache)
			vminVar = np.log10(vminCache)
		else:
			vmaxVar = np.nanpercentile(allSkyImageLog, maxPercentile)
			vminVar = np.nanpercentile(allSkyImageLog, minPercentile)

			vminCache.append(vminVar)
			vmaxCache.append(vmaxVar)

			vminVar = np.mean(vminCache)
			vmaxVar = np.mean(vmaxCache)

		pltIm = axImage.imshow(allSkyImageLog, cmap='jet', label = 'ax_image', interpolation = 'lanczos', vmax = vmaxVar, vmin = vminVar)
	else:
		if isinstance(vmaxCache, float):
			vmaxVar = vmaxCache
			vminVar = vminCache
		else:
			vmaxVar = np.nanpercentile(allSkyImage, maxPercentile)
			vminVar = np.nanpercentile(allSkyImage, minPercentile)

			vminCache.append(vminVar)
			vmaxCache.append(vmaxVar)

			vminVar = np.mean(vminCache)
			vmaxVar = np.mean(vmaxCache)

		pltIm = axImage.imshow(allSkyImage, cmap='jet', label = 'ax_image', interpolation = 'lanczos', vmax = vmaxVar, vmin = vminVar)

	axImage.axis('off')
	if colorBar:
		axColorBar = plt.subplot(gs[1])
		colorBarObj = plt.colorbar(pltIm, axColorBar)

		axColorBar.tick_params(which = 'minor', length = 2)
		axColorBar.tick_params(which = 'major', length = 4, width = 1)      
		axColorBar.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(10))

		if options['rfiMode']:
			cbCursor = axColorBar.plot( [ 0, 1 ], [ 0, 0 ], 'k-')
	   
	pltObj = fig.add_subplot(gs2[0], label = 'ax', polar = True)

	axImage.set_xlim((0, lVec.size))
	axImage.set_ylim((lVec.size, 0))
	pltObj.set_theta_zero_location("N")
	pltObj.set_theta_direction(1)

	if plotOptions['plotSkyObjects']:
		plotStatus = [[True, __plotSkyObject(axImage, skyObj, lVec.size, skyObjColor, knownSources[idx], fontSizeFactor = fontSizeFactor)] if (skyObj.alt.deg > 20.) or (options['rfiMode'] and skyObj.alt.deg > -10.) else [False, __plotSkyObject(axImage, skyObj, lVec.size, skyObjColor, knownSources[idx], offset = True, fontSizeFactor = fontSizeFactor)]  for idx, skyObj in enumerate(altAzObjects)]

		legend = axImage.legend( loc = 8, bbox_to_anchor = ( 0.5, -0.128 ), ncol = 4, framealpha = 0.0, fontsize = int(14 * fontSizeFactor), title = str(obsTime)[:-4])
		legend.get_title().set_fontsize(str(int(22 * fontSizeFactor)))

		for idx, skyText in enumerate(legend.get_texts()):
			if not plotStatus[idx][0]:
				plt.setp(skyText, color = 'red')
			elif options['rfiMode'] and altAzObjects[idx].alt.deg > -10. and altAzObjects[idx].alt.deg < 20.:
				plt.setp(skyText, color = 'orange')

	radii = []

	if plotOptions['graticule']:
		for radius in range(0, 90, 15): # r grid at 15 degree intervals
			radii.append(180 * np.cos(radius * np.pi/180)) # plot the radii so as to display as an orthographic grid
			pltObj.set_rgrids(radii)
		if radialLabelAngle: # you would not want to put y ticks on 0 anyhow as it would be messy
			yLabel = [ '', '15' + u'\xb0', '30' + u'\xb0', '45' + u'\xb0', '60' + u'\xb0', '75' + u'\xb0' ]
			pltObj.set_yticklabels(yLabel, color = foregroundColor)
			print(radialLabelAngle)
			pltObj.set_rlabel_position(radialLabelAngle)
		else:
			yLabel = []
			pltObj.set_yticklabels(yLabel)
	
		thetaticks = np.arange(0, 360, 45)
		pltObj.set_thetagrids(thetaticks, weight = 'bold', color = foregroundColor, fontsize = int(18 * fontSizeFactor))
		pltObj.tick_params('x', pad = 0, rotation = 'auto')
	
		pltObj.grid(False, 'both', color = backgroundColor, linewidth = gridThickness)
	
	pltObj.patch.set(alpha = 0.0)
	plt.sca(axImage)

	
	if options['rfiMode']:
		from rfiPlotter import onclick, hover, onaxesleave

		annot = pltObj.annotate( "", xy = ( 0, 0 ), xytext = ( 15, 15 ), textcoords = "offset points", bbox = dict( boxstyle = "round", fc = "w" ), arrowprops = dict( arrowstyle = "->" ) )
		annot.set_visible( False )

		global informationArr
		informationArr = {'pltObj': pltObj, 'annot': annot, 'pixels': options['imagingOptions']['pixelCount'][0], 'cbCursor': cbCursor, 'rawdata': allSkyImage, 'axColorBar': axColorBar}
		onclickLambda = lambda clickEvent: onclick(clickEvent, **informationArr)
		onaxesleaveLambda = lambda clickEvent: onaxesleave(clickEvent, **informationArr)
		hoverLambda = lambda event: hover(event, **informationArr)

		fig.canvas.mpl_connect("motion_notify_event", hoverLambda) 
		fig.canvas.mpl_connect( 'button_press_event', onclickLambda )
		fig.canvas.mpl_connect( 'axes_leave_event', onaxesleaveLambda )
		plt.ioff()
		plt.show()
		plt.pause(1)
	elif plotOptions['displayImages']:
		plt.ioff()
		plt.show()
		plt.pause(1)

	outputFolder += '/sb_{0}/'.format(subband)
	
	if not os.path.exists(outputFolder):
		os.makedirs(outputFolder)

	plotFilename = "{7}{0}_{1}_sb{2}_mode{3}{4}_{5}_{6}MHz.png".format(dateTime, obsSite, subband, rcuMode, polarity, options['imagingOptions']['method'].replace('/', '-'), int(frequency/1e6), outputFolder)
	plotFilename = plotFilename.replace(' ', '_').replace(':', '')
	
	print("Saving output to {0}".format(plotFilename))

	fig.savefig(plotFilename, facecolor=fig.get_facecolor(), edgecolor='none')

	plt.close(figNum)
	return plotFilename


def __plotSkyObject(axIm, skyObj, pixels, skyObjColor, sourceName, offset = False, fontSizeFactor = 1.):
	"""Summary
	
	Args:
	    axIm (TYPE): Description
	    skyObj (TYPE): Description
	    pixels (TYPE): Description
	    skyObjColor (TYPE): Description
	    sourceName (TYPE): Description
	    offset (bool, optional): Description
	    fontSizeFactor (float, optional): Description
	"""
	if not offset:
		if not skyObj.alt.deg > 20.:
			altRad = 0.
		else:
			altRad = skyObj.alt.rad

		rho = np.sin(np.pi / 2. - altRad )
		phi = skyObj.az.rad

		x, y = rho * np.sin(phi), rho * np.cos(phi)
		x, y = (pixels/2) - (pixels/2) * x, (pixels/2) - (pixels/2) * y
	else:
		x = 0. 
		y = 0.

	skyPlt = axIm.scatter(x, y, color = skyObjColor, marker = 'D', s = int(50 * fontSizeFactor), label = u"{0} - Az={1}\xb0, El={2}\xb0".format(sourceName, round(skyObj.az.deg, 1), round(skyObj.alt.deg, 1)), alpha = 1)
	
	if not offset:
			textObj = axIm.annotate(sourceName, xy = (x,y), xytext = (x+2,y+2), color = skyObjColor, fontsize = int(18 * fontSizeFactor))
			textObj.set_path_effects([mplPe.withStroke(linewidth = 5, foreground = 'w')])



def swhtPlot(hpMap, options, labelOptions, plottingFunc = healpy.newvisufunc.mollview, legendPlot = None, zenithArr = None, metaDataArr = None):
	"""Summary
	
	Args:
	    hpMap (TYPE): Description
	    options (TYPE): Description
	    labelOptions (TYPE): Description
	    plottingFunc (TYPE, optional): Description
	    legendPlot (None, optional): Description
	    zenithArr (None, optional): Description
	    metaDataArr (None, optional): Description
	
	Returns:
	    TYPE: Description
	"""
	if options['plottingOptions']['displayImages']:
		plt.switch_backend('TkAgg')
	else:
		plt.switch_backend('Agg')

	fileLoc = []
	figNum = np.random.randint(0, 65000) # If multiprocessed we don't want our figure IDs to overlap.
	print('Begining Plotting: The first frame will take longer than expected as we need to lookup the sky objects.')
	labelOptions[0] = [metaDataArr[0].split("'integrationMidpoint': '")[1].split("', ")[0], metaDataArr[-1].split("'integrationMidpoint': '")[1].split("', ")[0]]
	labelOptions[-1] = figNum
	figNum += 1
	hpMap = hpMap.filled(np.nan)


	fileName1 = swhtFullSkyImage(hpMap, options, labelOptions, plottingFunc, legendPlot, zenithArr)
	fileName2 = swhtFullSkyImage(hpMap, options, labelOptions, plottingFunc, legendPlot = None, zenithArr = None)

	namePattern = '-'.join(fileName2.split('-')[:-1])

	endStr = ''
	if options['plottingOptions']['logPlot']:
		endStr = '-log2'

	shutil.copyfile(fileName1, fileName1[:-4] + '-GIF' + '.png')
	shutil.copyfile(fileName2, fileName2[:-4] + '-GIF' + '.png')

	subprocess.call([options['plottingOptions']['ffmpegLoc'], '-y', '-framerate', '1/' + str(options['plottingOptions']['videoFramerate']), '-pattern_type', 'glob', '-i', '{0}*-GIF.png'.format(namePattern),  '-r', '15', options['plottingOptions']['outputFolder'] + 'swht-animation-{0}-rcu{1}{3}-{2:.1f}Mhz.gif'.format(labelOptions[0][0].split(' ')[0], str(labelOptions[1]) + labelOptions[4], labelOptions[3] / 1e6, endStr)])
	
	os.remove(fileName1[:-4] + '-GIF' + '.png')
	os.remove(fileName2[:-4] + '-GIF' + '.png')


def swhtFullSkyImage(hpMap, options, labelOptions, plottingFunc, legendPlot = None, zenithArr = None):
	"""Summary
	
	Args:
	    hpMap (TYPE): Description
	    options (TYPE): Description
	    labelOptions (TYPE): Description
	    plottingFunc (TYPE): Description
	    legendPlot (None, optional): Description
	    zenithArr (None, optional): Description
	
	Returns:
	    TYPE: Description
	"""
	origHpMap = hpMap.copy()
	plotOptions = options['plottingOptions']

	fileNameSuffix = ''
	
	logPlot, skyObjColor, backgroundColor, foregroundColor, colorBar, obsSite, outputFolder, maxPercentile, minPercentile, figureShape, fontSizeFactor, gridThickness = plotOptions['logPlot'], plotOptions['skyObjColor'], plotOptions['backgroundColor'], plotOptions['foregroundColor'], plotOptions['colorBar'], options['stationID'], plotOptions['outputFolder'], plotOptions['maxPercentile'], plotOptions['minPercentile'], plotOptions['figureShape'], plotOptions['fontSizeFactor'], plotOptions['gridThickness']
	dateTime, rcuMode, subband, frequency, polarity, figNum = labelOptions

	if fontSizeFactor is None:
		fontSizeFactor = figureShape[0] / 18.

	if logPlot:
		fileNameSuffix += '-log2'
		# Cleanup log plots
		hpMapNeg = hpMap < 0.
		hpMap[hpMapNeg] = 1.

		hpMap = np.log2(hpMap)
		if minPercentile != 0:
			minV = np.nanpercentile(hpMap[hpMap > 10.], minPercentile)
		else:
			minV = 1.

		maxV = np.nanpercentile(hpMap[hpMap > 10.], maxPercentile)
	else:
		if minPercentile != 0:
			minV = np.nanpercentile(hpMap, minPercentile)
		else:
			minV = 0.

		maxV = np.nanpercentile(hpMap, maxPercentile)

	hpMap[np.isnan(hpMap)] = minV - 0.5
	cmapVar = plt.get_cmap('jet')

	badCol = cmapVar(0.)
	cmapVar.set_under(badCol, alpha = 0.5)
	cmapVar.set_bad(badCol, alpha = 0.5)

	plt.rcParams["text.color"] = foregroundColor
	plt.rcParams["axes.labelcolor"] = foregroundColor
	plt.rcParams["xtick.color"] =  foregroundColor
	plt.rcParams["ytick.color"] = foregroundColor
	plt.rcParams['axes.edgecolor'] = foregroundColor

	graticuleBool = plotOptions['graticule']
	titleStr = 'LOFAR mode {0}{1} Full Sky Plot at {2}MHz (sb{3}) for {4}\nbetween {5} and {6}\n'.format(rcuMode, polarity, round(frequency / 1e6, 2), subband, obsSite, dateTime[0][:21], dateTime[1][:21])
	useGraticules = graticuleBool and ((legendPlot is not None) or (zenithArr is not None))
	pltArr = plottingFunc(hpMap, cmap = cmapVar, cbar = colorBar, min = minV, max = maxV, graticule = useGraticules, graticule_labels = useGraticules, shading = 'gouraud')
	
	if useGraticules:
		fileNameSuffix += '-graticules'
		plt.gca().set_longitude_grid(30.)
		plt.gca().set_latitude_grid(30.)
	else:
		plt.subplots_adjust(left=0.04, right=0.98, top=0.95, bottom=0.05)

	plt.title(titleStr,fontsize = int(28 * fontSizeFactor), color=foregroundColor)

	figObj = plt.gcf()
	figObj.set_size_inches(figureShape[0], figureShape[1]) # Setting the size before a healpy plot results in the figure size icnreasing, but a constant content size.
	plt.subplots_adjust(top = 1., hspace = 0.)
	figObj.patch.set_facecolor(backgroundColor)
	plt.gca().set_facecolor(backgroundColor)
	plt.gca().patch.set(alpha = 0.0) 

	# If we are plotting zeniths, add a buffer to plot status to account for it.
	plotStatus = {}

	if zenithArr and plotOptions['swhtZenithPointings']:
		fileNameSuffix = '-zenithpointings' + fileNameSuffix
		galCoordLon = lonToContinuousGalacticLon(np.array([skyCoord.l.rad for skyCoord in zenithArr]), True)
		galCoordLat = np.array([skyCoord.b.rad for skyCoord in zenithArr])
		lineVar = plt.plot(galCoordLon, galCoordLat, linestyle = '-', color = foregroundColor, label = None, alpha = 0.4)
		lineVar = plt.scatter(galCoordLon, galCoordLat, marker = 'D', c = foregroundColor, edgecolors = 'lime', label = 'Zenith Pointings', alpha = 0.4)

		plotStatus['Zenith Pointings'] =  True

	if legendPlot and plotOptions['plotSkyObjects']:
		fileNameSuffix = '-skyobjects' + fileNameSuffix
		obsTime, telescopeLoc = legendPlot

		print('Getting sky objects; if not cached this will take a while.')
		knownSources, referenceObject = getSkyObjects(options, obsTime, telescopeLoc)
		plotStatusList = [__plotFullSkyObject(origHpMap, skyObj, skyObjColor, knownSources[idx], fontSizeFactor = fontSizeFactor, lMax = options['imagingOptions']['swhtlMax']) for idx, skyObj in enumerate(referenceObject)]

		plotStatus.update({statusNameBool[0]: statusNameBool[1] for statusNameBool in plotStatusList})

		legend = plt.gca().legend(loc = 8, bbox_to_anchor = ( 0.5, -0.4 ), ncol = 4, title = u'Galactic Coordinates, Plotted L-> R, B -> T: l = [-180, 180]\xb0, b = [-90, 90]\xb0', framealpha = 0.0, fontsize = int(17 * fontSizeFactor))
		
		for idx, skyText in enumerate(legend.get_texts()):
			if plotStatus[skyText.get_text()]:
				plt.setp(skyText, color = 'red')

	outputFolder += '/sb_{0}/'.format(subband)
	if not os.path.exists(outputFolder):
		os.makedirs(outputFolder)

	fileName = '{0}swht-rcuMode{1}-sb{2}{3}.pdf'.format(outputFolder, str(rcuMode) + polarity, subband, fileNameSuffix)
	print(fileName)
	figObj.savefig(fileName, facecolor=figObj.get_facecolor(), edgecolor='none')
	
	if plotOptions['displayImages']:
		plt.show()

	hpMap = origHpMap.copy()
	
	return fileName

def __plotFullSkyObject(hpMap, skyObj, skyObjColor, sourceName, fontSizeFactor = 1., offset = False, lMax = 1):
	"""Summary
	
	Args:
	    hpMap (TYPE): Description
	    skyObj (TYPE): Description
	    skyObjColor (TYPE): Description
	    sourceName (TYPE): Description
	    fontSizeFactor (float, optional): Description
	    offset (bool, optional): Description
	
	Returns:
	    TYPE: Description
	"""
	alpha = 0.
	if isinstance(skyObj, list):
		lDegOrig = [skyObjEle.l.deg for skyObjEle in skyObj]
		lDeg = lonToHealpyLon(lDegOrig, False)
		bDeg = [skyObjEle.b.deg for skyObjEle in skyObj]
		if np.any(~np.isnan(hpMap[healpy.ang2pix(max(64, 2 * lMax), lDeg, bDeg, lonlat = True)])) and not offset:
			alpha = 1.
		else:
			offset = True

		lDeg = lonToContinuousGalacticLon(lDegOrig, False)

		lRad = lonToContinuousGalacticLon([skyObjEle.l.rad for skyObjEle in skyObj], True)
		bRad = [skyObjEle.b.rad for skyObjEle in skyObj]
		skyPlt = plt.plot(lRad, bRad, linestyle = '-', color = skyObjColor, label = None, alpha = alpha)

		labelText = u"{0} - l={1}\xb0, b={2}\xb0".format(sourceName, round(np.mean(lDeg), 1), round(np.mean(bDeg), 1))
		skyPlt_marker = plt.scatter(lRad, bRad, marker = 'X', c = 'w', s = int(50. * fontSizeFactor), edgecolors = skyObjColor, label = labelText, alpha = alpha)

	else:
		lDegOrig = skyObj.l.deg
		lDeg, bDeg = lonToHealpyLon(lDegOrig, False), skyObj.b.deg
		if not np.isnan(hpMap[healpy.ang2pix(max(64, 2 * lMax), lDeg, bDeg, lonlat = True)]) and not offset:
			alpha = 1.
		else:
			offset = True

		lDeg = lonToContinuousGalacticLon(lDegOrig, False)

		lRad, bRad = lonToContinuousGalacticLon(skyObj.l.rad, True), skyObj.b.rad

		labelText = u"{0} - l={1}\xb0, b={2}\xb0".format(sourceName, round(lDeg, 1), round(bDeg, 1))
		skyPlt = plt.scatter(lRad, bRad, marker = 'D', s = int(50 * fontSizeFactor), c = 'w', edgecolors = skyObjColor, label = labelText, alpha = alpha)

	if not offset:
		textObj = plt.gca().annotate(sourceName, xy = (np.mean(lRad),np.mean(bRad)), xytext = (np.mean(lRad),np.mean(bRad)), color = skyObjColor, fontsize = int(18 * fontSizeFactor))
		textObj.set_path_effects([mplPe.withStroke(linewidth = 5, foreground = 'w')])

	return [labelText, offset]


def getSkyObjects(options, obsTime, telescopeLoc):
	"""Summary
	
	Args:
	    options (TYPE): Description
	    obsTime (TYPE): Description
	    telescopeLoc (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	try:
		swhtBool = 'swht' in options['imagingOptions']['method']
		# Far celestial sources can be cached, bodies in the Solar System cannot.
		# Extend an empty lists a subdicts are annoying to handle.
		knownSources = []
		knownSources.extend(options['plottingOptions']['interstellarSources'])
		referenceObject = []
		referenceObject = [cachedSkyCoords(name, swhtBool) for name in knownSources]
	
		knownBodies = options['plottingOptions']['solarSystemSources']
		referenceObject.extend([astropyBodyHandler(sourceName, obsTime, telescopeLoc, swhtBool) for sourceName in knownBodies])

		if 'Polaris' in knownSources:
			knownSources[knownSources.index('Polaris')] = 'NCP (Polaris)' # Cleaner name for the legend
		knownSources.extend(knownBodies)

	except astropy.coordinates.name_resolve.NameResolveError as timeoutException:
		print("Unable to resolve all sources (likely timed out on database lookup), skipping plotting some sources.")
		print(timeoutException)

	return knownSources, referenceObject

def astropyBodyHandler(sourceName, obsTime, telescopeLoc, swhtBool):
	"""Summary
	
	Args:
	    sourceName (TYPE): Description
	    obsTime (TYPE): Description
	    telescopeLoc (TYPE): Description
	    swhtBool (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	if not swhtBool:
		return astropy.coordinates.get_body(sourceName, obsTime, telescopeLoc)
	
	return [astropy.coordinates.get_body(sourceName, obsTick, telescopeLoc).transform_to(astropy.coordinates.Galactic) for obsTick in obsTime]


@cache
def cachedSkyCoords(name, swhtBool = False):
	"""Summary
	
	Args:
	    name (TYPE): Description
	    swhtBool (bool, optional): Description
	
	Returns:
	    TYPE: Description
	"""
	if not swhtBool:
		return astropy.coordinates.SkyCoord.from_name(name)

	return astropy.coordinates.SkyCoord.from_name(name).transform_to(astropy.coordinates.Galactic)
