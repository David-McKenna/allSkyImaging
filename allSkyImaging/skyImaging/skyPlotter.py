"""Summary

Attributes:
	vmaxCache (TYPE): Description
	vminCache (TYPE): Description
"""
import numpy as np
import astropy.units as u
import astropy.coordinates
import datetime
import collections
import subprocess
import os
import healpy.newvisufunc
import ast

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
	"""

	#Ensure we are on the right backend
	if options['rfiMode'] or options['plottingOptions']['displayFigure']:
		plt.switch_backend('TkAgg')
	else:
		plt.switch_backend('Agg')

	global vmaxCache
	global vminCache

	print(options['plottingOptions']['colorBarLimits'].lower())
	if isinstance(options['plottingOptions']['colorBarLimits'], int):
		vmaxCache = collections.deque([], options['plottingOptions']['colorBarLimits'])
		vminCache = collections.deque([], options['plottingOptions']['colorBarLimits'])

		print('Colour bars will be averaged over {0} time steps ({1}s of video output if generated)'.format(options['plottingOptions']['colorBarLimits'], options['plottingOptions']['colorBarLimits'] / float(options['plottingOptions']['videoFramerate'])))
	else:
		if 'max' in options['plottingOptions']['colorBarLimits'].lower():
			vmaxCache = np.nanmax(allSkyIm)
			vminCache = np.nanmin(allSkyIm)

			print('Set the overall imaging limits (by max/min) to be between {0} and {1}'.format(vmaxCache, vminCache))
		elif 'percent' in options['plottingOptions']['colorBarLimits'].lower():
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

	if options['plottingOptions']['generateVideo'] and allSkyIm.shape[-1] > 20:
		dateTime, rcuMode, __, frequency, polarity, figNum = labelOptions
		filePrefix = dateTime[:2] # Century of observation as a prefix
		fileSuffix = "mode{0}{1}_{2}_{3}MHz.png".format(rcuMode, polarity, options['imagingOptions']['method'].replace('/', '-'), int(frequency/1e6))

		print("Exporting frames to video at " + "./{0}{1}.mpg".format(dateTime, fileSuffix))
		subprocess.call(['ffmpeg', '-y',  '-framerate',  str(options['plottingOptions']['videoFramerate']), '-pattern_type', 'glob',  '-i',  '{0}*{1}.png'.format(filePrefix, fileSuffix), '-r', str(max(options['plottingOptions']['videoFramerate'], 30)), "{0}{1}.mp4".format(dateTime, fileSuffix)])


def ftAllSkyImage(allSkyImage, options, labelOptions, stationLocation, lVec, mVec):
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

	knownSources, referenceObject = getSkyObjects(options, obsTime, telescopeLoc)

	altAzRef = astropy.coordinates.AltAz(obstime = obsTime, location = telescopeLoc)
	altAzObjects = [skyObj.transform_to(altAzRef) for skyObj in referenceObject]

	gs = mpl.gridspec.GridSpec(1, 2, width_ratios = [18, 1])
	gs2 = mpl.gridspec.GridSpec(1, 2, width_ratios = [18, 1])

	fig = plt.figure(figNum, figsize = (figureShape[0], figureShape[1]))
	fig.patch.set_facecolor(backgroundColor)
	plt.suptitle( 'LOFAR mode {0}{1} all sky plot at {2}MHz (sb{3}) for {4}\n'.format(rcuMode, polarity, round(frequency / 1e6, 2), subband, obsSite), fontsize = int(28 * fontSizeFactor), color=foregroundColor )#, va = 'top') 
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

	plotStatus = [[True, __plotSkyObject(axImage, skyObj, lVec.size, skyObjColor, knownSources[idx], fontSizeFactor = fontSizeFactor)] if (skyObj.alt.deg > 20.) or (options['rfiMode'] and skyObj.alt.deg > -10.) else [False, __plotSkyObject(axImage, skyObj, lVec.size, skyObjColor, knownSources[idx], offset = True, fontSizeFactor = fontSizeFactor)]  for idx, skyObj in enumerate(altAzObjects)]
	
	plt.rcParams["text.color"] = foregroundColor
	legend = axImage.legend( loc = 8, bbox_to_anchor = ( 0.5, -0.128 ), ncol = 4, framealpha = 0.0, fontsize = int(14 * fontSizeFactor), title = str(obsTime)[:-4])
	legend.get_title().set_fontsize(str(int(22 * fontSizeFactor)))

	for idx, skyText in enumerate(legend.get_texts()):
		if not plotStatus[idx][0]:
			plt.setp(skyText, color = 'red')
		elif options['rfiMode'] and altAzObjects[idx].alt.deg > -10. and altAzObjects[idx].alt.deg < 20.:
			plt.setp(skyText, color = 'orange')

	radii = []

	for radius in range(0, 90, 15): # r grid at 15 degree intervals
		radii.append(180 * np.cos(radius * np.pi/180)) # plot the radii so as to display as an orthographic grid
		pltObj.set_rgrids(radii)
	if radialLabelAngle: # you would not want to put y ticks on 0 anyhow as it would be messy
		yLabel = [ '', '15' + u'\xb0', '30' + u'\xb0', '45' + u'\xb0', '60' + u'\xb0', '75' + u'\xb0' ]
		pltObj.set_yticklabels(yLabel, color = skyObjColor)
		pltObj.set_rlabel_position(radialLabelAngle)
	else:
		yLabel = []
		pltObj.set_yticklabels(yLabel)

	thetaticks = np.arange(0, 360, 45)
	pltObj.set_thetagrids(thetaticks, weight = 'bold', color = skyObjColor, fontsize = int(18 * fontSizeFactor))
	pltObj.tick_params('x', pad = -35, rotation = 'auto')

	pltObj.grid(False, 'both', color = skyObjColor, linewidth = gridThickness)
	pltObj.patch.set(alpha = 0.0)
	plt.sca(axImage)

	if not os.path.exists(outputFolder):
		os.makedirs(outputFolder)

	plotFilename = "{7}{0}_{1}_sb{2}_mode{3}{4}_{5}_{6}MHz.png".format(dateTime, obsSite, subband, rcuMode, polarity, options['imagingOptions']['method'].replace('/', '-'), int(frequency/1e6), outputFolder)
	plotFilename = plotFilename.replace(' ', '_').replace(':', '')
	print("Saving output to {0}".format(plotFilename))
	
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

	fig.savefig(plotFilename, facecolor=fig.get_facecolor(), edgecolor='none')

	plt.close(figNum)
	return plotFilename


def __plotSkyObject(axIm, skyObj, pixels, skyObjColor, sourceName, offset = False, fontSizeFactor = 1.):
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



def swhtPlot(hpMap, options, labelOptions, plottingFunc = healpy.newvisufunc.mollview, legendPlot = None, zenithArr = None):
	if options['plottingOptions']['displayFigure']:
		plt.switch_backend('TkAgg')
	else:
		plt.switch_backend('Agg')

	hpMap = hpMap.filled(np.nan)

	swhtFullSkyImage(hpMap, options, labelOptions, plottingFunc, legendPlot, zenithArr)

	plotFileName = None

	return plotFileName


def swhtFullSkyImage(hpMap, options, labelOptions, plottingFunc, legendPlot = None, zenithArr = None):
	plotOptions = options['plottingOptions']
	
	logPlot, skyObjColor, backgroundColor, foregroundColor, colorBar, obsSite, outputFolder, maxPercentile, minPercentile, figureShape, fontSizeFactor, gridThickness = plotOptions['logPlot'], plotOptions['skyObjColor'], plotOptions['backgroundColor'], plotOptions['foregroundColor'], plotOptions['colorBar'], options['stationID'], plotOptions['outputFolder'], plotOptions['maxPercentile'], plotOptions['minPercentile'], plotOptions['figureShape'], plotOptions['fontSizeFactor'], plotOptions['gridThickness']
	dateTime, rcuMode, subband, frequency, polarity, figNum = labelOptions

	if fontSizeFactor is None:
		fontSizeFactor = figureShape[0] / 18.
	dateTime = ['', '']
	if logPlot:
		hpMapNeg = hpMap < 0.
		hpMap[hpMapNeg] = 1.

		#median = np.nanpercentile(hpMap, 80)

		#hpMap -= median

		hpMap = np.log2(hpMap)

		minV = np.nanpercentile(hpMap[hpMap > 10.], minPercentile)
		maxV = np.nanpercentile(hpMap[hpMap > 10.], maxPercentile)

	cmapVar = plt.get_cmap('jet')
	cmapVar.set_under(backgroundColor)
	cmapVar.set_bad(backgroundColor)

	#figRatio = [int(figureShape[0] * 0.8), int(figureShape[1] *0.2)]
	#gs = mpl.gridspec.GridSpec(2, 1, height_ratios = figRatio)
	#gs = mpl.gridspec.GridSpec(1, 1)

	plt.rcParams["text.color"] = foregroundColor
	plt.rcParams["axes.labelcolor"] = foregroundColor
	plt.rcParams["xtick.color"] =  foregroundColor
	plt.rcParams["ytick.color"] = foregroundColor
	plt.rcParams['axes.edgecolor'] = foregroundColor

	#axImage = figObj.add_subplot(gs[0], label = 'ax_image')
	#axImage.axis('off')

	#plt.sca(axImage) # Needed to force healpy to plot onto our figure
	graticuleBool = plotOptions['graticule']
	titleStr = 'LOFAR mode {0}{1} Full Sky Plot at {2}MHz (sb{3}) for {4} between {5} and {6}\n'.format(rcuMode, polarity, round(frequency / 1e6, 2), subband, obsSite, dateTime[0], dateTime[1])
	pltArr = plottingFunc(hpMap, cmap = cmapVar, cbar = colorBar, min = minV, max = maxV, graticule = graticuleBool, shading = 'gouraud')
	
	if graticuleBool:
		plt.gca().set_longitude_grid(30.)
		plt.gca().set_latitude_grid(30.)

		plt.gca().xaxis.set_tick_params(label1On=False)
		plt.gca().yaxis.set_tick_params(label1On=False)
	plt.title(titleStr)

	figObj = plt.gcf()
	figObj.set_size_inches(figureShape[0], figureShape[1])
	plt.subplots_adjust(top = 1., hspace = 0.)
	figObj.patch.set_facecolor(backgroundColor)
	plt.gca().set_facecolor(backgroundColor)

	# If we are plotting zeniths, add a buffer to plot status to account for it.
	plotStatus = []

	if legendPlot:
		if zenithArr:
			galCoordLon = lonToContinuousGalacticLon(np.array([skyCoord.l.rad for skyCoord in zenithArr]), True)
			galCoordLat = np.array([skyCoord.b.rad for skyCoord in zenithArr])
			lineVar = plt.plot(galCoordLon, galCoordLat, linestyle = '-', color = 'w', label = None)
			lineVar = plt.scatter(galCoordLon, galCoordLat, marker = 'D', c = 'w', edgecolors = 'lime', label = 'Zenith Pointings')

			plotStatus = plotStatus + [[False, True]]
		obsTime, telescopeLoc = legendPlot
		knownSources, referenceObject = getSkyObjects(options, obsTime, telescopeLoc)
		# Add an extra true if we have a zenithArr passed.
		print(referenceObject)
		print(referenceObject[0])
		plotStatus = plotStatus + [__plotFullSkyObject(hpMap, skyObj, skyObjColor, knownSources[idx], fontSizeFactor = fontSizeFactor) for idx, skyObj in enumerate(referenceObject)]

		legend = plt.gca().legend(loc = 8, bbox_to_anchor = ( 0.5, -0.4 ), ncol = 4, title = 'l = [-180, 180], b = [-90, 90]', framealpha = 0.0, fontsize = int(17 * fontSizeFactor))
		
		flatten = lambda listofLists: [subList for subListArr in listofLists for subList in subListArr]
		plotStatus = flatten(plotStatus)
		for idx, skyText in enumerate(legend.get_texts()):
			if plotStatus[idx]:
				plt.setp(skyText, color = 'red')
	

	figObj.savefig('./swht_debug.png', facecolor=figObj.get_facecolor(), edgecolor='none')
	if plotOptions['displayFigure']:
		plt.show()
	plotFileName = None, outputFolder

	return plotFileName

def __plotFullSkyObject(hpMap, skyObj, skyObjColor, sourceName, fontSizeFactor = 1., offset = False):
	alpha = 0.
	hpMap[hpMap < 1.] = np.nan
	if isinstance(skyObj, list):
		print(skyObj[0])
		lDegOrig = [skyObjEle.l.deg for skyObjEle in skyObj]
		lDeg = lonToHealpyLon(lDegOrig, False)
		bDeg = [skyObjEle.b.deg for skyObjEle in skyObj]
		print(lDeg, bDeg)
		if np.any(~np.isnan(hpMap[healpy.ang2pix(64, lDeg, bDeg, lonlat = True)])) and not offset:
			alpha = 1.
		else:
			offset = True

		lDeg = lonToContinuousGalacticLon(lDegOrig, False)

		lRad = lonToContinuousGalacticLon([skyObjEle.l.rad for skyObjEle in skyObj], True)
		bRad = [skyObjEle.b.rad for skyObjEle in skyObj]
		skyPlt = plt.plot(lRad, bRad, linestyle = '-', color = skyObjColor, label = None, alpha = alpha)
		skyPlt_marker = plt.scatter(lRad, bRad, marker = 'X', c = 'w', s = int(50. * fontSizeFactor), edgecolors = skyObjColor, label = u"{0} - l={1}\xb0, b={2}\xb0".format(sourceName, round(np.mean(lDeg), 1), round(np.mean(bDeg), 1)), alpha = alpha)

	else:
		lDegOrig = skyObj.l.deg
		lDeg, bDeg = lonToHealpyLon(lDegOrig, False), skyObj.b.deg
		if not np.isnan(hpMap[healpy.ang2pix(64, lDeg, bDeg, lonlat = True)]) and not offset:
			alpha = 1.
		else:
			offset = True

		lDeg = lonToContinuousGalacticLon(lDegOrig, False)

		lRad, bRad = lonToContinuousGalacticLon(skyObj.l.rad, True), skyObj.b.rad
		skyPlt = plt.scatter(lRad, bRad, marker = 'D', s = int(50 * fontSizeFactor), c = skyObjColor, label = u"{0} - l={1}\xb0, b={2}\xb0".format(sourceName, round(lDeg, 1), round(bDeg, 1)), alpha = alpha)

	if not offset:
		textObj = plt.gca().annotate(sourceName, xy = (np.mean(lRad),np.mean(bRad)), xytext = (np.mean(lRad),np.mean(bRad)), color = skyObjColor, fontsize = int(18 * fontSizeFactor))
		textObj.set_path_effects([mplPe.withStroke(linewidth = 5, foreground = 'w')])

	if isinstance(skyObj, list):
		return [offset, offset]

	return [offset]

def mapSWHTlValues(lVal, deg = False):
	if deg:
		a = 'a'
	return a

def getSkyObjects(options, obsTime, telescopeLoc):
	try:
		swhtBool = 'swht' in options['imagingOptions']['method']
		# Far celestial sources can be cached, bodies in the Solar System cannot.
		knownSources = options['plottingOptions']['interstellarSources']
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
	if not swhtBool:
		return astropy.coordinates.get_body(sourceName, obsTime, telescopeLoc)
	
	return [astropy.coordinates.get_body(sourceName, obsTick, telescopeLoc).transform_to(astropy.coordinates.Galactic) for obsTick in obsTime]


@cache
def cachedSkyCoords(name, swhtBool = False):
	"""Summary
	
	Args:
		name (TYPE): Description
	
	Returns:
		TYPE: Description
	"""
	if not swhtBool:
		return astropy.coordinates.SkyCoord.from_name(name)

	return astropy.coordinates.SkyCoord.from_name(name).transform_to(astropy.coordinates.Galactic)
