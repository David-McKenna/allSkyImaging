#allSkyImaging#

This is a module based around imaging LOFAR XST files through the use of either a Discrete Fourier Transform (DFT) or Spherical Wave Harmonic Transform (SWHT, Carrozi 2015). This allows for the creation of all sky (DFT) or full sky through combinations of integrations at several time intervals throughout the day (SWHT).

The only requirement for XSTs is their naming scheme: we expect a format of `date_time_sb<num>_xst.dat`.
##Usage##

In most cases, you'll be purely interested in using functions found within the `allSkyImager.py` script, where we provide two main functions of use, `image` and `monitor`. Both of these can be provided an options dictionary (described below), but work on different principles: `image` will take a folder and process any XST files within it, while `monitor` is pointed at a folder and monitors for changes. Once a threshold is reached, `monitor` will cann the same functions as image to process the data and repeat until the process is stopped.

Sample Commands:
```
allSkyImager.image('./sampleMode3/', options = optDictionary)

allSkyImager.monitor('./sampleMode3', options = optDictionary, processedOutput = '/processedSubfolder/', fileThreshold = 3, checkEvery = 10.)
```
****
* **Options (dict):** Controls the imaging and plotting processes, described in detail below
* **processedOutput (str):** Controls the name of the subfolder that processed XST files will be moved to by the monitor script
* **fileThreshold (int):** Sets the number of files required to be in a folder before monitor will start imaging the files
* **checkEvery (float/int):** Sets a delay between checks for files (in seconds)

##Options##

The options dictionary is a list of requested differences from the dictionary that can be found in `dataTools.defaultDict.py`. We modify the default dictionary by not providing a modified version, but a dictionary of requested changes. 

Since the default dictionary contains sub-dictionaries, we can change the sub-dictionary parameters by referencing the sub dictionary name, such as 'imagingOptions' and a value, such as 'logPlot' by using the key 'imagingOptions:logPlot'.

Changes can be made by referencing the key and the new value in a new dictionary. For example, if we wanted to change the RCU Mode to 5 and imaging mode to the SWHT, we would use the following dictionary,

```
optDictionary = {
		'rcuMode': 5,
		'imagingOptions:method': 'swht'
				}

```

This method can be used for any option within the dictionary, apart from the 'remoteURLLocations' subdictonary, which is used to download configurations from webpages online to perform imaging as needed.

The below options are described in the syntax

* **parameter (key for dictionary)**
    - **Default (type):** Sutable options
    - Description / Comment

###Overall Options###
These value do not need a subdictionary prefix to edit.


* **stationID **
    - **IE613 (str):** IE613, SE607, CS001...
    - Station code


* **rcuMode **
    - **3 (int):** 1-7
    - Setting RCU Observing Mode

* **activationPattern **
    - **Generic_Int_201512 (str):** Generic_Int_201512, Effelsberg_elements_20091110, Generic_International_Station_20091110, Generic_Core_201512, Generic_Remote_201512
    - HBA Activation Pattern if rcuMode > 4

* **h5GroupName**
    - **allSkyGroupName (str):** Anything
    - Name of the group of the imported / output data in the h5 file. If left at the default, it will me modified to include the first timestamp.

* **multiprocessing**
    - **True (bool):** True or False
    - Determines whether or not to enable multiprocessing. If enabled, we will spin up (n-1) processes to process the data, where n is the number of available CPU threads.

* **rfiMode**
    - **False (bool):** True or False
    - If enabled, this will force the imaging mode to DFT and enable some extra features to help analyse RFI presence in an XST observation. These include: making the intensity under the mouse on the color bar, giving the ra/dec of the location you click on, and the raw intentsity at the clicked location.

###File Locations###
All of these options are prefixed by 'fileLocations' as they are lcoated in a subdictionary. They control where different configuration and output files are saved.

* **fileLocations:antennaField**
    - **./config/{0}-AntennaField.conf (str):** folder + filename
    - Location for the LOFAR Blitz-style antenna field config. {0} is replaced with the station ID.

* **fileLocations:hbaDeltas**
    - **./config/{0}-iHBADeltas.conf (str):** folder + filename
    - Location for the LOFAR Blitz-style HBA antenna offsets config. {0} is replaced with the station ID.

* **fileLocations:lbaRotation**
    - **./config/stationrotations.txt (str):** folder + filename
    - Location for the lofarimaging text file containing the angular offsets of each LOFAR station.

* **fileLocations:outputH5Location**
    - **./{0}-mode{1}.h5 (str):** folder + filename
    - Location to put the output h5 data file. {0} is replaced with the station ID, {1} is replaced with the RCU mode.

* **fileLocations:calibrationLocation**
    - **./config/{0}-mode{1}.dat (str):** folder + filename
    - Location to save the calibration data file. {0} is replaced with the station ID, {1} is replaced with the RCU mode.

###Imaging Options
All of these options are prefixed by 'imagingOptions' as they are lcoated in a subdictionary. These options control the methodology for processing the antenna correlations.

* **imagingOptions:method**
    - **dft ():** dft, swht, (PLACEHOLDER, BROKEN) fft/kernelStr/opt-paramter(float)
    - Method for imaging the data: DFT implies the creation of an all sky map over several different frames, SWHT implies creating a full sky map out of whatever correlations are provided.

* **imagingOptions:calibrateData**
    - **True (bool):** True or False
    - If true, apply the calibration file calibration to the observation correlations.

* **imagingOptions:baselineLimits**
    - **[None, None] (list):** [float, float]
    - Maximum displacement between two antennas (in meters) for them to be processed. Otherwise, the correlations are set to 0.

* **imagingOptions:pixelCount**
    - **[255, 255] (list):** [int, int]
    - Size of the output image array. The output figure will differ in size, but this is the size of the raw data array that is created in the DFT method.

* **imagingOptions:fieldOfView**
    - **np.pi (float (radians)):** 0 - np.pi
    - The angular distance away from the zenith to image to. Possible issues: if not masked, we may have some stretching occurring.

* **imagingOptions:maskOutput**
    - **True (bool):** True or False
    - If true, mask any data outside the fieldOfView variable.

* **imagingOptions:correlationTypes**
    - **['I', 'YY'] (list):** Values of I, Q, U, V, XX, XY, YX, YY
    - Types of images generated from the correlation type. Can provides: Stokes Vectors or antenna correlations

* **imagingOptions:ftSubtractBackground**
    - **True (bool):** True or False
    - Whether or not to remove the minimum value in the output image from all values in the dataset.

* **imagingOptions:swhtlMax**
    - **32 (int):** Any value > 0
    - Controls the maximum l value that the SWHT processes up to. For every l value, we generate a m value l > m > 0

###Plotting Options###
All of these options are prefixed by 'plottingOptions' as they are lcoated in a subdictionary. These options control the parameters used for plotting the processed data.

* **plottingOptions:outputFolder**
    - **None (str):** : folder
    - Where output figures are saved. If kept as none, it will be placed in a subdirectory of the correlation folder all 'allSkyOutput'.

* **plottingOptions:plotImages**
    - **True (bool):** Tur eor False
    - If true, produce an output plot of the processed data.

* **plottingOptions:displayImages**
    - **False (bool):** True or False
    - If true, display images as they are plotted.

* **plottingOptions:generateVideo**
    - **True (bool):** True or False
    - If false, don't generate animations of the output. If true (DFT): if more than 10 frames are generated, a video will be generated for each correlation type requested, (SWHT) an animation of the clean / overplotted figures will be generated.

* **plottingOptions:videoFrameRate**
    - **8 (int):** int > 0
    - (DFT): Output video has a framerate of n frames per second. (SWHT): Output animation will change image every n seconds.

* **plottingOptions:figureShape**
    - **[18, 14] (list):** [int, int]
    - Size of the output figure (inches / units of 100 pixels).

* **plottingOptions:plotSkyObjects**
    - **True (bool):** True or False
    - If true, plot the locations of objects supplied in the next two parameters in the plots.

* **plottingOptions:interstellarSources**
    - **['Polaris', 'Cas A', 'Cyg A', 'Sgr A', 'Tau A', 'Vir A', 'Cen A', 'Vela'] (list):** [source names, comma separated]
    - Any interstellar sources (outside of the solar system) that you want to be plotted on the figure.

* **plottingOptions:solarSystemSources**
    - **['Sun', 'Jupiter', 'Moon', 'Uranus', 'Neptune'] (list):** [sources names, comma separated]
    - Any sources within the solar system that you want to be plotted on the figure. In the SWHT case, we will plot a point for the location of each source at every time step and put the name at the mean location.

* **plottingOptions:colorBarMemory**
    - **16 (int/str):** int, maxmin, percentile
    - Int: The last n frame's limits (by percentile options) will be averaged to set the colour bar limits. This can be considered to provide a sliding exposure window.
    - String containing 'max': The maximum / minimum value of all frames is the maximum / minimum for every frames plot.
    - String containing 'percentile': The percentiles (sampled from the below options) are samples for the entire array, and set for the entire array.

* **plottingOptions:maxPercentile**
    - **99 (int):** 0-100, > minPercentile
    - The sampled percentile for the upper value on the colour bar.

* **plottingOptions:minPercentile**
    - **33 (int):** 0-100 < maxPercentile
    - The sampled percentile for the lower value on the colour bar.

* **plottingOptions:logPlot**
    - **True (bool):** True or False
    - (DFT): Take the log10 of the data and plot it. (SWHT): Get the log2 of all values greater than 0.

* **plottingOptions:skyObjColor**
    - **black (str):** MatPlotLib Colour
    - Colour of the scatter point / text colour for sky objects.

* **plottingOptions:gridThickness**
    - **0.5 (float):** float (pixels)
    - (DFT Only): Thickness of the angular grid in pixels.

* **plottingOptions:backgroundColor**
    - **black (str):** MatPlotLib Colour
    - Background colour for the figures.

* **plottingOptions:forgroundColor**
    - **white (str):** MatPlotLib Colour
    - Foreground colours for most overlayed items.

* **plottingOptions:radialLabelAngle**
    - **0 (int):** 0-360
    - (DFT Only): Angular offset for graticule angle labels.

* **plottingOptions:fontSizeFactor**
    - **None (float):** None or float > 0
    - The scaling factor for all fonts on the figure. If left as None, this should handle itself to scale based on any change from the default Figure Shape.

* **plottingOptions:colorBar**
    - **True (bool):** True or False
    - Whether or not to provide a color bar with the plot.

* **plottingOptions:graticule**
    - **True (bool):** True or False
    - Whether or not to overlap the plot with the angular grid.

* **plottingOptions:swhtZenithPointings**
    - **True (bool):** True or False
    - If true, we will plot the locations of the zeniths on the full sky map.

##Sample Default Dictionaries##
Comments: 

* Anything above mode 5 should be plotted in log, not much can be seen in the raw data.
* SWHT Plots give clearer stucture in log in log images, but look prettier non-logged. 
* imagingOptions:ftSubtractBackground is useful for LBA images, but not HBA images.
```
optDict = {
	'method': 'dft',
	'correlationTypes': ['I', 'Q', 'U', 'V', 'YY'],
	
	'plottingOptions:videoFramerate': 30,
	'plottingOptions:minPercentile': 10,
	'plottingOptions:logPlot: True
}


optDict = {
	'method': 'swht',
	'correlationTypes': ['I', 'YY'],
	
	'plottingOptions:videoFramerate': 4,
	'plottingOptions:minPercentile': 1,
	'plottingOptions:logPlot': True
}


optDict = {
	'method': 'swht',
	'correlationTypes': ['I', 'YY'],
	
	'plottingOptions:videoFramerate': 4,
	'plottingOptions:minPercentile': 20,
	'plottingOptions:logPlot': False	
}

optDict = {
	'method': 'swht',
	'rcuMode': 5,
	'correlationTypes': ['I', 'YY'],

	'imagingOptions:ftSubtractBackground': False,
	
	'plottingOptions:videoFramerate': 4,
	'plottingOptions:minPercentile': 20,
	'plottingOptions:logPlot': False	
}

optDict = {
	'method': 'swht',
	'rcuMode': 5,
	'correlationTypes': ['I', 'YY'],
	
	'plottingOptions:videoFramerate': 4,
	'plottingOptions:minPercentile': 20,
	'plottingOptions:logPlot': False	
}
```
