"""RFI Interative plot handlers.

Originally ported to Python by Joe McCauley, slightly modified for this module. The source contained the following header:

@author: Joe McCauley (joe.mccauley@tcd.ie)
Written for Python 2.7
Based on a translated matlab script originally from ASTRON for processing 
xst data from an international LOFAR station.

"""
import numpy as np

from skyPlotter import informationArr

# Store some variables on a global level to help deal with the fact that these events are isolated.
global informationArr

def updateAnnot( xdata, ydata, pixels, annot, rawdata, **kwargs):
	"""Update the plotted annotation
	"""
	y, x = pol2cart( ydata/180, xdata, pixels )
	annot.xy = ( xdata, ydata )
	# Inconsistent wrapping; plot the right variable.
	if xdata < 0:
		xdata += 2 * np.pi
	text = 'Az=' + str( round( xdata * 180 / np.pi, 1 ) )+  ', El=' + str( round( np.arccos( ydata/180 ) * 180/np.pi, 1) ) + u'\xb0' + '\nInt.=' + '{:.3E}'.format((rawdata[int(y),int(x)]))
	annot.set_text( text )
	annot.get_bbox_patch().set_alpha( 0.66 )
	annot.set_color('black')

def onclick(event, annot, pltObj, pixels, rawdata, **kwargs):
	"""Handle the matplotlib click event
	"""
	vis = annot.get_visible()
	if event.inaxes == pltObj:
		if not vis:
			updateAnnot(event.xdata, event.ydata, pixels, annot, rawdata)
			annot.set_visible( True )
			event.canvas.draw()
		else:
			annot.set_visible( False )
			event.canvas.draw()
			
def hover(event, pltObj, pixels, rawdata, axColorBar, cbCursor, **kwargs):
	"""Handle cursor movement (for the colorbar line)
	"""
	if event.xdata:
	  	if event.inaxes == pltObj:
			y,x = pol2cart( event.ydata / 180, event.xdata, pixels )
			z=rawdata[ int( y ), int( x ) ]
			zline = ( z - np.nanmin( rawdata ) ) / np.nanmax( rawdata-np.nanmin( rawdata ) ) # calculate where to put the z line
			axColorBar = cleanCb(axColorBar)
			cbCursor = axColorBar.plot( [ 0, 1 ], [ zline, zline ], 'w-', linewidth = 4 ) #plot the new one
			event.canvas.draw()
	
			global informationArr
			informationArr['cbCursor'] = cbCursor
			#fig.canvas.draw_idle()

def onaxesleave(event, pltObj, axColorBar, cbCursor, **kwargs):
	"""Handle cursor leaving the plot
	"""
	cleanCb(axColorBar)
	cbCursor = axColorBar.plot([0, 1],[0, 0], 'k-')
	event.canvas.draw()

	annot = pltObj.annotate( "", xy = ( 0, 0 ), xytext = ( 15, 15 ), textcoords = "offset points", bbox = dict( boxstyle = "round", fc = "b" ), arrowprops = dict( arrowstyle = "->" ) )
	annot.set_visible( False )

	global informationArr
	informationArr['annot'] = annot
	informationArr['cbCursor'] = cbCursor

def cleanCb(axColorBar):
	"""Remove previous lines from the color bar
	"""
	for line in axColorBar.axes.lines:
		line.remove()

	return axColorBar

def pol2cart( rho, phi, pixels ):
    """Convert from polar coordinates to cartesian
    """
    x = rho * np.cos( phi )
    y = rho * np.sin( phi )
    x=( pixels/2 )-( pixels/2 )*x
    y=( pixels/2 )-( pixels/2 )*y

    return( x, y )
