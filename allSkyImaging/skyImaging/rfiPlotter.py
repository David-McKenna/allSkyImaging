"""Summary
"""
import numpy as np

from skyPlotter import informationArr
global informationArr

def update_annot( xdata, ydata, pixels, annot, rawdata, **kwargs):
	y, x = pol2cart( ydata/180, xdata, pixels )
	annot.xy = ( xdata, ydata )
	# Inconsistent wrapping; plot the right variable.
	if xdata < 0:
		xdata += 2 * np.pi
	text = 'Az=' + str( round( xdata * 180 / np.pi, 1 ) )+  ', El=' + str( round( np.arccos( ydata/180 ) * 180/np.pi, 1) ) + u'\xb0' + '\nInt.=' + '{:.3E}'.format((rawdata[int(y),int(x)]))
	annot.set_text( text )
	annot.get_bbox_patch().set_alpha( 0.66 )
	annot.set_color('black')

def onclick(event, annot, pltObj, pixels, rawdata, **kwargs): # Callbacks only work on the last plot made
	vis = annot.get_visible()
	if event.inaxes == pltObj:
		if not vis:
			update_annot(event.xdata, event.ydata, pixels, annot, rawdata)
			annot.set_visible( True )
			event.canvas.draw()
		else:
			annot.set_visible( False )
			event.canvas.draw()
			#fig.canvas.draw_idle()
			
def hover(event, pltObj, pixels, rawdata, axColorBar, cbCursor, **kwargs): # Callbacks only work on the last plot made
	if event.xdata:# if you have 2 plots on screen, for some reason the hover event is triggered if you hover over the first one giving errors
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

def onaxesleave(event, pltObj, axColorBar, cbCursor, **kwargs): # Callbacks only work on the last plot made
	cleanCb(axColorBar)
	cbCursor = axColorBar.plot([0, 1],[0, 0], 'k-')
	event.canvas.draw()

	annot = pltObj.annotate( "", xy = ( 0, 0 ), xytext = ( 15, 15 ), textcoords = "offset points", bbox = dict( boxstyle = "round", fc = "b" ), arrowprops = dict( arrowstyle = "->" ) )
	annot.set_visible( False )

	global informationArr
	informationArr['annot'] = annot
	informationArr['cbCursor'] = cbCursor

def cleanCb(axColorBar):
	for line in axColorBar.axes.lines:
		line.remove()

	return axColorBar

def pol2cart( rho, phi, pixels ):
    x = rho * np.cos( phi )
    y = rho * np.sin( phi )
    x=( pixels/2 )-( pixels/2 )*x
    y=( pixels/2 )-( pixels/2 )*y

    return( x, y )
