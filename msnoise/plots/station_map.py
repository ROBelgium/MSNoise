"""
This plots a very raw station map (needs improvement). This plot requires
cartopy !


.. include:: clickhelp/msnoise-plot-station_map.rst


Example:

``msnoise plot station_map`` :

.. image:: .static/station_map.png


It will also generate a HTML file showing the stations on the Leaflet Mapping
Service:

.. raw:: html

    <iframe src="_static/station_map.html" width=800 height=400></iframe>


.. versionadded:: 1.4 | Thanks to A. Mordret!

"""

import traceback
import folium
import os
import numpy as np
import sys

from ..api import *

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


def main(show=True, outfile=None):
    db = connect()
    stations  = get_stations(db, all=False)
    coords = [(sta.Y, sta.X) for sta in stations]
    coords = np.array(coords)

    sta_map = folium.Map(location=[np.mean(coords[:,0]), np.mean(coords[:,1])], zoom_start=3, tiles='OpenStreetMap')
    folium.RegularPolygonMarker(location=[np.mean(coords[:,1]), np.mean(coords[:,0])]).add_to(sta_map)
    for sta in stations:
        folium.RegularPolygonMarker(location=[sta.Y, sta.X],popup="%s_%s" % (sta.net, sta.sta),fill_color='red',
        number_of_sides=3,radius=12).add_to(sta_map)

    sta_map.add_child(folium.LatLngPopup())
    if outfile:
        tmp = outfile
        if outfile.startswith("?"):
            now = datetime.datetime.now()
            now = now.strftime('station map on %Y-%m-%d %H.%M.%S')
            tmp = outfile.replace('?', now)
        print( "output to:", tmp)
        sta_map.save('%s.html'%tmp)

    # plot topography/bathymetry as an image.
    bufferlat=(np.amax(coords[:,0])-np.amin(coords[:,0]))+.1
    bufferlon=(np.amax(coords[:,1])-np.amin(coords[:,1]))+.1
    m = Basemap(projection='mill',llcrnrlat=np.amin(coords[:,0])-bufferlat,urcrnrlat=np.amax(coords[:,0])+bufferlat,\
            llcrnrlon=np.amin(coords[:,1])-bufferlon,urcrnrlon=np.amax(coords[:,1])+bufferlon,resolution='i')

    # Draw station coordinates
    x, y = m(coords[:,1],coords[:,0])

    # create new figure, axes instance.
    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    # attach new axes image to existing Basemap instance.
    m.ax = ax
    #im = m.imshow(topodat,cm.GMT_haxby)
    try:
        m.shadedrelief()
    except:
        traceback.print_exc()
    m.scatter(x,y,50,marker='v',color='r')
    for sta in stations:
        xpt, ypt = m(sta.X,sta.Y)
        plt.text(xpt,ypt,"%s_%s" % (sta.net, sta.sta),fontsize=9,
                    ha='center',va='top',color='k',
                    bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.5)))
    # draw coastlines and political boundaries.
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    # draw parallels and meridians.
    # label on left and bottom of map.
    parallels = np.arange(-90,90,10.)
    m.drawparallels(parallels,labels=[1,0,0,1])
    meridians = np.arange(0.,360.,10.)
    m.drawmeridians(meridians,labels=[1,0,0,1])
    # add colorbar
    #cb = m.colorbar(im,"right", size="5%", pad='2%')
    ax.set_title('Station map')


    if outfile:
        if outfile.startswith("?"):
            now = datetime.datetime.now()
            now = now.strftime('station map on %Y-%m-%d %H.%M.%S')
            outfile = outfile.replace('?', now)
        print( "output to:", outfile)
        plt.savefig(outfile)
    if show:
        plt.show()

if __name__ == "__main__":
    main()
