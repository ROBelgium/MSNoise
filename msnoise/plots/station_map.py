"""
This plots a very raw station map (needs improvement). This plot requires
cartopy !


.. include:: clickhelp/msnoise-plot-station_map.rst


Example:

``msnoise plot station_map`` :

.. image:: .static/station_map.png

"""

import matplotlib.pyplot as plt

import matplotlib.gridspec as gridspec

import cartopy.crs as ccrs

import os
import numpy as np
import sys

from ..api import *


def main(show=True, outfile=None):
    db = connect()
    stations  = get_stations(db, all=False)
    
    plt.figure()
    ax = plt.subplot(111, projection=ccrs.PlateCarree())
    ax.coastlines()

    coords = [(sta.X, sta.Y) for sta in stations]
    coords = np.array(coords)
    plt.scatter(coords[:,0], coords[:,1], transform=ccrs.Geodetic())
    if outfile:
        if outfile.startswith("?"):
            now = datetime.datetime.now()
            now = now.strftime('station map on %Y-%m-%d %H.%M.%S')
            outfile = outfile.replace('?', now)
        print "output to:", outfile
        plt.savefig(outfile)
    if show:
        plt.show()

if __name__ == "__main__":
    main()