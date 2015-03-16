"""
This plots a very raw station map (needs improvement). This plot requires
cartopy !


.. code-block:: sh

    msnoise plot station_map --help

    Usage: msnoise-script.py plot station_map [OPTIONS]

      Plots the station map (very basic)

    Options:
      -s, --show BOOLEAN  Show interactively?
      --help              Show this message and exit.

Example:

``msnoise plot station_map`` :

.. image:: .static/station_map.png

"""
# plot interferogram
import matplotlib.pyplot as plt

import matplotlib.gridspec as gridspec

import cartopy.crs as ccrs

import os
import numpy as np
import sys

from ..api import *


def main(show=True):
    db = connect()
    stations  = get_stations(db, all=False)
    
    plt.figure()
    ax = plt.subplot(111, projection=ccrs.PlateCarree())
    ax.coastlines()

    coords = [(sta.X, sta.Y) for sta in stations]
    coords = np.array(coords)
    plt.scatter(coords[:,0], coords[:,1], transform=ccrs.Geodetic())

    plt.show()
            
        
                            

if __name__ == "__main__":
    main()