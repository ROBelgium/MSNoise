"""MSNoise is ...

Usage:
~~~~~~

.. code-block:: sh

    $ msnoise plot interferogram

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
    
    coords = [(sta.X, sta.Y) for sta in stations]
    print coords
    
    ax.coastlines()

    plt.show()
            
        
                            

if __name__ == "__main__":
    main()