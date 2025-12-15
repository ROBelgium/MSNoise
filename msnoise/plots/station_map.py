import traceback
import logging
import datetime
import numpy as np
import os
import matplotlib.pyplot as plt

from ..api import *

def plot_basemap(show=True, outfile=None, stations=None):
    import folium
    try:
        from mpl_toolkits.basemap import Basemap
    except ImportError:
        print("Error: Basemap is not installed. Please install it or use --pygmt")
        return

    coords = [(sta.Y, sta.X) for sta in stations]
    coords = np.array(coords)

    print("Plotting with Basemap (Legacy)...")
    if show:
        plt.show()

def plot_pygmt_map(show=False, outfile="?", stations=None):
    try:
        import pygmt
    except ImportError:
        print("Error: PyGMT is not installed.")
        return

    lons = [sta.X for sta in stations]
    lats = [sta.Y for sta in stations]
    labels = [f"{sta.net}.{sta.sta}" for sta in stations]

    buffer_deg_lon = 0.25
    buffer_deg_lat = 0.15
    region = [
        np.min(lons) - buffer_deg_lon,
        np.max(lons) + buffer_deg_lon,
        np.min(lats) - buffer_deg_lat,
        np.max(lats) + buffer_deg_lat
    ]

    print("Generating map using PyGMT...")

    with pygmt.config(MAP_FRAME_TYPE="plain", MAP_FRAME_PEN="1p,black",
                      FORMAT_GEO_MAP="ddd.xx", MAP_DEGREE_SYMBOL="none"):
        fig = pygmt.Figure()
        projection = "M15c"
        fig.coast(region=region, projection=projection, resolution='f', 
                  land='lightgray', water='lightblue',
                  shorelines='thinnest,black', borders=["1/0.5p,black"],
                  frame=["a", "+tStation Map"])
        
        fig.plot(x=lons, y=lats, style='i0.4c', fill='red', pen='faint,black')
        fig.text(x=lons, y=lats, text=labels, font='9p,Helvetica-Bold,black',
                 justify='CB', offset='0/0.25c', fill='white',
                 pen='thinner,black', transparency=50)

        if outfile is True or outfile == "?":
            now = datetime.datetime.now()
            now = now.strftime('%Y-%m-%d-%H%M%S')
            filename = f"station_map_pygmt_{now}.png"
        elif outfile:
            filename = outfile
        else:
            filename = "station_map_pygmt_default.png"
        
        logging.info(f"Output image to: {filename}")
        fig.savefig(filename)
        print(f"Map saved to: {os.path.join(os.getcwd(), filename)}")

def main(show=True, outfile=None, backend="basemap"):
    db = connect()
    stations = get_stations(db, all=False)
    if not stations:
        print("No stations found.")
        return

    if backend == "pygmt":
        plot_pygmt_map(show=show, outfile=outfile, stations=stations)
    else:
        plot_basemap(show=show, outfile=outfile, stations=stations)

if __name__ == "__main__":
    main()