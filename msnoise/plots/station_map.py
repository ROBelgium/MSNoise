"""
This plots a station map using PyGMT.

.. include:: /clickhelp/msnoise-plot-station_map.rst

Example:

``msnoise plot station_map`` :

.. image:: ../.static/station_map.png

It will also generate a HTML file showing the stations on the Leaflet Mapping
Service:

.. raw:: html

    <iframe src="_static/station_map.html" width=800 height=400></iframe>
"""

import datetime
import logging

import numpy as np

from ..core.db import connect
from ..core.stations import get_stations


def _outfile_name(outfile):
    if outfile and outfile.startswith("?"):
        now = datetime.datetime.now()
        now = now.strftime("station map on %Y-%m-%d %H.%M.%S")
        return outfile.replace("?", now)
    return outfile


def _plot_leaflet(stations, coords, outfile):
    import folium

    sta_map = folium.Map(
        location=[np.mean(coords[:, 0]), np.mean(coords[:, 1])],
        zoom_start=3,
        tiles="OpenStreetMap",
    )
    folium.RegularPolygonMarker(
        location=[np.mean(coords[:, 0]), np.mean(coords[:, 1])]
    ).add_to(sta_map)

    for sta in stations:
        folium.RegularPolygonMarker(
            location=[sta.Y, sta.X],
            popup="%s_%s" % (sta.net, sta.sta),
            fill_color="red",
            number_of_sides=3,
            radius=12,
        ).add_to(sta_map)

    sta_map.add_child(folium.LatLngPopup())

    if outfile:
        logging.debug("output to: %s.html" % outfile)
        sta_map.save("%s.html" % outfile)


def _plot_pygmt(stations, coords, show=True, outfile=None):
    try:
        import pygmt
    except ImportError:
        raise ImportError(
            "PyGMT is required to plot station maps. "
            "Please install pygmt to use `msnoise plot station_map`."
        )

    lats = coords[:, 0]
    lons = coords[:, 1]
    labels = ["%s_%s" % (sta.net, sta.sta) for sta in stations]

    bufferlat = (np.amax(lats) - np.amin(lats)) + 0.1
    bufferlon = (np.amax(lons) - np.amin(lons)) + 0.1
    region = [
        np.amin(lons) - bufferlon,
        np.amax(lons) + bufferlon,
        np.amin(lats) - bufferlat,
        np.amax(lats) + bufferlat,
    ]

    fig = pygmt.Figure()
    fig.coast(
        region=region,
        projection="M15c",
        resolution="i",
        land="lightgray",
        water="lightblue",
        shorelines="0.5p,black",
        borders=["1/0.5p,black"],
        frame=["a", "+tStation map"],
    )
    fig.plot(x=lons, y=lats, style="i0.4c", fill="red", pen="0.25p,black")
    fig.text(
        x=lons,
        y=lats,
        text=labels,
        font="9p,Helvetica,black",
        justify="CB",
        offset="0/0.25c",
        fill="white",
        pen="0.25p,black",
        transparency=30,
    )

    if outfile:
        logging.info("output to: %s" % outfile)
        fig.savefig(outfile)

    if show:
        fig.show()


def main(show=True, outfile=None):
    db = connect()
    stations = get_stations(db, all=False)
    coords = [(sta.Y, sta.X) for sta in stations]
    coords = np.array(coords)

    if not len(coords):
        logging.warning("No station found. Nothing to plot.")
        return

    outfile = _outfile_name(outfile)

    _plot_leaflet(stations, coords, outfile)
    _plot_pygmt(stations, coords, show=show, outfile=outfile)


if __name__ == "__main__":
    main()
