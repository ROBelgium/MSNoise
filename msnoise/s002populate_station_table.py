"""
This script is responsible for rapidly scanning the data archive and
identifying the Networks/Stations and insert them in the *stations* table in
the database.

The ``data_folder`` (as defined in the config) is scanned expecting the
``data_structure`` and possible values are defined in *data_structures.py*:

.. code-block:: python
    
    data_structure['SDS'] = "YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"
    data_structure['BUD'] = "NET/STA/STA.NET.LOC.CHAN.YEAR.DAY"
    data_structure['IDDS'] = "YEAR/NET/STA/CHAN.TYPE/DAY/NET.STA.LOC.CHAN.TYPE.YEAR.DAY.HOUR"
    data_structure['PDF'] = "YEAR/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"

If one's data structure is one of those, then the ``data_structure`` configuration
bit needs to be set to the acronym (SDS, BUD, IDDS or PDF).

More info on the recommended SDS ("SeisComP Data Structure") can be found here:
https://www.seiscomp3.org/wiki/doc/applications/slarchive/SDS
For other structures, one has to edit the data_structures.py file and define
the reader in this script.

By default, station coordinates are initialized at 0.

To run this script:

.. code-block:: sh

    $ msnoise populate

Custom data structure & station table population
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If one's data structure is not one of the pre-defined, it can be defined
directly in the ``data_structure`` configuration bit using forward slashes,
e.g.:

``data_structure`` = "NET/STA/YEAR/NET.STA.YEAR.DAY.MSEED"

MSNoise expects to find
a file named ``custom.py`` in the current folder. This python file will contain
a function called ``populate`` wich will accept one argument and return a list
of stations in the format ``NET_STA``:

.. code-block:: python

    def populate(data_folder):
        datalist = sorted(glob.glob(os.path.join(data_folder, "*", "*")))
        stations = []
        for di in datalist:
            tmp = os.path.split(di)
            sta = tmp[1]
            net = os.path.split(tmp[0])[1]
            stations.append("%s_%s" % (net, sta))
        return stations

"""

import glob
import sys
import traceback

from .api import *


def main():
    db = connect()
    print()
    print(">> Populating the Station table")
    print()
    data_folder = get_config(db, 'data_folder')
    data_structure = get_config(db, 'data_structure')
    if data_structure in ["SDS", "IDDS"]:
        datalist = sorted(glob.glob(os.path.join(data_folder, "*", "*", "*")))
        stations = []
        for di in datalist:
            tmp = os.path.split(di)
            sta = tmp[1]
            net = os.path.split(tmp[0])[1]
            stations.append("%s_%s" % (net, sta))
        del datalist
    elif data_structure in ["BUD", ]:
        datalist = sorted(glob.glob(os.path.join(data_folder, "*", "*",)))
        stations = []
        for di in datalist:
            tmp = os.path.split(di)
            sta = tmp[1]
            net = os.path.split(tmp[0])[1]
            stations.append("%s_%s" % (net, sta))
        del datalist
    elif data_structure in ["PDF", ]:
        datalist = sorted(glob.glob(os.path.join(data_folder, "*", "*",)))
        stations = []
        for di in datalist:
            tmp = os.path.split(di)
            sta = tmp[1]
            net = get_config(db, 'network')
            stations.append("%s_%s" % (net, sta))
        del datalist
    else:
        print("Can't parse the archive for format %s !" % data_structure)
        print("trying to import local parser (should return a station list)")
        print()
        try:
            sys.path.append(os.getcwd())
            from custom import populate
            stations = populate(data_folder)
        except:
            traceback.print_exc()
            print("No file named custom.py in the %s folder" % os.getcwd())
            return
    stations = np.unique(stations)
    
    db = connect()
    for station in stations:
        net, sta = station.split('_')
        print('Adding:', net, sta)
        X = 0.0
        Y = 0.0
        altitude = 0.0
        coordinates = 'UTM'
        instrument = 'N/A'
        update_station(db, net, sta, X, Y, altitude,
                       coordinates=coordinates, instrument=instrument)
    return True

if __name__ == "__main__":
    main()
