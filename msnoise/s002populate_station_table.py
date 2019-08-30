"""
This script is responsible for rapidly scanning the data archive,
identifying the Networks/Stations and inserting them in the *stations* table in
the database.

The ``data_folder`` (as defined in the config) is scanned following the
``data_structure``. Possible values for the data_structure are defined in
*data_structures.py*:

.. code-block:: python
    
    data_structure['SDS'] = "YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"
    data_structure['BUD'] = "NET/STA/STA.NET.LOC.CHAN.YEAR.DAY"
    data_structure['IDDS'] = "YEAR/NET/STA/CHAN.TYPE/DAY/NET.STA.LOC.CHAN.TYPE.YEAR.DAY.HOUR"
    data_structure['PDF'] = "YEAR/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"

If your data structure corresponds to one of these 4 structures, you need to
select the corresponding acronym (SDS, BUD, IDDS or PDF) for the data_structure
field.

More info on the recommended SDS ("SeisComP Data Structure") can be found here:
https://www.seiscomp3.org/wiki/doc/applications/slarchive/SDS
For other simple structures, one has to edit the `data_structure` configuration
(see below).

By default, station coordinates are initialized at 0.

To run this script:

.. code-block:: sh

    $ msnoise populate

Custom data structure & station table population
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If one's data structure does not belong to the pre-defined ones, it can be
defined directly in the ``data_structure`` configuration field using forward
slashes, e.g.:

``data_structure`` = "NET/STA/YEAR/NET.STA.YEAR.DAY.MSEED"

MSNoise expects to find a file named ``custom.py`` in the current folder.
This python file will contain a function called ``populate`` wich will accept
one argument and return a station dictionary with keys of the format ``NET_STA``
, and fields for the stations table in the database: Net,Sta,X,Y,Altitude,
Coordinates(UTM/DEG),Instrument.

.. code-block:: python
    
    import os, glob
    def populate(data_folder):
        datalist = sorted(glob.glob(os.path.join(data_folder, "*", "*")))
        stationdict = {}
        for di in datalist:
            tmp = os.path.split(di)
            sta = tmp[1]
            net = os.path.split(tmp[0])[1]
            stationdict[net+"_"+sta]=[net,sta,0.0,0.0,0.0,'UTM','N/A']
        return stationdict

.. _populate-expert:

Expert (lazy) mode:
~~~~~~~~~~~~~~~~~~~

If the `DataAvailability` has already been filled in by another process, for
example using the :ref:`"scan from path"<scan-archive-expert>` procedure, the
network/station names can be "populated" from the `DataAvailability` table
automatically. To do this, simply run:

.. code-block:: sh

    msnoise populate --fromDA

and MSNoise will insert the unique NET.STA in the `Stations` table.
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
        stationdict={}
        for di in datalist:
            tmp = os.path.split(di)
            sta = tmp[1]
            net = os.path.split(tmp[0])[1]
            stationdict[net+"_"+sta]=[net,sta,0.0,0.0,0.0,'UTM','N/A']
        del datalist
    elif data_structure in ["BUD", ]:
        datalist = sorted(glob.glob(os.path.join(data_folder, "*", "*",)))
        stationdict={}
        for di in datalist:
            tmp = os.path.split(di)
            sta = tmp[1]
            net = os.path.split(tmp[0])[1]
            stationdict[net+"_"+sta]=[net,sta,0.0,0.0,0.0,'UTM','N/A']
        del datalist
    elif data_structure in ["PDF", ]:
        datalist = sorted(glob.glob(os.path.join(data_folder, "*", "*",)))
        stationdict={}
        for di in datalist:
            tmp = os.path.split(di)
            sta = tmp[1]
            net = get_config(db, 'network')
            stationdict[net+"_"+sta]=[net,sta,0.0,0.0,0.0,'UTM','N/A']
        del datalist
    else:
        print("Can't parse the archive for format %s !" % data_structure)
        print("trying to import local parser (should return a station dictionary)")
        print("")
        try:
            sys.path.append(os.getcwd())
            from custom import populate
            stationdict = populate(data_folder)
        except:
            traceback.print_exc()
            print("No file named custom.py in the %s folder" % os.getcwd())
            return

    db = connect()
    for s in stationdict.keys() :
        net,sta,lon,lat,alt,coordinates,instype=stationdict[s]
        print('Adding:', net, sta)
        X = float(lon)
        Y = float(lat)
        altitude = float(alt)
        instrument = str(instype)
        update_station(db, net, sta, X, Y, altitude,
                       coordinates=coordinates, instrument=instrument)

    return True

if __name__ == "__main__":
    main()
