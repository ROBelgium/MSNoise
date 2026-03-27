"""MSNoise station and data-availability management."""
import datetime
import glob
import os

from .db import connect, get_logger
from .config import get_config
from .msnoise_table_def import Station, DataAvailability

def get_stations(session, all=False, net=None, format="raw"):
    """Get Stations from the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type all: bool
    :param all: Returns all stations from the database if True, or only
        stations where `used` = 1 if False (default)
    :type net: str
    :param net: if set, limits the stations returned to this network

    :rtype: list of :class:`msnoise.msnoise_table_def.declare_tables.Station`
    :returns: list of :class:`~msnoise.msnoise_table_def.declare_tables.Station`
    """
    q = session.query(Station)
    if all:
        if net is not None:
            stations = q.filter(Station.net == net).order_by(Station.net).\
                order_by(Station.sta)
        else:
            stations = q.order_by(Station.net).order_by(Station.sta).all()
    else:
        stations = q.filter(Station.used.is_(True)).order_by(Station.net).\
            order_by(Station.sta)
        if net is not None:
            stations = stations.filter(Station.net == net).\
                order_by(Station.net).order_by(Station.sta)
    if format == "raw":
        return stations
    if format == "seed_id":
        output = []
        for sta in stations:
            if sta.used_location_codes is None:
                location_codes = []
            else:
                location_codes = sta.used_location_codes.split(",")
            if sta.used_channel_names is None:
                channels = []
            else:
                channels = []
                for i, chan in enumerate(sta.used_channel_names.split(",")):
                    if chan.count("?"):
                        for comp in ["Z", "N", "E", "1", "2"]:
                            channels.append(chan.replace("?", comp))
                    else:
                        channels.append(chan)

            for loc in location_codes:
                for chan in channels:
                    output.append("%s.%s.%s.%s" % (sta.net, sta.sta, loc, chan))
        return output



def get_station(session, net, sta):
    """Get one Station from the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type net: str
    :param net: the network code
    :type sta: str
    :param sta: the station code

    :rtype: :class:`msnoise.msnoise_table_def.declare_tables.Station`
    :returns: a :class:`~msnoise.msnoise_table_def.declare_tables.Station` Object

    """
    station = session.query(Station).filter(Station.net == net).\
        filter(Station.sta == sta).first()
    return station



def update_station(session, net, sta, X, Y, altitude, coordinates='UTM',
                   instrument='N/A', used=1):
    """Updates or Insert a new Station in the database.

    .. seealso :: :class:`msnoise.msnoise_table_def.declare_tables.Station`

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type net: str
    :param net: The network code of the Station
    :type sta: str
    :param sta: The station code
    :type X: float
    :param X: The X coordinate of the station (Easting or Longitude)
    :type Y: float
    :param Y: The Y coordinate of the station (Northing or Latitude)
    :type altitude: float
    :param altitude: The altitude of the station
    :type coordinates: str
    :param coordinates: The coordinates system. "DEG" is WGS84 latitude/
        longitude in degrees. "UTM" is expressed in meters.
    :type instrument: str
    :param instrument: The instrument code, useful with PAZ correction
    :type used: bool
    :param used: Whether this station must be used in the computations.
    """

    if coordinates == "DEG" and (not -90 <= Y <= 90 or not -180 <= X <= 180):
        raise ValueError("Coordinates must be valid WGS84 latitude (%.4f) and longitude (%.4f). " % (Y, X))

    station = session.query(Station).filter(Station.net == net).\
        filter(Station.sta == sta).first()
    if station is None:
        station = Station(net, sta, X, Y, altitude, coordinates, used)
        session.add(station)
    else:
        station.X = X
        station.Y = Y
        station.altitude = altitude
        station.coordinates = coordinates
        station.used = used
    session.commit()
    return True



def get_station_pairs(session, used=None, net=None):
    """Returns an iterator over all possible station pairs.
    If auto-correlation is configured in the database, returns N*N pairs,
    otherwise returns N*(N-1)/2 pairs.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type used: bool, int
    :param used: Select only stations marked used if False (default) or all
        stations present in the database if True
    :type net: str
    :param net: Network code to filter for the pairs.

    :rtype: iterable
    :returns: An iterable of :class:`~msnoise.msnoise_table_def.declare_tables.Station` object
        pairs
    """
    stations = get_stations(session, all=False, net=net)
    if len(get_config(session, name="components_to_compute_single_station")):
        return itertools.combinations_with_replacement(stations, 2)
    else:
        return itertools.combinations(stations, 2)



def check_stations_uniqueness(session, station):
    """

    :param session:
    :param station:
    :return:
    """
    # if the station is net.sta.loc, nothing to do
    if station.count(".") == 2:
        return station

    logging.info("It seems you're voluntarily missing the location code for"
                 " \"%s\". We'll handle this automatically, if there are no "
                 "conflicts." % station)
    net, sta = station.split(".")
    locs = get_station(session, net, sta).locs()
    if len(locs) != 1:
        logging.info("There are more than 1 location codes for this station: "
                     "%s" % locs)
        return station
    station += ".%s" % locs[0]
    logging.info("Found %s to be the unique solution for this station" % station)
    return station



def get_interstation_distance(station1, station2, coordinates="DEG"):
    """Returns the distance in km between `station1` and `station2`.

    .. warning:: Currently the stations coordinates system have to be the same!

    :type station1: :class:`~msnoise.msnoise_table_def.declare_tables.Station`
    :param station1: A Station object
    :type station2: :class:`~msnoise.msnoise_table_def.declare_tables.Station`
    :param station2: A Station object
    :type coordinates: str
    :param coordinates: The coordinates system. "DEG" is WGS84 latitude/
        longitude in degrees. "UTM" is expressed in meters.

    :rtype: float
    :returns: The interstation distance in km
    """
    from obspy.geodetics import gps2dist_azimuth
    if coordinates == "DEG":
        dist, azim, bazim = gps2dist_azimuth(station1.Y, station1.X,
                                             station2.Y, station2.X)
        return dist / 1.e3
    else:
        dist = np.hypot(float(station1.X - station2.X),
                        float(station1.Y - station2.Y)) / 1.e3
        return dist

# ============================================================


def update_data_availability(session, net, sta, loc, chan, path, file, starttime,
                             endtime, data_duration, gaps_duration,
                             samplerate):
    """
    Updates a DataAvailability object in the database

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type net: str
    :param net: The network code of the Station
    :type sta: str
    :param sta: The station code
    :type chan: str
    :param chan: The component (channel)
    :type path: str
    :param path: The full path to the folder containing the file
    :type file: str
    :param file: The name of the file
    :type starttime: datetime.datetime
    :param starttime: Start time of the file
    :type endtime: datetime.datetime
    :param endtime: End time of the file
    :type data_duration: float
    :param data_duration: Cumulative duration of available data in the file
    :type gaps_duration: float
    :param gaps_duration: Cumulative duration of gaps in the file
    :type samplerate: float
    :param samplerate: Sample rate of the data in the file (in Hz)
    """

    data = session.query(DataAvailability).\
        filter(DataAvailability.path == path). \
        filter(DataAvailability.file == file).\
        filter(DataAvailability.net == net).\
        filter(DataAvailability.sta == sta). \
        filter(DataAvailability.loc == loc). \
        filter(DataAvailability.chan == chan).first()
    if data is None:
        flag = "N"
        data = DataAvailability(net, sta, loc, chan, path, file, starttime, endtime,
                                data_duration, gaps_duration, samplerate, flag)
        session.add(data)
        toreturn = 1
    else:
        modified = False
        for item in ['net', 'sta', 'loc', 'chan', 'path', 'starttime',
                     'endtime', 'data_duration', 'gaps_duration', 'samplerate']:
            if eval("data.%s != %s" % (item, item)):
                modified = True
                break
        if modified:
            data.net = net
            data.sta = sta
            data.loc = loc
            data.chan = chan
            data.path = path
            data.starttime = starttime
            data.endtime = endtime
            data.data_duration = data_duration
            data.gaps_duration = gaps_duration
            data.samplerate = samplerate
            data.flag = "M"
            toreturn = -1
        else:
            toreturn = 0
    session.commit()
    return toreturn



def get_new_files(session):
    """
    Returns the files marked "N"ew or "M"odified in the database

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: list
    :returns: list of :class:`~msnoise.msnoise_table_def.declare_tables.DataAvailability`
    """

    files = session.query(DataAvailability).\
        filter(DataAvailability.flag != 'A').\
        order_by(DataAvailability.starttime).all()
    return files



def get_data_availability(session, net=None, sta=None, loc=None, chan=None,
                          starttime=None, endtime=None):
    """
    Returns the :class:`~msnoise.msnoise_table_def.declare_tables.DataAvailability` objects
    for specific `net`, `sta`, `starttime` or `endtime`

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type net: str
    :param net: Network code
    :type sta: str
    :param sta: Station code
    :type starttime: datetime.datetime, datetime.date
    :param starttime: Start time of the search
    :type endtime: datetime.datetime, datetime.date
    :param endtime: End time of the search

    :rtype: list
    :returns: list of :class:`~msnoise.msnoise_table_def.declare_tables.DataAvailability`
    """
    from sqlalchemy.sql.expression import func
    if not starttime:
        data = session.query(DataAvailability).\
            filter(DataAvailability.net == net).\
            filter(DataAvailability.sta == sta). \
            filter(DataAvailability.loc == loc). \
            filter(DataAvailability.chan == chan).all()
    elif not net:
        data = session.query(DataAvailability).\
            filter(DataAvailability.starttime <= endtime).\
            filter(DataAvailability.endtime >= starttime).all()
    else:
        data = session.query(DataAvailability).\
            filter(DataAvailability.net == net).\
            filter(DataAvailability.sta == sta). \
            filter(DataAvailability.loc == loc). \
            filter(func.DATE(DataAvailability.starttime) <= endtime.date()).\
            filter(func.DATE(DataAvailability.endtime) >= starttime.date()).all()
        if not len(data):
            data = session.query(DataAvailability). \
                filter(DataAvailability.sta == "MULTIPLEX"). \
                filter(func.DATE(DataAvailability.starttime) <= endtime.date()).\
                filter(
                func.DATE(DataAvailability.endtime) >= starttime.date()).all()
    return data



def mark_data_availability(session, net, sta, flag):
    """
    Updates the flag of all
    :class:`~msnoise.msnoise_table_def.declare_tables.DataAvailability` objects matching
    `net.sta` in the database

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type net: str
    :param net: Network code
    :type sta: str
    :param sta: Station code
    :type flag: str
    :param flag: Status of the DataAvailability object: New, Modified or
        Archive. Values accepted are {'N', 'M', 'A'}
    """
    logging.debug("Updating: %s %s to flag=%s" %(net, sta, flag))
    da = DataAvailability.__table__
    stmt = da.update().where(da.c.sta==sta).where(da.c.net==net).\
        values(flag=flag)
    session.execute(stmt)
    session.commit()
