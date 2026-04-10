"""MSNoise station and data-availability management."""

__all__ = [
    "ARCHIVE_STRUCTURES",
    "add_data_source",
    "check_stations_uniqueness",
    "get_data_availability",
    "get_data_source",
    "get_default_data_source",
    "get_interstation_distance",
    "get_new_files",
    "get_station",
    "get_station_pairs",
    "get_stations",
    "get_waveform_path",
    "import_stationxml",
    "list_data_sources",
    "mark_data_availability",
    "read_waveforms_from_availability",
    "resolve_data_source",
    "set_all_stations_source",
    "set_network_source",
    "set_station_source",
    "to_sds",
    "update_data_availability",
    "update_data_source",
    "update_station",
    "populate_stations",
]

import itertools
import logging

import numpy as np

from ..msnoise_table_def import Station, DataAvailability

logger = logging.getLogger('msnoise.stations')

#: Archive path templates keyed by structure name.
#: Used by :func:`~msnoise.s01_scan_archive` and :func:`populate_stations`.
ARCHIVE_STRUCTURES = {
    "SDS":  "YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY",
    "BUD":  "NET/STA/STA.NET.LOC.CHAN.YEAR.DAY",
    "IDDS": "YEAR/NET/STA/CHAN.TYPE/DAY/NET.STA.LOC.CHAN.TYPE.YEAR.DAY.HOUR",
    "PDF":  "YEAR/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY",
}


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



def get_station_pairs(session, used=None, net=None, include_single_station=False):
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
    if include_single_station:
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
                             samplerate, data_source_id=None):
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
    :param path: Path to the folder containing the file, **relative** to
        ``DataSource.uri``.
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
    :type data_source_id: int or None
    :param data_source_id: FK to DataSource. ``None`` → project default.
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
                                data_duration, gaps_duration, samplerate, flag,
                                data_source_id=data_source_id)
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
            if data_source_id is not None:
                data.data_source_id = data_source_id
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


# ============================================================
# DataSource management
# ============================================================

def add_data_source(session, name, uri="", data_structure="SDS",
                    auth_env="MSNOISE"):
    """Add a new :class:`~msnoise.msnoise_table_def.DataSource` to the database.

    :param session: SQLAlchemy session.
    :param name: Human label e.g. ``"local"``, ``"IRIS"``, ``"EIDA"``.
    :param uri: Data location URI. Schemes:
        - bare path or ``sds:///path`` → local SDS archive
        - ``fdsn://http://...``         → FDSN web service
        - ``eida://http://...``         → EIDA routing client
    :param data_structure: SDS sub-path format for local sources (default ``"SDS"``).
        Ignored for FDSN/EIDA.
    :param auth_env: Environment variable prefix for credentials
        (default ``"MSNOISE"``). Worker looks up
        ``{auth_env}_FDSN_USER``, ``{auth_env}_FDSN_PASSWORD``,
        ``{auth_env}_FDSN_TOKEN``.
    :returns: The created :class:`~msnoise.msnoise_table_def.DataSource` object.
    """
    from ..msnoise_table_def import DataSource
    ds = DataSource(name=name, uri=uri, data_structure=data_structure,
                    auth_env=auth_env)
    session.add(ds)
    session.commit()
    session.refresh(ds)
    logger.info(f"Added DataSource {name!r} (id={ds.ref}, uri={uri!r})")
    return ds


def get_data_source(session, id=None, name=None):
    """Retrieve a :class:`~msnoise.msnoise_table_def.DataSource` by id or name.

    :param session: SQLAlchemy session.
    :param id: Primary key of the DataSource.
    :param name: Human label of the DataSource.
    :returns: :class:`~msnoise.msnoise_table_def.DataSource` or ``None``.
    :raises ValueError: If neither *id* nor *name* is provided.
    """
    from ..msnoise_table_def import DataSource
    if id is not None:
        return session.query(DataSource).filter(DataSource.ref == id).first()
    if name is not None:
        return session.query(DataSource).filter(DataSource.name == name).first()
    raise ValueError("Provide either id or name")


def get_default_data_source(session):
    """Return the project default :class:`~msnoise.msnoise_table_def.DataSource`.

    The default is the ``DataSource`` with ``ref=1`` (created by the installer).

    :param session: SQLAlchemy session.
    :returns: :class:`~msnoise.msnoise_table_def.DataSource`.
    :raises RuntimeError: If no DataSource exists (installer not run).
    """
    from ..msnoise_table_def import DataSource
    ds = session.query(DataSource).filter(DataSource.ref == 1).first()
    if ds is None:
        raise RuntimeError(
            "No default DataSource found (ref=1). Has the installer been run?"
        )
    return ds


def list_data_sources(session):
    """Return all :class:`~msnoise.msnoise_table_def.DataSource` rows.

    :param session: SQLAlchemy session.
    :returns: List of :class:`~msnoise.msnoise_table_def.DataSource`.
    """
    from ..msnoise_table_def import DataSource
    return session.query(DataSource).order_by(DataSource.ref).all()


def update_data_source(session, id, **kwargs):
    """Update fields on an existing :class:`~msnoise.msnoise_table_def.DataSource`.

    :param session: SQLAlchemy session.
    :param id: Primary key of the DataSource to update.
    :param kwargs: Fields to update: ``name``, ``uri``, ``data_structure``,
        ``auth_env``.
    :returns: Updated :class:`~msnoise.msnoise_table_def.DataSource`.
    :raises ValueError: If the DataSource does not exist.
    """
    ds = get_data_source(session, id=id)
    if ds is None:
        raise ValueError(f"DataSource with id={id} not found")
    allowed = {"name", "uri", "data_structure", "auth_env", "archive_format", "network_code", "channels"}
    for key, val in kwargs.items():
        if key not in allowed:
            raise ValueError(f"Unknown DataSource field: {key!r}")
        setattr(ds, key, val)
    session.commit()
    logger.info(f"Updated DataSource id={id}: {kwargs}")
    return ds


def resolve_data_source(session, station):
    """Return the effective :class:`~msnoise.msnoise_table_def.DataSource` for
    a station.

    If the station has ``data_source_id`` set, that DataSource is returned.
    Otherwise the project default (``ref=1``) is used.

    :param session: SQLAlchemy session.
    :param station: :class:`~msnoise.msnoise_table_def.Station` ORM object.
    :returns: :class:`~msnoise.msnoise_table_def.DataSource`.
    """
    if station.data_source_id is not None:
        return get_data_source(session, id=station.data_source_id)
    return get_default_data_source(session)


# ============================================================
# Station ↔ DataSource assignment
# ============================================================

def set_station_source(session, net, sta, data_source_id):
    """Assign a specific :class:`~msnoise.msnoise_table_def.DataSource` to a station.

    :param session: SQLAlchemy session.
    :param net: Network code.
    :param sta: Station code.
    :param data_source_id: Primary key of the DataSource to assign.
        Pass ``None`` to revert to the project default.
    :raises ValueError: If the station or DataSource does not exist.
    """
    from ..msnoise_table_def import DataSource
    station = get_station(session, net, sta)
    if station is None:
        raise ValueError(f"Station {net}.{sta} not found")
    if data_source_id is not None:
        ds = session.query(DataSource).filter(
            DataSource.ref == data_source_id).first()
        if ds is None:
            raise ValueError(f"DataSource with id={data_source_id} not found")
    station.data_source_id = data_source_id
    session.commit()
    logger.info(f"Set {net}.{sta} → DataSource id={data_source_id}")


def set_network_source(session, net, data_source_id):
    """Assign a DataSource to all stations of a network.

    :param session: SQLAlchemy session.
    :param net: Network code.
    :param data_source_id: Primary key of the DataSource. ``None`` reverts to default.
    """
    stations = session.query(Station).filter(Station.net == net).all()
    if not stations:
        raise ValueError(f"No stations found for network {net!r}")
    for sta in stations:
        sta.data_source_id = data_source_id
    session.commit()
    logger.info(
        f"Set all {len(stations)} stations of network {net!r} "
        f"→ DataSource id={data_source_id}"
    )


def set_all_stations_source(session, data_source_id):
    """Assign a DataSource to every station in the database.

    :param session: SQLAlchemy session.
    :param data_source_id: Primary key of the DataSource. ``None`` reverts all to default.
    """
    n = session.query(Station).update({"data_source_id": data_source_id})
    session.commit()
    logger.info(f"Set all {n} stations → DataSource id={data_source_id}")


# ============================================================
# StationXML import
# ============================================================

def import_stationxml(session, path_or_url, data_source_id=None,
                      save_to_response_path=True):
    """Import stations from a StationXML file or URL.

    Parses the inventory with ObsPy and creates or updates
    :class:`~msnoise.msnoise_table_def.Station` rows.  Also populates
    ``used_location_codes`` and ``used_channel_names`` from the channel-level
    inventory (requires ``level=channel`` in FDSN queries).  Empty location
    codes are stored as ``\"--\"`` (MSNoise convention).

    When *save_to_response_path* is ``True`` (default), the parsed inventory is
    written as a StationXML file into the project's ``response_path`` directory
    (read from global config).  This makes the instrument responses immediately
    available to the preprocessing step without any manual file copying.  The
    file is named ``<NET>.<STA>.xml`` for single-network inventories, or
    ``inventory_<timestamp>.xml`` when the inventory spans multiple networks.
    Saving failures are logged as warnings but never abort the import.

    :param session: SQLAlchemy session.
    :param path_or_url: Path to a local StationXML file, or a URL
        (e.g. an FDSN station web service query URL with ``level=channel``).
    :param data_source_id: DataSource to assign to imported stations.
        ``None`` → project default (``DataSource.ref=1``).
    :param save_to_response_path: Write the inventory to ``response_path``
        after a successful import (default ``True``).
    :returns: Tuple ``(created, updated, saved_path)`` where *saved_path* is
        the absolute path of the written file, or ``None`` if not saved.
    """
    import os
    import datetime as _dt
    from obspy import read_inventory
    logger.info(f"Importing StationXML from {path_or_url!r}")
    inv = read_inventory(path_or_url)

    created = 0
    updated = 0
    for network in inv:
        for station in network:
            net = network.code
            sta = station.code
            X = station.longitude
            Y = station.latitude
            alt = station.elevation

            # Collect unique location codes and channel names from channel list.
            # Empty location code → "--" (MSNoise convention).
            locs  = sorted({ch.location_code or "--" for ch in station.channels})
            chans = sorted({ch.code for ch in station.channels})
            loc_str  = ",".join(locs)  if locs  else None
            chan_str = ",".join(chans) if chans else None

            existing = get_station(session, net, sta)
            if existing is None:
                new_sta = Station(
                    net=net, sta=sta,
                    X=X, Y=Y, altitude=alt,
                    coordinates="DEG",
                    used=True,
                    data_source_id=data_source_id,
                )
                if loc_str is not None:
                    new_sta.used_location_codes = loc_str
                if chan_str is not None:
                    new_sta.used_channel_names = chan_str
                session.add(new_sta)
                created += 1
                logger.debug(
                    f"  Created station {net}.{sta}"
                    + (f" locs={loc_str}" if loc_str else "")
                    + (f" chans={chan_str}" if chan_str else "")
                )
            else:
                existing.X = X
                existing.Y = Y
                existing.altitude = alt
                existing.coordinates = "DEG"
                if data_source_id is not None:
                    existing.data_source_id = data_source_id
                if loc_str is not None:
                    existing.used_location_codes = loc_str
                if chan_str is not None:
                    existing.used_channel_names = chan_str
                updated += 1
                logger.debug(
                    f"  Updated station {net}.{sta}"
                    + (f" locs={loc_str}" if loc_str else "")
                    + (f" chans={chan_str}" if chan_str else "")
                )
    session.commit()
    logger.info(f"StationXML import done: {created} created, {updated} updated")

    # ── Optionally save inventory to response_path ────────────────────────────
    saved_path = None
    if save_to_response_path:
        try:
            from .config import get_config
            response_path = get_config(session, "response_path",
                                        category="global")
            if not response_path:
                logger.warning(
                    "response_path not configured — StationXML not saved to disk."
                )
            else:
                os.makedirs(response_path, exist_ok=True)
                # Build a deterministic filename from the inventory contents
                nets = list({net.code for net in inv})
                stas = list({sta.code for net in inv for sta in net})
                if len(nets) == 1 and len(stas) == 1:
                    fname = f"{nets[0]}.{stas[0]}.xml"
                elif len(nets) == 1:
                    fname = f"{nets[0]}.xml"
                else:
                    ts = _dt.datetime.now(tz=_dt.timezone.utc).strftime("%Y%m%dT%H%M%S")
                    fname = f"inventory_{ts}.xml"
                dest = os.path.join(response_path, fname)
                inv.write(dest, format="STATIONXML")
                saved_path = dest
                logger.info(f"Inventory saved to {dest!r}")
        except Exception as exc:
            logger.warning(
                f"Could not save inventory to response_path: {exc}"
            )

    return created, updated, saved_path
    session.commit()
    logger.info(f"StationXML import done: {created} created, {updated} updated")
    return created, updated

def get_waveform_path(session, da):
    """Reconstruct the full absolute path for a DataAvailability record.

    Joins ``DataSource.uri`` with ``da.path`` and ``da.file``.  If the DA
    record has no ``data_source_id`` (legacy or NULL), falls back to treating
    ``da.path`` as an absolute path (backward-compatible behaviour).

    :param session: SQLAlchemy session.
    :param da: :class:`~msnoise.msnoise_table_def.DataAvailability` ORM object.
    :returns: Absolute path string to the waveform file.
    """
    import os
    if da.data_source_id is not None:
        ds = get_data_source(session, id=da.data_source_id)
    else:
        ds = get_default_data_source(session)
    root = ds.uri if ds and ds.uri else ""
    return os.path.join(root, da.path, da.file) if root else os.path.join(da.path, da.file)


def read_waveforms_from_availability(session, da_records, t_start, t_end, logger=None):
    """Read waveforms from a list of DataAvailability records into a Stream.

    Builds the full path for each record via :func:`get_waveform_path`
    (honouring the DataSource root), deduplicates paths, and reads each file
    sliced to ``[t_start, t_end]`` using ObsPy.  Files that cannot be read
    are silently skipped with a debug log.

    This is the canonical "DA records → waveform Stream" helper that any
    worker (PSD, future steps) should use instead of hand-rolling
    ``os.path.join(f.path, f.file)``.

    :param session: SQLAlchemy session.
    :param da_records: Iterable of DataAvailability ORM objects.
    :param t_start: :class:`~obspy.core.utcdatetime.UTCDateTime` start.
    :param t_end: :class:`~obspy.core.utcdatetime.UTCDateTime` end.
    :param logger: Optional :class:`logging.Logger`; uses module logger if None.
    :returns: :class:`~obspy.core.stream.Stream` (may be empty).
    """
    import logging as _logging
    import obspy as _obspy

    _log = logger or _logging.getLogger("msnoise.stations")

    paths = list(dict.fromkeys(get_waveform_path(session, da) for da in da_records))
    st = _obspy.Stream()
    for fpath in paths:
        try:
            st += _obspy.read(fpath, starttime=t_start, endtime=t_end)
        except Exception as exc:
            _log.debug(f"Could not read {fpath}: {exc}")
    return st



def to_sds(stats, year, jday):
    """Build an SDS-format relative file path from ObsPy trace stats.

    Returns a path string of the form::

        YYYY/NET/STA/CHAN.D/NET.STA.LOC.CHAN.D.YYYY.DDD

    :param stats: :class:`obspy.core.trace.Stats` object.
    :param year: 4-digit year integer.
    :param jday: Julian day-of-year integer (1-366).
    :returns: Relative SDS path string.
    """
    SDS = "YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.JDAY"
    f = SDS.replace("YEAR", "%04i" % year)
    f = f.replace("NET",  stats.network)
    f = f.replace("STA",  stats.station)
    f = f.replace("LOC",  stats.location)
    f = f.replace("CHAN", stats.channel)
    f = f.replace("JDAY", "%03i" % jday)
    f = f.replace("TYPE", "D")
    return f


def populate_stations(db, loglevel="INFO"):
    """Scan the default DataSource archive and populate the Station table.

    Walks the archive directory according to the configured ``data_structure``
    (SDS, BUD, IDDS, PDF, or a custom format).  Station coordinates are
    initialised to zero — use :func:`import_stationxml` afterwards to set
    real coordinates.

    For custom data structures, a ``custom.py`` file in the current working
    directory must export a ``populate(data_folder)`` function that returns a
    ``{NET_STA: [net, sta, lon, lat, alt, coordinates]}`` dictionary.

    :param db: SQLAlchemy session.
    :param loglevel: Logging verbosity (default ``"INFO"``).
    :returns: ``True`` on success.
    """
    import glob
    import os
    import sys
    import traceback
    import logging as _logging
    from .db import get_logger

    logger = get_logger("msnoise.db_populate", loglevel, with_pid=True)
    logger.info("Populating the Station table")

    _ds = get_default_data_source(db)
    data_folder = _ds.uri or ""
    data_structure = _ds.data_structure or "SDS"
    network_override = _ds.network_code or "*"

    stationdict = {}
    if data_structure in ("SDS", "IDDS"):
        for di in sorted(glob.glob(os.path.join(data_folder, "*", "*", "*"))):
            sta = os.path.basename(di)
            net = os.path.basename(os.path.dirname(di))
            stationdict[f"{net}_{sta}"] = [net, sta, 0.0, 0.0, 0.0, "UTM"]
    elif data_structure == "BUD":
        for di in sorted(glob.glob(os.path.join(data_folder, "*", "*"))):
            sta = os.path.basename(di)
            net = os.path.basename(os.path.dirname(di))
            stationdict[f"{net}_{sta}"] = [net, sta, 0.0, 0.0, 0.0, "UTM"]
    elif data_structure == "PDF":
        for di in sorted(glob.glob(os.path.join(data_folder, "*", "*"))):
            sta = os.path.basename(di)
            stationdict[f"{network_override}_{sta}"] = [
                network_override, sta, 0.0, 0.0, 0.0, "UTM"]
    else:
        logger.warning(f"Unknown data_structure {data_structure!r} — "
                       "trying local custom.py parser.")
        try:
            sys.path.insert(0, os.getcwd())
            from custom import populate  # noqa: PLC0415
            stationdict = populate(data_folder)
        except Exception:
            traceback.print_exc()
            _logging.error("No custom.py found in %s" % os.getcwd())
            return False
        finally:
            if os.getcwd() in sys.path:
                sys.path.remove(os.getcwd())

    for key, (net, sta, lon, lat, alt, coordinates) in stationdict.items():
        logger.info(f"Adding: {net}.{sta}")
        update_station(db, net=net, sta=sta,
                       X=float(lon), Y=float(lat), altitude=float(alt),
                       coordinates=coordinates)
    return True
