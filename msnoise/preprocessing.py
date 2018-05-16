import sys
import traceback

from obspy.core import UTCDateTime, Stream, read

try:
    from scikits.samplerate import resample
except:
    pass

from .api import *


def preprocess(db, stations, comps, goal_day, params, responses=None):
    """
    Fetches data for each ``stations`` and each ``comps`` using the
    data_availability table in the database.

    To correct for instrument responses, make sure to set ``remove_response``
    to "Y" in the config and to provide the ``responses`` DataFrame.

    :Example:
    >>> from msnoise.api import connect, get_params, preload_instrument_responses
    >>> from msnoise.preprocessing import preprocess
    >>> db = connect()
    >>> params = get_params(db)
    >>> responses = preload_instrument_responses(db)
    >>> st = preprocess(db, ["YA.UV06","YA.UV10"], ["Z",], "2010-09-01", params, responses)
    >>> st
     2 Trace(s) in Stream:
    YA.UV06.00.HHZ | 2010-09-01T00:00:00.000000Z - 2010-09-01T23:59:59.950000Z | 20.0 Hz, 1728000 samples
    YA.UV10.00.HHZ | 2010-09-01T00:00:00.000000Z - 2010-09-01T23:59:59.950000Z | 20.0 Hz, 1728000 samples

    :type db: :class:`sqlalchemy.orm.session.Session`
    :param db: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`msnoise.api.connect`.
    :type stations: list of str
    :param stations: a list of station names, in the format NET.STA.
    :type comps: list of str
    :param comps: a list of component names, in Z,N,E,1,2.
    :type goal_day: str
    :param goal_day: the day of data to load, ISO 8601 format: e.g. 2016-12-31.
    :type params: class
    :param params: an object containing the config parameters, as obtained by
        :func:`msnoise.api.get_params`.
    :type responses: :class:`pandas.DataFrame`
    :param responses: a DataFrame containing the instrument responses, as
        obtained by :func:`msnoise.api.preload_instrument_responses`.
    :rtype: :class:`obspy.core.stream.Stream`
    :return: A Stream object containing all traces.
    """
    datafiles = {}
    output = Stream()
    for station in stations:
        datafiles[station] = {}
        net, sta = station.split('.')
        gd = datetime.datetime.strptime(goal_day, '%Y-%m-%d')
        files = get_data_availability(
            db, net=net, sta=sta, starttime=gd, endtime=gd)
        for comp in comps:
            datafiles[station][comp] = []
        for file in files:
            if file.comp[-1] not in comps:
                continue
            fullpath = os.path.join(file.path, file.file)
            datafiles[station][file.comp[-1]].append(fullpath)

    for istation, station in enumerate(stations):
        net, sta = station.split(".")
        for comp in comps:
            files = eval("datafiles['%s']['%s']" % (station, comp))
            if len(files) != 0:
                logging.debug("%s.%s Reading %i Files" %
                              (station, comp, len(files)))
                stream = Stream()
                for file in sorted(files):
                    try:
                        st = read(file, dytpe=np.float,
                              starttime=UTCDateTime(gd),
                              endtime=UTCDateTime(gd)+86400)
                    except:
                        logging.debug("ERROR reading file %s"%file)
                        continue
                    tmp = st.select(network=net, station=sta, component=comp)
                    if not len(tmp):
                        for tr in st:
                            tr.stats.network = net
                        st = st.select(network=net, station=sta, component=comp)
                    else:
                        st = tmp
                    for tr in st:
                        tr.data = tr.data.astype(np.float)
                        tr.stats.network = tr.stats.network.upper()
                        tr.stats.station = tr.stats.station.upper()
                        tr.stats.channel = tr.stats.channel.upper()

                    stream += st
                    del st
                stream.sort()
                try:
                    # HACK not super clean... should find a way to prevent the
                    # same trace id with different sps to occur
                    stream.merge(method=1, interpolation_samples=3, fill_value=None)
                except:
                    logging.debug("Error while merging...")
                    traceback.print_exc()
                    continue
                stream = stream.split()

                logging.debug("%s Checking sample alignment" % stream[0].id)
                for i, trace in enumerate(stream):
                    stream[i] = check_and_phase_shift(trace)

                logging.debug("%s Checking Gaps" % stream[0].id)
                if len(getGaps(stream)) > 0:
                    max_gap = 10
                    only_too_long = False
                    while getGaps(stream) and not only_too_long:
                        too_long = 0
                        gaps = getGaps(stream)
                        for gap in gaps:
                            if int(gap[-1]) <= max_gap:
                                try:
                                    stream[gap[0]] = stream[gap[0]].__add__(stream[gap[1]], method=1,
                                                                        fill_value="interpolate")
                                    stream.remove(stream[gap[1]])
                                except:
                                    stream.remove(stream[gap[1]])

                                break
                            else:
                                too_long += 1
                        if too_long == len(gaps):
                            only_too_long = True

                stream = stream.split()
                for tr in stream:
                    if tr.stats.sampling_rate < (params.goal_sampling_rate-1):
                        stream.remove(tr)
                taper_length = 20.0  # seconds
                for trace in stream:
                    if trace.stats.npts < 4 * taper_length * trace.stats.sampling_rate:
                        stream.remove(trace)
                    else:
                        trace.detrend(type="demean")
                        trace.detrend(type="linear")
                        trace.taper(max_percentage=None, max_length=1.0)

                if not len(stream):
                    logging.debug(" has only too small traces, skipping...")
                    continue

                for trace in stream:
                    logging.debug(
                        "%s.%s Highpass at %.2f Hz" % (station, comp, params.preprocess_highpass))
                    trace.filter("highpass", freq=params.preprocess_highpass, zerophase=True)

                    if trace.stats.sampling_rate != params.goal_sampling_rate:
                        logging.debug(
                            "%s.%s Lowpass at %.2f Hz" % (station, comp, params.preprocess_lowpass))
                        trace.filter("lowpass", freq=params.preprocess_lowpass, zerophase=True, corners=8)

                        if params.resampling_method == "Resample":
                            logging.debug("%s.%s Downsample to %.1f Hz" %
                                          (station, comp, params.goal_sampling_rate))
                            trace.data = resample(
                                trace.data, params.goal_sampling_rate / trace.stats.sampling_rate, 'sinc_fastest')

                        elif params.resampling_method == "Decimate":
                            decimation_factor = trace.stats.sampling_rate / params.goal_sampling_rate
                            if not int(decimation_factor) == decimation_factor:
                                logging.warning("%s.%s CANNOT be decimated by an integer factor, consider using Resample or Lanczos methods"
                                                " Trace sampling rate = %i ; Desired CC sampling rate = %i" %
                                                (station, comp, trace.stats.sampling_rate, params.goal_sampling_rate))
                                sys.stdout.flush()
                                sys.exit()
                            logging.debug("%s.%s Decimate by a factor of %i" %
                                          (station, comp, decimation_factor))
                            trace.data = trace.data[::int(decimation_factor)]

                        elif params.resampling_method == "Lanczos":
                            logging.debug("%s.%s Downsample to %.1f Hz" %
                                          (station, comp, params.goal_sampling_rate))
                            trace.data = np.array(trace.data)
                            trace.interpolate(method="lanczos", sampling_rate=params.goal_sampling_rate, a=1.0)

                        trace.stats.sampling_rate = params.goal_sampling_rate

                if params.remove_response:
                    logging.debug('%s Removing instrument response'%stream[0].id)

                    response = responses[responses["channel_id"] == stream[0].id]
                    if len(response) > 1:
                        response = response[response["start_date"] <= UTCDateTime(gd)]
                    if len(response) > 1:
                        response = response[response["end_date"] >= UTCDateTime(gd)]
                    elif len(response) == 0:
                        logging.info("No instrument response information "
                                     "for %s, skipping" % stream[0].id)
                        continue
                    try:
                        datalesspz = response["paz"].values[0]
                    except:
                        logging.error("Bad instrument response information "
                                      "for %s, skipping" % stream[0].id)
                        continue
                    stream.simulate(paz_remove=datalesspz,
                                    remove_sensitivity=True,
                                    pre_filt=params.response_prefilt,
                                    paz_simulate=None, )
                for tr in stream:
                    tr.data = tr.data.astype(np.float32)
                output += stream
                del stream
            del files
    clean_scipy_cache()
    return output
