import sys
import traceback

from obspy.core import UTCDateTime, Stream, read

try:
    from scikits.samplerate import resample
except:
    pass

from .api import *
import io
import logbook
logger = logbook.Logger(__name__)

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
    MULTIPLEX = False
    MULTIPLEX_files = {}
    for station in stations:
        datafiles[station] = {}
        net, sta = station.split('.')
        gd = datetime.datetime.strptime(goal_day, '%Y-%m-%d')
        files = get_data_availability(
            db, net=net, sta=sta, starttime=gd, endtime=gd)
        for comp in comps:
            datafiles[station][comp] = []
        for file in files:
            if file.sta != "MULTIPLEX":
                if file.comp[-1] not in comps:
                    continue
                fullpath = os.path.join(file.path, file.file)
                datafiles[station][file.comp[-1]].append(fullpath)
            else:
                MULTIPLEX = True
                print("Mutliplex mode, reading the files")
                fullpath = os.path.join(file.path, file.file)
                multiplexed = sorted(glob.glob(fullpath))
                for comp in comps:
                    for fn in multiplexed:
                        if fn in MULTIPLEX_files:
                            _ = MULTIPLEX_files[fn]
                        else:
                            # print("Reading %s" % fn)
                            _ = read(fn, format=params.archive_format or None)
                            traces = []
                            for tr in _:
                                if "%s.%s" % (tr.stats.network, tr.stats.station) in stations and tr.stats.channel[-1] in comps:
                                    traces.append(tr)
                            del _
                            _ = Stream(traces=traces)
                            MULTIPLEX_files[fn] = _
                        datafiles[station][comp].append(_)

    for istation, station in enumerate(stations):
        net, sta = station.split(".")
        for comp in comps:
            files = eval("datafiles['%s']['%s']" % (station, comp))
            if len(files) != 0:
                logger.debug("%s.%s Reading %i Files" %
                              (station, comp, len(files)))
                traces = []
                for file in files:
                    if isinstance(file, Stream):
                        st = file.select(network=net, station=sta, component=comp).copy()
                    else:
                        try:
                            # print("Reading %s" % file)
                            # t=  time.time()
                            st = read(file, dytpe=np.float,
                                      starttime=UTCDateTime(gd),
                                      endtime=UTCDateTime(gd)+86400,
                                      station=sta,
                                      format=params.archive_format or None)
                            # print("done in", time.time()-t)
                        except:
                            logger.debug("ERROR reading file %s" % file)
                            continue
                    for tr in st:
                        if len(tr.stats.channel) == 2:
                            tr.stats.channel += tr.stats.location
                            tr.stats.location = "00"
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

                        traces.append(tr)
                    del st
                stream = Stream(traces=traces)
                if not(len(stream)):
                    continue
                f = io.BytesIO()
                stream.write(f, format='MSEED')
                f.seek(0)
                stream = read(f, format="MSEED")

                stream.sort()
                # try:
                #     # HACK not super clean... should find a way to prevent the
                #     # same trace id with different sps to occur
                #     stream.merge(method=1, interpolation_samples=3, fill_value=None)
                # except:
                #     logger.debug("Error while merging...")
                #     traceback.print_exc()
                #     continue
                # stream = stream.split()
                if not len(stream):
                    continue
                logger.debug("%s Checking sample alignment" % stream[0].id)
                for i, trace in enumerate(stream):
                    stream[i] = check_and_phase_shift(trace)

                logger.debug("%s Checking Gaps" % stream[0].id)
                if len(getGaps(stream)) > 0:
                    max_gap = params.preprocess_max_gap*stream[0].stats.sampling_rate

                    gaps = getGaps(stream)
                    while len(gaps):
                        too_long = 0
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
                            break
                        gaps = getGaps(stream)
                    del gaps

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
                    logger.debug(" has only too small traces, skipping...")
                    continue

                for trace in stream:
                    logger.debug(
                        "%s Highpass at %.2f Hz" % (trace.id, params.preprocess_highpass))
                    trace.filter("highpass", freq=params.preprocess_highpass, zerophase=True, corners=4)

                    if trace.stats.sampling_rate != params.goal_sampling_rate:
                        logger.debug(
                            "%s Lowpass at %.2f Hz" % (trace.id, params.preprocess_lowpass))
                        trace.filter("lowpass", freq=params.preprocess_lowpass, zerophase=True, corners=8)

                        if params.resampling_method == "Resample":
                            logger.debug("%s Downsample to %.1f Hz" %
                                          (trace.id, params.goal_sampling_rate))
                            trace.data = resample(
                                trace.data, params.goal_sampling_rate / trace.stats.sampling_rate, 'sinc_fastest')

                        elif params.resampling_method == "Decimate":
                            decimation_factor = trace.stats.sampling_rate / params.goal_sampling_rate
                            if not int(decimation_factor) == decimation_factor:
                                logger.warning("%s CANNOT be decimated by an integer factor, consider using Resample or Lanczos methods"
                                                " Trace sampling rate = %i ; Desired CC sampling rate = %i" %
                                                (trace.id, trace.stats.sampling_rate, params.goal_sampling_rate))
                                sys.stdout.flush()
                                sys.exit()
                            logger.debug("%s Decimate by a factor of %i" %
                                          (trace.id, decimation_factor))
                            trace.data = trace.data[::int(decimation_factor)]

                        elif params.resampling_method == "Lanczos":
                            logger.debug("%s Downsample to %.1f Hz" %
                                          (trace.id, params.goal_sampling_rate))
                            trace.data = np.array(trace.data)
                            trace.interpolate(method="lanczos", sampling_rate=params.goal_sampling_rate, a=1.0)

                        trace.stats.sampling_rate = params.goal_sampling_rate
                    del trace

                if params.remove_response:
                    logger.debug('%s Removing instrument response'%stream[0].id)

                    response = responses[responses["channel_id"] == stream[0].id]
                    if len(response) > 1:
                        response = response[response["start_date"] <= UTCDateTime(gd)]
                    if len(response) > 1:
                        response = response[response["end_date"] >= UTCDateTime(gd)]
                    elif len(response) == 0:
                        logger.info("No instrument response information "
                                     "for %s, skipping" % stream[0].id)
                        continue
                    try:
                        datalesspz = response["paz"].values[0]
                    except:
                        logger.error("Bad instrument response information "
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
    del MULTIPLEX_files
    return output
