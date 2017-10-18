import calendar
import glob
import sys
import time
import traceback


from obspy.core import UTCDateTime, Stream, read

try:
    from scikits.samplerate import resample
except:
    pass

from .api import *


def preprocess(db, stations, comps, goal_day, params, responses=None):

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
                    st = read(file, dytpe=np.float,
                              starttime=UTCDateTime(gd),
                              endtime=UTCDateTime(gd)+86400)
                    tmp = st.select(network=net, station=sta, component=comp)
                    if not len(tmp):
                        for tr in st:
                            tr.stats.network = net
                        st = st.select(network=net, station=sta, component=comp)
                    else:
                        st = tmp
                    for tr in st:
                        tr.data = tr.data.astype(np.float)
                    stream += st
                    del st
                stream.sort()
                stream.merge(method=1, interpolation_samples=3, fill_value=None)
                stream = stream.split()

                logging.debug("Checking sample alignment")
                for i, trace in enumerate(stream):
                    stream[i] = check_and_phase_shift(trace)

                logging.debug("Checking Gaps")
                if len(getGaps(stream)) > 0:
                    max_gap = 10
                    only_too_long = False
                    while getGaps(stream) and not only_too_long:
                        too_long = 0
                        gaps = getGaps(stream)
                        for gap in gaps:
                            if int(gap[-1]) <= max_gap:
                                stream[gap[0]] = stream[gap[0]].__add__(stream[gap[1]], method=1,
                                                                        fill_value="interpolate")
                                stream.remove(stream[gap[1]])
                                break
                            else:
                                too_long += 1
                        if too_long == len(gaps):
                            only_too_long = True
                stream = stream.split()
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

                if get_config(db, 'remove_response', isbool=True):
                    logging.debug('%s Removing instrument response'%stream[0].id)
                    response_prefilt = eval(get_config(db, 'response_prefilt'))

                    response = responses[responses["channel_id"] == stream[0].id]
                    if len(response) > 1:
                        response = response[response["start_date"]<UTCDateTime(gd)]
                        response = response[response["end_date"]>UTCDateTime(gd)]
                    elif len(response) == 0:
                        logging.info("No instrument response information "
                                     "for %s, exiting" % stream[0].id)
                        sys.exit()
                    datalesspz = response["paz"].values[0]
                    stream.simulate(paz_remove=datalesspz,
                                    remove_sensitivity=True,
                                    pre_filt=response_prefilt,
                                    paz_simulate=None, )
                for tr in stream:
                    tr.data = tr.data.astype(np.float32)
                output += stream
                del stream
            del files
    clean_scipy_cache()
    return 0, output
