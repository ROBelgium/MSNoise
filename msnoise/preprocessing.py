import sys
import calendar
import glob
import time
import traceback

from obspy import read_inventory
from obspy.core import utcdatetime, UTCDateTime
from obspy.io.xseed import Parser

try:
    from scikits.samplerate import resample
except:
    pass

from .api import *


def preprocess(db, stations, comps, goal_day, params, tramef_Z, tramef_E=np.array([]), tramef_N=np.array([])):
    datafilesZ = {}
    datafilesE = {}
    datafilesN = {}

    for station in stations:
        datafilesZ[station] = []
        datafilesE[station] = []
        datafilesN[station] = []
        net, sta = station.split('.')
        gd = datetime.datetime.strptime(goal_day, '%Y-%m-%d')
        files = get_data_availability(
            db, net=net, sta=sta, starttime=gd, endtime=gd)
        for file in files:
            comp = file.comp
            fullpath = os.path.join(file.path, file.file)
            if comp[-1] == 'Z':
                datafilesZ[station].append(fullpath)
            elif comp[-1] == 'E':
                datafilesE[station].append(fullpath)
            elif comp[-1] == 'N':
                datafilesN[station].append(fullpath)

    j = 0
    for istation, station in enumerate(stations):
        for comp in comps:
            files = eval("datafiles%s['%s']" % (comp, station))
            if len(files) != 0:
                logging.debug("%s.%s Reading %i Files" %
                              (station, comp, len(files)))
                stream = Stream()
                for file in sorted(files):
                    st = read(file, dytpe=np.float,
                              starttime=UTCDateTime(gd),
                              endtime=UTCDateTime(gd) + 86400)
                    for tr in st:
                        tr.data = tr.data.astype(np.float)
                    stream += st
                    del st

                logging.debug("Checking sample alignment")
                for i, trace in enumerate(stream):
                    stream[i] = check_and_phase_shift(trace)

                stream.sort()
                logging.debug("Checking Gaps")
                if len(getGaps(stream)) > 0:
                    max_gap = 10
                    only_too_long = False
                    while getGaps(stream) and not only_too_long:
                        too_long = 0
                        gaps = getGaps(stream)
                        for gap in gaps:
                            if int(gap[-1]) <= max_gap:
                                stream[gap[0]] = stream[gap[0]].__add__(stream[gap[1]], method=0,
                                                                        fill_value="interpolate")
                                stream.remove(stream[gap[1]])
                                break
                            else:
                                too_long += 1
                        if too_long == len(gaps):
                            only_too_long = True

                taper_length = 20.0  # seconds
                for trace in stream:
                    if trace.stats.npts < 4 * taper_length * trace.stats.sampling_rate:
                        trace.data = np.zeros(trace.stats.npts)
                    else:
                        trace.detrend(type="demean")
                        trace.detrend(type="linear")
                        taper_1s = taper_length * float(trace.stats.sampling_rate) / trace.stats.npts
                        cp = cosine_taper(trace.stats.npts, taper_1s)
                        trace.data *= cp
                try:
                    stream.merge(method=0, fill_value=0.0)
                except:
                    continue

                logging.debug("%s.%s Slicing Stream to %s:%s" % (station, comp, utcdatetime.UTCDateTime(
                    goal_day.replace('-', '')), utcdatetime.UTCDateTime(
                    goal_day.replace('-', '')) + params.goal_duration - stream[0].stats.delta))
                stream[0].trim(utcdatetime.UTCDateTime(goal_day.replace('-', '')), utcdatetime.UTCDateTime(
                    goal_day.replace('-', '')) + params.goal_duration - stream[0].stats.delta, pad=True, fill_value=0.0,
                               nearest_sample=False)

                if get_config(db, 'remove_response', isbool=True):
                    logging.debug('Removing instrument response')
                    response_format = get_config(db, 'response_format')
                    response_prefilt = eval(get_config(db, 'response_prefilt'))
                    files = glob.glob(os.path.join(get_config(db,
                                                              'response_path'),
                                                   "*"))
                    if response_format == "inventory":
                        firstinv = True
                        inventory = None
                        for file in files:
                            try:
                                inv = read_inventory(file)
                                if firstinv:
                                    inventory = inv
                                    firstinv = False
                                else:
                                    inventory += inv
                            except:
                                traceback.print_exc()
                                pass
                        if inventory:
                            stream.attach_response(inventory)
                            stream.remove_response(output='VEL',
                                                   pre_filt=response_prefilt)
                    elif response_format == "dataless":
                        for file in files:
                            p = Parser(file)
                            try:
                                p.getPAZ(stream[0].id,
                                         datetime=UTCDateTime(gd))
                                break
                            except:
                                traceback.print_exc()
                                del p
                                continue
                        stream.simulate(seedresp={'filename': p, "units": "VEL"},
                                        pre_filt=response_prefilt,
                                        paz_remove=None,
                                        paz_simulate=None, )
                    elif response_format == "paz":
                        msg = "Unexpected type for `response_format`: %s" % \
                              response_format
                        raise TypeError(msg)
                    elif response_format == "resp":
                        msg = "Unexpected type for `response_format`: %s" % \
                              response_format
                        raise TypeError(msg)
                    else:
                        msg = "Unexpected type for `response_format`: %s" % \
                              response_format
                        raise TypeError(msg)
                trace = stream[0]

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
                        trace.data = trace.data[::decimation_factor]

                    elif params.resampling_method == "Lanczos":
                        logging.debug("%s.%s Downsample to %.1f Hz" %
                                      (station, comp, params.goal_sampling_rate))
                        trace.data = np.array(trace.data)
                        trace.interpolate(method="lanczos", sampling_rate=params.goal_sampling_rate, a=1.0)

                    trace.stats.sampling_rate = params.goal_sampling_rate


                year, month, day, hourf, minf, secf, wday, yday, isdst = trace.stats.starttime.utctimetuple()

                if j == 0:
                    t = time.strptime("%04i:%02i:%02i:%02i:%02i:%02i" %
                                      (year, month, day, hourf, minf, secf), "%Y:%m:%d:%H:%M:%S")
                    basetime = calendar.timegm(t)

                if len(trace.data) % 2 != 0:
                    trace.data = np.append(trace.data, 0.)
                if len(trace.data) != len(tramef_Z[istation]):
                    missing = len(tramef_Z[istation]) - len(trace.data)
                    for i in range(missing):
                        trace.data = np.append(trace.data, 0.)
                if comp == "Z":
                    tramef_Z[istation] = trace.data
                elif comp == "E":
                    tramef_E[istation] = trace.data
                elif comp == "N":
                    tramef_N[istation] = trace.data

                del trace, stream
    if len(tramef_E) != 0:
        return basetime, tramef_Z, tramef_E, tramef_N
    else:
        return basetime, tramef_Z
