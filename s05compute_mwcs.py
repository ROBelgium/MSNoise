"""
.. warning:: if using only ``mov_stack`` = 1, no DTT jobs is inserted in the
    database and consequently, no MWCS calculation will be done! FIX!
    
Following Clarke et al (2011), we apply the :ref:`mwcs`
to study the relative dephasing between Moving-Window stacks ("Current") and a
Reference using Moving-Window Cross-Spectral analysis. The *jobs* "T"o do have
been inserted in the datavase during the stack procedure.

In short, both time series are sliced in several overlapping windows and
preprocessed. The similarity of the two time-series is assessed using the
cross-coherence between energy densities in the frequency domain. The time
delay between the two cross correlations is found in the unwrapped phase of
the cross spectrum and is linearly proportional to frequency. This "Delay" for
each window between two signals is the slope of a weighted linear regression
(WLS) of the samples within the frequency band of interest.

For each filter, the frequency band can be configured using ``mwcs_low``
and ``mwcs_high``, and the window and overlap lengths using ``mwcs_wlen`` and
``mwcs_overlap``.

The output of this process is a table of delays measured at each window in the
functions. The following is an example for lag times between -115 and -90.
In this case, the window length was 10 seconds with an overlap of 5 seconds.

.. code-block:: python

            LAG_TIME                    DELAY                   ERROR                  MEAN COHERENCE
    -1.150000000000000000e+02 -1.478114638393185076e-01 5.372711913598642725e-02 2.758524391177089585e-01
    -1.100000000000000000e+02 -6.820752699231921734e-02 2.054664431127230587e-02 3.162099935286825647e-01
    -1.050000000000000000e+02 -1.033702957794639388e-01 8.664515540289926751e-03 4.243926988068138506e-01
    -1.000000000000000000e+02 -2.866877569687806965e-02 6.252221598805602840e-03 5.715984952889150428e-01
    -9.500000000000000000e+01  4.180394100884333997e-02 1.510228578922533961e-02 4.123855778998536947e-01
    -9.000000000000000000e+01  4.813940023363388887e-02 3.270065701850408124e-02 3.058618779249135389e-01

This process is job-based, so it is possible to run several instances in
parallel.

To run this script:

.. code-block:: sh

    python s05compute_mwcs.py

"""


from obspy.core import read

from database_tools import *
from MWCS import mwcs
import logging


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG,
                        filename="./compute_mwcs.log",
                        format='%(asctime)s [%(levelname)s] %(message)s',
                        filemode='w')
    
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)
    
    logging.info('*** Starting: Compute MWCS ***')
    
    db = connect()
    components_to_compute = get_components_to_compute(db)
    
    mov_stack = get_config(db, "mov_stack")
    if mov_stack.count(',') == 0:
        mov_stacks = [int(mov_stack), ]
    else:
        mov_stacks = [int(mi) for mi in mov_stack.split(',')]
    
    goal_sampling_rate = float(get_config(db, "cc_sampling_rate"))
    maxlag = float(get_config(db, "maxlag"))
    start, end, datelist = build_movstack_datelist(db)
    # allow_large_concats(db)
    
    # First we reset all DTT jobs to "T"odo if the REF is new for a given pair
    
    for station1, station2 in get_station_pairs(db, used=True):
        sta1 = "%s.%s" % (station1.net, station1.sta)
        sta2 = "%s.%s" % (station2.net, station2.sta)
        pair = "%s:%s" % (sta1, sta2)
        if is_dtt_next_job(db, type='DTT', ref=pair):
            logging.info(
                "We will recompute all MWCS based on the new REF for %s" % pair)
            reset_dtt_jobs(db, pair)
            update_job(db, "REF", pair, type='DTT', flag='D')
    
    # Then we compute the jobs
    while is_dtt_next_job(db, flag='T', type='DTT'):
        pair, days, refs = get_dtt_next_job(db, flag='T', type='DTT')
        logging.info(
            "There are MWCS jobs for some days to recompute for %s" % pair)
        for day in days:
            for f in get_filters(db, all=False):
                filterid = int(f.ref)
                for components in components_to_compute:
                    ref_name = pair.replace('.', '_').replace(':', '_')
                    rf = os.path.join("STACKS", "%02i" %
                                      filterid, "REF", components, ref_name + ".MSEED")
                    ref = read(rf)[0].data
                    for mov_stack in mov_stacks:
                        df = os.path.join(
                            "STACKS", "%02i" % filterid, "%03i_DAYS" %
                            mov_stack, components, ref_name, str(day) + ".MSEED")
                        if os.path.isfile(df):
                            cur = read(df)[0].data
                            logging.debug(
                                'Processing MWCS for: %s.%s.%02i - %s - %02i days' %
                                (ref_name, components, filterid, day, mov_stack))
                            output = mwcs(
                                cur, ref, f.mwcs_low, f.mwcs_high, goal_sampling_rate, -maxlag, f.mwcs_wlen, f.mwcs_step)
                            outfolder = os.path.join(
                                'MWCS', "%02i" % filterid, "%03i_DAYS" % mov_stack, components, ref_name)
                            if not os.path.isdir(outfolder):
                                os.makedirs(outfolder)
                            np.savetxt(
                                os.path.join(outfolder, "%s.txt" % str(day)), output)
                            del output
            update_job(db, day, pair, type='DTT', flag='D')
    
    logging.info('*** Finished: Compute MWCS ***')
