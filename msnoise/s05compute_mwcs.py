"""
.. warning:: if using only ``mov_stack`` = 1, no DTT jobs is inserted in the
    database and consequently, no MWCS calculation will be done! FIX!
    
Following Clarke et al (2011), we apply the :ref:`mwcs`
to study the relative dephasing between Moving-Window stacks ("Current") and a
Reference using Moving-Window Cross-Spectral analysis. The *jobs* "T"o do have
been inserted in the datavase during the stack procedure.


Filter Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``mwcs_low``: The lower frequency bound of the linear regression done in
  MWCS (in Hz)
* ``mwcs_high``: The upper frequency bound of the linear regression done in
  MWCS (in Hz)
* ``mwcs_wlen``: Window length (in seconds) to perform MWCS
* ``mwcs_step``: Step (in seconds) of the windowing procedure in MWCS


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

          LAG_TIME          DELAY           ERROR         MEAN COHERENCE
    -1.1500000000e+02 -1.4781146383e-01 5.3727119135e-02 2.7585243911e-01
    -1.1000000000e+02 -6.8207526992e-02 2.0546644311e-02 3.1620999352e-01
    -1.0500000000e+02 -1.0337029577e-01 8.6645155402e-03 4.2439269880e-01
    -1.0000000000e+02 -2.8668775696e-02 6.2522215988e-03 5.7159849528e-01
    -9.5000000000e+01  4.1803941008e-02 1.5102285789e-02 4.1238557789e-01
    -9.0000000000e+01  4.8139400233e-02 3.2700657018e-02 3.0586187792e-01

This process is job-based, so it is possible to run several instances in
parallel.

To run this step:

.. code-block:: sh

    $ msnoise compute_mwcs

This step also supports parallel processing/threading:

.. code-block:: sh

    $ msnoise -t 4 compute_mwcs

will start 4 instances of the code (after 1 second delay to avoid database
conflicts). This works both with SQLite and MySQL but be aware problems
could occur with SQLite.

.. versionadded:: 1.4
    Parallel Processing
"""

import logging
from obspy.core import read

from .api import *
from .MWCS import mwcs



def main():
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s [%(levelname)s] %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    logging.info('*** Starting: Compute MWCS ***')
    
    db = connect()
    components_to_compute = get_components_to_compute(db)
    export_format = get_config(db, 'export_format')
    if export_format == "BOTH":
        extension = ".MSEED"
    else:
        extension = "."+export_format
    mov_stack = get_config(db, "mov_stack")
    if mov_stack.count(',') == 0:
        mov_stacks = [int(mov_stack), ]
    else:
        mov_stacks = [int(mi) for mi in mov_stack.split(',')]
    
    goal_sampling_rate = float(get_config(db, "cc_sampling_rate"))
    maxlag = float(get_config(db, "maxlag"))

    # First we reset all DTT jobs to "T"odo if the REF is new for a given pair
    for station1, station2 in get_station_pairs(db, used=True):
        sta1 = "%s.%s" % (station1.net, station1.sta)
        sta2 = "%s.%s" % (station2.net, station2.sta)
        pair = "%s:%s" % (sta1, sta2)
        if is_dtt_next_job(db, jobtype='DTT', ref=pair):
            logging.info(
                "We will recompute all MWCS based on the new REF for %s" % pair)
            reset_dtt_jobs(db, pair)
            update_job(db, "REF", pair, jobtype='DTT', flag='D')
    
    logging.debug('Ready to compute')
    # Then we compute the jobs
    outfolders = []
    while is_dtt_next_job(db, flag='T', jobtype='DTT'):
        pair, days, refs = get_dtt_next_job(db, flag='T', jobtype='DTT')
        logging.info(
            "There are MWCS jobs for some days to recompute for %s" % pair)
        for f in get_filters(db, all=False):
            filterid = int(f.ref)
            for components in components_to_compute:
                ref_name = pair.replace('.', '_').replace(':', '_')
                rf = os.path.join("STACKS", "%02i" %
                                  filterid, "REF", components,
                                  ref_name + extension)
                ref = read(rf)[0].data
                for day in days:
                    for mov_stack in mov_stacks:
                        df = os.path.join(
                            "STACKS", "%02i" % filterid, "%03i_DAYS" %
                            mov_stack, components, ref_name,
                            str(day) + extension)
                        if os.path.isfile(df):
                            cur = read(df)[0].data
                            logging.debug(
                                'Processing MWCS for: %s.%s.%02i - %s - %02i days' %
                                (ref_name, components, filterid, day, mov_stack))
                            output = mwcs(
                                cur, ref, f.mwcs_low, f.mwcs_high, goal_sampling_rate, -maxlag, f.mwcs_wlen, f.mwcs_step)
                            outfolder = os.path.join(
                                'MWCS', "%02i" % filterid, "%03i_DAYS" % mov_stack, components, ref_name)
                            if outfolder not in outfolders:
                                if not os.path.isdir(outfolder):
                                    os.makedirs(outfolder)
                                outfolders.append(outfolder)
                            np.savetxt(os.path.join(outfolder, "%s.txt" % str(day)), output)
                            del output, cur
        for day in days:
            update_job(db, day, pair, jobtype='DTT', flag='D')
    
    logging.info('*** Finished: Compute MWCS ***')

if __name__ == "__main__":
    main()