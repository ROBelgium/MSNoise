"""MSNoise is capable of using a reference function defined by absolute or
relative dates span. For example, an absolute range could be "from 1 January
2010 to 31 December 2011" and a relative range could be "the last 200 days".
In the latter case, the REF will need to be exported at every run, meaning the
following steps (MWCS and DTT) will be executed on the whole configured period.
If the REF is defined between absolute dates, excluding "today", the MWCS and
DTT will only be calculated for new data (e.g. "yesterday" and "today").
The corresponding configuration bits are ``ref_begin`` and ``ref_end``. In the
future, we plan on allowing multiple references to be defined.

Only data for new/modified dates need to be exported. If any CC-job has been
marked "Done" within the last day and triggered the creation of STACK jobs,
the stacks will be calculated and a new MWCS job will be inserted in the
database. For dates in the period of interest, the moving-window stack will
only be exported if new/modified CCF is available.
The export directory are "REF/" and "DAY%03i/" where %03i will be replaced by
the number of days stacked together (DAYS_005 for a 5-days stack, e.g.).

Please note that within MSNoise, stacks are always *inclusive* of the time/day
mentioned. For example, a 5-days stack on January 10, will contain
cross-correlation functions computed for January 6, 7, 8, 9 AND 10!
The graphical representation centered on a "January 10" tick might then display
changes in the CCF that occurred *on* the 10th !

Moving-window stacks are configured using the ``mov_stack`` parameter in
``msnoise admin``.

If ``stack_method`` is 'linear', then a simple mean CFF of all daily is saved
as the mov or ref CCF. On the other hand, if ``stack_method`` is 'pws', then
all the Phase Weighted Stack (PWS) is computed and saved as the mov or ref CCF.
The PWS is done in two steps: first the mean coherence between the instantaneous
phases of all windows is calculated, and eventually serves a weighting factor
on the mean. The smoothness of this weighting array is defined using the
``pws_timegate`` parameter in the configuration. The weighting array is the
power of the mean coherence array. If ``pws_power`` is equal to 0, a linear
stack is done (then it's faster to do set ``stack_method`` = 'linear'). Usual
value is 2.

.. warning:: PWS is largely untested, not cross-validated. It looks good, but
    that doesn't mean a lot, does it? Use with Caution! And if you
    cross-validate it, please let us know!!

    Schimmel, M. and Paulssen H., "Noise reduction and detection
    of weak, coherent signals through phase-weighted stacks". Geophysical Journal
    International 130, 2 (1997): 497-505.

Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

* |ref_begin|
* |ref_end|
* |mov_stack|
* |stack_method| | *new in 1.4*
* |pws_timegate| | *new in 1.4*
* |pws_power| | *new in 1.4*
* |hpc| | *new in 1.6*


Once done, each job is marked "D"one in the database and, unless ``hpc`` is
``Y``, MWCS jobs are inserted/updated in the database.

Usage:
~~~~~~

.. include:: ../clickhelp/msnoise-stack.rst


For most users, the REF stack will need to be computed only once for specific
dates and then, on routine basis, only compute the MOV stacks:

.. warning With MSNoise 1.6, we have split the two actions, and the REF
    stacks need to be computed first ! This process will put the corresponding
    STACK jobs "I"n progress and you will need to reset them before running the
    MOV stacks.

.. code-block:: sh

    $ msnoise stack -r
    $ msnoise reset STACK
    $ msnoise stack -m

as for all other steps, this procedure can be run in parallel:

.. code-block:: sh

    $ msnoise -t 4 stack -r
    $ msnoise reset STACK
    $ msnoise -t 4 stack -m


.. versionadded:: 1.4
    The Phase Weighted Stack.

.. versionadded:: 1.6
    The ``hpc`` parameter that can prevent the automatic creation of MWCS jobs.
    The REF and MOV stacks have been separated and need to be run independently.
"""

import argparse

import scipy.signal

from .api import *

import logbook


def main(stype, interval=1.0, loglevel="INFO"):
    """Computes the REF/MOV stacks.

    Parameters
    ----------
    stype : {'mov', 'ref'}
        Defines which of the REF or Moving-window stacks must be exported
    interval : float, optional
        Number of days before now to search for modified CC jobs

    """
    logger = logbook.Logger(__name__)
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.stack_child', loglevel,
                        with_pid=True)
    logger.debug('Starting the %s stack' % stype)
    db = connect()

    start, end, datelist = build_movstack_datelist(db)

    params = get_params(db)
    mov_stack = params.mov_stack
    if mov_stack.count(',') == 0:
        mov_stacks = [int(mov_stack), ]
    else:
        mov_stacks = [int(mi) for mi in mov_stack.split(',')]
    if 1 in mov_stacks:
        mov_stacks.remove(1)  # remove 1 day stack, it should exist already

    filters = get_filters(db, all=False)
    while is_dtt_next_job(db, flag='T', jobtype='STACK'):
        jobs = get_dtt_next_job(db, flag='T', jobtype='STACK')

        if not len(jobs):
            # edge case, should only occur when is_next returns true, but
            # get_next receives no jobs (heavily parallelised calls).
            time.sleep(np.random.random())
            continue
        pair = jobs[0].pair
        refs, days = zip(*[[job.ref, job.day] for job in jobs])

        print(
            "There are STACKS jobs for some days to recompute for %s" % pair)
        sta1, sta2 = pair.split(':')
        for f in filters:
            filterid = int(f.ref)
            for components in params.all_components:
                pair = "%s:%s" % (sta1, sta2)
                sta1 = sta1
                sta2 = sta2
                print('Processing %s-%s-%i' %
                      (pair, components, filterid))
                c = get_results(db, sta1, sta2, filterid, components, datelist,
                                mov_stack=1, format="dataframe")
                fn = r"STACKS2\%02i\001_DAYS\%s\%s_%s.h5" % (
                    filterid, components, sta1, sta2)
                if not os.path.isdir(os.path.split(fn)[0]):
                    os.makedirs(os.path.split(fn)[0])
                c.to_hdf(fn, key="CCF")
                for mov_stack in mov_stacks:
                    tmp = c.rolling("%iD" % mov_stack).mean()
                    fn = r"STACKS2\%02i\%03i_DAYS\%s\%s_%s.h5" % (
                        filterid, mov_stack, components, sta1, sta2)
                    if not os.path.isdir(os.path.split(fn)[0]):
                        os.makedirs(os.path.split(fn)[0])
                    tmp.to_hdf(fn, key="CCF")

