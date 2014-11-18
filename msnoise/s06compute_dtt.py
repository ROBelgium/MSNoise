"""
This code is responsible for the calculation of dt/t using the result of the
MWCS calculations.

.. warning:: Currently, all pairs are analysed using the same parameters, which
    are hard-coded in this file. This should change in the future, as different
    inter-station distances should/could be treated with different parameters.

The dt/t is determined as the slope of the delays vs time lags. The slope is
calculated a weighted linear regression (WLS) through selected points. The selection
of points is based on criteria:

* ``lMlag``: Max time lag on the Left side (ex: -30)
* ``lmlag``: Min time lag on the left side (ex: -10)
* ``rmlag``: Min time lag on the right side (ex: 10)
* ``rMlag``: Max time lag on the right side (ex: 30)
* ``minCoh``: Minimum coherence (ex: 0.5)
* ``maxErr``: Maximum error on the delay value (ex:0.1)
* ``maxDt``: Maximum absolute delay value (ex: 0.5)

Using example values above, we chose to use only 10-30 s coda part of the
signal, neglecting direct waves in the 0-10 seconds range. Graphically, this
data selection looks like:

.. image:: .static/Figure04_dttmatrix_01_005DAYS_ZZ-2010-10-12_cmyk.png

Each of the 4 left subplot of this figure shows a colormapper matrix of which each row
corresponds to the data of 1 station pair and each column corresponds to
different time lags. The cells are then colored using, from left to right:
Delays, Errors, Phase Coherence and Data Selection. 

Once data (cells) have been selected, they are analyzed two times: first using
a WLS that is forced to pass the origin (0,0) and second when a constant is
added to allow for the WLS to be offset from the origin. For each value, the
error is computed and stored. M0 and EM0 are the slope and its error for the
first WLS, and M, EM together with A and EA are the slope, its error, the
constant and its error for the second WLS. The output of this calculation
is a table, with one row for each station pair.


.. code-block:: python

    Date,               A,                  EA,                 EM,                     EM0,                M,                  M0,                 Pairs
    2013-01-06,-0.16837287183494098,0.05266069549188657,0.00208377243783003,0.000965214393762639, 0.0068202122517553850, 0.000377577217283868,BE_GES_BE_HOU
    2013-01-06,-0.00804644347723505,0.05779368290438131,0.00291327704047938,0.000972986173494356,-0.0022691002608228083,-0.002643541217653975,BE_GES_BE_MEM
    2013-01-06, 0.10074727166431066,0.01446482303784344,0.00179566184493113,0.004541720342080022,-0.0014573842288369784, 0.007414785762726537,BE_GES_BE_RCHB
    2013-01-06,-0.05568112085310946,0.00989268830665443,0.00057839595252212,0.001081021896160096,-0.0032896527406604500,-0.001360750598930064,BE_GES_BE_SKQ
    2013-01-06, 0.01508666636244094,0.02022437208406759,0.00096543261879024,0.000898329875688285, 0.0008371422712720330, 0.001045072240756187,BE_GES_BE_STI
    2013-01-06, 0.02683099320453874,0.03289970421407011,0.00153137352293737,0.001502611018541102, 0.0030233174200087437, 0.003024510068343402,BE_GES_BE_UCC
    2013-01-06,-0.01212934374358136,0.00433513329889883,0.00039019542460264,0.000413471992144575, 0.0002583635832057903,-0.000427090531110851,BE_HOU_BE_MEM
    2013-01-06, 0.10762476137331067,0.01886622600766513,0.00076824344993628,0.002163833958287798,-0.0003079139649344223, 0.001126925437222496,BE_HOU_BE_RCHB
    2013-01-06,-0.04684857637025102,0.01944924279920320,0.00069968472500258,0.000782078143903881,-0.0006613328679275146, 0.000271023281627172,BE_HOU_BE_SKQ
    2013-01-06, 0.02030579417237702,0.01613165714957967,0.00131522430756707,0.001311822061778742, 0.0005162635594570485,-3.103066113367655e-0,BE_HOU_BE_STI
    ...
    2013-01-06,-0.00225882255299991,0.00371411997707409,0.00010340935420880,9.19916693442618e-05, 0.0007363569140943659, 0.000762389133912101,ALL


To run this script:

.. code-block:: sh

    python s06compute_dtt.py 


Grouping Station Pairs
~~~~~~~~~~~~~~~~~~~~~~~
Although not clearly visible on the figure above, the very last row of the
matrix doesn't contain information about one station pair, but contains
a weighted mean of all delays (from all pairs) for each time lag. For each time
lag, delays from each pair is taken into account if it satisfies the same
criteria as for the individual data selection. Once the last row (the ALL line)
has been calculated, it goes through the normal process of the double WLS and
is saved to the output file, as visible above. In the future, MSNoise will be
able to treat as many groups as the user want, allowing, e.g. a "crater" and
a "slopes" groups.

Mean of All Pairs vs Mean Pair
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The dt/t calculated using the mean pair (ALL, in red on subplots 4 and 5)
and by calculating the weighted mean of the dt/t of all pairs (in green)
don't show a significant difference. The standard deviation around the latter
is more spread than on the former, but this has to be investigated.

Forcing vs No Forcing through Origin
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The reason for allowing the WLS to cross the axis elsewhere than on (0,0) is,
for example, to study the potential clock drifts or noise source position
variations.

"""

import numpy as np
import os

from database_tools import *
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import statsmodels.api as sm

import logging


logging.basicConfig(level=logging.DEBUG,
                    filename="/dev/null",
                    format='%(asctime)s [%(levelname)s] %(message)s',
                    filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)




def wavg_wstd(data, errors):
    d = data
    errors[errors == 0] = 1e-6
    w = 1. / errors
    wavg = (d * w).sum() / w.sum()
    N = len(np.nonzero(w)[0])
    wstd = np.sqrt(np.sum(w * (d - wavg) ** 2) / ((N - 1) * np.sum(w) / N))
    return wavg, wstd


def main():
    
    logging.info('*** Starting: Compute DT/T ***')
    db = connect()
    
    #~ PLOT= True
    PLOT = False
    
    dtt_lag = get_config(db, "dtt_lag")
    dtt_v = float(get_config(db, "dtt_v"))
    dtt_minlag = float(get_config(db, "dtt_minlag"))
    dtt_width = float(get_config(db, "dtt_width"))
    dtt_sides = get_config(db, "dtt_sides")
    minCoh = get_config(db, "dtt_mincoh")
    maxErr = get_config(db, "dtt_maxerr")
    maxDt = get_config(db, "dtt_maxdt")
    
    # lMlag = -100  # HARD
    # lmlag = -10  # HARD
    # rmlag = 10  # HARD
    # rMlag = 100  # HARD
    # minCoh = 0.5  # HARD
    # maxErr = 0.1  # HARD
    # maxDt = 1.0  # HARD
    # lagmethod="rel"
    
    start, end, datelist = build_movstack_datelist(db)
    
    mov_stack = get_config(db, "mov_stack")
    if mov_stack.count(',') == 0:
        mov_stacks = [int(mov_stack), ]
    else:
        mov_stacks = [int(mi) for mi in mov_stack.split(',')]
    
    components_to_compute = get_components_to_compute(db)
    updated_dtt = updated_days_for_dates(
        db, start, end, '%', type='DTT', returndays=True,interval=datetime.timedelta(days=0.1))
    
    for f in get_filters(db, all=False):
        filterid = int(f.ref)
        for components in components_to_compute:
            for mov_stack in mov_stacks:
                logging.info('Loading mov=%i days for filter=%02i' %
                             (mov_stack, filterid))
                for current in updated_dtt:
                    if current > datetime.date.today():
                        break
                    logging.debug("Processing %s - %02i - %02i mov" % (current, filterid, mov_stack))
                    first = True
                    for station1, station2 in get_station_pairs(db, used=True):
                        sta1 = "%s_%s" % (station1.net, station1.sta)
                        sta2 = "%s_%s" % (station2.net, station2.sta)
                        pair = "%s_%s" % (sta1, sta2)
                        day = os.path.join('MWCS', "%02i" % filterid, "%03i_DAYS" %
                                           mov_stack, components, pair, '%s.txt' % current)
                        if os.path.isfile(day):
                            #~ print day
                            df = pd.read_csv(
                                day, delimiter=' ', header=None, index_col=0, names=['t', 'dt', 'err', 'coh'])
                            tArray = df.index.values
                            if dtt_lag == "dynamic":
                                tindex = np.where(((tArray >= lMlag) & (tArray <= lmlag)) | (
                                (tArray >= rmlag) & (tArray <= rMlag)))[0]
                            else:
                                # TLQ: artificially set error to very high value where lag time < dtt_v propagation time
                                dist = get_interstation_distance(station1, station2)
                                lmlag = -dist * dtt_v
                                rmlag =  dist * dtt_v
                                lMlag = lmlag - dtt_width
                                rMlag = rmlag + dtt_width
                                if dtt_sides == "both":
                                    tindex = np.where(((tArray >= lMlag) & (tArray <= lmlag)) | ((tArray >= rmlag) & (tArray <= rMlag)))[0]
                                elif dtt_sides == "left":
                                    tindex = np.where((tArray >= lMlag) & (tArray <= lmlag))[0]
                                else:
                                    tindex = np.where((tArray >= rmlag) & (tArray <= rMlag))[0]
                            
                            tmp = np.setdiff1d(np.arange(len(tArray)),tindex)
                            df['err'][tmp] = 1.0
                            df['coh'][tmp] = 0.0
                            
                            #END TLQ
                            if first:
                                tArray = df.index.values
                                dtArray = df['dt']
                                errArray = df['err']
                                cohArray = df['coh']
                                pairArray = [pair, ]
                                first = False
                            else:
                                dtArray = np.vstack((dtArray, df['dt']))
                                errArray = np.vstack((errArray, df['err']))
                                cohArray = np.vstack((cohArray, df['coh']))
                                pairArray.append(pair)
                            del df
                        del day
    
                    if not first:
                        #~ tindex = np.tindwhere(((tArray >= lMlag) & (tArray <= lmlag)) | (
                            #~ (tArray >= rmlag) & (tArray <= rMlag)))[0]
    
                        Dates = []
                        Pairs = []
                        M = []
                        EM = []
                        A = []
                        EA = []
                        M0 = []
                        EM0 = []
                        if len(pairArray) != 1:
                            # first stack all pairs to a ALL mean pair, using
                            # indexes of selected values:
                            new_dtArray = np.zeros(len(tArray))
                            new_errArray = np.zeros(len(tArray)) + 9999
                            new_cohArray = np.zeros(len(tArray))
                            for i in range(len(tArray)):
                                #~ if i in tindex:
                                if 1:
                                    cohindex = np.where(
                                        cohArray[:, i] >= minCoh)[0]
                                    errindex = np.where(
                                        errArray[:, i] <= maxErr)[0]
                                    dtindex = np.where(
                                        np.abs(dtArray[:, i]) <= maxDt)[0]
    
                                    index = np.intersect1d(cohindex, errindex)
                                    index = np.intersect1d(index, dtindex)
    
                                    wavg, wstd = wavg_wstd(
                                        dtArray[:, i][index], errArray[:, i][index])
                                    new_dtArray[i] = wavg
                                    new_errArray[i] = wstd
                                    new_cohArray[i] = 1.0
    
                            dtArray = np.vstack((dtArray, new_dtArray))
                            errArray = np.vstack((errArray, new_errArray))
                            cohArray = np.vstack((cohArray, new_cohArray))
                            pairArray.append("ALL")
                            del new_cohArray, new_dtArray, new_errArray, cohindex, errindex, dtindex, wavg, wstd
                            
                            # then stack selected pais to GROUPS:
                            groups = {}
                            # groups['CRATER'] = ["CSS","FJS","BON","SNE","FLR","RVL","FOR",]
                            # groups['GPENTES'] = ["CRA","GPN","GPS","HDL","GBS"]
                            # groups['VOLCAN'] = groups['CRATER'] + groups['GPENTES'] + ['HIM','VIL']
                            # groups['NORD'] = ['PRO','MVL','RSL','MAT','PCR','OBS','CIL','CAM']
                            npairs = len(pairArray)-1
                            for group in groups.keys():
                                new_dtArray = np.zeros(len(tArray))
                                new_errArray = np.zeros(len(tArray)) + 9999
                                new_cohArray = np.zeros(len(tArray))
                                pairindex = []
                                for j, pair in enumerate(pairArray[:npairs]):
                                    net1, sta1, net2, sta2 = pair.split('_')
                                    if sta1 in groups[group] and sta2 in groups[group]:
                                        pairindex.append(j)
                                pairindex = np.array(pairindex)
                                #~ print "found", len(pairindex), " pairs for group", group
                                for i in range(len(tArray)):
                                    #~ if i in tindex:
                                    if 1:
                                        cohindex = np.where(
                                            cohArray[:, i] >= minCoh)[0]
                                        errindex = np.where(
                                            errArray[:, i] <= maxErr)[0]
                                        dtindex = np.where(
                                            np.abs(dtArray[:, i]) <= maxDt)[0]
        
                                        index = np.intersect1d(cohindex, errindex)
                                        index = np.intersect1d(index, dtindex)
                                        index = np.intersect1d(index, pairindex)
                                        
        
                                        wavg, wstd = wavg_wstd(
                                            dtArray[:, i][index], errArray[:, i][index])
                                        new_dtArray[i] = wavg
                                        new_errArray[i] = wstd
                                        new_cohArray[i] = 1.0
        
                                dtArray = np.vstack((dtArray, new_dtArray))
                                errArray = np.vstack((errArray, new_errArray))
                                cohArray = np.vstack((cohArray, new_cohArray))
                                pairArray.append(group)
                                del new_cohArray, new_dtArray, new_errArray, cohindex, errindex, dtindex, wavg, wstd
                                # END OF GROUP HANDLING
    
                        if PLOT:
                            plt.figure(figsize=(15, 10))
                            gs = gridspec.GridSpec(
                                1, 6, width_ratios=[3, 3, 3, 3, 3, 3])
                            plt.subplots_adjust(
                                hspace=0.005, wspace=0.15, bottom=0.05, top=0.92, left=0.07, right=0.95)
                            plt.subplot(gs[0])
                            ymax = dtArray.shape[0]
                            extent = (-117.5, 117.5, ymax + 0.5, -0.5)  # HARD
    
                            plt.imshow(
                                dtArray, interpolation='none', aspect='auto',
                                vmin=-maxDt, vmax=maxDt, cmap='bwr', extent=extent)
                            for item in [lMlag, lmlag, rmlag, rMlag]:
                                plt.axvline(item, ls='--', zorder=15, c='k', lw=2)
                            plt.title('Delays')
                            cb = plt.colorbar(
                                orientation='horizontal', pad=0.08, shrink=0.8)
                            cb.set_label('Delay (seconds)')
                            cb.set_ticks([-maxDt, 0, maxDt])
                            plt.xlabel('Time Lag')
                            plt.ylabel('Pair Number')
                            plt.subplot(gs[1])
                            plt.title('Errors')
                            plt.imshow(
                                errArray, interpolation='none', aspect='auto', vmax=0.1, extent=extent, cmap='hot')
                            for item in [lMlag, lmlag, rmlag, rMlag]:
                                plt.axvline(item, ls='--', zorder=15, c='k', lw=2)
                            plt.yticks([], [])
                            cb = plt.colorbar(
                                orientation='horizontal', pad=0.08, shrink=0.8)
                            cb.set_label('Error (seconds)')
                            cb.set_ticks([-0.1, 0, 0.1])
                            plt.xlabel('Time Lag')
                            plt.subplot(gs[2])
                            plt.title('Phase Coherences')
                            plt.imshow(
                                cohArray, interpolation='none', aspect='auto', extent=extent, cmap='hot')
                            for item in [lMlag, lmlag, rmlag, rMlag]:
                                plt.axvline(item, ls='--', zorder=15, c='k', lw=2)
                            cb = plt.colorbar(
                                orientation='horizontal', pad=0.08, shrink=0.8)
                            cb.set_label('Coherence Coefficient')
                            cb.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
                            plt.yticks([], [])
                            plt.xlabel('Time Lag')
    
                        # then process all pairs + the ALL
                        if len(dtArray.shape) == 1:  # if there is only one pair:
                            dtArray = dtArray.reshape((1, dtArray.shape[0]))
                            cohArray = cohArray.reshape((1, cohArray.shape[0]))
                            errArray = errArray.reshape((1, errArray.shape[0]))
    
                        used = np.zeros(dtArray.shape)
    
                        for i, pair in enumerate(pairArray):
                            cohindex = np.where(cohArray[i] >= minCoh)[0]
                            errindex = np.where(errArray[i] <= maxErr)[0]
                            dtindex = np.where(np.abs(dtArray[i]) <= maxDt)[0]
    
                            #~ index = np.intersect1d(tindex, cohindex)
                            index = np.intersect1d(cohindex, errindex)
                            index = np.intersect1d(index, dtindex)
    
                            used[i][index] = 1.0
    
                            w = 1.0 / (errArray[i][index] ** 2)
    
                            VecXfilt = tArray[index]
                            VecYfilt = dtArray[i][index]
                            if len(VecYfilt) >= 2:
                                B = sm.tools.tools.add_constant(
                                    VecXfilt, prepend=False)
                                res = sm.regression.linear_model.WLS(
                                    VecYfilt, B, w).fit()
                                res0 = sm.regression.linear_model.WLS(
                                    VecYfilt, VecXfilt, w).fit()
                                if res.df_resid > 0:
                                    m, a = res.params
                                    em, ea = res.bse
    
                                    m0 = res0.params[0]
                                    em0 = res0.bse[0]
    
                                    M.append(m)
                                    EM.append(em)
                                    A.append(a)
                                    EA.append(ea)
    
                                    M0.append(m0)
                                    EM0.append(em0)
    
                                    Dates.append(current)
                                    Pairs.append(pair)
    
                                    del m, a, em, ea, m0, em0
                                del res, res0, B
                            del VecXfilt, VecYfilt, w
                            del index, cohindex, errindex, dtindex
    
                        if PLOT:
                            plt.subplot(gs[3])
                            plt.title('Data Selection')
                            plt.imshow(
                                used, extent=extent, cmap='binary', interpolation='none', aspect='auto')
                            cb = plt.colorbar(
                                orientation='horizontal', pad=0.08, shrink=0.8)
                            cb.set_label('Selected')
                            cb.set_ticks([0.0, 1.0],)
                            cb.set_ticklabels(['No', 'Yes'])
                            plt.yticks([], [])
                            plt.xlabel('Time Lag')
    
                            plt.subplot(gs[4])
                            plt.title('Delay Time Variations\nNo Forcing')
                            Ma = np.array(M)
                            EMa = np.array(EM)
                            plt.errorbar(
                                100 * Ma, np.arange(len(Ma)), xerr=EMa * 100, c='k', lw=0, elinewidth=2.0)
                            plt.scatter(100 * Ma, np.arange(len(Ma)), c=Ma * 100,
                                        edgecolor='k', vmin=-0.2, vmax=0.2, zorder=100, s=75)
    
                            cb = plt.colorbar(
                                orientation='horizontal', pad=0.08, shrink=0.8,)
                            cb.set_ticks([-0.2, 0.2],)
                            cb.set_label('dt/t in %')
                            plt.axvline(0, c='k')
                            plt.axvline(
                                100 * Ma[-1], c='r', lw=2, label="ALL: %.2e" % Ma[-1])
                            plt.axvspan(100 * Ma[-1] - 100 * EM[-1], 100 * Ma[
                                        -1] + 100 * EM[-1], color='r', alpha=0.2)
                            # print "M,EM", 100*M[-1],100*EM[-1]
    
                            wdtt, werr = wavg_wstd(
                                np.array(Ma[:-1]), np.array(EMa[:-1]))
                            plt.axvline(
                                100 * wdtt, c='g', lw=3, label="Weighted $\overline{x}$: %.2e" % wdtt)
                            plt.axvspan(
                                100 * wdtt - 100 * werr, 100 * wdtt + 100 * werr, color='g', alpha=0.2)
                            # print 'wdtt, werr', 100*wdtt, 100*werr
    
                            plt.xlim(-0.2, 0.2)
                            plt.xticks([-0.2, 0.0, 0.2])
                            plt.xlabel('dt/t in %')
                            plt.ylim(len(Ma), 0)
                            plt.yticks([], [])
    
                            plt.subplot(gs[5])
                            plt.title('Delay Time Variations\nForcing Origin')
                            Ma = np.array(M0)
                            EMa = np.array(EM0)
                            plt.errorbar(
                                100 * Ma, np.arange(len(Ma)), xerr=EMa * 100, c='k', lw=0, elinewidth=2.0)
                            plt.scatter(100 * Ma, np.arange(len(Ma)), c=Ma * 100,
                                        edgecolor='k', vmin=-0.2, vmax=0.2, zorder=100, s=75)
    
                            cb = plt.colorbar(
                                orientation='horizontal', pad=0.08, shrink=0.8,)
                            cb.set_ticks([-0.2, 0.2],)
                            cb.set_label('dt/t in %')
                            plt.axvline(0, c='k')
                            plt.axvline(
                                100 * Ma[-1], c='r', lw=2, label="ALL: %.2e" % Ma[-1])
                            plt.axvspan(100 * Ma[-1] - 100 * EM0[-1], 100 * Ma[
                                        -1] + 100 * EM0[-1], color='r', alpha=0.2)
                            # print "M0,EM0", 100*M0[-1],100*EM0[-1]
    
                            wdtt, werr = wavg_wstd(
                                np.array(Ma[:-1]), np.array(EMa[:-1]))
                            plt.axvline(
                                100 * wdtt, c='g', lw=3, label="Weighted $\overline{x}$: %.2e" % wdtt)
                            plt.axvspan(
                                100 * wdtt - 100 * werr, 100 * wdtt + 100 * werr, color='g', alpha=0.2)
                            # print 'wdtt0, werr0', 100*wdtt, 100*werr
    
                            plt.xlim(-0.2, 0.2)
                            plt.xticks([-0.2, 0.0, 0.2])
                            plt.xlabel('dt/t in %')
                            plt.ylim(len(Ma), 0)
                            plt.yticks([], [])
    
                            plt.savefig('dttmatrix_%02i_%03iDAYS_%s-%s.png' %
                                        (filterid, mov_stack, components, current), dpi=300)
    
                            #~ stations = get_stations(db, all=False).all()
                            #~ nstations = len(stations)
                            #~ dttmatrix = np.zeros((nstations, nstations))
                            #~ for station1, station2 in get_station_pairs(db, used=True):
                                #~ sta1 = "%s_%s" % (station1.net, station1.sta)
                                #~ sta2 = "%s_%s" % (station2.net, station2.sta)
                                #~ pair = "%s_%s" % (sta1, sta2)
                                #~ if pair in pairArray:
                                    #~ id = pairArray.index(pair)
                                    #~ dttmatrix[i, j] = M[id] * 100
                            #~ plt.figure(figsize=(8, 8))
                            #~ plt.subplots_adjust(
                                #~ hspace=0.005, wspace=0.1, bottom=0.05, top=0.90, left=0.1, right=0.90)
                            #~ plt.imshow(
                                #~ dttmatrix, interpolation='none', vmin=-0.2, vmax=0.2, cmap='bwr')
                            #~ cb = plt.colorbar(shrink=0.8, pad=0.05)
                            #~ cb.set_label('dt/t in %')
                            #~ plt.yticks(np.arange(nstations), [
                                       #~ "%s.%s" % (st.net, st.sta) for st in stations])
                            #~ plt.xticks(np.arange(nstations), ["%s.%s" % (st.net, st.sta)
                                       #~ for st in stations], rotation=90, verticalalignment='bottom')
                            #~ for tick in plt.gca().xaxis.iter_ticks():
                                #~ tick[0].label2On = True
                                #~ tick[0].label1On = False
                                #~ tick[0].label2.set_rotation('vertical')
                            #~ plt.grid()
                            #~ plt.savefig('dttmatrix2_%02i_%03iDAYS_%s-%s.png' %
                                        #~ (filterid, mov_stack, components, current), dpi=300)
                            plt.show()
    
                        logging.debug(
                            "%s: exporting: %i pairs" % (current, len(pairArray)))
                        df = pd.DataFrame(
                            {'Pairs': Pairs, 'M': M, 'EM': EM, 'A': A, 'EA': EA, 'M0': M0, 'EM0': EM0}, index=pd.DatetimeIndex(Dates))
                        # Needs to be changed !
                        output = os.path.join(
                            'DTT', "%02i" % filterid, "%03i_DAYS" % mov_stack, components)
                        if not os.path.isdir(output):
                            os.makedirs(output)
                        df.to_csv(
                            os.path.join(output, '%s.txt' % current), index_label='Date')
                        del df, M, EM, A, EA, M0, EM0, Pairs, Dates, used
                        del tArray, dtArray, errArray, cohArray, pairArray
                        del output
    
    
    logging.info('*** Finished: Compute DT/T ***')


if __name__ == "__main__":
    main()