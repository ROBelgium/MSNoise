'''
method - MWCS or STRETCH
mov_stack_idx (str or int)
components = 'ZZ, ZE, ..'
filterid (str or int)
bystation : median, mean or None (None plot all the pairs)
keep_nopair : True or False (remove the unused pairs)
show : True or False
outfile : Export the figure (somewhere) - Specify path

'''
    
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np
import pandas as pd
from scipy.signal import medfilt
from msnoise.api import *
from msnoise._version import *
import os
import warnings
from matplotlib.dates import DateFormatter
from datetime import datetime

def main(method='MWCS', mov_stack_idx=1, components='ZZ', filterid='01', bystation=None, keep_nopair=False, wiggles=False, show=True, outfile=None, loglevel="INFO") :
    logger = get_logger('msnoise.cc_dvv_plot_aurelogram', loglevel,with_pid=True)
    db = connect()
    params = get_params(db)
    start, end, datelist = build_movstack_datelist(db)
    mov_stacks = params.mov_stack 

    stations = get_stations(db)
    
    ### verifying the inputs
    if mov_stack_idx > len(mov_stacks) :
        print('Index of mov. stack is too high! Change the value')
        return
    Testpath = os.path.join("STR2","%02i" % int(filterid),"%s_%s" %(mov_stacks[int(mov_stack_idx)-1][0], mov_stacks[int(mov_stack_idx)-1][1])
                                , components)
    if os.path.isdir(Testpath) == False :
        print('The path : %s does not exist' %Testpath)
        return

    list_sta1 = np.array([], dtype='str')
    for sta in stations :
        list_sta1 = np.append(list_sta1, '%s.%s.%s' %(sta.net, sta.sta,sta.used_location_codes)) 

    pair_stations = get_station_pairs(db)
    list_pair = np.array([], dtype='str')

    ## Listing pairs and first stations ##
    for _, pair in enumerate(pair_stations) :
        list_pair = np.append(list_pair, '%s.%s.%s_%s.%s.%s' %(pair[0].net, pair[0].sta,pair[0].used_location_codes, pair[1].net, pair[1].sta,pair[1].used_location_codes))


    dts = pd.Timedelta(mov_stacks[int(mov_stack_idx)-1][1]).total_seconds()
    datearray = np.arange(np.datetime64(start,'s'), np.datetime64(end,'D')+1, int(dts))

    dvv_grid = np.full((list_pair.shape[0], datearray.shape[0]), np.nan, dtype='float') # creating the grid based on dt and numbers of pairs


    for i, pair in enumerate(list_pair):

        if method == 'STRETCH' :
            #path = './STR2/%02i/%s_%s/%s/%s.csv' %(int(filterid),str(mov_stacks[int(mov_stack_idx)-1][0]),str(mov_stacks[int(mov_stack_idx)-1][1]), str(components), pair)
            path = os.path.join("STR2","%02i" % int(filterid),"%s_%s" %(mov_stacks[int(mov_stack_idx)-1][0], mov_stacks[int(mov_stack_idx)-1][1])
                                , components, pair+'.csv')
            if os.path.isfile(path) == False : 
                continue
            data = pd.read_csv(path, index_col=0, parse_dates=True).iloc[:,0]
            timing = np.array(data.index, dtype='datetime64')
            s = (data-1.)*100

        if method == 'MWCS' :
            stt1, stt2 = pair.split('_')[0],pair.split('_')[1]
            try : data = xr_get_dtt(stt1, stt2, components,int(filterid),mov_stacks[int(mov_stack_idx)-1])
            except : 
                warnings.warn('No dtt data for the component : %s, filterid : %02i and mov_stack : %s' %(components,int(filterid),str(mov_stacks[int(mov_stack_idx)-1][0])))
                continue
            timing = np.array(data.index, dtype='datetime64')
            s = data.m * -100

        startidx = np.argmin(np.abs(datearray - timing[0])) 
        s = np.array(s, dtype='float')
        dvv_grid[i, startidx:startidx+s.shape[0]] = s

    if bystation == 'median' :
        dvv_by_station = np.full((list_sta1.shape[0], dvv_grid.shape[1]), np.nan)
        for j, sta in enumerate(list_sta1) :
            cond = np.char.startswith(list_pair, sta) + np.char.endswith(list_pair, sta)
            dvv_per_sta = np.where(cond.reshape((len(cond),1))==True, dvv_grid,np.nan)
            dvv_by_station[j,:] = np.nanmedian(dvv_per_sta, axis=0)

        dvv_grid, list_pair = dvv_by_station, list_sta1

    if bystation == 'mean' :
        dvv_by_station = np.full((list_sta1.shape[0], dvv_grid.shape[1]), np.nan)
        for j, sta in enumerate(list_sta1) :
            cond = np.char.startswith(list_pair, sta) + np.char.endswith(list_pair, sta)
            dvv_per_sta = np.where(cond.reshape((len(cond),1))==True, dvv_grid,np.nan)
            dvv_by_station[j,:] = np.nanmean(dvv_per_sta, axis=0)

        dvv_grid, list_pair = dvv_by_station, list_sta1

    elif keep_nopair == False :
        cond = np.sum(np.isfinite(dvv_grid), axis=1)
        condidx = np.extract(np.where(cond == 0, 0,1), np.arange(len(list_pair)))
        dvv_grid_nopair = np.full((condidx.shape[0], dvv_grid.shape[1]), np.nan, dtype='float')
        list_pair_nopair = []
        for j, index in enumerate(condidx) :
            dvv_grid_nopair[j,:] = dvv_grid[index,:]
            list_pair_nopair.append(list_pair[index])
        dvv_grid, list_pair = dvv_grid_nopair, list_pair_nopair

    ### plotting ###

    fig, ax = plt.subplots(figsize=(10,10), dpi=300)
    
    if wiggles == True :
        axb = ax.twinx()
        for z, dvvVal in enumerate(dvv_grid) :
            color='indianred'
            if (z//2 - z/2) == 0: 
                color = 'cornflowerblue'
                linewidth = .5
            if len(list_pair) > 30 :
                color = 'dimgrey'
                linewidth = .3
            ax.plot(datearray, z + (dvvVal), c=color, linewidth=linewidth)
            axb.plot(datearray, z +(dvvVal), c='b', linewidth=.3, ls='--', visible = False)
        
        ds = 3
        poslabel = np.arange(0,len(list_pair), ds)
        if len(list_pair) <= 30 :
            defticks = np.arange(0,len(list_pair), 3)
            posticks = np.array([defticks-.5,defticks,defticks+.5]).T.ravel()
            poslabel = [-.5,0,.5]*len(defticks)
            axb.set_yticks(posticks, labels=poslabel)
            axb.grid(ls='--')
            axb.set_ylabel('dv/v = -dt/t (%)')
        else :
            axb.set_yticks([])
        ax.set_ylim(-1, len(list_pair)+.5), axb.set_ylim(-1, len(list_pair)+.5)

    else : 
        norm = clr.Normalize(vmin=-np.nanmax(np.abs(dvv_grid)), vmax=np.nanmax(np.abs(dvv_grid)))
        TIMEPLT, PAIRIDX = np.meshgrid(datearray, np.arange(len(list_pair)))
        pclr = ax.pcolormesh(TIMEPLT, PAIRIDX, dvv_grid, norm=norm, cmap='coolwarm')

        for g in np.arange(len(list_pair))+1 :
            ax.axhline(y=g-.5, ls='--', c='grey', linewidth=.7)
        clb = plt.colorbar(pclr, fraction=0.03)
        clb.set_label('dv/v = -dt/t (%)')
        
    ax.set_yticks(np.arange(len(list_pair)), labels=list_pair)
    plt.gca().xaxis.set_major_formatter(DateFormatter("%Y-%m-%d\n%H:%M"))
    date_plot = datetime.now()
    plt.gcf().text(0.95, 0.045, 'MSNoise - Lecocq et al. Aurelogram plot\nv%i.%i_%s' %(date_plot.year, date_plot.month, get_git_version())
                   , fontsize=8, horizontalalignment='right', c='grey')
    #ax.set_xlim(datearray[0],datearray[-1])
    if bystation in ('median', 'mean') :
        typeaur = 'by stations (%s)' %bystation
    else :
         typeaur = 'by pair of stations'
    
    plt.title('Aurelogram %s - Method %s\nFilter id : %s, Component : %s, Mov. Stack : %s %s' %(typeaur, method, filterid,components,str(mov_stacks[int(mov_stack_idx)-1][0]), str(mov_stacks[int(mov_stack_idx)-1][1])))
    if outfile :
        plt.savefig(outfile, bbox_inches='tight')
    if show :
        plt.show()