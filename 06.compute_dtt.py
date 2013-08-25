import numpy as np
import os

from database_tools import *
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import statsmodels.api as sm

import logging
logging.basicConfig(level=logging.DEBUG,
                    filename="./compute_dtt.log",
                    format='%(asctime)s [%(levelname)s] %(message)s',
                    filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

logging.info('*** Starting: Compute DT/T ***')


def wavg_wstd(data, errors):
    d = data
    errors[errors==0] = 1e-6
    w = 1./errors
    wavg = (d * w).sum() / w.sum()
    N = len(np.nonzero(w)[0])
    wstd = np.sqrt( np.sum(w*(d-wavg)**2) / ( (N-1) * np.sum(w) /N) )
    return wavg, wstd



db = connect()

# PLOT= True
PLOT= False

lMlag = -50 #HARD
lmlag = -5 #HARD
rmlag = 5 #HARD
rMlag = 50 #HARD
minCoh = 0.5 #HARD
maxErr = 0.1 #HARD
maxDt = 0.5 #HARD

start, end, datelist = build_movstack_datelist(db)

mov_stack = get_config(db,"mov_stack")
if mov_stack.count(',') == 0:
    mov_stacks = [int(mov_stack),]
else:
    mov_stacks = [int(mi) for mi in mov_stack.split(',')]

components_to_compute = get_components_to_compute(db)
updated_dtt = updated_days_for_dates(db,start,end,'%',type='DTT',returndays=True)

for f in get_filters(db,all=False):
    filterid = int(f.ref)
    for components in components_to_compute:
        for mov_stack in mov_stacks:
            logging.info('Loading mov=%i days for filter=%02i'%(mov_stack,filterid))
            for current in updated_dtt:
                logging.debug("Processing %s" % current)
                first = True
                for station1, station2 in get_station_pairs(db,used=True):
                    sta1 = "%s_%s"%(station1.net, station1.sta)
                    sta2 = "%s_%s"%(station2.net, station2.sta)
                    pair = "%s_%s"%(sta1,sta2)
                    day = os.path.join('MWCS',"%02i"%filterid,"%03i_DAYS"%mov_stack,components,pair,'%s.txt'%current)
                    if os.path.isfile(day):
                        df = pd.read_csv(day,delimiter=' ',header=None,index_col=0,names=['t','dt','err','coh'])
                        if first:
                            tArray = df.index.values
                            dtArray = df['dt']
                            errArray = df['err']
                            cohArray = df['coh']
                            pairArray = [pair,]
                            first = False
                        else:
                            dtArray = np.vstack((dtArray,df['dt']))
                            errArray = np.vstack((errArray,df['err']))
                            cohArray = np.vstack((cohArray,df['coh']))
                            pairArray.append(pair)
                        del df
                    del day
                
                if not first:
                    tindex = np.where(((tArray >= lMlag) & (tArray<=lmlag)) | ((tArray>=rmlag) & (tArray<=rMlag)))[0]
                    
                    Dates = []
                    Pairs = []
                    M = []
                    EM = []
                    A = []
                    EA = []
                    M0 = []
                    EM0 = []
                    if len(pairArray) != 1:
                        # first stack all pairs to a ALL mean pair, using indexes of selected values:
                        new_dtArray = np.zeros(len(tArray))
                        new_errArray = np.zeros(len(tArray)) + 9999
                        new_cohArray = np.zeros(len(tArray))
                        for i in range(len(tArray)):
                            if i in tindex:
                                cohindex = np.where(cohArray[:,i] >= minCoh)[0]
                                errindex = np.where(errArray[:,i] <= maxErr)[0]
                                dtindex = np.where(np.abs(dtArray[:,i]) <= maxDt)[0]
                            
                                index = np.intersect1d(cohindex,errindex)
                                index = np.intersect1d(index,dtindex)
                                
                                wavg, wstd = wavg_wstd(dtArray[:,i][index], errArray[:,i][index])
                                new_dtArray[i] = wavg
                                new_errArray[i] = wstd
                                new_cohArray[i] = 1.0
                        
                        dtArray = np.vstack((dtArray,new_dtArray))
                        errArray = np.vstack((errArray,new_errArray))
                        cohArray = np.vstack((cohArray,new_cohArray))
                        pairArray.append("ALL")
                        del new_cohArray,new_dtArray,new_errArray,cohindex,errindex,dtindex,wavg,wstd
                    
                    if PLOT:
                        plt.figure(figsize=(15,10))
                        gs = gridspec.GridSpec(1, 6,width_ratios=[3, 3, 3, 3, 3, 3]) 
                        plt.subplots_adjust(hspace=0.005,wspace=0.15,bottom=0.05,top=0.92,left=0.07,right=0.95)
                        plt.subplot(gs[0])
                        ymax = dtArray.shape[0] 
                        extent = (-117.5,117.5,ymax+0.5,-0.5) #HARD
        
                        plt.imshow(dtArray,interpolation='none',aspect='auto',vmin=-maxDt,vmax=maxDt,cmap='bwr',extent=extent)
                        for item in [lMlag,lmlag,rmlag,rMlag]:
                            plt.axvline(item,ls='--',zorder=15,c='k',lw=2)
                        plt.title('Delays')
                        cb = plt.colorbar(orientation='horizontal',pad=0.08,shrink=0.8)
                        cb.set_label('Delay (seconds)')
                        cb.set_ticks([-maxDt,0,maxDt])
                        plt.xlabel('Time Lag')
                        plt.ylabel('Pair Number')
                        plt.subplot(gs[1])
                        plt.title('Errors')
                        plt.imshow(errArray,interpolation='none',aspect='auto',vmax=0.1,extent=extent,cmap='hot')
                        for item in [lMlag,lmlag,rmlag,rMlag]:
                            plt.axvline(item,ls='--',zorder=15,c='k',lw=2)
                        plt.yticks([],[])
                        cb = plt.colorbar(orientation='horizontal',pad=0.08,shrink=0.8)
                        cb.set_label('Error (seconds)')
                        cb.set_ticks([-0.1,0,0.1])
                        plt.xlabel('Time Lag')
                        plt.subplot(gs[2])
                        plt.title('Phase Coherences')
                        plt.imshow(cohArray,interpolation='none',aspect='auto',extent=extent,cmap='hot')
                        for item in [lMlag,lmlag,rmlag,rMlag]:
                            plt.axvline(item,ls='--',zorder=15,c='k',lw=2)
                        cb = plt.colorbar(orientation='horizontal',pad=0.08,shrink=0.8)
                        cb.set_label('Coherence Coefficient')
                        cb.set_ticks([0.0,0.25,0.5,0.75,1.0])
                        plt.yticks([],[])
                        plt.xlabel('Time Lag')
                        
                    # then process all pairs + the ALL
                    if len(dtArray.shape) == 1: #if there is only one pair:
                        dtArray = dtArray.reshape((1,dtArray.shape[0]))
                        cohArray = cohArray.reshape((1,cohArray.shape[0]))
                        errArray = errArray.reshape((1,errArray.shape[0]))

                    
                    used = np.zeros(dtArray.shape)
                    
                    for i, pair in enumerate(pairArray):
                        cohindex = np.where(cohArray[i] >= minCoh)[0]
                        errindex = np.where(errArray[i] <= maxErr)[0]
                        dtindex = np.where(np.abs(dtArray[i]) <= maxDt)[0]
                        
                        index = np.intersect1d(tindex,cohindex)
                        index = np.intersect1d(index,errindex)
                        index = np.intersect1d(index,dtindex)
                        
                        used[i][index] = 1.0
                        
                        w = 1.0 / (errArray[i][index] ** 2)
                        
                        VecXfilt = tArray[index]
                        VecYfilt = dtArray[i][index]
                        if len(VecYfilt) >= 2:
                            B = sm.tools.tools.add_constant(VecXfilt,prepend=False)
                            res = sm.regression.linear_model.WLS(VecYfilt,B,w).fit()
                            res0 = sm.regression.linear_model.WLS(VecYfilt,VecXfilt,w).fit()
                            if res.df_resid >0:
                                m, a = res.params
                                em,ea = res.bse
                                
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
        
                                del  m, a, em, ea, m0, em0
                            del res, res0, B
                        del VecXfilt, VecYfilt, w
                        del index, cohindex, errindex, dtindex
                        
                    
                    if PLOT:
                        plt.subplot(gs[3])
                        plt.title('Data Selection')
                        plt.imshow(used,extent=extent,cmap='binary',interpolation='none',aspect='auto')
                        cb = plt.colorbar(orientation='horizontal',pad=0.08,shrink=0.8)
                        cb.set_label('Selected')
                        cb.set_ticks([0.0,1.0],)
                        cb.set_ticklabels(['No','Yes'])
                        plt.yticks([],[])
                        plt.xlabel('Time Lag')
        
        
                        plt.subplot(gs[4])
                        plt.title('Delay Time Variations\nNo Forcing')
                        Ma = np.array(M)
                        EMa = np.array(EM)
                        plt.errorbar(100*Ma,np.arange(len(Ma)),xerr=EMa*100,c='k',lw=0,elinewidth=2.0)
                        plt.scatter(100*Ma, np.arange(len(Ma)), c=Ma*100,edgecolor='k',vmin=-0.2,vmax=0.2,zorder=100,s=75)
                        
                        cb = plt.colorbar(orientation='horizontal',pad=0.08,shrink=0.8,)
                        cb.set_ticks([-0.2,0.2],)
                        cb.set_label('dt/t in %')
                        plt.axvline(0,c='k')
                        plt.axvline(100*Ma[-1],c='r',lw=2,label="ALL: %.2e"%Ma[-1])
                        plt.axvspan(100*Ma[-1]-100*EM[-1],100*Ma[-1]+100*EM[-1],color='r',alpha=0.2)
                        # print "M,EM", 100*M[-1],100*EM[-1]
                        
                        wdtt,werr = wavg_wstd(np.array(Ma[:-1]),np.array(EMa[:-1]))
                        plt.axvline(100*wdtt,c='g',lw=3,label="Weighted $\overline{x}$: %.2e"%wdtt)
                        plt.axvspan(100*wdtt-100*werr,100*wdtt+100*werr,color='g',alpha=0.2)
                        # print 'wdtt, werr', 100*wdtt, 100*werr
                        
                        plt.xlim(-0.2,0.2)
                        plt.xticks([-0.2,0.0,0.2])
                        plt.xlabel('dt/t in %')
                        plt.ylim(len(Ma),0)
                        plt.yticks([],[])
                        
                        plt.subplot(gs[5])
                        plt.title('Delay Time Variations\nForcing Origin')
                        Ma = np.array(M0)
                        EMa = np.array(EM0)
                        plt.errorbar(100*Ma,np.arange(len(Ma)),xerr=EMa*100,c='k',lw=0,elinewidth=2.0)
                        plt.scatter(100*Ma, np.arange(len(Ma)), c=Ma*100,edgecolor='k',vmin=-0.2,vmax=0.2,zorder=100,s=75)
                        
                        cb = plt.colorbar(orientation='horizontal',pad=0.08,shrink=0.8,)
                        cb.set_ticks([-0.2,0.2],)
                        cb.set_label('dt/t in %')
                        plt.axvline(0,c='k')
                        plt.axvline(100*Ma[-1],c='r',lw=2,label="ALL: %.2e"%Ma[-1])
                        plt.axvspan(100*Ma[-1]-100*EM0[-1],100*Ma[-1]+100*EM0[-1],color='r',alpha=0.2)
                        # print "M0,EM0", 100*M0[-1],100*EM0[-1]
                        
                        wdtt,werr = wavg_wstd(np.array(Ma[:-1]),np.array(EMa[:-1]))
                        plt.axvline(100*wdtt,c='g',lw=3,label="Weighted $\overline{x}$: %.2e"%wdtt)
                        plt.axvspan(100*wdtt-100*werr,100*wdtt+100*werr,color='g',alpha=0.2)
                        # print 'wdtt0, werr0', 100*wdtt, 100*werr
                        
                        plt.xlim(-0.2,0.2)
                        plt.xticks([-0.2,0.0,0.2])
                        plt.xlabel('dt/t in %')
                        plt.ylim(len(Ma),0)
                        plt.yticks([],[])            
                        
                        plt.savefig('dttmatrix_%02i_%03iDAYS_%s-%s.png'%(filterid,mov_stack,components,current),dpi=300)
        
                        stations = get_stations(db,used=True)
                        nstations = len(stations)
                        dttmatrix = np.zeros((nstations,nstations))
                        for station1, station2 in get_station_pairs(db,used=True):
                            sta1 = "%s_%s"%(station1.net, station1.sta)
                            sta2 = "%s_%s"%(station2.net, station2.sta)
                            pair = "%s_%s"%(sta1,sta2)
                            if pair in pairArray:
                                id = pairArray.index(pair)
                                dttmatrix[i,j] = M[id]*100
                        plt.figure(figsize=(8,8))
                        plt.subplots_adjust(hspace=0.005,wspace=0.1,bottom=0.05,top=0.90,left=0.1,right=0.90)
                        plt.imshow(dttmatrix,interpolation='none',vmin=-0.2,vmax=0.2,cmap='bwr')
                        cb = plt.colorbar(shrink=0.8,pad=0.05)
                        cb.set_label('dt/t in %')
                        plt.yticks(np.arange(nstations),["%s.%s"%(st.net, st.sta) for st in stations])
                        plt.xticks(np.arange(nstations),["%s.%s"%(st.net, st.sta) for st in stations],rotation=90,verticalalignment='bottom')
                        for tick in plt.gca().xaxis.iter_ticks():
                            tick[0].label2On = True
                            tick[0].label1On = False
                            tick[0].label2.set_rotation('vertical')
                        plt.grid()
                        plt.savefig('dttmatrix2_%02i_%03iDAYS_%s-%s.png'%(filterid,mov_stack,components,current),dpi=300)
                        plt.show()
                        
                    logging.debug("%s: exporting: %i pairs"%(current, len(pairArray)))
                    df = pd.DataFrame({'Pairs':Pairs,'M':M,'EM':EM,'A':A,'EA':EA,'M0':M0,'EM0':EM0},index=pd.DatetimeIndex(Dates))
                    #Needs to be changed !
                    output = os.path.join('DTT',"%02i"%filterid,"%03i_DAYS"%mov_stack,components)
                    if not os.path.isdir(output):
                        os.makedirs(output)
                    df.to_csv(os.path.join(output,'%s.txt'%current),index_label='Date')
                    del df, M, EM, A, EA,M0,EM0, Pairs, Dates,used
                    del tArray ,dtArray,errArray,cohArray,pairArray
                    del output
                

logging.info('*** Finished: Compute DT/T ***')