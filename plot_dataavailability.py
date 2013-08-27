from tmp_database_tools import *
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
import matplotlib.dates
import numpy as np
import matplotlib.gridspec as gridspec

db = connect()
start, end, datelist = build_movstack_datelist(db)
print dir(start)
dates = []
stations = []
for day in datelist:
    # print day
    data = get_data_availability(db,starttime=day, endtime=day)
    for di in data:
        net, sta, comp, starttime, endtime, data_duration, gaps_duration, samplerate, flag = di
        stations.append("%s.%s"%(net,sta))
        dates.append(starttime)

data = pd.DataFrame({"stations":stations},index=dates)
data = data.groupby('stations')

llen = (end-start).days +1
ngroups = len(data.groups.keys())
matrix = np.zeros((ngroups,llen))
start = datetime.datetime.combine(start, datetime.time(0,0,0))

for i, group in enumerate(sorted(data.groups.keys())):
    print group
    new  = True
    for di in data.groups[group]:
        if new:
            print group, di
            new=  False
        dt= (di-start).days
        matrix[i,dt] = 1

gs = gridspec.GridSpec(2, 1,height_ratios=[4, 1]) 

plt.figure(figsize=(11.8,8.4))
ax = plt.subplot(gs[0])
plt.imshow(matrix,interpolation="none",aspect='auto',cmap='bwr',vmin=-1,vmax=1,extent=(date2num(start),date2num(end),0,ngroups),origin='lower')

plt.yticks(np.arange(ngroups)+0.5,sorted(data.groups.keys()))
ax.xaxis.set_major_locator(
    matplotlib.dates.MonthLocator())

ax.xaxis.set_major_formatter(
    matplotlib.dates.DateFormatter('%Y-%m-%d')
)
plt.gcf().autofmt_xdate()
plt.grid()


ax = plt.subplot(gs[1])
plt.plot(datelist,np.sum(matrix,axis=0))
plt.ylabel('N stations')
plt.gcf().autofmt_xdate()
plt.grid()

plt.show()

