import glob, os
import numpy as np
from database_tools import *

db = connect()
print
print ">> Populating the Station table"
print
data_folder = get_config(db,'data_folder')
data_structure = get_config(db,'data_structure')
if data_structure in ["SDS","IDDS"]:
    datalist = sorted(glob.glob( os.path.join(data_folder, "*","*","*")))
    stations = []
    for di in datalist:
        tmp = os.path.split(di)
        sta = tmp[1]
        net = os.path.split(tmp[0])[1]
        stations.append("%s_%s"%(net,sta))
    del datalist
elif data_structure in ["BUD",]:
    datalist = sorted(glob.glob( os.path.join(data_folder, "*","*",)))
    stations = []
    for di in datalist:
        tmp = os.path.split(di)
        sta = tmp[1]
        net = os.path.split(tmp[0])[1]
        stations.append("%s_%s"%(net,sta))
    del datalist
elif data_structure in ["PDF",]:
    datalist = sorted(glob.glob( os.path.join(data_folder, "*","*",)))
    stations = []
    for di in datalist:
        tmp = os.path.split(di)
        sta = tmp[1]
        net = get_config(db,'network')
        stations.append("%s_%s"%(net,sta))
    del datalist

stations = np.unique(stations)

db = connect()
for station in stations:
    net, sta = station.split('_')
    print 'Adding:', net, sta
    X = 0.0
    Y = 0.0
    altitude = 0.0
    coordinates = 'UTM'
    instrument = 'N/A'
    update_station(db, net, sta, X, Y, altitude, coordinates = coordinates, instrument=instrument)