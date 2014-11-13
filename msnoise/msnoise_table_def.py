# table_def.py
from sqlalchemy import create_engine
from sqlalchemy import Column, Date, Integer, String, Float, Boolean, DateTime, text, TIMESTAMP, Enum
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
import datetime

Base = declarative_base()
########################################################################
class Filter(Base):
    """"""
    __tablename__ = "filters"
 
    ref = Column(Integer, primary_key=True)
    low = Column(Float)  
    mwcs_low = Column(Float)  
    high = Column(Float)  
    mwcs_high = Column(Float)  
    rms_threshold = Column(Float)  
    mwcs_wlen = Column(Float)  
    mwcs_step = Column(Float)  
    used = Column(Boolean)  
 
    #----------------------------------------------------------------------
    def __init__(self, low, mwcs_low, high, mwcs_high, rms_threshold, mwcs_wlen, mwcs_step, used):
        """"""
        self.low = low   
        self.mwcs_low = mwcs_low   
        self.high = high   
        self.mwcs_high = mwcs_high   
        self.rms_threshold = rms_threshold   
        self.mwcs_wlen = mwcs_wlen   
        self.mwcs_step = mwcs_step   
        self.used = used   
 
########################################################################

class Job(Base):
    """"""
    __tablename__ = "jobs"
 
    ref = Column(Integer, primary_key=True)
    day = Column(String(10))  
    pair = Column(String(20))  
    type = Column(String(10))  
    flag = Column(String(1))
    lastmod = Column(TIMESTAMP, server_onupdate=text('CURRENT_TIMESTAMP'))
 
    #----------------------------------------------------------------------
    def __init__(self, day, pair, type,flag,lastmod=None):
        """"""
        self.day = day   
        self.pair = pair   
        self.type = type   
        self.flag = flag
        if lastmod is None:
            self.lastmod=datetime.datetime.utcnow()

 
########################################################################

class Station(Base):
    """"""
    __tablename__ = "stations"
    ref = Column(Integer, primary_key=True)
    net = Column(String(10))  
    sta = Column(String(10))  
    X = Column(Float)  
    Y = Column(Float)  
    altitude = Column(Float)  
    coordinates = Column(Enum('DEG','UTM'))
    instrument = Column(String(20))  
    used = Column(Boolean)  

 
    #----------------------------------------------------------------------
    def __init__(self, net,sta,X,Y,altitude,coordinates,instrument,used):
        """"""
        self.net = net   
        self.sta = sta   
        self.X = X   
        self.Y = Y   
        self.altitude = altitude   
        self.coordinates = coordinates   
        self.instrument = instrument   
        self.used = used   

########################################################################

class Config(Base):
    """"""
    __tablename__ = "config"
    name = Column(String(255), primary_key=True)  
    value = Column(String(255))  
 
    #----------------------------------------------------------------------
    def __init__(self, name, value):
        """"""
        self.name = name
        self.value = value 

########################################################################

class DataAvailability(Base):
    """"""
    __tablename__ = "data_availability"
    ref = Column(Integer, primary_key=True,autoincrement=True)
    net = Column(String(10))  
    sta = Column(String(10))  
    comp = Column(String(20))  
    path = Column(String(255))  
    file = Column(String(255))  
    starttime = Column(DateTime)  
    endtime = Column(DateTime)  
    data_duration = Column(Float)
    gaps_duration = Column(Float)
    samplerate = Column(Float)
    flag = Column(String(1))

 
    #----------------------------------------------------------------------
    def __init__(self, net,sta,comp,path,file,starttime,endtime,data_duration,gaps_duration,samplerate,flag):
        """"""
        self.net = net
        self.sta = sta
        self.comp = comp
        self.path = path
        self.file = file
        self.starttime = starttime
        self.endtime = endtime
        self.data_duration = data_duration
        self.gaps_duration = gaps_duration
        self.samplerate = samplerate
        self.flag = flag

########################################################################
 
