import matplotlib.pyplot as plt
import numpy as np
import datetime
import pandas as pd
import random
import obspy

catalogpath='/Users/timlin/Documents/Project/EQrequest/Hawaii_ALL_M3.dat' #EQ information
A=pd.read_csv(catalogpath,header=None,sep=',',names=['time','eqlat','eqlon','eqdep','eqmag','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x'],skiprows=0)
Cat_Date=pd.to_datetime(A['time'],format='%Y-%m-%dT%H:%M:%S.%fZ')#convert time to Datetime format

INfile=open('pairs_BP0.8-2_wind30s_new_rob.seq','r')

nseq=0
sav_seq=[]
plt.figure(1)
sav_dist=[]
for line in INfile.readlines():
    #random color
    r_col=[random.random(), random.random(), random.random()]
    eqs = line.split()
    sav_time=[]
    sav_seq_lon=[]
    sav_seq_lat=[]
    pre_lon=[]
    pre_lat=[]
    for eq in eqs:
        EQtime=datetime.datetime.strptime(eq,'%Y%m%d%H%M%S')
        sav_time.append(EQtime)
        findCat=Cat_Date[(Cat_Date>EQtime-datetime.timedelta(seconds=30)) & (Cat_Date<EQtime+datetime.timedelta(seconds=30))]
        if len(findCat.keys())>1:
            continue #the two earthquake too close may contaminant each other!
        idx_cat=findCat.keys()[0]#now looking for this event, catalog may have duplicated, use only the first one!
        eqlon=float(A['eqlon'][idx_cat])
        eqlat=float(A['eqlat'][idx_cat])
        eqdep=float(A['eqdep'][idx_cat])
        if pre_lon==[]:
            pre_lon=eqlon
            pre_lat=eqlat
        dist_deg=obspy.geodetics.locations2degrees(lat1=pre_lat,long1=pre_lon,lat2=eqlat,long2=eqlon)
        if dist_deg>0.1:
            plt.figure(1)
            plt.plot(eqlon,eqlat,'r*',markersize=15)
            plt.text(eqlon,eqlat,eq)
        sav_dist.append(dist_deg)
        sav_seq_lon.append(eqlon)
        sav_seq_lat.append(eqlat)
    plt.figure(1)
    plt.plot(sav_seq_lon,sav_seq_lat,'o-',color=r_col)
    plt.figure(2)
    plt.plot(sav_time,[nseq]*len(sav_time),'o-',color=r_col)
    nseq += 1
#    if nseq==10000:
#        break

plt.show()

plt.plot(sav_dist)
plt.show()

INfile.close()
