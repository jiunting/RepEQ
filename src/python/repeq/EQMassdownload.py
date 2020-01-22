#Download waveforms based on catalog
"""
# -*- coding: utf-8 -*-
Created on Mon Feb 12 12:40:50 2018
@author: timlin
"""
import obspy
from obspy.clients.fdsn.mass_downloader import CircularDomain,RectangularDomain,Restrictions, MassDownloader
import datetime
import glob
import pandas as pd
import numpy as np
import time

def cattime2normal(timestr):
    '''
    convert '2000-01-01T06:58:39.780Z' to
    [2000,1,1,6,58,39,780000]
    '''
    tmp1,tmp2=timestr.split('T')
    yyyy,mm,dd=tmp1.split('-')
    yyyy=int(yyyy);mm=int(mm);dd=int(dd)
    HH,MM,SS=tmp2.split(':')
    HH=int(HH);MM=int(MM);SS=float(SS[:-1]);NS=int((SS-int(SS))*1E6);SS=int(SS)
    evtime=[yyyy,mm,dd,HH,MM,SS,NS]
    return(evtime)


def EQfilter(catalog,BT_time=['20170607020500','20191210000000'],BT_lon=[132,133],BT_lat=[30,31],BT_dep=[0,50],BT_mag=[5.0,9.0]):
    '''
    IN:
     Catalog file from USGS API (run: USGS_catalogAPI.py)
    OUT:
     EQfilter return a list of EQ time string
     EQ=['2000-01-01T06:55:51.000Z','2018-02-11T18:06:39.000Z']
    EXAMPLE:
     T,lon,lat,dep,mag=EQfilter('Catalog2000.dat',BT_time=[t1,t2],
           BT_lon=[119,123],BT_lat=[21,26],BT_dep=[0,100],BT_mag=[5.0,9.0])
    ''' 
    #input adjustment
    time1=BT_time[0]
    time2=BT_time[1]
    if not( type(time1) is datetime.datetime ):
        yyyy=int(time1[0:4])
        mm=int(time1[4:6])
        dd=int(time1[6:8])
        HH=int(time1[8:10])
        MM=int(time1[10:12])
        SS=int(time1[12:14])
        time1=datetime.datetime(yyyy,mm,dd,HH,MM,SS)
    
    if not( type(time2) is datetime.datetime ):
        yyyy=int(time2[0:4])
        mm=int(time2[4:6])
        dd=int(time2[6:8])
        HH=int(time2[8:10])
        MM=int(time2[10:12])
        SS=int(time2[12:14])
        time2=datetime.datetime(yyyy,mm,dd,HH,MM,SS)

    with open(catalog,'r') as IN1:
        T=[];lat=[];lon=[];dep=[];mag=[];#magTp=[];nst=[];gap=[]
        #dmin=[];rms=[];net=[];evid=[];updated=[];place=[];evtype=[];
        #horizontalError=[];depthError=[];magError=[];magNst=[];status=[];locationSource=[];magSource=[]
        for line in IN1.readlines():
            #print(line) #only for checking, dangerous!!! do not run this   
            elems=line.split(',')
            try:
                tmp_T=elems[0] #e.g. tmp_T='2000-01-01T06:55:51.000Z'
                tmpyyyy,tmpmm,tmpdd,tmpHH,tmpMM,tmpSS,tmpNS=cattime2normal(tmp_T)
                tmp2_T=datetime.datetime(tmpyyyy,tmpmm,tmpdd,tmpHH,tmpMM,tmpSS) #convert to datetime
                tmp_lat=float(elems[1])
                tmp_lon=float(elems[2])
                tmp_dep=float(elems[3])
                tmp_mag=float(elems[4])
                #tmp_magTp=elems[5]
                #filter
                if (time1<=tmp2_T<=time2) and (BT_lon[0]<=tmp_lon<=BT_lon[1]) and (BT_lat[0]<=tmp_lat<=BT_lat[1]) and (BT_dep[0]<=tmp_dep<=BT_dep[1]) and (BT_mag[0]<=tmp_mag<=BT_mag[1]):
                    T.append(tmp_T)
                    lat.append(tmp_lat)
                    lon.append(tmp_lon)
                    dep.append(tmp_dep)
                    mag.append(tmp_mag)
                    #magTp.append()
                else:
                    continue
                
            except:
                print("formating wrong\n",line)
                continue #something wring
    return(T,lon,lat,dep,mag)
                  

def download(evtime,sec_before,sec_after,lon,lat,minrad,maxrad,provider=["IRIS"],OUT='./'):
    '''
    #example input
    evtime='2000-01-01T06:58:39.780Z'
    sec_before=120
    sec_after=600
    lon=120
    lat=24
    minrad=0
    maxrad=10
    provider=["IRIS"]
    OUT='./TESTDL'
    download(evtime,sec_before,sec_after,lon,lat,minrad,maxrad,provider,OUT)
    '''
    yyyy,mm,dd,HH,MM,SS,NS=cattime2normal(evtime) 
    origin_time = obspy.UTCDateTime(yyyy,mm,dd,HH,MM,SS,NS)
    domain = CircularDomain(latitude=lat, longitude=lon,
                        minradius=minrad, maxradius=maxrad)
    
    
    #domain = RectangularDomain(minlatitude=34.452, maxlatitude=38.72, minlongitude=-123.201, maxlongitude=-118.015)
    
    
    restrictions = Restrictions(
    # Get data from 5 minutes before the event to one hour after the
    # event. This defines the temporal bounds of the waveform data.
    starttime=origin_time - sec_before,
    endtime=origin_time + sec_after,
    # You might not want to deal with gaps in the data. If this setting is
    # True, any trace with a gap/overlap will be discarded.
    reject_channels_with_gaps=True,
    # And you might only want waveforms that have data for at least 95 % of
    # the requested time span. Any trace that is shorter than 95 % of the
    # desired total duration will be discarded.
    minimum_length=0.9,
    # No two stations should be closer than 10 km to each other. This is
    # useful to for example filter out stations that are part of different
    # networks but at the same physical station. Settings this option to
    # zero or None will disable that filtering.
    minimum_interstation_distance_in_m=10E1,
    # Only HH or BH channels. If a station has HH channels, those will be
    # downloaded, otherwise the BH. Nothing will be downloaded if it has
    # neither. You can add more/less patterns if you like.
    #channel_priorities=["HH[ZNE]", "BH[ZNE]"],
    #channel_priorities=["BH[ZNE]"],
    #channel_priorities=["BH[ZNE]"],
    #channel_priorities=["BHZ","HNZ"],
    #channel_priorities=["BHZ"],
    #channel_priorities=["HNZ"],
    channel_priorities=["BHZ","HHZ"],
    #channel_priorities=["HN[ZNE]"],
    # Location codes are arbitrary and there is no rule as to which
    # location is best. Same logic as for the previous setting.
    location_priorities=["", "00", "10","100"])
    
    # No specified providers will result in all known ones being queried.
    #mdl = MassDownloader(providers=provider)
    mdl = MassDownloader(providers=["IRIS"])
    # The data will be downloaded to the ``./waveforms/`` and ``./stations/``
    # folders with automatically chosen file names.
    outstr=evtime.split('T')[0].replace('-','')+str(HH).zfill(2)+str(MM).zfill(2)+str(SS).zfill(2) #save dir as evid
    outmsdir=OUT+'/'+outstr+"/waveforms"
    outstadir=OUT+'/'+outstr+"/stations"
    
    mdl.download(domain, restrictions,threads_per_client=20, mseed_storage=outmsdir,stationxml_storage=outstadir)
    return(outmsdir,outstadir)


    
def rm_response(mseedpath,stapath,setlon,setlat,stainfo_path=None):
    mseeds=glob.glob(mseedpath+'/*mseed')
    pre_filt = (0.001, 0.002, 40, 50)
    
    if stainfo_path != None:
        stainfo=pd.read_csv(stainfo_path,sep='|',header=0)
    
    for long_seed in mseeds: #seedpath
        seed=long_seed.split('/')[-1]
        net=seed.split('.')[0]
        sta=seed.split('.')[1]
        #find corresponding sta file
        respf=glob.glob(stapath+'/'+net+'.'+sta+'.xml')
        if respf:
            D=obspy.read(long_seed)
            D[0].write(long_seed,format='SAC') #convert to SAC
            D.clear
            D=obspy.read(long_seed)
            
            inv=obspy.read_inventory(respf[0])
            D.detrend(type='linear')
            D.taper(max_percentage=0.05)
            D.remove_response(inventory=inv, output='vel', pre_filt=pre_filt)
            
            outname=long_seed.replace('.mseed','.SAC')
            #print('outname=',outname)
            
            #D=obspy.read(outname) #reopen the sac file and write EQloc and staloc if exist
            D[0].stats.sac['evlo']=setlon
            D[0].stats.sac['evla']=setlat

            if stainfo_path != None:
                idx=np.where(stainfo['Station'] == sta)[0]
                if len(idx)!=0: 
                    D[0].stats.sac['stlo']=stainfo['Longitude'][idx[0]]
                    D[0].stats.sac['stla']=stainfo['Latitude'][idx[0]]
            D[0].stats.sac['lcalda']=True
            D[0].write(outname,format='SAC')

            D.clear
            
        else:
            print('Cannot find station file for',seed)



'''
#################start download waveforms##########################
T,lon,lat,dep,mag=EQfilter('Catalog_Sumatra_N1.dat',BT_time=['20000101000000','20191201000000'],
                           BT_lon=[-180,180],BT_lat=[-90,90],BT_dep=[0,600],BT_mag=[3.0,9.0])

print('number of EQ:%s found:'%(len(T)))
for i,evtime in enumerate(T):
    print('%d out of %d'%(i,len(T)))
    print(evtime)
    wave_path,resp_path=download(evtime,sec_before=30,sec_after=1200,lon=lon[i],lat=lat[i],minrad=0,maxrad=6.0,provider=['IRIS'],OUT='./Sumatra_EQ')
    time.sleep(5)


##################download END######################################
'''
