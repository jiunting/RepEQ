#Download catalog or waveforms data
import requests
import datetime
import obspy


def catalog_USGS(time1,time2,area=[-156.357,-154.061,18.407,20.437],magnitude=3.0,outname='testout.cat'):
    #download the earthquake catalog from USGS earthquake API
    import requests
    import datetime
    import sys

    #formatting check
    if type(time1)==str:
        if len(time1)<4:
            print('Time must be at least [YYYY]MMDD  or <datetime> type')
            sys.exit(2)
        elif len(time1)==4:
            time1=datetime.datetime(float(time1[:4]),1,1)
        elif len(time1)==6:
            time1=datetime.datetime(float(time1[:4]),float(time1[4:6]),1)
        elif len(time1)==8:
            time1=datetime.datetime(float(time1[:4]),float(time1[4:6]),float(time1[6:8]))
        else:
            print('Please make sure time1 format should be [YYYY]MMDD or in <datetime> type')
    if type(time2)==str:
        if len(time2)<4:
            print('Time must be at least [YYYY]MMDD  or <datetime> type')
            sys.exit(2)
        elif len(time2)==4:
            time2=datetime.datetime(float(time2[:4]),1,1)
        elif len(time2)==6:
            time2=datetime.datetime(float(time2[:4]),float(time2[4:6]),1)
        elif len(time2)==8:
            time2=datetime.datetime(float(time2[:4]),float(time2[4:6]),float(time2[6:8]))
        else:
            print('Please make sure time1 format should be [YYYY]MMDD or in <datetime> type')
    
    if type(time1)!=datetime.datetime or type(time2)!=datetime.datetime:
        print('Please make sure time1/time2 format should be in [YYYY]MMDD or in <datetime> type')
        sys.exit(2)

    #time1=datetime.datetime(2000,1,1,0,0,0),time2=datetime.datetime(2019,12,31)
    #header:
    #time,latitude,longitude,depth,mag,magType,nst,gap,dmin,rms,net,id,updated,place,type,horizontalError,depthError,magError,magNst,status,locationSource,magSource'
    #Time you're interested
    time0=time1
    timemax=time2
    OUT1=open(outname,'w') #output catalog file
    #For Hawaii All
    #Define the area where you're interested
    minlongitude=area[0]
    maxlongitude=area[1]
    minlatitude=area[2]
    maxlatitude=area[3]
    #the minimum M (don't want to go too small)
    minmag=magnitude
    deltadays=30
    if minmag<2:
        deltadays=10
    elif minmag>=2 and minmag<5:
        deltadays=30
    elif minmag>=5 and minmag<6:
        deltadays=50
    else:
        deltadays=60
    time1=time0
    while True:
        time2=time1+datetime.timedelta(days=deltadays) #if too many earthquakes that exceed the API limit, change 60 days to smaller
        time1str=time1.strftime('%Y-%m-%d')
        time2str=time2.strftime('%Y-%m-%d')
        url='https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime='+time1str+'&endtime='+time2str+'&minmagnitude='+str(minmag)+'&maxlatitude='+str(maxlatitude)+'&minlatitude='+str(minlatitude)+'&maxlongitude='+str(maxlongitude)+'&minlongitude='+str(minlongitude)
        data=requests.get(url)
        #write file
        lines=data.text.split('\n')[::-1][1:-1]
        for line in lines:
            OUT1.write('%s\n'%(line))
        time1=time2
        if time1>timemax:
            break

    OUT1.close()



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
                  

def download_waves(evtime,sec_bef_aft=[120,600],Ftype='circ',lon_lat=[120,24],range_rad=[0,6],provider=["IRIS"],OUT='./'):
    import obspy
    from obspy.clients.fdsn.mass_downloader import CircularDomain,RectangularDomain,Restrictions, MassDownloader
    import datetime
    import glob
    import pandas as pd
    import numpy as np
    import time
    '''
    print('evtime:2000-01-01T06:58:39.780')
    print('Seconds before and after the evtime, sec_bef_aft=[120,600]')
    print('Restriction type: FType=="circ" or "rect"')
    print('Center/boundary of the filtering reference ex:[120,24] for "circ" or [119,123,21,26] for "rect"')
    print('range_rad=[0,6]')
    '''
    
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
    assert len(sec_bef_aft)==2, 'Give seconds before and after the origin time! Ex:[60,1200]'
    if Ftype=='circ':
        assert len(lon_lat)==2, 'length of lon_lat should be 2. Ex:[123,32.1]'
        assert len(range_rad)==2, 'length of range_rad should be 2'
        domain = CircularDomain(latitude=lon_lat[1], longitude=lon_lat[0],minradius=range_rad[0], maxradius=range_rad[1])
    elif Ftype=='rect':
        assert len(lon_lat)==4, 'length of lon_lat should be 4. Ex:[119,123,21,26]'
        #assert len(range_rad)==2, 'length of range_rad should be 2'
        domain = RectangularDomain(minlatitude=lon_lat[2], maxlatitude=lon_lat[3], minlongitude=lon_lat[0], maxlongitude=lon_lat[1])
    
    
    restrictions = Restrictions(
    # Get data from 5 minutes before the event to one hour after the
    # event. This defines the temporal bounds of the waveform data.
    starttime=origin_time - sec_bef_aft[0],
    endtime=origin_time + sec_bef_aft[1],
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
    mdl = MassDownloader(providers=provider)
#    mdl = MassDownloader(providers=["IRIS"])
    # The data will be downloaded to the ``./waveforms/`` and ``./stations/``
    # folders with automatically chosen file names.
    outstr=evtime.split('T')[0].replace('-','')+str(HH).zfill(2)+str(MM).zfill(2)+str(SS).zfill(2) #save dir as evid
    outmsdir=OUT+'/'+outstr+"/waveforms"
    outstadir=OUT+'/'+outstr+"/stations"
    
    mdl.download(domain, restrictions,threads_per_client=20, mseed_storage=outmsdir,stationxml_storage=outstadir)
    return(outmsdir,outstadir)


    
def rm_response(mseedpath,stapath,setlon,setlat,stainfo_path=None):
    import obspy
    from obspy.clients.fdsn.mass_downloader import CircularDomain,RectangularDomain,Restrictions, MassDownloader
    import datetime
    import glob
    import pandas as pd
    import numpy as np
    import time
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
