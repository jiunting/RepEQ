#Download catalog or waveforms data
import requests
import datetime
import obspy
from obspy import UTCDateTime
import numpy as np

def catalog_USGS(times=['2000','2020'],area=[-156.357,-154.061,18.407,20.437],magnitude=[3.0,6.5],outname='testout.cat'):
    #download the earthquake catalog from USGS earthquake API
    import requests
    import datetime
    import sys
    assert len(times)==2, 'Please provide a range for time [time1,time2] the format can be [YYYY]MMDD  or <datetime> type'
    assert len(area)==4, 'Please provide valid range for searching area [lon1,lon2,lat1,lat2]'
    assert len(magnitude)==2, 'Length of magnitude should be 2, EX:[0.5,6.0]'
    time1=times[0]
    time2=times[1]
    #formatting check
    if type(time1)==str:
        if len(time1)<4:
            print('Time must be at least [YYYY]MMDD  or <datetime> type')
            sys.exit(2)
        elif len(time1)==4:
            time1=datetime.datetime(int(time1[:4]),1,1)
        elif len(time1)==6:
            time1=datetime.datetime(int(time1[:4]),int(time1[4:6]),1)
        elif len(time1)==8:
            time1=datetime.datetime(int(time1[:4]),int(time1[4:6]),int(time1[6:8]))
        else:
            print('Please make sure time1 format should be [YYYY]MMDD or in <datetime> type')
    if type(time2)==str:
        if len(time2)<4:
            print('Time must be at least [YYYY]MMDD  or <datetime> type')
            sys.exit(2)
        elif len(time2)==4:
            time2=datetime.datetime(int(time2[:4]),1,1)
        elif len(time2)==6:
            time2=datetime.datetime(int(time2[:4]),int(time2[4:6]),1)
        elif len(time2)==8:
            time2=datetime.datetime(int(time2[:4]),int(time2[4:6]),int(time2[6:8]))
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
    minmag=magnitude[0]
    maxmag=magnitude[1]
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
        print('Now at:%s %s'%(time2.strftime('%Y/%m/%d'), timemax.strftime('%Y/%m/%d')) )
        time2=time1+datetime.timedelta(days=deltadays) #if too many earthquakes that exceed the API limit, change 60 days to smaller
        time1str=time1.strftime('%Y-%m-%d')
        time2str=time2.strftime('%Y-%m-%d')
        url='https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime='+time1str+'&endtime='+time2str+'&minmagnitude='+str(minmag)+'&maxmagnitude='+str(maxmag)+'&maxlatitude='+str(maxlatitude)+'&minlatitude='+str(minlatitude)+'&maxlongitude='+str(maxlongitude)+'&minlongitude='+str(minlongitude)
        data=requests.get(url)
        #write file
        lines=data.text.split('\n')[::-1][1:-1]
        for line in lines:
            OUT1.write('%s\n'%(line))
        time1=time2
        if time1>timemax:
            break

    OUT1.close()



def make_catalog(times=['2000','2020'],dt=86400,lon_lat=[120,24],outname='testout.cat'):
    #make fake catalog for downloading continuous data
    assert len(times)==2, 'Please provide a range for time [time1,time2] the format can be [YYYY]MMDD  or <datetime> type'
    time1=times[0]
    time2=times[1]
    #formatting check
    if type(time1)==str:
        if len(time1)<4:
            print('Time must be at least [YYYY]MMDD  or <datetime> type')
            sys.exit(2)
        elif len(time1)==4:
            time1=datetime.datetime(int(time1[:4]),1,1)
        elif len(time1)==6:
            time1=datetime.datetime(int(time1[:4]),int(time1[4:6]),1)
        elif len(time1)==8:
            time1=datetime.datetime(int(time1[:4]),int(time1[4:6]),int(time1[6:8]))
        else:
            print('Please make sure time1 format should be [YYYY]MMDD or in <datetime> type')
    if type(time2)==str:
        if len(time2)<4:
            print('Time must be at least [YYYY]MMDD  or <datetime> type')
            sys.exit(2)
        elif len(time2)==4:
            time2=datetime.datetime(int(time2[:4]),1,1)
        elif len(time2)==6:
            time2=datetime.datetime(int(time2[:4]),int(time2[4:6]),1)
        elif len(time2)==8:
            time2=datetime.datetime(int(time2[:4]),int(time2[4:6]),int(time2[6:8]))
        else:
            print('Please make sure time1 format should be [YYYY]MMDD or in <datetime> type')
    time0=time1
    timemax=time2
    OUT1=open(outname,'w') #output catalog file
    while time0<=timemax:
        format_time=time0.strftime('%Y-%m-%dT%H:%M:%S.%fZ')
        OUT1.write('%s,%f,%f,%s\n'%(format_time,lon_lat[1],lon_lat[0],'''6.81,2.95,ml,55,195,0.07168,0.11,us,us12345678,2000-01-01T00:00:00.000Z,"Fake Catalog",earthquake,0.0,0.0,0.0,0,reviewed,us,us'''))
        time0+=datetime.timedelta(seconds=dt)
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
        #time1 is string
        #determine input length
        L_inp=len(time1)
        assert L_inp>=8 #sould be at least yyyymmdd
        st_slice=(0,4,6,8,10,12)
        ed_slice=(4,6,8,10,12,14)
        max_idx=6
        T=[]
        for i in range(max_idx):
            try:
                T.append(int(time1[st_slice[i]:ed_slice[i]]))
            except:
                T.append(0)

        #yyyy=int(time1[0:4])
        #mm=int(time1[4:6])
        #dd=int(time1[6:8])
        #HH=int(time1[8:10])
        #MM=int(time1[10:12])
        #SS=int(time1[12:14])
        time1=datetime.datetime(T[0],T[1],T[2],T[3],T[4],T[5])
    
    if not( type(time2) is datetime.datetime ):
        L_inp=len(time2)
        assert L_inp>=8 #sould be at least yyyymmdd
        st_slice=(0,4,6,8,10,12)
        ed_slice=(4,6,8,10,12,14)
        max_idx=6
        T=[]
        for i in range(max_idx):
            try:
                T.append(int(time2[st_slice[i]:ed_slice[i]]))
            except:
                T.append(0)
    
        #yyyy=int(time1[0:4])
        #yyyy=int(time2[0:4])
        #mm=int(time2[4:6])
        #dd=int(time2[6:8])
        #HH=int(time2[8:10])
        #MM=int(time2[10:12])
        #SS=int(time2[12:14])
        #time2=datetime.datetime(yyyy,mm,dd,HH,MM,SS)
        time2=datetime.datetime(T[0],T[1],T[2],T[3],T[4],T[5])

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
                  

def download_waves(evtime,sec_bef_aft=[120,600],Ftype='circ',lon_lat=[120,24],range_rad=[0,6],channel=["BHZ","HHZ"],provider=["IRIS"],OUT='.'):
    import obspy
    #from obspy.clients.fdsn.mass_downloader import CircularDomain,RectangularDomain,Restrictions, MassDownloader
    from My_mass_downloader.mass_downloader import CircularDomain,RectangularDomain,Restrictions, MassDownloader
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
    #reject_channels_with_gaps=True,
    reject_channels_with_gaps=False,
    # And you might only want waveforms that have data for at least 95 % of
    # the requested time span. Any trace that is shorter than 95 % of the
    # desired total duration will be discarded.
    minimum_length=0.8,
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
    channel_priorities=channel,
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


def download_waves_catalog(cata_name,cata_filters,sec_bef_aft=[120,600],range_rad=[0,6],channel=["BHZ","HHZ"],provider=["IRIS"],waveforms_outdir='.'):
    #T,lon,lat,dep,mag=EQfilter(cata_name,BT_time=['20170607020500','20191210000000'],BT_lon=[132,133],BT_lat=[30,31],BT_dep=[0,50],BT_mag=[5.0,9.0])
    T,lon,lat,dep,mag=EQfilter(cata_name,BT_time=cata_filters['filt_times'],BT_lon=cata_filters['filt_lon'],BT_lat=cata_filters['filt_lat'],
                               BT_dep=cata_filters['filt_dep'],BT_mag=cata_filters['filt_m'])
    import time
    for i,eqT in enumerate(T):
        if i==0:
            print('--------------start downloading------------')
        if i%10==0:
            print('Now at #%d / %d'%(i,len(T)))
        
        download_waves(eqT,sec_bef_aft=sec_bef_aft,Ftype='circ',lon_lat=[lon[i],lat[i]],range_rad=range_rad,channel=channel,provider=provider,OUT=waveforms_outdir)
        time.sleep(10)


def get_waveforms(net,sta,comp,chn,t1,t2):
    #simple waveform query if net,sta,comp,chn are known
    from obspy.clients.fdsn import Client
    client = Client("IRIS")
    #t1 = UTCDateTime("2010-02-27T06:30:00.000")
    #t2 = t1 + 5
    st = client.get_waveforms(net, sta, chn, comp, t1, t2)
    return st


def get_stations(net,sta,comp="BH*,HH*",t1=UTCDateTime("2015-01-01"),t2=UTCDateTime("2015-01-02"),longitude=-155.27,latitude=19.34,minradius=0,maxradius=1):
    from obspy.clients.fdsn import Client
    client = Client("IRIS")
    #starttime = UTCDateTime("2001-01-01")
    #endtime = UTCDateTime("2001-01-02")
    inventory = client.get_stations(network=net,station=sta,channel=comp,starttime=t1,endtime=t2,longitude=longitude,latitude=latitude,minradius=minradius,maxradius=maxradius)
    sav_net=[]
    sav_sta=[]
    sav_lon=[]
    sav_lat=[]
    sav_heigh=[]
    for i in inventory:
        tmp_net = i.code
        #stations below are all in the same network
        for ii in i:
            tmp_sta = ii.code
            sav_net.append(tmp_net)
            sav_sta.append(tmp_sta)
            sav_lon.append(ii.longitude)
            sav_lat.append(ii.latitude)
            sav_heigh.append(ii.elevation)

    return sav_net,sav_sta,sav_lon,sav_lat,sav_heigh




def download_continuous_cent(net,sta,comp,chn,sampl,t1,t2,lon,lat,r1=0,r2=1):
    #download continuous data based on a center point
    sav_net,sav_sta,sav_lon,sav_lat,sav_heigh = get_stations(net=net,sta=sta,comp=comp,
                                                             t1=t1,t2=t2,
                                                             longitude=lon,latitude=lat,minradius=r1,maxradius=r2)
    #all_tr=[]
    staInfo_pre={
    'net':sav_net,
    'sta':sav_sta,
    'lon':sav_lon,
    'lat':sav_lat,
    'heigh':sav_heigh,
    }
    staInfo_fail = staInfo_pre.copy()
    staInfo_post = {
    'net':[],
    'sta':[],
    'lon':[],
    'lat':[],
    'heigh':[],
    'comp':[],
    'chn':[],
    }
    print('Total %d stations to be downloaded'%(len(sav_net)))
    for i,net in enumerate(sav_net):
        try:
            tr = get_waveforms(net,sav_sta[i],comp=comp,chn=chn,t1=t1,t2=t2)
        except:
            staInfo_fail['net'].append(net)
            staInfo_fail['sta'].append(sav_sta[i])
            staInfo_fail['lon'].append(sav_lon[i])
            staInfo_fail['lat'].append(sav_lat[i])
            staInfo_fail['heigh'].append(sav_heigh[i])
            continue #data unavailable
        print(net,sav_sta[i],comp,chn,t1,t2)
        tr.detrend()
        tr.merge() #now tr can be same station but different component
        for itr in range(len(tr)):
            if isinstance(tr[itr].data, np.ma.masked_array):
                tr[itr].data.fill_value = 0
                tr[itr].data = tr[itr].data.filled()
        tr.detrend()
        tr.filter("bandpass",freqmin=2,freqmax=7)
        #tr.filter("bandpass",freqmin=1,freqmax=8)
        tr.trim(starttime=t1-2, endtime=t2+2, nearest_sample=True, pad=True, fill_value=0)
        tr.interpolate(sampling_rate=sampl, starttime=t1)
        tr.trim(starttime=t1, endtime=t2, nearest_sample=True, pad=True, fill_value=0)
        for itr in tr:
            staInfo_post['net'].append(itr.stats.network)
            staInfo_post['sta'].append(itr.stats.station)
            staInfo_post['comp'].append(itr.stats.channel)
            staInfo_post['chn'].append(itr.stats.location)
            staInfo_post['lon'].append(sav_lon[i])
            staInfo_post['lat'].append(sav_lat[i])
            staInfo_post['heigh'].append(sav_heigh[i])
        try:
            all_tr += tr
        except:
            all_tr = tr
    return all_tr,staInfo_pre,staInfo_post,staInfo_fail


def bulk_download_continuous_cent(home,project_name,download_params,n_cores=1,waveforms_outdir='.'):
    import os
    import time
    #download all continuous data given the time range t1-t2 and lon/lat range
    net = download_params['net']
    sta = download_params['sta']
    comp = download_params['comp']
    chn = download_params['chn']
    sampl = download_params['sampl']
    cent_lon = download_params['cent_lon']
    cent_lat = download_params['cent_lat']
    min_radius = download_params['min_radius']
    max_radius = download_params['max_radius']
    t1 = download_params['t1']
    t2 = download_params['t2']
    #loop through t1-t2
    N = int((t2 - t1) / 86400) #number of days to be dealt with
    def run_loop(n):
        st = t1 + 86400 * n
        ed = st + 86400
        print('start downloading:',st,'-',ed)
        #time.sleep(5)
        #return True
        all_tr,staInfo_pre,staInfo_post,staInfo_fail = download_continuous_cent(net,sta,comp,chn,sampl,st,ed,lon=cent_lon,lat=cent_lat,r1=min_radius,r2=max_radius)
        #write the waveform file
        outdir = st.strftime('%Y%m%d%H%M%S')
        if not(os.path.exists(waveforms_outdir+'/'+outdir)):
            os.makedirs(waveforms_outdir+'/'+outdir)
        all_tr.write(waveforms_outdir+'/'+outdir+'/merged.ms',format='MSEED')
        #write station npy files in dictionary format
        np.save(waveforms_outdir+'/'+outdir+'/staInfo_pre.npy',staInfo_pre)
        np.save(waveforms_outdir+'/'+outdir+'/staInfo_post.npy',staInfo_post)
        np.save(waveforms_outdir+'/'+outdir+'/staInfo_fail.npy',staInfo_fail)
        
    #run loop by one core or multiple cores
    if n_cores==1:
        for i in range(N):
            run_loop(i)
    else:
        from joblib import Parallel, delayed
        results = Parallel(n_jobs=n_cores,prefer="threads")(delayed(run_loop)(i) for i in range(N))
        print('Parallel finished')



def make_template(df,sampling_rate,filter=[0.2,8],tcs_length=[1,9]):
    '''
        download template by event ID
        original function by: Amanda Thomas
        modified by: Tim Lin
         -merge data with interpolate data point 2020.09
         -add multiple attemps when download has unexpected issue 2020.09
         
    '''
    from obspy.clients.fdsn import Client
    from libcomcat.search import get_event_by_id
    from libcomcat.dataframes import get_phase_dataframe
    from obspy import Stream
    import time
    client = Client("IRIS")
    # make templates
    regional=df['Regional']
    eventid = regional+str(df['ID'])
    detail = get_event_by_id(eventid, includesuperseded=True)
    OT = UTCDateTime(detail.time) #event origin time
    phases = get_phase_dataframe(detail, catalog=regional)
    phases = phases[phases['Status'] == 'manual']
    phases=phases[~phases.duplicated(keep='first',subset=['Channel','Phase'])]
    st=Stream()
    tr=Stream()
    sav_net_sta_comp = []
    sav_phase = []
    sav_arr = []
    sav_travel = []
    All_info = {}
    for ii in range(len(phases)):
        elems = phases.iloc[ii]['Channel'].split('.')
        net = elems[0]
        sta = elems[1]
        comp = elems[2]
        location = elems[3]
        if location=='--':
            location = '' #sometime the location use -- instead of empty str
        Phase = phases.iloc[ii]['Phase'] #P or S wave, save this info for further relocation
        #phase = phases.iloc[ii]['Phase']
        arr = UTCDateTime(phases.iloc[ii]['Arrival Time'])
        #print(int(np.round(arr.microsecond/(1/sampling_rate*10**6))*1/sampling_rate*10**6)==1000000)
        if int(np.round(arr.microsecond/(1/sampling_rate*10**6))*1/sampling_rate*10**6)==1000000:
            arr.microsecond=0
            #arr.second=arr.second+1
            arr=arr+1
        #print('arrival Time after arrangement',arr)
        else:
            arr.microsecond=int(np.round(arr.microsecond/(1/sampling_rate*10**6))*1/sampling_rate*10**6)
        t1 = arr-tcs_length[0]
        t2 = arr+tcs_length[1]
        #try:
        i_attempt = 0 #reset number of attempt
        tr_exist = False
        while i_attempt<5:
            try:
                #tr = client.get_waveforms(net, sta, "*", comp, t1-2, t2+2) #Be careful, this may query more than 1 trace channel i.e. '00','01'...
                tr = client.get_waveforms(net, sta, location, comp, t1-2, t2+2) #Be careful, this may query more than 1 trace channel i.e. '00','01'...
                tr_exist = True
                break
            except:
                #print('Attempts %d'%(i_attempt),net, sta, "*", comp, t1-2, t2+2)
                time.sleep(2) #wait 2 sec and try again later...
                i_attempt += 1
                continue
        
        if not(tr_exist):
            #some unexpected error when attemp to client.get_waveforms
            #print("No data for "+net+" "+sta+" "+comp+" "+str(t1)+" "+str(t2))
            continue
        else:
            #print("Data available:",net, sta, location, comp, t1-2, t2+2)
            tr.merge(method=1,interpolation_samples=-1,fill_value='interpolate')
            tr.detrend()
            tr.trim(starttime=t1-2, endtime=t2+2, nearest_sample=1, pad=1, fill_value=0)
            if filter:
                tr.filter("bandpass",freqmin=filter[0],freqmax=filter[1])
            tr.interpolate(sampling_rate=sampling_rate, starttime=t1)
            tr.trim(starttime=t1, endtime=t2, nearest_sample=1, pad=1, fill_value=0)
            assert len(tr)==1, 'Unexpecting error when query:%s %s %s %s %s %s'%(net, sta, location, comp, t1-2, t2+2)
            #tr[0].stats.location = tr[0].stats.location+'.'+Phase
            st += tr.copy()
            #save name, time and "phase" info for later relocation
            #.ms only gives starttime (know arrival time) but not P or S wave
            sav_net_sta_comp.append(net+'.'+sta+'.'+comp+'.'+location)
            sav_phase.append(Phase)
            #tmp_arrT = arr.strftime('%Y-%m-%dT%H:%M:%S.%f') #arrival time in isoformat i.e. 2018-05-05T17:44:18 or 2018-05-05T17:44:23.960000
            tmp_arrT = arr.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-4] #accuracy to 0.01 sec
            #if len(tmp_arrT)==26:
            #    tmp_arrT = tmp_arrT[:-4]
            #print('Time:',tmp_arrT)
            sav_arr.append(tmp_arrT)
            sav_travel.append(arr-OT) #travel time (relative to the origin)
    #======================================IMPORTANT NOTE!!! the order of these array are NOT necessarily the same as st=============================
    #======================================When there are both P and S data in the same net.sta.chn.loc, the order will be messed up=================
    #======================================Use these info as extra caution===========================================================================
    All_info['net_sta_comp'] = np.array(sav_net_sta_comp)
    All_info['phase'] = np.array(sav_phase)
    All_info['arrival'] = np.array(sav_arr) #absolute arrival time
    All_info['travel'] = np.array(sav_travel) #phase travel time (sec)
    All_info['OT_template'] = OT.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-4]
    All_info['tcs_length'] = tcs_length
    All_info['filter'] = filter
    All_info['eqloc'] = [df['Lon'],df['Lat'],df['Depth'] ]
    All_info['eqmag'] = df['Magnitude']
    All_info['evid'] = eventid
    return st,All_info

'''
def get_order(st_out,All_info):
    from obspy import UTCDateTime
    #reorder the All_info based on the st_out
    assert len(st_out)==len(All_info['net_sta_comp']), "length of Stream and pick_info are not the same"
    sav_order = [] #new order based on the .ms file
    for n in range(len(st_out)):
        NET = st_out[n].stats.network
        STA = st_out[n].stats.station
        CHN = st_out[n].stats.channel
        LOC = st_out[n].stats.location
        temp_name = '.'.join([NET,STA,CHN,LOC]) #name from st
        #Check #1
        tmpidx = np.where( All_info['net_sta_comp']==temp_name )[0]
        #print('search %s ,found '%(temp_name),All_info['net_sta_comp'][tmpidx],All_info['phase'][tmpidx])
        if All_info['net_sta_comp'][n]!=temp_name:
            print('Order in .ms has been changed, change the pick_info file')
        assert len(tmpidx)>0, 'cannot find any name: %s'%(temp_name)
        if len(tmpidx)==1:
            sav_order.append(tmpidx[0])
        elif len(tmpidx)==2:
            info_PS = All_info['phase'][tmpidx] #one is P the other is S
            #find which one, P or S, is the nth st
            selected_tr = st_out.select(network=NET,station=STA,channel=CHN,location=LOC)
            selected_tr = selected_tr.copy()
            assert len(selected_tr)==len(tmpidx), 'number of data not consistent'
            #find T_n, and the other T1 (P or S)
            T1 = selected_tr[0].stats.starttime
            T2 = selected_tr[1].stats.starttime
            T_n = st_out[n].stats.starttime
            if T_n==T1:
                if T1<T2:
                    PS_n = 'P'
                else:
                    PS_n = 'S'
            else:
                #T_n==T2
                if T1<T2:
                    PS_n = 'S'
                else:
                    PS_n = 'P'
            #know the phase in st[n] is PS_n, now the len(tmpidx) should only be 1
            phase_list = np.array([iph.capitalize()[0] for iph in All_info['phase']])
            tmpidx = np.where( (All_info['net_sta_comp']==temp_name) & (All_info['phase']==PS_n) )[0]
            #print('tmpidx=',tmpidx); print('PS=',PS_n)
            assert len(tmpidx)==1, 'find same net_sta_comp and same phase, impossible!'
            sav_order.append(tmpidx[0])
    sav_order = np.array(sav_order)
    new_All_info = All_info.copy()
    new_All_info['net_sta_comp'] = All_info['net_sta_comp'][sav_order]
    new_All_info['phase'] = All_info['phase'][sav_order]
    new_All_info['arrival'] = All_info['arrival'][sav_order]
    new_All_info['travel'] = All_info['travel'][sav_order]
    return new_All_info
'''


def All_info_reorder(st_out,All_info):
    #reorder the info file in .ms file, this is basically the same as get_order but better method
    #reorder the All_info based on the st_out
    assert len(st_out)==len(All_info['net_sta_comp']), "length of Stream and pick_info are not the same"
    sav_order = [] #new order based on the .ms file
    info_start = np.array([UTCDateTime(it)-All_info['tcs_length'][0] for it in All_info['arrival']])
    for n in range(len(st_out)):
        NET = st_out[n].stats.network
        STA = st_out[n].stats.station
        CHN = st_out[n].stats.channel
        LOC = st_out[n].stats.location
        temp_name = '.'.join([NET,STA,CHN,LOC]) #name from st
        ST_start = st_out[n].stats.starttime #starttime
        print('Now in',temp_name,ST_start)
        #Check #1
        tmpidx = np.where( (All_info['net_sta_comp']==temp_name) & (np.abs(info_start-ST_start)==np.min(np.abs(info_start-ST_start))) )[0]
        print('Found',All_info['net_sta_comp'][tmpidx],info_start[tmpidx])
        assert len(tmpidx)==1, "find same net_sta_comp and same starttime, impossible!"
        sav_order.append(tmpidx[0]) #this should be the only one
    sav_order = np.array(sav_order)
    #do final check before return info file
    assert len(sav_order)==len(list(set(sav_order))), "Repeat index is impossible in this case!"
    new_All_info = All_info.copy()
    new_All_info['net_sta_comp'] = All_info['net_sta_comp'][sav_order]
    new_All_info['phase'] = All_info['phase'][sav_order]
    new_All_info['arrival'] = All_info['arrival'][sav_order]
    new_All_info['travel'] = All_info['travel'][sav_order]
    return new_All_info







def bulk_make_template(home,project_name,dfs,sampling_rate,filter=[0.2,8],tcs_length=[1,9]):
    from obspy import UTCDateTime
    #looping through all df, download templates, and output .txt file for checking
    OUT1 = open(home+'/'+project_name+'/waveforms_template/'+'template_summary.txt','w') #clean the previous summary file(if exist)
    OUT1.write('%s %s %s %s %s %s %s %s %s\n'%('#idx','Template_path','Nphases','Date','EVID','Lon','Lat','Dep','Mag')) #write header
    OUT1.close()
    OUT2 = open(home+'/'+project_name+'/waveforms_template/'+'template_fail.txt','w')
    OUT2.close()
    for idf in range(len(dfs)):
        DT = dfs.iloc[idf].Date+'T'+dfs.iloc[idf].Time #DateTime
        EVID = dfs.iloc[idf].Regional+dfs.iloc[idf].ID #eventID
        st,All_info = make_template(dfs.iloc[idf],sampling_rate,filter,tcs_length) #assume this always work!
        if len(st)==0:
            print('Template cannot be downloaded, probably no manual picks')
            OUT2 = open(home+'/'+project_name+'/waveforms_template/'+'template_fail.txt','a')
            OUT2.write('%d %s %s %.5f %.5f %.3f %.3f\n'%(idf,DT,EVID,dfs.iloc[idf].Lon,dfs.iloc[idf].Lat,dfs.iloc[idf].Depth,dfs.iloc[idf].Magnitude))
            OUT2.close()
            continue
        else:
            #template successfully downloaded
            outfile = home+'/'+project_name+'/waveforms_template/'+'template_%05d.ms'%(idf)
            st.write(outfile,format='MSEED') #after saving, the order may be messed up
            
            #read st data out and check the sorting
            st_out = obspy.read(outfile)
            #All_info_new = get_order(st_out,All_info)
            All_info_new = All_info_reorder(st_out,All_info)
            
            #write phase information
            outfile_info = home+'/'+project_name+'/waveforms_template/'+'template_%05d.npy'%(idf)
            np.save(outfile_info,All_info_new)
            
            #write log file
            OUT1 = open(home+'/'+project_name+'/waveforms_template/'+'template_summary.txt','a')
            OUT1.write('%d %s %d %s %s  %.5f %.5f %.3f %.3f\n'%(idf,outfile,len(st),DT,EVID,dfs.iloc[idf].Lon,dfs.iloc[idf].Lat,dfs.iloc[idf].Depth,dfs.iloc[idf].Magnitude))
            OUT1.close()



def check_order(home,project_name):
    #this is just to check if data in template.ms and template.npy has the same order
    import glob
    files_ms = glob.glob(home+'/'+project_name+'/waveforms_template/template_*.ms')
    for file_ms in files_ms:
        print('Now in:',file_ms)
        #read ms and info file
        tr = obspy.read(file_ms)
        pick_info = np.load(file_ms.replace('.ms','.npy'),allow_pickle=True)
        pick_info = pick_info.item()
        for n in range(len(tr)):
            NET = tr[n].stats.network
            STA = tr[n].stats.station
            CHN = tr[n].stats.channel
            LOC = tr[n].stats.location
            temp_name = '.'.join([NET,STA,CHN,LOC])
            #Check #1
            assert temp_name==pick_info['net_sta_comp'][n], 'order not the same! %s'%(file_ms)
            #Check #2 (if there are more than 1 name due to P,S in the same data)
            idx = np.where(pick_info['net_sta_comp']==temp_name)[0]
            if len(idx)>1:
                #check P and S, find another pair fron tr
                selected_tr = tr.select(network=NET,station=STA,channel=CHN,location=LOC)
                selected_tr = selected_tr.copy()
                assert len(selected_tr)==len(idx), 'number of data not consistent'
                if selected_tr[0].stats.starttime==tr[n].stats.starttime:
                    selected_tr1 = selected_tr[1].stats.starttime
                else:
                    selected_tr1 = selected_tr[0].stats.starttime
                #get P or S from pick_info
                PS = pick_info['phase'][n][0] #only get the first char
                if PS=='P':
                    assert tr[n].stats.starttime<selected_tr1, 'P wave NOT arrives earlier'
                else:
                    assert tr[n].stats.starttime>selected_tr1, 'S wave NOT arrives later'
    return True






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
