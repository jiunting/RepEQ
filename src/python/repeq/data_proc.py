#some data (pre/post)-processing tools

import glob
import obspy
import numpy as np
from obspy import UTCDateTime
from scipy import signal

def cat2pd(cata_path):
    '''
        convert USGS catalog format to pandas
    '''
    import pandas as pd
    import os
    #load catalog in pandas
    #read the full path
    cat = np.genfromtxt(cata_path, delimiter=',', skip_header=0,usecols=(0,1,2,3,4,10,11), dtype=("|U22",float,float,float,float,"|U2","|U22"))
    df = pd.DataFrame(cat, columns=['ID','Time','Magnitude','Lat','Lon','Depth','Regional'])
    for ii in range(len(cat)):
        df = df.append({'ID': cat[ii][6][2:],'Time': cat[ii][0][11:], 'Magnitude': cat[ii][4], 'Date': cat[ii][0][:10],'Lat': cat[ii][1], 'Lon': cat[ii][2], 'Depth': cat[ii][3], 'Regional': cat[ii][5]}, ignore_index=True)
    return df


def cal_CCC(data1,data2):
    #calculate normalize CCF, find max CC, and its lag idx
    tmpccf=signal.correlate(data1,data2,'full')
    auto1=signal.correlate(data1,data1,'full')
    auto2=signal.correlate(data2,data2,'full')
    tmpccf=tmpccf/np.sqrt(np.max(auto1)*np.max(auto2))
    maxCCC=np.max(tmpccf)
    lag=tmpccf.argmax()
    return maxCCC,lag


def merge_daily(home,project_name,sampling_rate,filter=[0.2,8],pattern='*000000',fmt='MS'):
    '''
    merge daily/continuous data after running download_tools.download_waves_catalog
    the script loops into the directory in home/project_name/waveforms/{pattern}
    fmt= 'MS' or 'npy'
    example:
     from repeq import pre_proc
     pre_proc.merge_daily(home,project_name,sampling_rate,filter=filter,pattern='20180211*')
    Output:
     merged.ms
    '''
    Ds = glob.glob(home+'/'+project_name+'/waveforms/'+pattern)
    Ds.sort()
    print('Number of dirs=%d'%(len(Ds)))
    if fmt.upper()=='MS':
        for D in Ds:
            print('In dir:',D)
            t1 = UTCDateTime(D.split('/')[-1])
            t2 = t1 + 86400 #this is daily so exactly +86400 sec
            st = obspy.read(D+'/waveforms/*.mseed')
            st.merge(method=1,interpolation_samples=-1,fill_value='interpolate')
            st.detrend()
            if filter:
                st.filter("bandpass",freqmin=filter[0],freqmax=filter[1])
            st.trim(starttime=t1-2, endtime=t2+2, nearest_sample=1, pad=1, fill_value=0)
            st.interpolate(sampling_rate=sampling_rate, starttime=t1,method='linear')
            st.trim(starttime=t1, endtime=t2, nearest_sample=1, pad=1, fill_value=0)
            st.write(D+'/waveforms/merged.ms',format="MSEED")
    elif fmt.upper()=='NPY':
        for D in Ds:
            print('In dir:',D)
            t1 = UTCDateTime(D.split('/')[-1])
            t2 = t1 + 86400 #this is daily so exactly +86400 sec
            st = obspy.read(D+'/waveforms/*.mseed')
            st.merge(method=1,interpolation_samples=-1,fill_value='interpolate')
            st.detrend()
            if filter:
                st.filter("bandpass",freqmin=filter[0],freqmax=filter[1])
            st.trim(starttime=t1-2, endtime=t2+2, nearest_sample=1, pad=1, fill_value=0)
            st.interpolate(sampling_rate=sampling_rate, starttime=t1,method='linear')
            st.trim(starttime=t1, endtime=t2, nearest_sample=1, pad=1, fill_value=0)
            np.save(D+'/waveforms/merged.npy',st)


def read_obspy(filename):
    #obspy.read cannot read mseed larger than 2GB, use this function instead
    import io
    reclen = 512
    chunksize = 100000 * reclen # Around 50 MB
    with io.open(filename, "rb") as fh:
        while True:
            with io.BytesIO() as buf:
                c = fh.read(chunksize)
                if not c:
                    break
                buf.write(c)
                buf.seek(0, 0)
                st = obspy.read(buf)
                try:
                    allst += st
                except:
                    allst = st.copy()
    allst.merge()
    return allst





def read_detections(home,project_name,filter_params={'diff_t':60,'min_sta':5,'min_CC':0.3},fmt=1):
    '''
    fmt = 1, format in:
     #OriginTime meanCC nSTA templateIDX
     2018-04-27T04:22:16.840000 1.000 41 template_00001
     2018-04-27T05:21:44.280000 0.329 41 template_00001
    fmt = 2
     some other formats will be added in soon
    '''
    files=glob.glob(home+'/'+project_name+'/output/Template_match/Detections/Detected_tmp_*.txt')
    files.sort()
    
    All_eqs={} #initial dictionary
    clean_All_eqs={} #clean version of All_eqs

    for file in files:
        print('In:',file)
        IN1 = open(file,'r')
        templateID = int(file.split('/')[-1].split('_')[-1].split('.')[0]) #this will be the key for All_eqs
        sav_gps = {} #groups in a template detection
        Ngps = 0 #number of group
        tmp_DT = UTCDateTime("1900-01-01") #initial tmp_DT
        for line in IN1.readlines():
            if line[0] == '#':
                continue #header

            #1.filter by Nstations
            if int(line.split()[2]) < filter_params['min_sta']:
                continue
            #2.filter by CC value
            if float(line.split()[1])<filter_params['min_CC']:
                continue

            #2018-05-01T04:25:52.525000 0.260 6 template_393
            Date_Time = line.split()[0]
            DT = UTCDateTime(Date_Time)
            #templateID = int(line.split()[-1].split('_')[1]) #this will be the key for All_eqs

            if (DT-tmp_DT)>=filter_params['diff_t']:
                Ngps += 1
                sav_gps[Ngps] = {'DT':DT.datetime,'CC':float(line.split()[1])} #add a new group
                tmp_DT = DT
            else:
                if sav_gps[Ngps]['CC'] < float(line.split()[1]):
                    sav_gps[Ngps]['DT'] = DT.datetime #new replace the old
                    sav_gps[Ngps]['CC'] = float(line.split()[1])
                tmp_DT=DT
        All_eqs[templateID] = sav_gps
        clean_All_eqs[templateID] = np.array([j['DT'] for i,j in sav_gps.items()]) #a clean version of All_eqs
        IN1.close()
    np.save(home+'/'+project_name+'/output/Template_match/Detections/summaryEQ_detail.npy',All_eqs)
    np.save(home+'/'+project_name+'/output/Template_match/Detections/summaryEQ_clean.npy',clean_All_eqs)
            
    return All_eqs,clean_All_eqs



def clean_detc(detc,filter_detc):
    '''
        clean detailed repeating EQs detections by filter_detc
        example:
        filter_detc = {
        'min_stan':5, #number of non-zero CC measurements
        'min_CC':0.2, #min mean(CC) value
        'diff_t':60, #time difference between events should larger than this
        }
    '''
    new_detc = {}
    Ngps = -1 #initialize group index
    sav_gps = {} #create tmp dictionary for grouping
    sav_k = []
    #input dictionary detc and output a clean dictionary based on the filter
    tmp_DT = UTCDateTime("1900-01-01")
    tmp_CC = 0
    keys = list(detc.keys())
    keys.sort()
    for k in keys:
        #print('in k',k)
        CC = detc[k]['CC']
        CC = np.array(CC)
        #1.filter by Nstations that are not zero CC (due to data missing)
        #if int(len(CC))<filter_detc['min_stan']: #the old filter only consider number of stations but not considering missing data (zeros data and thus, zero CC)
        if len(np.where(CC!=0)[0])<filter_detc['min_stan']:
            #print('number of stations=',len(CC))
            continue
        #2.filter by meanCC value
        if np.mean(CC)<filter_detc['min_CC']:
            #print('mean CC=',np.mean(CC))
            continue
        #dealing with time
        DT = UTCDateTime(k)
        if np.abs(DT-tmp_DT)>=filter_detc['diff_t']:
            #found new group
            Ngps += 1
            #print('find new group',k)
            #new_detc[k] = detc[k]
            #sav_gps[Ngps] = {'DT':DT,'CC':np.mean(CC)} #add a new group
            sav_gps[Ngps] = detc[k] #add a new group
            sav_k.append(k) #append a templory k
            tmp_DT = DT
            tmp_CC = np.mean(CC)
        else:
            #print('replace old group',k)
            if tmp_CC<np.mean(CC): #previous CC lower than new CC
                sav_gps[Ngps] = detc[k]
                sav_k[Ngps] = k #replace the previous k
            #sav_gps[Ngps]['DT'] = DT #new replace the old
            #sav_gps[Ngps]['CC'] = np.mean(CC)
            ####removed--2020.12.28 check if these are correct####
            #tmp_DT = DT
            #tmp_CC = np.mean(CC)
            ################
    for i_gp in range(len(sav_k)):
        new_detc[sav_k[i_gp]] = sav_gps[i_gp]
    return new_detc



def clean_data_cut(D,filter_detc):
    '''
    clean data_cut from the EQs detections
    example:
    filter_detc = {
    'min_stan':5, #number of non-zero CC measurements
    'min_CC':0.2, #min mean(CC) value
    'diff_t':60, #time difference between events should larger than this
    }
    '''
    #1. get all detection tcs has CC>min_CC
    sav_ot = [] #all the key that pass the CC threshold
    for k in D['meanCC'].keys():
        if np.mean(D['meanCC'][k]) >= filter_detc['min_CC']:
            sav_ot.append(k)

    #2. dealing with Nsta
    sav_ot2 = []
    for k in sav_ot:
        if len(D['phase'][k]) >= filter_detc['min_stan']:
            sav_ot2.append(k)

    new_detc_tcs = {} #new D['detc_tcs']
    new_phase = {}
    new_meanCC = {}
    for k in sav_ot2:
        new_detc_tcs[k] = D['detc_tcs'][k]
        new_phase[k] = D['phase'][k]
        new_meanCC[k] = D['meanCC'][k]

    new_D = {}
    new_D['detc_tcs'] = new_detc_tcs
    new_D['phase'] = new_phase
    new_D['meanCC'] = new_meanCC
    new_D['OT_template'] = D['OT_template']
    return new_D

    #3. do not deal with diff time since it has been dealt before
    '''
    tmp_DT = UTCDateTime("1900-01-01")
    flag_prev_exist = 0 #previous data within dt threshold, comparison is necessarily
    tmp_CC = 0
    for k in sav_ot2:
        if np.abs(UTCDateTime(k)-tmp_DT)>=filter_detc['diff_t']:
            #find new group
            tmp_DT = UTCDateTime(k)
            tmp_CC = D['meanCC'][k]
            #check if previous data in the same group
        else:
            #start compare
            if D['meanCC'][k]>tmp_CC:
                #new replace old
                tmp_CC = D['meanCC'][k]
    '''



def clean_events_time(time,min_time=5):
    #filter event by assuming a minimum interevent time
    #time: array of UTCDateTime
    #min_time: any two events should > min_time(sec)
    sav_time = [] #new time array
    tmp_t = UTCDateTime('19000101')
    n_drop = 0
    for t in time:
        if np.abs(t-tmp_t)>min_time:
            #new event found
            sav_time.append(t)
            tmp_t = t
        else:
            n_drop += 1
            print('T1,T2 within the min_time, considering the same:',tmp_t,t)
            continue
    print('Number of event dropped:',n_drop)
    sav_time = np.array(sav_time)
    return sav_time





def make_sta_table(home,project_name,pattern='*000000'):
    #make station table in Lon Lat Name
    import glob
    from bs4 import BeautifulSoup

    #Read stations file and get stlon,stlat
    def get_staloc(date_path):
        xml_files = glob.glob(date_path+'/'+'stations/'+'*.xml')
        sav_stlon = []
        sav_stlat = []
        sav_stelev = []
        sav_stname = []
        for xml_file in xml_files:
            net_sta = xml_file.split('/')[-1].split('.')
            net_sta = net_sta[0]+'.'+net_sta[1]
            tmpIN1 = open(xml_file,'r').read()
            soup = BeautifulSoup(tmpIN1)
            stlon = float(soup.find_all('longitude' or 'Longitude')[0].text)
            stlat = float(soup.find_all('latitude' or 'Latitude')[0].text)
            stelev = float(soup.find_all('elevation' or 'Elevation')[0].text)
            sav_stlon.append(stlon)
            sav_stlat.append(stlat)
            sav_stelev.append(stelev)
            sav_stname.append(net_sta)
        sav_stlon = np.array(sav_stlon)
        sav_stlat = np.array(sav_stlat)
        sav_stelev = np.array(sav_stelev)
        sav_stname = np.array(sav_stname)
        return sav_stlon,sav_stlat,sav_stelev,sav_stname

    #search all the date in waveforms directory
    Ds = glob.glob(home+'/'+project_name+'/waveforms/'+pattern)
    Ds.sort()
    #create list without repeatly count
    all_stlon = []
    all_stlat = []
    all_stelev = []
    all_stname = []
    for D in Ds:
        sav_stlon,sav_stlat,sav_stelev,sav_stname = get_staloc(D)
        for i,i_name in enumerate(sav_stname):
            if not (i_name in all_stname):
                all_stlon.append(sav_stlon[i])
                all_stlat.append(sav_stlat[i])
                all_stelev.append(sav_stelev[i])
                all_stname.append(sav_stname[i])

    #search all the date in waveforms_template
    Ds = glob.glob(home+'/'+project_name+'/waveforms_template/'+pattern)
    Ds.sort()
    for D in Ds:
        sav_stlon,sav_stlat,sav_stelev,sav_stname = get_staloc(D)
        for i,i_name in enumerate(sav_stname):
            if not (i_name in all_stname):
                all_stlon.append(sav_stlon[i])
                all_stlat.append(sav_stlat[i])
                all_stelev.append(sav_stelev[i])
                all_stname.append(sav_stname[i])

    #write table txt file in project_name/stations
    OUT1 = open(home+'/'+project_name+'/stations/'+'stations.txt','w')
    for i in range(len(all_stname)):
        OUT1.write('%f %f %f %s\n'%(all_stlon[i],all_stlat[i],all_stelev[i],all_stname[i]))

    OUT1.close()
    return all_stlon,all_stlat,all_stelev,all_stname




def cal_moving(time,data,typ,window_pts,mov_pts):
    #calculate moving [std,avg,downsamp] of a daily data
    #typ = 'avg', 'std', or samp
    sav_time = []
    sav_data = []
    total_pts = len(data)
    #initial st,ed
    st = 0
    ed = st+window_pts
    if typ=='std':
        while ed<=total_pts:
            sav_data.append(np.std(data[st:ed]))
            #sav_time.append((T0+np.mean(time[st:ed])).datetime) #datetime format
            sav_time.append(np.mean(time[st:ed])  ) #second format
            st += mov_pts
            ed += mov_pts
    elif typ=='avg':
        while ed<=total_pts:
            sav_data.append(np.mean(data[st:ed]))
            #sav_time.append((T0+np.mean(time[st:ed])).datetime) #datetime format
            sav_time.append(np.mean(time[st:ed])  ) #second format
            st += mov_pts
            ed += mov_pts
    elif typ=='samp':
        while ed<=total_pts:
            sav_data.append(data[ int((st+ed)/2) ]) #just a sample in the middle of window (downsampling)
            #sav_time.append((T0+np.mean(time[st:ed])).datetime) #datetime format
            sav_time.append(np.mean(time[st:ed])  ) #second format
            st += mov_pts
            ed += mov_pts
    else:
        print('please specify typ=[std,avg,samp]')

    return sav_time,sav_data



def cal_moving_all(home,project_name,pattern='20*',typ='std',window_pts=45000,mov_pts=4500):
    '''
        get moving std for all the data in waveforms dir
    '''
    daily_dirs = glob.glob(home+'/'+project_name+'/waveforms/'+pattern)
    daily_dirs.sort()
    sav_all_movstd = {}
    avail_date = []
    for daily_dir in daily_dirs:
        try:
            D = obspy.read(daily_dir+'/waveforms/merged.ms')
        except:
            D = read_obspy(daily_dir+'/waveforms/merged.ms')
        T0 = D[0].stats.starttime.strftime('%Y-%m-%d')
        avail_date.append(T0)
        for ist in range(len(D)):
            #get net.sta.comp.loc
            net = D[ist].stats.network
            sta = D[ist].stats.station
            comp = D[ist].stats.channel
            loc = D[ist].stats.location
            name = '.'.join([net,sta,comp,loc])
            #calculate STD
            mov_T,mov_std = cal_moving(D[ist].times(),D[ist].data,typ,window_pts,mov_pts)
            if name in sav_all_movstd:
                sav_all_movstd[name][T0] = np.array(mov_std)
            else:
                sav_all_movstd[name] = {}
                sav_all_movstd[name][T0] = np.array(mov_std)
    #save the result in project_name/waveforms/
    np.save(home+'/'+project_name+'/waveforms/'+'cal_moving_%s.npy'%(typ),sav_all_movstd)
    np.save(home+'/'+project_name+'/waveforms/'+'avail_date.npy',avail_date)
    np.save(home+'/'+project_name+'/waveforms/'+'moving_T_%s.npy'%(typ),mov_T)
    return sav_all_movstd,avail_date,mov_T


#class waveform_QC():
#    def __init__(self,home,project_name,pattern='20*'):
#    '''
#        home: path of the home working directory <str>
#        project_name: project name <str>
#        pattern: search name pattern in the home/project_name/waveforms/ <str>
#    '''
#        self.home = home
#        self.project_name = project_name
#        self.pattern = pattern
#        self.typ = typ
#        self.window_pts = window_pts
#        self.mov_pts = mov_pts
#
#    def cal_moving(self,typ='avg',window_pts=45000,mov_pts=4500):
#    '''
#        typ: tpye for calculation 'avg' or 'std' <str>
#        window_pts: points per window <int>
#        mov_pts: point for each move <int>
#    '''
#        import glob
#        import numpy as np
#        home = self.home
#        project_name = self.project_name



def cut_dailydata(home,project_name,detc_file,filter_detc,cut_window=[5,20]):
    '''
        detc_file: detailed detection file in home/project_name/output/Template_match/Detections (either just name or full path)
        filter_detc: filter used by clean_detc
        example
        detc_file = 'Detected_tmp_00000.npy'
        filter_detc = {
        'min_stan':5, #number of stations
        'min_CC':0.2, #min CC value
        'diff_t':60, #time difference between events should larger than this
        }
        
    '''
    #load detailed Detected_tmp_xxxxxx.npy file
    try:
        #detc_file in pure name
        detc = np.load(home+'/'+project_name+'/output/Template_match/Detections/'+detc_file,allow_pickle=True)
    except:
        #detc_file in full path name
        detc = np.load(detc_file,allow_pickle=True)
    
    detc = detc.item()
    detc = clean_detc(detc,filter_detc) #detc={'net_sta_comp':['HV.JOKA.HHZ.', 'HV.KNHD.EHZ.00'...],'phase':['P','S'...],'CCC':[0.99,0.98...],'CC':[0.93,0.9...],'shift':[0.0,0.04...]}

    #get eqid and its corresponding phase info
    eqid = detc_file.split('_')[-1].split('.')[0]
    phase_info_file = home+'/'+project_name+'/waveforms_template/template_'+eqid+'.npy'
    phase_info = np.load(phase_info_file,allow_pickle=True)
    phase_info = phase_info.item() #phase_info= {'net_sta_comp':['HV.PHOD.HNZ.', 'HV.PHOD.HNE.',...] , 'phase':['P','S',..],'arrival':['2018-05-02T12:54:15.68','2018-05-02T12:54:16.44',...],'travel':[1.2,2.231,3.21...]}
    OT_template = phase_info['OT_template']

    #load template waveforms
    temp_file = home+'/'+project_name+'/waveforms_template/template_'+eqid+'.ms'
    temp = obspy.read(temp_file)


    def get_travel(phase_info,net_sta_comp,PS):
        #get travel time(sec) from phase_info file by specifing a net_sta_comp (and loc) e.g. 'HV.PHOD.HNZ.'
        #print('phaseinfo=',phase_info['net_sta_comp'])
        #print('looking for',net_sta_comp)
        #print('PS list=',phase_info['phase'])
        #print('looking for',PS)
        idx = np.where((phase_info['net_sta_comp']==net_sta_comp) & (phase_info['phase']==PS) )[0]
        assert len(idx)==1, 'only one data matches the searching net_sta_comp and phase name'
        return phase_info['travel'][idx][0]
    
    #cut_window = [1,9] #window for daily data, sec prior arrival and after arrival
    sampling_rate = temp[0].stats.sampling_rate #all the sampling rate should be same
    #---loop every detection---
    sum_tcs_phase = {} #info of tcs and PS phase
    sav_tcs = {}
    sav_CC = {}
    all_sav_PS = {} #record P or S wave info for all detections
    prev_dir = '' #record dir that previous read on. An optimal way to prevent reading D over and over again
    for i_eq_time,eq_time in enumerate(detc.keys()):
        print('in:%d / %d'%(i_eq_time,len(detc.keys())))
        #find which daily data it is
        YMD = eq_time.split('T')[0].replace('-','')
        dir = glob.glob(home+'/'+project_name+'/waveforms/'+YMD+'*')[0]
        #an optimal way to recycle the D so that it won't read D over and over again
        if dir!=prev_dir:
            #read merged daily data
            if prev_dir!='':
                D.clear() #clear the previous D to prevent memory blows up
            #start loading new daily data
            try:
                D = obspy.read(dir+'/waveforms/merged.ms')
            except:
                D = read_obspy(dir+'/waveforms/merged.ms')
        prev_dir = dir
        #get the CC value
        meanCC = np.mean(detc[eq_time]['CC'])
        sav_CC[eq_time] = meanCC
        #select net_sta_comp
        St = obspy.Stream()
        sav_PS = [] #record P or S wave info
        for ista,net_sta_comp in enumerate(detc[eq_time]['net_sta_comp']):
            #P or S?
            PS = detc[eq_time]['phase'][ista]
            #get travel time(sec) for this net_sta_comp (and loc)
            #print('Travel Time for=',net_sta_comp,PS)
            travel_time = get_travel(phase_info,net_sta_comp,PS)
            #travel_time = detc[eq_time]['travel'][ista]
            #print('    ',travel_time)
            elems = net_sta_comp.split('.')
            selected_D = D.select(network=elems[0],station=elems[1],channel=elems[2],location=elems[3])
            selected_D = selected_D.copy()
            assert len(selected_D)==1, 'selected multiple daily data, something wrong' #template might have multiple data (P/S phases) but not for daily data
            #cut data
            t1 = UTCDateTime(eq_time)+travel_time-cut_window[0]
            t2 = UTCDateTime(eq_time)+travel_time+cut_window[1]
            #print('    cut from ',t1.isoformat(),t2.isoformat())
            selected_D.trim(starttime=t1-2, endtime=t2+2, nearest_sample=True, pad=True, fill_value=0)
            selected_D.interpolate(sampling_rate=sampling_rate, starttime=t1,method='linear')
            selected_D.trim(starttime=t1, endtime=t2, nearest_sample=1, pad=1, fill_value=0)
            St += selected_D[0]
            sav_PS.append(PS)
        #return St #test the script
        sav_tcs[eq_time] = St
        all_sav_PS[eq_time] = sav_PS
    sum_tcs_phase['detc_tcs'] = sav_tcs
    sum_tcs_phase['phase'] = all_sav_PS
    sum_tcs_phase['OT_template'] = OT_template   #origin time for template
    sum_tcs_phase['meanCC'] = sav_CC
    #return sav_tcs,all_sav_PS
    if len(detc.keys())==0:
        return False
    else:
        return sum_tcs_phase


def bulk_cut_dailydata(home,project_name,filter_detc,cut_window=[5,20],overwrite=False):
    '''
        bulk cut daily data from Detection results
    '''
    import glob,os
    #get all data path
    detc_files = glob.glob(home+'/'+project_name+'/output/Template_match/Detections/Detected_tmp_*.npy')
    detc_files.sort()

    #loop every detection
    for detc_file in detc_files:
        eqid = int(detc_file.split('/')[-1].split('.')[0].split('_')[-1])
        if not overwrite:
            if os.path.exists(home+'/'+project_name+'/output/Template_match/Data_detection_cut/Detected_data_%05d.npy'%(eqid)):
                continue
        #start cut the data
        print('Start cutting based on:',detc_file)
        sum_tcs_phase = cut_dailydata(home,project_name,detc_file,filter_detc,cut_window)
        #save the results
        if sum_tcs_phase:
            np.save(home+'/'+project_name+'/output/Template_match/Data_detection_cut/Detected_data_%05d.npy'%(eqid),sum_tcs_phase)




def cal_lag(template,daily_cut,tcs_length_temp,tcs_length_daily,align_wind,measure_params):
    #--------------------------------
    #align the timeseries and calculate moving shift
    #template<obspy Trace>: template waveform (single Trace)
    #daily_cut<obspy Trace>: target waveform (single Trace) Must match the template net.sta.chan.loc and phase info
    #Note that it can only check net.sta.chan.loc is the same, make sure matchs the right phase before using this function
    #tcs_length_temp<list or np.array>:tcs_length_temp=[t1,t2], t1 sec before the arrival time(pick), t2 sec after the arrival time
    #tcs_length_daily<list or np.array>:tcs_length_daily=[tt1,tt2], tt1 sec before the arrival time, tt2 sec after the arrival time
    #align_wind<list or array>: window for alignment (could be multiple/fine alignment)
    '''
    measure_params={
    'wind':[1,1],
    'mov':0.01,
    'interp':0.01,
    'taper':0.05,     #taper percentage
    }
    
    measure_params={
    'wind':[1,9],
    'mov':0.01,
    'interp':False,
    'taper':False,     #taper percentage
    }
    '''
    #--------------------------------
    #check net.sta.chan.loc name
    import matplotlib.pyplot as plt
    NET_temp = template.stats.network
    STA_temp = template.stats.station
    CHN_temp = template.stats.channel
    LOC_temp = template.stats.location
    NET_daily = daily_cut.stats.network
    STA_daily = daily_cut.stats.station
    CHN_daily = daily_cut.stats.channel
    LOC_daily = daily_cut.stats.location
    assert '.'.join([NET_temp,STA_temp,CHN_temp,LOC_temp])=='.'.join([NET_daily,STA_daily,CHN_daily,LOC_daily]), "data does not match!"
    temp_OT = template.stats.starttime
    daily_OT = daily_cut.stats.starttime
    #initial measure_params
    if not 'interp' in measure_params:
        measure_params['interp'] = False
    if not 'taper' in measure_params:
        measure_params['taper'] = False
    
    #align the phase
    align_wind = np.array(align_wind)
    if measure_params['interp']:
        delta = measure_params['interp'] #interpolate to this time interval value
    else:
        #measure_params['interp']==False
        delta = template.stats.delta #1/sampling rate
    #----------dealing with phase alignment----------
    if np.ndim(align_wind)==2:
        #multiple stage alignment
        for i in range(len(align_wind)):
            phs_wind = align_wind[i]
            #Note that phase arrival time for template = temp_OT+tcs_length_temp[0]
            temp_t1 = temp_OT+tcs_length_temp[0]-phs_wind[0]
            temp_t2 = temp_OT+tcs_length_temp[0]+phs_wind[1]
            #Note phase arrival time for daily = daily_OT+tcs_length_daily[0]
            daily_t1 = daily_OT+tcs_length_daily[0]-phs_wind[0]
            daily_t2 = daily_OT+tcs_length_daily[0]+phs_wind[1]
            #cut data
            D_temp = template.slice(starttime=temp_t1,endtime=temp_t2)
            D_temp = D_temp.data
            D_daily = daily_cut.slice(starttime=daily_t1,endtime=daily_t2)
            D_daily = D_daily.data
            #plt.plot(D_daily,'k')
            #plt.plot(D_temp,'r')
            #plt.show()
            maxCCC,lag = cal_CCC(D_temp,D_daily)
            midd = len(D_daily)-1  #length of b, at this idx, refdata align with target data
            shft = (lag-midd)*delta #convert to second (dt correction of P)
            #print('In %d iter-shift:%s sec, CC=%f'%(i,shft,maxCCC))
            #print('phase arr=',daily_OT+tcs_length_daily[0])
            #if shft is positive, daily_cut is earlier than template. So the OT needs to be earlier,and vice versa
            daily_OT -= shft
    else:
        #only align once
        phs_wind = align_wind
        #Note that phase arrival time for template = temp_OT+tcs_length_temp[0]
        temp_t1 = temp_OT+tcs_length_temp[0]-phs_wind[0]
        temp_t2 = temp_OT+tcs_length_temp[0]+phs_wind[1]
        #Note phase arrival time for daily = daily_OT+tcs_length_daily[0]
        daily_t1 = daily_OT+tcs_length_daily[0]-phs_wind[0]
        daily_t2 = daily_OT+tcs_length_daily[0]+phs_wind[1]
        #cut data
        D_temp = template.slice(starttime=temp_t1,endtime=temp_t2)
        D_temp = D_temp.data
        D_daily = daily_cut.slice(starttime=daily_t1,endtime=daily_t2)
        D_daily = D_daily.data
        maxCCC,lag = cal_CCC(D_temp,D_daily)
        midd = len(D_daily)-1  #length of b, at this idx, refdata align with target data
        shft = (lag-midd)*delta #convert to second (dt correction of P)
        #print('In %d iter-shift:%s sec, CC=%f'%(i,shft,maxCCC))
        daily_OT -= shft #if shft is positive, daily_cut is earlier than template, vice versa
        print('shift=',shft)
    #----------dealing with phase alignment END and already got the phase arr for daily_cut----------
    temp_arr = temp_OT+tcs_length_temp[0]
    #print('#################################')
    #print('template OT-arr=',temp_OT,temp_arr)
    daily_arr = daily_OT+tcs_length_daily[0]
    wind = measure_params['wind']
    mov = measure_params['mov']
    #starting time of measured window for template
    t_st_temp = temp_arr-wind[0]
    t_ed_temp = temp_arr+wind[1]
    #print('template st-ed=',t_st_temp,t_ed_temp)
    #print('#################################')
    #for dailydata
    t_st_daily = daily_arr-wind[0]
    t_ed_daily = daily_arr+wind[1]
    sav_t = [] #relative time(sec) to arrival time
    sav_shft = []
    sav_CCC = []
    #sav_temp = []
    #sav_daily = []
    #while (t_ed_temp+wind[1]<=template.stats.endtime) and (t_ed_daily+wind[1]<=daily_cut.stats.endtime):
    while (t_ed_temp+mov<=template.stats.endtime) and (t_ed_daily+mov<=daily_cut.stats.endtime):
        #cut the data
        #print(template.stats)
        #print('cut template:',t_st_temp,t_ed_temp)
        #print(daily_cut.stats)
        #print('cut daily:',t_st_daily,t_ed_daily)
        D_temp = template.copy()
        #print('--------------------------------')
        #print('orig D_temp=',D_temp)
        D_temp.interpolate(sampling_rate=(1.0/delta),method='linear')
        #print('--interp D_temp=',D_temp)
        #D_temp.trim(starttime=t_st_temp-wind[0]-1,endtime=t_ed_temp+wind[1]+1,nearest_sample=1, pad=1, fill_value=0)
        D_temp.trim(starttime=t_st_temp-1,endtime=t_ed_temp+1,nearest_sample=1, pad=1, fill_value=0)
        #print('----set st-ed=',t_st_temp-1,t_ed_temp+1)
        #print('----After trim D_temp=',D_temp)
        #print('--------Attemp to interp st=',t_st_temp)
        #interpolate data (either new sampling or original sampling)
        D_temp.interpolate(sampling_rate=(1.0/delta),starttime=t_st_temp,method='linear') #force the starttime to be "exactly" st(no 0.0001 difference)
        D_temp.trim(starttime=t_st_temp,endtime=t_ed_temp,nearest_sample=1, pad=1, fill_value=0)
        if measure_params['taper']:
            D_temp.taper(measure_params['taper'])
    
        D_temp = D_temp.data
        
        #--dealing with daily cut---
        D_daily = daily_cut.copy()
        D_daily.interpolate(sampling_rate=(1.0/delta),method='linear')
#        print('daily data from:',D_daily.stats.starttime,D_daily.stats.endtime)
#        print('Trim st=',t_st_daily-2)
        D_daily.trim(starttime=t_st_daily-1,endtime=t_ed_daily+1,nearest_sample=1, pad=1, fill_value=0)
        #D_daily.trim(starttime=t_st_daily-wind[0],endtime=t_ed_daily+wind[1],nearest_sample=1, pad=1, fill_value=0)
        #D_daily = D_daily.slice(starttime=t_st_daily-2,endtime=t_ed_daily+2)
#        print('After trim, daily data from:',D_daily.stats.starttime,D_daily.stats.endtime)
        #interpolate data
#        print('orig_data=',daily_cut)
#        print('cut_data= ',D_daily)
#        print('interp from st=',t_st_daily)
#        sav_debug = {}
#        sav_debug['cut_data'] = D_daily
#        sav_debug['interp'] = t_st_daily
#        np.save('sav_debug.npy',sav_debug)
#        D_daily = np.load('sav_debug.npy',allow_pickle=True)
#        D_daily = D_daily.item()
#        D_daily = D_daily['cut_data']
        D_daily.interpolate(sampling_rate=(1.0/delta),starttime=t_st_daily,method='linear')
#        print('interp success!',D_daily)
        D_daily.trim(starttime=t_st_daily,endtime=t_ed_daily,nearest_sample=1, pad=1, fill_value=0)
        if measure_params['taper']:
            D_daily.taper(measure_params['taper'])
        #D_daily = daily_cut.slice(starttime=t_st_daily,endtime=t_ed_daily)
        D_daily = D_daily.data

        #measure lag
        #sav_temp.append(D_temp)
        #sav_daily.append(D_daily)
        maxCCC,lag = cal_CCC(D_temp,D_daily)
        midd = len(D_daily)-1  #length of b, at this idx, refdata align with target data
        shft = (lag-midd)*delta #convert to second (dt correction of P)
        windt = ((t_ed_daily-daily_arr)+(t_st_daily-daily_arr))/2 #current window time (WRS to arrival)
        sav_t.append(windt)
        sav_shft.append(shft)
        sav_CCC.append(maxCCC)
        #add moving time to the next iter
        t_st_temp += mov
        t_ed_temp += mov
        t_st_daily += mov
        t_ed_daily += mov
    sav_t = np.array(sav_t)
    sav_shft = np.array(sav_shft)
    sav_CCC = np.array(sav_CCC)
    return sav_t,sav_shft,sav_CCC   #,sav_temp,sav_daily


#sav_t,sav_shft,sav_CCC,sav_temp,sav_daily=cal_lag(template,daily_cut,tcs_length_temp,tcs_length_daily,align_wind,measure_params)
#sav_t,sav_shft,sav_CCC=cal_lag(template,daily_cut,tcs_length_temp,tcs_length_daily,align_wind,measure_params)
#plt.scatter(sav_t,sav_shft,c=sav_CCC,cmap='jet')
#for i in range(len(sav_temp)):
#    plt.plot(sav_temp[i]/np.max(sav_temp[i])+i,'r')
#    plt.plot(sav_daily[i]/np.max(sav_daily[i])+i,'k')

def bulk_cal_lag(home,project_name,tcs_length_temp,tcs_length_daily,align_wind,measure_params,overwrite=False,n_jobs=1,i_par=0):
    #make lag calculation from all the daily_cut data in the home/project_name/output/Template_match/Data_detection_cut
    #tcs_length_temp [t1,t2]: the same value of tcs_length used by template.Template
    #tcs_length_daily [t1,t2]: the same value of cut_window used by data_proc.bulk_cut_dailydata
    '''
    measure_params={
    'wind':[1,1],
    'mov':0.01,
    'interp':0.01,
    'taper':0.05,     #taper percentage
    }
    '''
    #overwrite=False
    #n_jobs <int>: number of multiprocessing
    #i_par <int>: the i-th partitioning. e.g. from 0 to n_jobs-1. If n_jobs=8, i_par can be 0~7 (i.e. range(8))
    
    import glob,os
    #load all daily_cut data
    daily_cuts = glob.glob(home+'/'+project_name+'/output/Template_match/Data_detection_cut/Detected_data_*.npy')
    daily_cuts.sort()
    #start data partition
    all_idx = np.arange(len(daily_cuts))
    partition_idx = np.where(all_idx%n_jobs==i_par)[0]
    daily_cuts = np.array(daily_cuts)
    daily_cuts = daily_cuts[partition_idx] #new daily_cuts based on the partitioning
    
    #initial dictionary that save all info
    #lag_measure = {} #with templateID as key
    for daily_cut in daily_cuts:
        print('-------------------------------------')
        print('Now in',daily_cut)
        tempID = int(daily_cut.split('_')[-1].split('.')[0])
        #if data already exist, skip
        if os.path.exists(home+'/'+project_name+'/output/Template_match/Measure_lag/measure_lag_temp_%05d.npy'%(tempID)):
            print('--The data %s already exist, set overwrite to skip/overwrite it.'%(home+'/'+project_name+'/output/Template_match/Measure_lag/measure_lag_temp_%05d.npy'%(tempID)))
            if not(overwrite):
                continue #pass calculation if not overwritting data
        daily_cut = np.load(daily_cut,allow_pickle=True)
        daily_cut = daily_cut.item()
        #find the template .ms file
        template = obspy.read(home+'/'+project_name+'/waveforms_template/template_%05d.ms'%(tempID))
        template_info = np.load(home+'/'+project_name+'/waveforms_template/template_%05d.npy'%(tempID),allow_pickle=True)
        template_info = template_info.item()
        template_phases = template_info['phase']
        #add key in lag_measure
        #lag_measure[tempID] = {'template_OT':template_info['OT_template'],'detc_OT':{}} #initial keys in lag_measure[tempID]
        lag_measure_sub = {'template_OT':template_info['OT_template'],'detc_OT':{}} #initial keys in lag_measure_sub
        for ik in daily_cut['detc_tcs']:
            print('  det=',ik)
            #the ith detection e.g. ik='2018-04-22T16:24:34.44'
            daily_data = daily_cut['detc_tcs'][ik]
            daily_phases = daily_cut['phase'][ik]
            #lag_measure[tempID]['detc_OT'][ik] = {}   # initial (net_sta_comp.phase) key in lag_measure[tempID]['detc_OT'][ik]
            lag_measure_sub['detc_OT'][ik] = {}
            for i_cut in range(len(daily_data)):
                D_daily = daily_data[i_cut]
                PS_daily = daily_phases[i_cut] #.capitalize()[0]. Both PS_daily and template_info['phase'] can be Pg,Sg... no problem
                #match the corresponding template data
                NET = D_daily.stats.network
                STA = D_daily.stats.station
                CHN = D_daily.stats.channel
                LOC = D_daily.stats.location
                daily_net_sta_comp = '.'.join([NET,STA,CHN,LOC])
                #lag_measure[tempID]['detc_OT'][ik][daily_net_sta_comp+'.'+PS_daily] = {} #to save measurements
                lag_measure_sub['detc_OT'][ik][daily_net_sta_comp+'.'+PS_daily] = {}
                #template selection
                #Method #1, assume order of daily cut_data, net_sta_comp, and phase is the same
                selected_idx = np.where((template_info['net_sta_comp']==daily_net_sta_comp) & (template_info['phase']==PS_daily))[0][0]
                #Method #2, If there're two data selected, assume they are one P and one S, and P is always faster.....
                selected_temp = template.select(network=NET,station=STA,channel=CHN,location=LOC)
                selected_temp = selected_temp.copy()
                #most of the case should return only 1 data, but if there's P and S in 1 station...
                if len(selected_temp)!=1:
                    daily_phase = daily_phases[i_cut]
                    t1 = selected_temp[0].stats.starttime
                    t2 = selected_temp[1].stats.starttime
                    if t2-t1>0:
                        if PS_daily.capitalize()[0]=='P':
                            selected_temp = obspy.Stream(selected_temp[0])
                            print('return first one')
                        elif PS_daily.capitalize()[0]=='S':
                            selected_temp = obspy.Stream(selected_temp[1])
                            print('return second one')
                    else:
                        if PS_daily.capitalize()[0]=='P':
                            selected_temp = obspy.Stream(selected_temp[1])
                        elif PS_daily.capitalize()[0]=='S':
                            selected_temp = obspy.Stream(selected_temp[0])
                #---End of Method#2 selection---
                #print('starttime:',template[selected_idx].stats.starttime,selected_temp[0].stats.starttime)
                assert template[selected_idx].stats.starttime==selected_temp[0].stats.starttime, 'Selection inconsistent! check the Method #1 & #2'
                #if the assert always work, delete the Method2 and only use the method1
                sav_t,sav_shft,sav_CCC = cal_lag(template[selected_idx],D_daily,tcs_length_temp,tcs_length_daily,align_wind,measure_params) #select data based on Method#1
                #sav_t,sav_shft,sav_CCC = cal_lag(selected_temp[0],D_daily,tcs_length_temp,tcs_length_daily,align_wind,measure_params) #based on Method#2
                #lag_measure[tempID]['detc_OT'][ik][daily_net_sta_comp+'.'+PS_daily]['time'] = sav_t
                #lag_measure[tempID]['detc_OT'][ik][daily_net_sta_comp+'.'+PS_daily]['shift'] = sav_shft
                #lag_measure[tempID]['detc_OT'][ik][daily_net_sta_comp+'.'+PS_daily]['CCC'] = sav_CCC
                lag_measure_sub['detc_OT'][ik][daily_net_sta_comp+'.'+PS_daily]['time'] = sav_t
                lag_measure_sub['detc_OT'][ik][daily_net_sta_comp+'.'+PS_daily]['shift'] = sav_shft
                lag_measure_sub['detc_OT'][ik][daily_net_sta_comp+'.'+PS_daily]['CCC'] = sav_CCC
        np.save(home+'/'+project_name+'/output/Template_match/Measure_lag/measure_lag_temp_%05d.npy'%(tempID),lag_measure_sub)
    if n_jobs!=1:
        return #if using parallel
    #merge all the measure_lag files into a large file
    all_files = glob.glob(home+'/'+project_name+'/output/Template_match/Measure_lag/measure_lag_temp_*.npy')
    all_files.sort()
    lag_measure = {} #with templateID as key
    for lag_file in all_files:
        #get the template ID
        tmpID = lag_file.split('/')[-1].split('_')[-1].split('.')[0]
        lag_measure_sub = np.load(lag_file,allow_pickle=True)
        lag_measure_sub = lag_measure_sub.item()
        lag_measure[tmpID] = lag_measure_sub
    np.save(home+'/'+project_name+'/output/Template_match/Measure_lag/measure_lag_all.npy',lag_measure)


def bulk_cal_lag_parallel(home,project_name,tcs_length_temp,tcs_length_daily,align_wind,measure_params,overwrite=False,n_jobs=4):
    #parallelized version of data_proc.bulk_cal_lag
    #***always use this function then data_proc.bulk_cal_lag when n_job != 1
    from joblib import Parallel, delayed
    results = Parallel(n_jobs=n_jobs,verbose=10,backend='multiprocessing')(delayed(bulk_cal_lag)(home,project_name,tcs_length_temp,tcs_length_daily,align_wind,measure_params,overwrite=overwrite,n_jobs=n_jobs,i_par=i_par) for i_par in range(n_jobs)  )
    print(results)
    #-----after finish all the lag calculation--------
    if n_jobs==1:
        return #merged file done in the bulk_cal_lag
    #merge all the measure_lag files into a large file
    all_files = glob.glob(home+'/'+project_name+'/output/Template_match/Measure_lag/measure_lag_temp_*.npy')
    all_files.sort()
    lag_measure = {} #with templateID as key
    for lag_file in all_files:
        #get the template ID
        tmpID = lag_file.split('/')[-1].split('_')[-1].split('.')[0]
        lag_measure_sub = np.load(lag_file,allow_pickle=True)
        lag_measure_sub = lag_measure_sub.item()
        lag_measure[tmpID] = lag_measure_sub
    np.save(home+'/'+project_name+'/output/Template_match/Measure_lag/measure_lag_all.npy',lag_measure)
    
    '''
    #normal way to do the job without multiprocessing
    for i_par in range(n_job):
        bulk_cal_lag(home,project_name,tcs_length_temp,tcs_length_daily,align_wind,measure_params,overwrite=False,n_jobs=n_job,i_par=i_par)
    '''






def GMD_solve(G,D):
    #calculate inversion here
    GT=np.transpose(G)
    try:
        M=np.dot(np.dot(np.linalg.inv(np.dot(GT,G)),np.transpose(G)),D)
    except:
        M=np.dot(np.dot(np.linalg.pinv(np.dot(GT,G)),np.transpose(G)),D) #pseudo inv
    return M

def cal_slope(t,y):
    #input time and lag measurement
    #output slope
    assert len(t)==len(y),'len of time and measurement should be the same'
    G = np.ones([len(t),1])
    G = np.hstack([G,t.reshape(-1,1)])
    M = GMD_solve(G,y)
    return M


def cal_accum(template_time,time1,time2,dt=3600):
    #calculate accumulated num within a time range. Called by data_visual.plot_accNumber
    time1 = UTCDateTime(time1)
    time2 = UTCDateTime(time2)
    #get the data within time range
    idx = np.where((template_time>=time1) & (template_time<=time2) )[0]
    template_time = template_time[idx]
    #moving window calculate accumulated num
    tmp1 = time1
    sav_num = []
    sav_t = []
    while tmp1+dt<time2:
        sav_num.append(len(np.where(template_time<=tmp1)[0]))
        sav_t.append(tmp1)
        tmp1 += dt
    return sav_t,sav_num


#get staChn,phs from data_cut (inputs for data_visual.plot_reptcs(home,project_name,tempID,staChn,phs,cut_window))
def get_cut_info(home,project_name,tempID):
    import numpy as np
    import os
    '''
        input tempID (e.g. '00010' stands for file Detected_data_00836.npy in home/project_name/output/Template_match/Data_detection_cut/)
        output all the available station full name (net.station.channel.loc), and phase
    '''
    tcs_cut = home+'/'+project_name+'/'+'output/Template_match/Data_detection_cut/Detected_data_'+tempID+'.npy'
    if not(os.path.exists(tcs_cut)):
        print('file:%s do not exist!'%(tcs_cut))
        return

    D = np.load(tcs_cut,allow_pickle=True)
    D = D.item()
    all_StaChnPhs = []
    #loop all the ik and save station if not exist
    for ik in D['detc_tcs'].keys():
        for i in range(len(D['detc_tcs'][ik])):
            StaChnPhs = '.'.join([D['detc_tcs'][ik][i].stats.network,D['detc_tcs'][ik][i].stats.station,D['detc_tcs'][ik][i].stats.channel,D['detc_tcs'][ik][i].stats.location])
            StaChnPhs = StaChnPhs+'_'+D['phase'][ik][i]
            if not StaChnPhs in all_StaChnPhs:
                all_StaChnPhs.append(StaChnPhs)

    #split the all_StaChnPhs into array
    A = np.array([i.split('_') for i in all_StaChnPhs])
    return A[:,0],A[:,1]

















