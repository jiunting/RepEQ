#some data (pre/post)-processing tools

import glob
import obspy
import numpy as np
from obspy import UTCDateTime


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


def merge_daily(home,project_name,sampling_rate,filter=[0.2,8],pattern='*000000'):
    '''
    merge daily/continuous data after running download_tools.download_waves_catalog
    the script loops into the directory in home/project_name/waveforms/{pattern}
    exp:
     from repeq import pre_proc
     pre_proc.merge_daily(home,project_name,sampling_rate,filter=filter,pattern='20180211*')
    Output:
     merged.ms
    '''
    Ds = glob.glob(home+'/'+project_name+'/waveforms/'+pattern)
    Ds.sort()
    print('Number of dirs=%d'%(len(Ds)))
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
        st.interpolate(sampling_rate=sampling_rate, starttime=t1)
        st.trim(starttime=t1, endtime=t2, nearest_sample=1, pad=1, fill_value=0)
        st.write(D+'/waveforms/merged.ms',format="MSEED")


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
        'min_stan':5, #number of stations
        'min_CC':0.2, #min CC value
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
        #1.filter by Nstations
        if int(len(CC))<filter_detc['min_stan']:
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
            tmp_DT = DT
            tmp_CC = np.mean(CC)
    for i_gp in range(len(sav_k)):
        new_detc[sav_k[i_gp]] = sav_gps[i_gp]
    return new_detc





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
    #calculate moving std,avg,downsamp of a daily data
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
        D = obspy.read(daily_dir+'/waveforms/merged.ms')
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
    all_sav_PS = {} #record P or S wave info for all detections
    for eq_time in detc.keys():
        #find which daily data it is
        YMD = eq_time.split('T')[0].replace('-','')
        dir = glob.glob(home+'/'+project_name+'/waveforms/'+YMD+'*')[0]
        #read merged daily data
        D = obspy.read(dir+'/waveforms/merged.ms')
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
            assert len(selected_D)==1, 'selected multiple data, something wrong'
            #cut data
            t1 = UTCDateTime(eq_time)+travel_time-cut_window[0]
            t2 = UTCDateTime(eq_time)+travel_time+cut_window[1]
            #print('    cut from ',t1.isoformat(),t2.isoformat())
            selected_D.trim(starttime=t1-2, endtime=t2+2, nearest_sample=True, pad=True, fill_value=0)
            selected_D.interpolate(sampling_rate=sampling_rate, starttime=t1)
            selected_D.trim(starttime=t1, endtime=t2, nearest_sample=1, pad=1, fill_value=0)
            St += selected_D[0]
            sav_PS.append(PS)
        #return St #test the script
        sav_tcs[eq_time] = St
        all_sav_PS[eq_time] = sav_PS
    sum_tcs_phase['detc_tcs'] = sav_tcs
    sum_tcs_phase['phase'] = all_sav_PS
    sum_tcs_phase['OT_template'] = OT_template   #origin time for template
    #return sav_tcs,all_sav_PS
    if len(detc.keys())==0:
        return False
    else:
        return sum_tcs_phase


def bulk_cut_dailydata(home,project_name,filter_detc,cut_window=[5,20]):
    '''
        bulk cut daily data from Detection results
    '''
    import glob
    #get all data path
    detc_files = glob.glob(home+'/'+project_name+'/output/Template_match/Detections/Detected_tmp_*.npy')
    detc_files.sort()

    #loop every detection
    for detc_file in detc_files:
        sum_tcs_phase = cut_dailydata(home,project_name,detc_file,filter_detc,cut_window)
        #save the results
        eqid = int(detc_file.split('/')[-1].split('.')[0].split('_')[-1])
        if sum_tcs_phase:
            np.save(home+'/'+project_name+'/output/Template_match/Data_detection_cut/Detected_data_%05d.npy'%(eqid),sum_tcs_phase)




def cal_lag(template,daily_cut,tcs_length_temp,tcs_length_daily,phase_wind):
    #--------------------------------
    #align the timeseries and calculate moving shift
    #template<obspy Trace>: template waveform (single Trace)
    #daily_cut<obspy Trace>: target waveform (single Trace) Must match the template net.sta.chan.loc and phase info
    #Note that it can only check net.sta.chan.loc is the same, make sure matchs the right phase before using this function
    #tcs_length_temp<list or np.array>:tcs_length_temp=[t1,t2], t1 sec before the arrival time(pick), t2 sec after the arrival time
    #tcs_length_daily<list or np.array>:tcs_length_daily=[tt1,tt2], tt1 sec before the arrival time, tt2 sec after the arrival time
    #phase_wind<list or array>: window for alignment (could be multiple/fine alignment)
    #--------------------------------
    #check net.sta.chan.loc name
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
    #align the phase
    phase_wind = np.array(phase_wind)
    if phase_wind.ndim==2:
        for i in range(len(phase_wind)):
            phs_wind = phase_wind[i]
            temp_t1 =

    else:
        #only align once






