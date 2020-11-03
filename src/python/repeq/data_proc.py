#some data processing tools

import glob
import obspy
import numpy as np
from obspy import UTCDateTime


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











