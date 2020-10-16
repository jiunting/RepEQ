#some data processing tools

import glob
import obspy
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


def read_detections(home,project_name,filter_params={'diff_t':60,'min_sta':5,'min_CC':0.5},fmt=1):
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
        sav_gps = {} #groups in a template detection
        Ngps = 0 #number of group
        tmp_DT = UTCDateTime("1900-01-01") #initial tmp_DT
        for line in IN1.readlines():
            #print(line)
            if line[0] == '#':
                continue #header

            #1.filter by Nstations
            if int(line.split()[2]) < filter_param['min_sta']:
                continue
            #2.filter by CC value
            if float(line.split()[1])<filter_param['min_CC']:
                continue
            
            #2018-05-01T04:25:52.525000 0.260 6 template_393
            Date_Time = line.split()[0]
            DT = UTCDateTime(Date_Time)
            templateID = int(line.split()[-1].split('_')[1]) #this will be the key for All_eqs

            if (DT-tmp_DT).total_seconds()>=filter_param['diff_t']:
                Ngps += 1
                sav_gps[Ngps] = {'DT':DT,'CC':float(line.split()[1])} #add a new group
                tmp_DT = DT
            else:
                if sav_gps[Ngps]['CC'] < float(line.split()[1]):
                    sav_gps[Ngps]['DT'] = DT#new replace the old
                    sav_gps[Ngps]['CC'] = float(line.split()[1])
                tmp_DT=DT

        All_eqs[template] = sav_gps
        clean_All_eqs[template] = np.array([j['DT'] for i,j in sav_gps.items()]) #a clean version of All_eqs
            
    return All_eqs, clean_All_eqs












