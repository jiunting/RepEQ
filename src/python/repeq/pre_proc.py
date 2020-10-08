#some data pre-processing tools

import glob
import obspy
from obspy import UTCDateTime


def merge_daily(home,project_name,sampling_rate,filter=[0.2,8]):
    Ds = glob.glob(home+'/'+project_name+'/waveforms/'+'*000000')
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
        st.interpolate(sampling_rate=sr, starttime=t1)
        st.trim(starttime=t1, endtime=t2, nearest_sample=1, pad=1, fill_value=0)

