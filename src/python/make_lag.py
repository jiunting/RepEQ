import glob
import matplotlib.pyplot as plt
import datetime
import obspy
from bs4 import BeautifulSoup
import pandas as pd
from obspy.taup import TauPyModel
import numpy as np
from scipy import signal
import os
import shutil
import time

def get_staloc(net_sta_key,n_date):
    xml_file=glob.glob(n_date+'/'+'stations/'+net_sta_key+'.xml')[0]
    tmpIN1=open(xml_file,'r').read()
    soup=BeautifulSoup(tmpIN1)
    stlon=float(soup.find_all('longitude' or 'Longitude')[0].text)
    stlat=float(soup.find_all('latitude' or 'Latitude')[0].text)
    return(stlon,stlat)


def cal_CCF(data1,data2):
    #calculate normalize CCF, find max CCC, and lag idx
    tmpccf=signal.correlate(data1,data2,'full')
    auto1=signal.correlate(data1,data1,'full')
    auto2=signal.correlate(data2,data2,'full')
    tmpccf=tmpccf/np.sqrt(np.max(auto1)*np.max(auto2))
    maxCCC=np.max(tmpccf)
    lag=tmpccf.argmax()
    return(maxCCC,lag)


def cal_CCCscore(ndata,sav_ij_date,sav_CCC):
    CCCscore=np.zeros(ndata)
    for i in range(len(sav_ij_date)):
        CCCscore[sav_ij_date[i][0]]=CCCscore[sav_ij_date[i][0]]+sav_CCC[i]
        CCCscore[sav_ij_date[i][1]]=CCCscore[sav_ij_date[i][1]]+sav_CCC[i]
    return(CCCscore)

def make_legid(longname):
    sav_legend=[]
    for i in longname:
        sav_legend.append(i.split('/')[-1][:12])
    return(sav_legend)


def evloc(Cat_Date,evdate):
    '''
    Cat_Date:catalog date
    evdate:event datetime that you want to search from the catalog
    '''
    tmp_dt_all=Cat_Date-evdate# see which event it is from the catalog
    junk_time=[] #find which one has the min time difference
    for tmp_dt in tmp_dt_all:
        junk_time.append(np.abs(tmp_dt.total_seconds()))
    idx_cat=np.where(junk_time==np.min(junk_time))[0] #now looking for this event
    eqlon=float(A['eqlon'][idx_cat])
    eqlat=float(A['eqlat'][idx_cat])
    eqdep=float(A['eqdep'][idx_cat])
    return(eqlon,eqlat,eqdep)

#-------------------------------------------------------------------------------#
#pairsf='pairs_BP0.8-2_wind30s.out' #pairs file from read_log.py
#pairsf='test_pairs.out' #pairs file from read_log.py
pairsf='seq12.inp'
#pairsf='pairs_BP0.8-2_wind30s_one.out' #pairs file from read_log.py
eqpath='/Users/timlin/Documents/Project/EQrequest/Hawaii/Hawaii_ALL/' #where you put your waveform data (EQ folders)
catalogpath='/Users/timlin/Documents/Project/EQrequest/Hawaii_ALL_M3.dat' #EQ information
A=pd.read_csv(catalogpath,header=None,sep=',',names=['time','eqlat','eqlon','eqdep','eqmag','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x'],skiprows=0)
Cat_Date=pd.to_datetime(A['time'],format='%Y-%m-%dT%H:%M:%S.%fZ')#convert time to Datetime format

#filt_freq_HR=(0.8,2) #Yu & Wen, (2012)
#filt_freq_HR=(0.5,2) #correct P arrival
#filt_freq_HR=(0.5,2) #
#filt_freq_HR=(1,4)
#filt_freq_HR=(0.03,0.1)
#------------parameters for correcting P arrival-----------------#
filt_freq_HR=[(0.5,2),(0.5,2)] #set n-step Pwave corrections
p_wind=[(5,15),(2,4)]#window for Pwave correction. Seconds before(positive!) and after theoritical P arrival
CCC_thres=0.9 #threshold for repEQ from log file
CCsta_thres=0.9 #threshold for individual station
min_num=1 #at least n stations got this threshold
#-----------parameters for lag measurement after correcting P arrival----------------#
L_wind=(20,150) #Total(large window) data to be measured. Seconds before, after corrected P arrival
filt_L_wind=(0.5,2) #filter for the large window
S_wind=6 # n seconds for S(small window) of measurement each time
mov=0.2 # moving seconds
sampt=0.005 #interpolate to this interval
Write_out=True #write measured lag?
#-------------------------------------------------------------------------------#
EQfolders=glob.glob(eqpath+'*')
EQfolders.sort()

IN1=open(pairsf,'r')

if not(os.path.exists('./lag_INFO')):
    os.makedirs('lag_INFO')

for line in IN1.readlines():
    tmpelems=line.split()
    p1_str=tmpelems[0].split('-')[0]
    p2_str=tmpelems[0].split('-')[1]
    #convert to datetime
    p1=datetime.datetime.strptime(p1_str,'%Y%m%d%H%M%S')
    p2=datetime.datetime.strptime(p2_str,'%Y%m%d%H%M%S')
    count_CCC=0 #number of station has CCC above the CCC_thres
    for n_cal in range( int((len(tmpelems)-1)/2) ):
        #number of calculations available for p1-p2 pair in the pairsf
        #-1 because the first is p1-p2 id, /2 because 2 elems output/calculation
        sta=tmpelems[2*n_cal+1]
        CCC=float(tmpelems[2*n_cal+2])
        if CCC>=CCC_thres:
            count_CCC+=1
    if count_CCC<min_num:
        continue
    #else, p1-p2 it is a repeating EQ pair
    print(line.strip(),'is repEQ')
    for n_cal in range( int((len(tmpelems)-1)/2) ):
        sta=tmpelems[2*n_cal+1]
        CCsta=float(tmpelems[2*n_cal+2])
        if CCsta<CCsta_thres:
            continue #CCsta too small, skip
        #evlo1,evla1,evdp1=evloc(Cat_Date,p1) #information should all be there in the sac file. If not, use this
        #evlo2,evla2,evdp2=evloc(Cat_Date,p2)
        #get_staloc(net_sta_key,n_date)
        p1_D=0;p2_D=0
        #p1_D=[];p2_D=[]
        p1_D=obspy.read(eqpath+p1_str+'/waveforms/'+sta+'.sac')
        p2_D=obspy.read(eqpath+p2_str+'/waveforms/'+sta+'.sac')
        p1_D[0].interpolate(1/sampt,method='linear')
        p2_D[0].interpolate(1/sampt,method='linear')
        #tP is the time WRS to b not origin(O).
        #tP1=p1_D[0].stats.sac['user1'] #Pwave arrival time for pair#1 in the sac file
        #tP2=p2_D[0].stats.sac['user1'] #Pwave arrival time for pair#2 in the sac file
        tP1=p1_D[0].stats.sac['t1']
        tP2=p2_D[0].stats.sac['t1']
        Otime1=p1_D[0].stats.sac['o']
        Otime2=p2_D[0].stats.sac['o']
        tP1=tP1-Otime1
        tP2=tP2-Otime2
        #set a larger window#1 around P arrival to first order correct
        p1_D_w1=p1_D.copy() #copy from the original file since the original data will be used again later. w1:window1
        p2_D_w1=p2_D.copy()
        p1_D_w1[0].detrend('linear') #use linear detrend otherwise trend is just start point to end point
        p2_D_w1[0].detrend('linear')
        p1_D_w1[0].taper(max_percentage=0.05)
        p2_D_w1[0].taper(max_percentage=0.05)
        #bandpass filter for w#1
        p1_D_w1[0].filter('bandpass',freqmin=filt_freq_HR[0][0],freqmax=filt_freq_HR[0][1],corners=4,zerophase=True)
        p2_D_w1[0].filter('bandpass',freqmin=filt_freq_HR[0][0],freqmax=filt_freq_HR[0][1],corners=4,zerophase=True)
        #-------P arrival correction #1--------#
        #cut window for pair#1, in this case, tcs start time=folderID, p1=evid in datetime format
        p1_wind1_st=p1+datetime.timedelta(seconds=tP1-p_wind[0][0]) #for pair#1, window#1(large),starttime
        p1_wind1_ed=p1+datetime.timedelta(seconds=tP1+p_wind[0][1]) #endtime
        #print('p1=',p1)
        #print('tP1-wind=',tP1-p_wind[0][0])
        #print('1.cutting window:',p1_wind1_st,p1_wind1_ed)
        if not( p1_D_w1[0].stats.starttime<=p1_wind1_st and p1_D_w1[0].stats.endtime >= p1_wind1_ed ):
            print('window outside the timeseries:W1-p1')
            continue
        p1_D_w1_slice=p1_D_w1.slice(starttime=obspy.UTCDateTime(p1_wind1_st),endtime=obspy.UTCDateTime(p1_wind1_ed))
        #p1_D_w1_slice.plot()
        #cut window for pair#2
        p2_wind1_st=p2+datetime.timedelta(seconds=tP2-p_wind[0][0]) #for pair#2, window#1(large),starttime
        p2_wind1_ed=p2+datetime.timedelta(seconds=tP2+p_wind[0][1])
        if not( p2_D_w1[0].stats.starttime<=p2_wind1_st and p2_D_w1[0].stats.endtime >= p2_wind1_ed ):
            print('window outside the timeseries:W1-p2')
            continue
        p2_D_w1_slice=p2_D_w1.slice(starttime=obspy.UTCDateTime(p2_wind1_st),endtime=obspy.UTCDateTime(p2_wind1_ed))
        #p2_D_w1_slice.plot()
        #detrend them
        #p1_D_w1_slice.detrend()
        #p2_D_w1_slice.detrend()
        '''
        #This is a bug in the interpolate() method, use interpolate(rate,method='linear') instead
        print('maxval:',np.max(np.abs(p1_D_w1_slice[0].data)),np.max(np.abs(p2_D_w1_slice[0].data)))
        if np.max(np.abs(p1_D_w1_slice[0].data))>1e10 or np.max(np.abs(p2_D_w1_slice[0].data))>1e10:
            #p1_D=obspy.read(eqpath+p1_str+'/waveforms/'+sta+'.sac')
            #p2_D=obspy.read(eqpath+p2_str+'/waveforms/'+sta+'.sac')
            #p1_D[0].interpolate(1/sampt)
            #p2_D[0].interpolate(1/sampt)
            tmpD1=p1_D.copy()
            tmpD2=p2_D.copy()
            print(tmpD1[0].stats)
            tmpD1.plot()
            print(tmpD2[0].stats)
            tmpD2.plot()
            time.sleep(1)
        '''
        CCC,lag=cal_CCF(p1_D_w1_slice[0].data,p2_D_w1_slice[0].data)
        midd=(p2_D_w1_slice[0].stats.npts)-1  #length of b?? at this idx, refdata align with target data
        dt=p2_D_w1_slice[0].stats.delta
        shP=(lag-midd)*(dt) #convert to second (dt correction of P)
        print('shift:%s sec'%(shP))
        #if np.abs(shP)>2:
        #    plt.plot(p1_D_w1_slice[0].data,'k')
        #    plt.plot(p2_D_w1_slice[0].data,'r')
        #    plt.show()
        tP2_cor=tP2-shP #tP2_cor=P-shP
        p2_wind1_st=p2+datetime.timedelta(seconds=tP2_cor-p_wind[0][0]) #for pair#2, window#1(large),starttime
        p2_wind1_ed=p2+datetime.timedelta(seconds=tP2_cor+p_wind[0][1])
        if not( p2_D_w1[0].stats.starttime<=p2_wind1_st and p2_D_w1[0].stats.endtime >= p2_wind1_ed ):
            print('window outside the timeseries:W1-p2cor')
            continue
        p2_D_w1_slice=p2_D_w1.slice(starttime=obspy.UTCDateTime(p2_wind1_st),endtime=obspy.UTCDateTime(p2_wind1_ed)) #new slice for corrected p arrival
        if np.abs(shP)>10:
            print('large shift:',shP)
            #plt.figure()
            #plt.plot(p1_D_w1_slice[0].data,'k')
            #plt.plot(p2_D_w1_slice[0].data,'r')
            #plt.title('correction of shP='+str(shP))
            #plt.show()
        #plt.plot(p1_D_w1_slice[0].times(),p1_D_w1_slice[0].data,'k--')
        #plt.plot(p1_D_w1_slice[0].times(),p2_D_w1_slice[0].data,'r--')
        #plt.show()
        #----------------P arrival correction #1 finished, start #2 correction ------------#
        p1_D_w2=p1_D.copy() #copy from the original file since the original data will be used again later. w1:window1
        p2_D_w2=p2_D.copy()
        p1_D_w2[0].detrend('linear')
        p2_D_w2[0].detrend('linear')
        p1_D_w2[0].taper(max_percentage=0.05)
        p2_D_w2[0].taper(max_percentage=0.05)
        #bandpass filter for w#1
        p1_D_w2[0].filter('bandpass',freqmin=filt_freq_HR[0][0],freqmax=filt_freq_HR[0][1],corners=4,zerophase=True)
        p2_D_w2[0].filter('bandpass',freqmin=filt_freq_HR[0][0],freqmax=filt_freq_HR[0][1],corners=4,zerophase=True)
        #cut window for pair#1, in this case, tcs start time=folder name!!! p1=evid in datetime format
        p1_wind2_st=p1+datetime.timedelta(seconds=tP1-p_wind[1][0]) #for pair#1, window#1(large),starttime
        p1_wind2_ed=p1+datetime.timedelta(seconds=tP1+p_wind[1][1]) #endtime
        #print('2.cutting window for p1:',p1_wind2_st,p1_wind2_ed)
        if not( p1_D_w2[0].stats.starttime<=p1_wind2_st and p1_D_w2[0].stats.endtime >= p1_wind2_ed ):
            print('window outside the timeseries:W2-p1')
            continue
        p1_D_w2_slice=p1_D_w2.slice(starttime=obspy.UTCDateTime(p1_wind2_st),endtime=obspy.UTCDateTime(p1_wind2_ed))
        #p1_D_w1_slice.plot()
        #cut window for pair#2
        p2_wind2_st=p2+datetime.timedelta(seconds=tP2_cor-p_wind[1][0]) #for pair#2, window#1(large),starttime
        p2_wind2_ed=p2+datetime.timedelta(seconds=tP2_cor+p_wind[1][1])
        if not( p2_D_w2[0].stats.starttime<=p2_wind2_st and p2_D_w2[0].stats.endtime >= p2_wind2_ed ):
            print('window outside the timeseries:W2-p2')
            continue
        p2_D_w2_slice=p2_D_w2.slice(starttime=obspy.UTCDateTime(p2_wind2_st),endtime=obspy.UTCDateTime(p2_wind2_ed))
        #p2_D_w1_slice.plot()
        #detrend them
        #p1_D_w1_slice.detrend()
        #p2_D_w1_slice.detrend()
        CCC,lag=cal_CCF(p1_D_w2_slice[0].data,p2_D_w2_slice[0].data)
        midd=(p2_D_w2_slice[0].stats.npts)-1  #length of b?? at this idx, refdata align with target data
        dt=p2_D_w2_slice[0].stats.delta
        shP=(lag-midd)*(dt) #convert to second (dt correction of P)
        tP2_cor2=tP2_cor-shP #tP2_cor=P-shP
        p2_wind2_st=p2+datetime.timedelta(seconds=tP2_cor2-p_wind[1][0]) #for pair#2, window#1(large),starttime
        p2_wind2_ed=p2+datetime.timedelta(seconds=tP2_cor2+p_wind[1][1])
        if not( p2_D_w2[0].stats.starttime<=p2_wind2_st and p2_D_w2[0].stats.endtime >= p2_wind2_ed ):
            print('window outside the timeseries:W2-p2cor')
            continue
        p2_D_w2_slice=p2_D_w2.slice(starttime=obspy.UTCDateTime(p2_wind2_st),endtime=obspy.UTCDateTime(p2_wind2_ed)) #new slice for corrected p arrival
        #plt.plot(p1_D_w2_slice[0].times(),p1_D_w2_slice[0].data,'k*')
        #plt.plot(p1_D_w2_slice[0].times(),p2_D_w2_slice[0].data,'r*')
        #plt.close()
        #plt.show()
        #----------------P arrival correction #2 finished------------#
        #----------------make moving window CC-----------------------#
        L_data1=p1_D.copy()
        L_data2=p2_D.copy()
        #clean the large(L) window data
        L_data1[0].detrend('linear')
        L_data2[0].detrend('linear')
        L_data1[0].taper(max_percentage=0.05)
        L_data2[0].taper(max_percentage=0.05)
        L_data1[0].filter('bandpass',freqmin=filt_L_wind[0],freqmax=filt_L_wind[1],corners=4,zerophase=True)
        L_data2[0].filter('bandpass',freqmin=filt_L_wind[0],freqmax=filt_L_wind[1],corners=4,zerophase=True)
        #cut data by P info
        p1_L_wind_st=p1+datetime.timedelta(seconds=tP1-L_wind[0])
        p1_L_wind_ed=p1+datetime.timedelta(seconds=tP1+L_wind[1])
        #print('3.cutting window for calculation',p1_L_wind_st,p1_L_wind_ed)
        if not( L_data1[0].stats.starttime<=p1_L_wind_st and L_data1[0].stats.endtime >= p1_L_wind_ed ):
            print('window outside the timeseries:cut-Lp1')
            continue
        L_data1=L_data1.slice(starttime=obspy.UTCDateTime(p1_L_wind_st),endtime=obspy.UTCDateTime(p1_L_wind_ed))
        p2_L_wind_st=p2+datetime.timedelta(seconds=tP2_cor2-L_wind[0])
        p2_L_wind_ed=p2+datetime.timedelta(seconds=tP2_cor2+L_wind[1])
        if not( L_data2[0].stats.starttime<=p2_L_wind_st and L_data2[0].stats.endtime >= p2_L_wind_ed ):
            print('window outside the timeseries:cut-Lp2')
            continue
        L_data2=L_data2.slice(starttime=obspy.UTCDateTime(p2_L_wind_st),endtime=obspy.UTCDateTime(p2_L_wind_ed))
        #L_data1.plot()
        #L_data2.plot()
        #small window measurements
        #print(L_data1[0].stats.npts,L_data2[0].stats.npts)
        sampl=int(L_data1[0].stats.sampling_rate)
        movpts=int(mov*sampl)
        st=0
        ed=st+S_wind*sampl
        sav_lags=[]
        sav_st=[]
        while ed<len(L_data1[0].data):
            S_D1=L_data1[0].data[st:ed].copy() #small window, data 1
            S_D2=L_data2[0].data[st:ed].copy()
            S_D1=signal.detrend(S_D1,type='linear')
            S_D2=signal.detrend(S_D2,type='linear')
            #put array into obspy.Trace and use the taper function
            #very detail cal
            S_D1=obspy.Trace(S_D1)
            S_D2=obspy.Trace(S_D2)
            S_D1.taper(max_percentage=0.05)
            S_D2.taper(max_percentage=0.05)
            S_D1=S_D1.data.copy()
            S_D2=S_D2.data.copy()
            S_D1=S_D1/(np.max(np.abs(S_D1)))
            S_D2=S_D2/(np.max(np.abs(S_D2)))
            CCC,lag=cal_CCF(S_D1,S_D2)
            midd=(len(S_D2))-1  #length of b?? at this idx, refdata align with target data
            lags=(lag-midd)*(1/sampl) #lag seconds
            sav_lags.append(lags)
            sav_st.append( ((st+ed)*0.5)/sampl-L_wind[0] )
            st+=movpts
            ed+=movpts
            '''
            plt.figure()
            plt.plot(S_D1,'k')
            plt.plot(S_D2,'r')
            plt.show()
            '''
        
        plt.figure(1)
        plt.subplot(2,1,1)
        plt.plot(L_data1[0].times(),L_data1[0].data/np.max(np.min(L_data1[0].data)),'k')
        plt.plot(L_data1[0].times(),L_data2[0].data/np.max(np.min(L_data2[0].data)),'r') #L_data1[0].times()&L_data2[0].times() are the same
        if Write_out:
            OUT_lagf=open('./lag_INFO/lg_'+sta+'_'+p1_str+'_'+p2_str+'.txt','w')
            OUT_lagf.write('sec(after_alligned_P) lag(sec)\n') #header
            for nline in range(len(sav_st)):
                OUT_lagf.write('%f %f\n'%(sav_st[nline],sav_st[nline]))
            OUT_lagf.close()
        plt.title('%s CC=%3.2f  GC=%3.2f, %3.2f'%(sta,CCsta,p1_D[0].stats.sac.gcarc,p2_D[0].stats.sac.gcarc))
        plt.legend([p1_str,p2_str])
        #plt.xlim([0,50])
        plt.subplot(2,1,2)
        plt.plot(sav_st,sav_lags)
        plt.plot([sav_st[0],sav_st[-1]],[0,0],'k--')
        plt.ylim([-0.2,0.2])
        #plt.xlim([0,50])
        print('saving figure to:','./lag_INFO/lg_'+sta+'_'+p1_str+'_'+p2_str+'.png')
        plt.savefig('./lag_INFO/lg_'+sta+'_'+p1_str+'_'+p2_str+'.png')
        plt.clf()
        plt.close('all')
        #plt.show()

        #----------------moming windog CC finished-------------------#
        L_data1.clear()
        L_data2.clear()

        p1_D.clear()
        p2_D.clear()
    #break

IN1.close()










