# -*- coding: utf-8 -*-
"""
Created on Wed May 08 18:19:05 2019
    
@author: TimLin
"""


def build_TauPyModel(home,project_name,vel_mod_file,background_model='PREM'):
    #function modified from mudpy
    #vel_mod_file: the full path of .mod file
    #return a .npz file that can be read by TauPyModel
    '''
    This function will take the structure from the .mod file
    and paste it on top of a pre computed mantle structure such as PREM.
    This assumes that the .mod file provided by the user ends with a 0 thickness 
    layer on the MANTLE side of the Moho
    '''
    from numpy import genfromtxt
    from os import environ,path
    from obspy.taup import taup_create, TauPyModel
    import obspy
    import shutil
    #mudpy source folder

    #load user specified .mod infromation
    structure = genfromtxt(vel_mod_file)
    shutil.copy(vel_mod_file,home+'/'+project_name+'/structure/'+vel_mod_file.split('/')[-1])
    #load background velocity structure
    if background_model=='PREM':
        #get the background file from obspy
        bg_model_file=obspy.__path__[0]+'/taup/data/'+'prem.nd'
        #Q values
        Qkappa=1300
        Qmu=600
        #Write new _nd file one line at a time
        nd_name=path.basename(vel_mod_file).split('.')[0]
        nd_name=nd_name+'.nd'
        f=open(home+'/'+project_name+'/structure/'+nd_name,'w')
        #initalize
        ztop=0
        for k in range(len(structure)-1):
            #Variables for writing to file
            zbot=ztop+structure[k,0]
            vp=structure[k,2]
            vs=structure[k,1]
            rho=structure[k,3]

            # Write to the file
            line1=('%8.2f\t%8.5f   %7.5f   %7.5f\t%6.1f     %5.1f\n' % (ztop,vp,vs,rho,Qkappa,Qmu))
            line2=('%8.2f\t%8.5f   %7.5f   %7.5f\t%6.1f     %5.1f\n' % (zbot,vp,vs,rho,Qkappa,Qmu))
            f.write(line1)
            f.write(line2)

            #update
            ztop=zbot

        #now read PREM file libe by libne and find appropriate depth tos tart isnerting
        fprem=open(bg_model_file,'r')
        found_depth=False

        while True:

            line=fprem.readline()

            if line=='': #End of file
                break

            if found_depth==False:
                #Check that it's not a keyword line like 'mantle'
                if len(line.split())>1:

                    #not a keyword, waht depth are we at?
                    current_depth=float(line.split()[0])

                    if current_depth > zbot: #PREM depth alrger than .mod file last line
                        found_depth=True
                        f.write('mantle\n')

            #Ok you have found the depth write it to file
            if found_depth == True:
                f.write(line)

        fprem.close()
        f.close()
        # make TauPy npz
        taup_in=home+'/'+project_name+'/structure/'+nd_name
        taup_out=home+'/'+project_name+'/structure/'
        taup_create.build_taup_model(taup_in,output_folder=taup_out)
    else: #To be done later (ha)
        print('ERROR: That background velocity model does not exist')

    return home+'/'+project_name+'/structure/'+nd_name.replace('nd','npz')



def searchRepEQ(home,project_name,vel_model,cata_name,data_filters,startover=False,make_fig_CC=2,QC=True,save_note=True):
    '''
        startover=True: re-run everything, or False: check the files that already exist to see if they need to be updated
        make_fig_CC: plot the figure when CC value >= this number. set a number >1 if you dont want plot
        QC: run simple check to see if the data are just zeros 
    '''
    #import matplotlib.pyplot as plt
    import obspy
    from bs4 import BeautifulSoup
    import pandas as pd
    from obspy.taup import TauPyModel
    import numpy as np
    from scipy import signal
    import os,sys,shutil,glob,datetime
    
    def get_staloc(net_sta_key,n_date):
        xml_file=glob.glob(n_date+'/'+'stations/'+net_sta_key+'.xml')[0]
        tmpIN1=open(xml_file,'r').read()
        soup=BeautifulSoup(tmpIN1)
        stlon=float(soup.find_all('longitude' or 'Longitude')[0].text)
        stlat=float(soup.find_all('latitude' or 'Latitude')[0].text)
        return stlon,stlat

    def get_traveltime(stlon,stlat,eqlon,eqlat,eqdep,model_name='iasp91'):
        dist_degree=obspy.geodetics.locations2degrees(lat1=eqlat,long1=eqlon,lat2=stlat,long2=stlon)
        model = TauPyModel(model=model_name)
        P=model.get_travel_times(source_depth_in_km=eqdep, distance_in_degree=dist_degree, phase_list=('P','p'), receiver_depth_in_km=0)
        S=model.get_travel_times(source_depth_in_km=eqdep, distance_in_degree=dist_degree, phase_list=('S','s'), receiver_depth_in_km=0)
        return P[0].time,S[0].time,dist_degree

    def cal_CCF(data1,data2):
        #calculate normalize CCF, find max CCC, and lag idx
        tmpccf=signal.correlate(data1,data2,'full')
        auto1=signal.correlate(data1,data1,'full')
        auto2=signal.correlate(data2,data2,'full')
        tmpccf=tmpccf/np.sqrt(np.max(auto1)*np.max(auto2))
        maxCCC=np.max(tmpccf)
        lag=tmpccf.argmax()
        return maxCCC,lag

    def cal_CCCscore(ndata,sav_ij_date,sav_CCC):
        CCCscore=np.zeros(ndata)
        for i in range(len(sav_ij_date)):
            CCCscore[sav_ij_date[i][0]]=CCCscore[sav_ij_date[i][0]]+sav_CCC[i]
            CCCscore[sav_ij_date[i][1]]=CCCscore[sav_ij_date[i][1]]+sav_CCC[i]
        return CCCscore

    #make .npz file
    TauPy_name=build_TauPyModel(home,project_name,vel_model) #make .npz file

    #-------------------------------------------------------------------------------#
    eqpath=home+'/'+project_name+'/waveforms/' #where you put your waveform data (EQ folders)
    catalogpath=home+'/'+project_name+'/catalog/'+cata_name #EQ information

    repeq_dir=home+'/'+project_name+'/output/logs' #This is the output directory for .log files
    A=pd.read_csv(catalogpath,header=None,sep=',',names=['time','eqlat','eqlon','eqdep','eqmag','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x'],skiprows=0)
    Cat_Date=pd.to_datetime(A['time'],format='%Y-%m-%dT%H:%M:%S.%fZ')#convert time to Datetime format

    #filt_freq_HR=(0.5,2)
    #filt_freq_HR=(0.8,2) #Yu & Wen, (2012)
    #filt_freq_HR=(1,4)
    #filt_freq_HR=(0.03,0.1)
    #p_wind=(0,30)#seconds before and after theoritical P arrival

    filt_freq_HR=data_filters['freq']
    p_wind=data_filters['window']
    max_sepr=data_filters['max_sepr']
    #-------------------------------------------------------------------------------#
    EQfolders=glob.glob(eqpath+'*')
    EQfolders.sort()

    #make calculation note
    if save_note:
        OUT2=open(repeq_dir+'/'+project_name+'_searchRepEQ.note','w')

    ###Search through all the event directories and find all the stations (NET.CODE) ###
    same_net_sta=[]
    sav_info={} #dictionary that save the net_station time
    for EQfolder in EQfolders:
        print('In:',EQfolder)
        sacs=glob.glob(EQfolder+'/waveforms/*.mseed')
        for sac in sacs:
            net=sac.split('/')[-1].split('.')[0]
            sta=sac.split('/')[-1].split('.')[1]
            print(net+'_'+sta)
            net_sta=net+'.'+sta
            try:
                sav_info[net_sta].append(EQfolder)
            except:
                sav_info[net_sta]=[EQfolder]
            if not(net_sta in same_net_sta):
                same_net_sta.append(net_sta)

    same_net_sta.sort()
    print(same_net_sta) #all the different stations in folder(eqpath)

    # This is a new calculation, re-calculation, or overwrite everything
    if not os.path.isdir(repeq_dir):
        #repeq_dir has been deleted
        os.mkdir(repeq_dir)

    #-----looping all the stations,data and make CC calculations---------------
    for net_sta_key in same_net_sta:
        #for net_sta_key in ['HV.TOUO']:
        #check_file_flag=0 #if the log file exist, and wanna check into the information
        if os.path.exists(repeq_dir+'/'+net_sta_key+'.log'):
            if not startover:
                continue #just skip the file if it's already exist
                '''
                #but not just skip the file, read the file and check if all the calculations are there
                UPD1=open(repeq_dir+'/'+net_sta_key+'.log','r') #open the pre-existing file and check
                Chk_lines=UPD1.read()
                UPD1.close()
                OUT1=open(repeq_dir+'/'+net_sta_key+'.log','a')
                check_file_flag=1 #check into the file and see if theres any missing
                '''
            else:
                #overwrite the existing .log file
                OUT1=open(repeq_dir+'/'+net_sta_key+'.log','w')
        else:
            #log file not exist, make a new file
            OUT1=open(repeq_dir+'/'+net_sta_key+'.log','w')

        if save_note:
            OUT2.write('--------Starting Sta:%s --------\n'%(net_sta_ley))

        sav_Date=[]
        n=0
        pre_sampr=[]      #make sure sampling rate for a same station always keeps the same. If not, don't use it
        pre_date=[]       #EQID for the reference sampling rate (i.e. eqid for pre_sampr) 
        sav_data=[]       #save all the cutted P waves (same station, different eqs), loop them by Cn,2  later to see if it is repeat EQ
        sav_t=[]          #saved time corresponding to sav_data, so plt.plot(sav_t[0],sav_data[0]) should make sense
        sav_OT=[]         #event origin time with respect to KZtime
        sav_Parrivl=[]    #saved P travel time
        sav_evlon=[]
        sav_evlat=[]
        sav_legend=[] #
        for n_date in sav_info[net_sta_key]:
            #for the same net_sta, different date
            YMDHms=n_date.split('/')[-1]
            Y=int(YMDHms[:4])
            M=int(YMDHms[4:6])
            D=int(YMDHms[6:8])
            H=int(YMDHms[8:10])
            m=int(YMDHms[10:12])
            s=int(YMDHms[12:14])
            print(Y,M,D,H,m,s)
            Date=datetime.datetime(Y,M,D,H,m,s) #event start time
            tmp_dt_all=Cat_Date-Date# see which event it is from the catalog
            junk_time=[] #find which one has the min time difference
            for tmp_dt in tmp_dt_all:
                junk_time.append(np.abs(tmp_dt.total_seconds()))
            idx_cat=np.where(junk_time==np.min(junk_time))[0][0] #now looking for this event, catalog may have duplicated, use only the first one!
            eqlon=float(A['eqlon'][idx_cat])
            eqlat=float(A['eqlat'][idx_cat])
            eqdep=float(A['eqdep'][idx_cat])
            if eqdep<0.0:
                if save_note:
                    OUT2.write(' -dropping negative eqdep: %s\n'%(n_date))
                continue
            #read the data path by obspy
            data=obspy.read(n_date+'/waveforms/'+net_sta_key+'*.mseed')
            if QC:
                #quick quality check 
                dlen=len(data[0].data)
                small_wind_std=np.std(data[0].data[:dlen//10]) #first 10% of the whole window
                big_wind_std=np.std(data[0].data) #whole window
                if small_wind_std<1e-9 or big_wind_std<1e-9:
                    if save_note:
                        OUT2.write(' -dropping no data from QC: %s\n'%(n_date))
                    continue #no data in small or big window 
                if (small_wind_std/big_wind_std)>2.0 :
                    if save_note:
                        OUT2.write(' -dropping no P/S waveforms from QC: %s\n'%(n_date))
                    continue #small window has larger std than big window, probably no waveforms
            tmp_sampr=data[0].stats.sampling_rate #check sampling rate
            if pre_sampr==[]:
                pre_sampr=tmp_sampr #set the sampr for first station
                pre_date=n_date
            if tmp_sampr!=pre_sampr:
                if save_note:
                    OUT2.write(' -dropping sampling rate inconsistent: %s=%f; %s=%f\n'%(pre_date,pre_sampr,n_date,tmp_sampr))
                #print('Sampling rate inconsistent')
                continue
            #-------------get info------------
            stlon,stlat=get_staloc(net_sta_key,n_date)
            #get travel time
            #tP,tS,GCARC=get_traveltime(stlon,stlat,eqlon,eqlat,eqdep,model_name='iasp91')
            tP,tS,GCARC=get_traveltime(stlon,stlat,eqlon,eqlat,eqdep,model_name=TauPy_name)
            #---------hang on, write the information in the sac file------
            #make sac a dictionary
            Otime=obspy.UTCDateTime(Date)-data[0].stats.starttime #event origin time with respect to KZTime/Date
            tP=tP+Otime #P arrival with respect to KZTime/Date
            tS=tS+Otime #S arrival
            '''
                user1: P travel time from source to station
                user2: S travel time from origin to station
            '''
            data[0].stats.update({'sac':{'t1':tP,'t2':tS,'o':Otime,'user1':tP-Otime,'user2':tS-Otime,'stlo':stlon,'stla':stlat,'evlo':eqlon,'evla':eqlat,'evdp':eqdep,'gcarc':GCARC}})
            data.write(n_date+'/waveforms/'+net_sta_key+'.sac',format='SAC')
            #--------ok, continue to process the data--------
            data.detrend('linear')
            data.taper(max_percentage=0.05)
            t=data[0].times()
            #bandpass filter
            data.filter('bandpass',freqmin=filt_freq_HR[0],freqmax=filt_freq_HR[1],corners=4,zerophase=True)
            y=data[0].data
            if np.isnan(y).any():
                #print('detect Nan value in the data, continue to the next')
                if save_note:
                    OUT2.write(' -dropping Nan value in the data: %s\n'%(n_date))
                continue
            #cut data by P-arrival
            idxt=np.where( (tP-p_wind[0]<=t) & (t<=tP+p_wind[1]) )[0] #find the index of P-waves
            Pwave_t=t[idxt]
            Pwave_y=y[idxt]/np.max(np.abs(y[idxt]))
            sav_Date.append(Date)
            sav_data.append(Pwave_y) #sav_data: cutted P wave
            sav_t.append(Pwave_t)    #sav_t: save Pwave time. 0 at the KZtime.  plt.plot(sav_t[0],sav_data[0]) should make sense. t=0 at kztime
            sav_Parrivl.append(tP)   #P wave travel time
            sav_OT.append(Otime)
            #sav_data_long.append(y)
            #sav_t_long.append(t)
            #save earthquake location, for filtering
            sav_evlon.append(eqlon)
            sav_evlat.append(eqlat)
            sav_legend.append(n_date.split('/')[-1])
            #        plt.plot(Pwave_t,Pwave_y)
            #        plt.plot(stlon,stlat,'^')
            #        plt.plot(eqlon,eqlat,'r*')
            #        plt.text(stlon,stlat+0.08*n,net_sta_key+'_'+n_date.split('/')[-1][:4],rotation=0,fontsize=6)
            #clear the memory
            data.clear()
            n+=1
        #make cross correlation for any two date
        print('making CC for:%s, Total data:%s'%(net_sta_key,len(sav_data)))
        sav_ij_date=[]
        sav_CCC=[] #save the CCC for all ij_date
        for idate in range(len(sav_data)-1):
            for jdate in range(len(sav_data)):
                #print('CC for',idate,jdate)
                '''
                if check_file_flag:
                    #log file exist and you must have already read it into a large string, right?
                    if '%s-%s'%(sav_legend[idate],sav_legend[jdate]) in Chk_lines:
                        continue
                '''
                eqdist_degree=obspy.geodetics.locations2degrees(lat1=sav_evlat[idate],long1=sav_evlon[idate],lat2=sav_evlat[jdate],long2=sav_evlon[jdate])
                if (eqdist_degree>max_sepr):
                    if save_note:
                        OUT2.write(' -dropping events too far : %s - %s\n'%(sav_Date[idate],sav_Date[jdate]))
                    continue #events are too far
                if jdate<=idate:
                    #if save_note:
                    #    OUT2.write(' -dropping redundent calculation : %s - %s\n'%(sav_Date[idate],sav_Date[jdate]))
                    continue #skip the repeat calculation
                CCC,lag=cal_CCF(sav_data[idate],sav_data[jdate])
                if np.isnan(CCC):
                    if save_note:
                        OUT2.write(' -dropping unknow reason causing CCC = nan : %s - %s\n'%(sav_Date[idate],sav_Date[jdate]))
                    continue #some wried(e.g. nan) value in data
                sav_CCC.append(np.max(CCC))
                sav_ij_date.append((idate,jdate))
                #Output as a file
                OUT1.write('%s-%s %s %f\n'%(sav_legend[idate],sav_legend[jdate],net_sta_key,CCC))
                #Output figure, make figures for debug
                if CCC>=make_fig_CC:
                    import matplotlib.pyplot as plt
                    plt.figure()
                    plt.plot(sav_t[idate]-sav_OT[idate],sav_data[idate],'k',linewidth=1)
                    plt.plot(sav_t[jdate]-sav_OT[jdate],sav_data[jdate],linewidth=1)
                    plt.xlabel('Origin time(s)',fontsize=15)
                    plt.title(net_sta_key+' CC=%3.2f (unshifted)'%(CCC))
                    plt.legend([sav_legend[idate],sav_legend[jdate]])
                    plt.savefig(repeq_dir+'/'+net_sta_key+'.'+sav_legend[idate]+'_'+sav_legend[jdate]+'.png')
                    plt.close() #close figure, don't want to show
                    #plt.legend(use_legend)
        OUT1.close()
    #-----looping all the stations,data and make CC calculations END---------------



def read_logs(home,project_name):
    import glob
    log_dir=home+'/'+project_name+'/output/logs'
    logs=glob.glob(log_dir+'/'+'*.log')
#    outname='%s_freq%.3f-%.3f_wind%d-%d.summary'%(project_name,data_filters['freq'][0],data_filters['freq'][1],data_filters['window'][0],data_filters['window'][1])
    outname='%s.summary'%(project_name)
#    if outdir:
#        try:
#            if outdir[-1]=='/':
#                outdir=outdir[:-1]
#            OUT1=open(outdir+'/'+outname,'w')
#        except:
#            print('Output directory/name:%s do not exist!'%(outdir))
#    else:
    OUT1=open(home+'/'+project_name+'/output/logs/'+outname,'w')

    #####save logs into a very large dictionary#####
    Large_D={}
    #nsta=0
    for logname in logs:
        net_sta=logname.split('/')[-1].split('.')[0]+'.'+logname.split('/')[-1].split('.')[1]
        print(net_sta)
        Large_D[net_sta]=[] #create new key
        IN1=open(logname,'r')
        for line in IN1.readlines():
            Large_D[net_sta].append(line.strip())
        IN1.close()

    #all the keys with net_staname
    ref_keys=list(Large_D.keys())
    sav_D={}
    for nkey in range(len(ref_keys)):
        print('-Dealing with',ref_keys[nkey])
        total_lines=len(Large_D[ref_keys[nkey]])
        for nline,line in enumerate(Large_D[ref_keys[nkey]]):
            if nline%100==0:
                print('n_line=',nline,'out of',total_lines)
            p1p2=line.split()[0] # p1-p2
            CCC=float(line.split()[-1])
            if p1p2 in sav_D:
                sav_D[p1p2].append(ref_keys[nkey])
                sav_D[p1p2].append(CCC)
            else:
                sav_D[p1p2]=[]
                sav_D[p1p2].append(ref_keys[nkey])
                sav_D[p1p2].append(CCC)

    #sort the order
    all_keys=list(sav_D.keys())
    all_keys.sort()

    #Write sav_D to output
    for all_key in all_keys:
        OUT1.write('%s '%(all_key))
        for i_elem in range( int(len(sav_D[all_key])/2) ):
            OUT1.write('%s %4.2f '%(sav_D[all_key][2*i_elem],sav_D[all_key][2*i_elem+1] ))
        OUT1.write('\n')

    OUT1.close()


def sequence(home,project_name,seq_filters):
    #make sequence file based on parameters setting
    import datetime
    import numpy as np
    import glob
    summary_name='%s.summary'%(project_name)
    OUTfile=open(home+'/'+project_name+'/output/logs/'+'%s.sequence'%(project_name),'w')
    INfile=open(home+'/'+project_name+'/output/logs/'+summary_name,'r')
    
    min_nsta_HiCC=seq_filters['min_nsta']  #at least x stations has CC greater than the below threshold
    min_CC=seq_filters['min_CC']
    time_sep=seq_filters['time_sep'] #p1-p2 needs to be separated by at least x sec
    
    def det_repeq(line,min_nsta_HiCC,min_CC,time_sep)->'Boolean':
        #determine whether p1-p2 is a pair
        elems=line.split()
        p1p2=elems[0]
        t1=datetime.datetime.strptime(p1p2.split('-')[0],'%Y%m%d%H%M%S')
        t2=datetime.datetime.strptime(p1p2.split('-')[1],'%Y%m%d%H%M%S')
        if (t2-t1).total_seconds()<time_sep:
            return False,None,None
        count_CC=0
        count_sta=0
        for n_measu in range(int((len(elems) - 1) / 2)):
            sta=elems[1 + n_measu * 2]
            CC=float(elems[2 + n_measu * 2])
            count_sta += 1
            if CC>=min_CC:
                count_CC += 1
        if (count_CC>=min_nsta_HiCC):
            return True,p1p2.split('-')[0],p1p2.split('-')[1]
        else:
            return False,None,None

    EQseq=[]
    nseq=0 #number of sequences
    for line in INfile.readlines():
        isrepEQ,p1,p2 = det_repeq(line,min_nsta_HiCC,min_CC,time_sep)
        if isrepEQ:
            nseq += 1
            comp_p1p2={allEQ:i for i,subset in enumerate(EQseq) for allEQ in subset} #dict of {'eqid':num of seq}
            #p1 exist, p2 not
            if (p1 in comp_p1p2):
                #check if p2 exist, if not, append p2 also
                if not(p2 in comp_p1p2):
                    EQseq[comp_p1p2[p1]].append(p2)
                    continue
                else:
                    #both p1,p2 exist
                    continue
    
            #p2 exist, p1 not
            if (p2 in comp_p1p2):
                #check if p1 exist, if not, append p1 also
                if not(p1 in comp_p1p2):
                    EQseq[comp_p1p2[p2]].append(p1)
                    continue

            #p1, p2 does not exist
            if (not(p1 in comp_p1p2)) and (not(p2 in comp_p1p2)):
                EQseq.append([p1,p2])
                continue

    INfile.close()
    #write EQ sequences to file
    for nseq in EQseq:
        nseq.sort()
        for neq in nseq:
            OUTfile.write('%s '%(neq))
        OUTfile.write('\n')

    OUTfile.close()



def measure_lag(home,project_name,lag_params,sequence_file,cata_name):
    import glob
    import datetime
    import obspy
    from bs4 import BeautifulSoup
    import pandas as pd
    import numpy as np
    from scipy import signal
    import os
    import shutil
    import time
    import matplotlib
    matplotlib.use('pdf') #instead using interactive backend
    import matplotlib.pyplot as plt


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

    def read_summary(summary_path):
        #make summary into a dictionary
        IN1=open(summary_path,'r')
        pairs={}
        for line in IN1.readlines():
            elems=line.split()
            pairs[elems[0]]=line.strip()
        IN1.close()
        return pairs

    #-------------------------------------------------------------------------------#
    #pairsf='pairs_BP0.8-2_wind30s.out' #pairs file from read_log.py
    #pairsf='test_pairs.out' #pairs file from read_log.py
    #pairsf='seq12.inp'
    # pairsf='/Users/timlin/Documents/Project/TestREPEQ/QQQ/output/logs/QQQ.sequence'
    pairsf=home+'/'+project_name+'/'+'output'+'/'+'logs'+'/'+sequence_file
    #pairsf='pairs_BP0.8-2_wind30s_one.out' #pairs file from read_log.py
    #eqpath='/Users/timlin/Documents/Project/EQrequest/Hawaii/Hawaii_ALL/' #where you put your waveform data (EQ folders)
    # eqpath='/Users/timlin/Documents/Project/TestREPEQ/QQQ/waveforms/' #where you put your waveform data (EQ folders)
    eqpath=home+'/'+project_name+'/'+'waveforms/'
    # summary_path='/Users/timlin/Documents/Project/TestREPEQ/QQQ/output/logs/QQQ.summary'
    summary_path=home+'/'+project_name+'/'+'output'+'/'+'logs'+'/'+project_name+'.summary'
    #catalogpath='/Users/timlin/Documents/Project/EQrequest/Hawaii_ALL_M3.dat' #EQ information
    # catalogpath='/Users/timlin/Documents/Project/TestREPEQ/QQQ/catalog/area1.cat' #EQ information
    catalogpath=home+'/'+project_name+'/'+'catalog'+'/'+cata_name
    A=pd.read_csv(catalogpath,header=None,sep=',',names=['time','eqlat','eqlon','eqdep','eqmag','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x'],skiprows=0)
    Cat_Date=pd.to_datetime(A['time'],format='%Y-%m-%dT%H:%M:%S.%fZ')#convert time to Datetime format

    #filt_freq_HR=(0.8,2) #Yu & Wen, (2012)
    #filt_freq_HR=(0.5,2) #correct P arrival
    #filt_freq_HR=(0.5,2) #
    #filt_freq_HR=(1,4)
    #filt_freq_HR=(0.03,0.1)
    #------------parameters for correcting P arrival-----------------#
    filt_freq_HR=lag_params['filt_freq_HR'] #set n-step Pwave corrections
    p_wind=lag_params['p_wind']  #window for Pwave correction. Seconds before(positive!) and after theoritical P arrival
    CCC_thres=lag_params['CCC_thres'] #threshold for repEQ from log file
    CCsta_thres=lag_params['CCsta_thres'] #threshold for individual station
    min_num=lag_params['min_num']  #at least n stations got this threshold
    #-----------parameters for lag measurement after correcting P arrival----------------#
    L_wind=lag_params['L_wind']  #Total(large window) data to be measured. Seconds before, after corrected P arrival
    filt_L_wind=lag_params['filt_L_wind'] #filter for the large window
    S_wind=lag_params['S_wind'] # n seconds for S(small window) of measurement each time
    mov=lag_params['mov'] # moving seconds
    sampt=lag_params['sampt'] #interpolate to this interval
    Write_out=lag_params['Write_out'] #write measured lag?
    #-------------------------------------------------------------------------------#

    EQfolders=glob.glob(eqpath+'*')
    EQfolders.sort()
    IN1=open(pairsf,'r')

    #load summary file into dictionary for later check
    Suminfo=read_summary(summary_path) #Suminfo is a dictionary with p1-p2 keys and p1-p2 measurements content

    for line in IN1.readlines():
        print('----------Starting New Line--------------')
        print(line)
        tmpelems_seq=line.split()
        #measure the lag between each other p1 p2 p3 p4......
        for i in range(len(tmpelems_seq)-1):
            for j in range(i+1,len(tmpelems_seq)):
                print('Now dealing with i,j',i,j)
                p1_str=tmpelems_seq[i]
                p2_str=tmpelems_seq[j]
                #convert to datetime
                p1=datetime.datetime.strptime(p1_str,'%Y%m%d%H%M%S')
                p2=datetime.datetime.strptime(p2_str,'%Y%m%d%H%M%S')
                tmp_pair=p1_str+'-'+p2_str # find '20180618133904-20180621043156', for example, in the .summary file
                if not tmp_pair in  Suminfo:
                    continue
                pair_info=Suminfo[tmp_pair]
                tmpelems=pair_info.split() #tmpelems is, for example ['20180618133904-20180621043156', 'HV.ERZ4', '0.96', 'HV.ERZ2', '0.67']
                '''
                #To check again (or add filter) if they are a repeating EQ pair
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
                print(pair_info,'is repEQ')
                '''
                for n_cal in range( int((len(tmpelems)-1)/2) ):
                    sta=tmpelems[2*n_cal+1]
                    CCsta=float(tmpelems[2*n_cal+2])
                    if CCsta<CCsta_thres:
                        continue #low CC value at this station (maybe bad data), skip this station
                    print(' Sta:',sta)
                    #evlo1,evla1,evdp1=evloc(Cat_Date,p1) #information should all be there in the sac file. If not, use this
                    #evlo2,evla2,evdp2=evloc(Cat_Date,p2)
                    #get_staloc(net_sta_key,n_date)
                    p1_D=0;p2_D=0 #reset obspy data stream, or use p1_D=[];p2_D=[]
                    p1_D=obspy.read(eqpath+p1_str+'/waveforms/'+sta+'.sac') #p1_D: data stream for p1
                    p2_D=obspy.read(eqpath+p2_str+'/waveforms/'+sta+'.sac')
                    if p1_D[0].stats.delta!=p2_D[0].stats.delta:
                        continue #sampling rate are different, did the station change?
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
                    CCC_correct1,lag=cal_CCF(p1_D_w1_slice[0].data,p2_D_w1_slice[0].data)
                    midd=(p2_D_w1_slice[0].stats.npts)-1  #length of b?? at this idx, refdata align with target data
                    dt=p2_D_w1_slice[0].stats.delta
                    shP=(lag-midd)*(dt) #convert to second (dt correction of P)
                    print(' -1st correction shift:%s sec'%(shP))
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
                    CCC_correct2,lag=cal_CCF(p1_D_w2_slice[0].data,p2_D_w2_slice[0].data)
                    midd=(p2_D_w2_slice[0].stats.npts)-1  #length of b?? at this idx, refdata align with target data
                    dt=p2_D_w2_slice[0].stats.delta
                    shP=(lag-midd)*(dt) #convert to second (dt correction of P)
                    print(' -2nd correction shift:%s sec'%(shP))
                    tP2_cor2=tP2_cor-shP #tP2_cor=P-shP , second correction
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
                    sav_st=[] #sec for measured window (P at 0sec)
                    sav_windCCC=[] #save CCC for measured windows
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
                        sav_windCCC.append(CCC)
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
                    plt.plot(L_data1[0].times()-L_wind[0],L_data1[0].data/np.max(np.min(L_data1[0].data)),'k')#now 0 sec is the Parrival
                    plt.plot(L_data1[0].times()-L_wind[0],L_data2[0].data/np.max(np.min(L_data2[0].data)),'r') #L_data1[0].times()&L_data2[0].times() are the same
                    if Write_out:
                        # OUT_lagf=open('./lag_INFO/lg_'+sta+'_'+p1_str+'_'+p2_str+'.txt','w')
                        OUT_lagf=open(home+'/'+project_name+'/'+'output'+'/'+'lags'+'/'+'lg_'+sta+'_'+p1_str+'_'+p2_str+'.txt','w')
                        OUT_lagf.write('#sec(after_alligned_P) lag(sec) windowCCC CCC_corr1 CCC_corr2\n') #header, P = 0 sec
                        for nline in range(len(sav_st)):
                            OUT_lagf.write('%f %f %f %f %f\n'%(sav_st[nline],sav_lags[nline],sav_windCCC[nline],CCC_correct1,CCC_correct2))
                        OUT_lagf.close()
                    plt.yticks([],[])
                    ax=plt.gca()
                    ax.tick_params(labelbottom=False)
                    plt.xlim([-1*L_wind[0],L_wind[1]])
                    plt.title('%s CC=%3.2f  GC=%3.2f, %3.2f'%(sta,CCsta,p1_D[0].stats.sac.gcarc,p2_D[0].stats.sac.gcarc))
                    #plt.legend([p1_str,p2_str])
                    plt.legend([p1.strftime('%Y-%m-%d %H:%M:%S'),p2.strftime('%Y-%m-%d %H:%M:%S')])
                    #plt.xlim([0,50])
                    plt.subplot(2,1,2)
                    plt.plot(sav_st,sav_lags,'b')
                    plt.plot([sav_st[0],sav_st[-1]],[0,0],'k--')
                    plt.ylim([-0.2,0.2])
                    plt.xlim([-1*L_wind[0],L_wind[1]])
                    plt.xlabel('Seconds since p-wave arrival',fontsize=14)
                    plt.ylabel('Delay(sec)',fontsize=14)
                    #plt.xlim([0,50])
                    print('  saving figure to:',home+'/'+project_name+'/'+'output'+'/'+'lags'+'/'+'lg_'+sta+'_'+p1_str+'_'+p2_str+'.png')
                    # plt.savefig('./lag_INFO/lg_'+sta+'_'+p1_str+'_'+p2_str+'.png')
                    plt.savefig(home+'/'+project_name+'/'+'output'+'/'+'lags'+'/'+'lg_'+sta+'_'+p1_str+'_'+p2_str+'.png')
                    plt.clf()
                    plt.close('all')
                    plt.show()

                    #----------------moming windog CC finished-------------------#
                    L_data1.clear()
                    L_data2.clear()

                    p1_D.clear()
                    p2_D.clear()
                    #break #break at first station
                #break #break at first pair
            #break #break at first pair seq(e.g. p1-p2,p1-p3)
        #break #break at first line in seq file
    IN1.close()
                    






