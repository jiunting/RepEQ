


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



def searchRepEQ(home,project_name,vel_model,cata_name,data_filters,startover=False,make_fig_CC=2):
    '''
        startover=True: re-run everything, or False: check the files that already exist to see if they need to be updated
        make_fig_CC: plot the figure when CC value >= this number. set a number >1 if you dont want plot
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
    #-------------------------------------------------------------------------------#
    EQfolders=glob.glob(eqpath+'*')
    EQfolders.sort()

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
        check_file_flag=0 #if the log file exist, and wanna check into the information
        if os.path.exists(repeq_dir+'/'+net_sta_key+'.log'):
            if not startover:
                #but not just skip the file, read the file and check if all the calculations are there
                UPD1=open(repeq_dir+'/'+net_sta_key+'.log','r') #open the pre-existing file and check
                Chk_lines=UPD1.read()
                UPD1.close()
                OUT1=open(repeq_dir+'/'+net_sta_key+'.log','a')
                check_file_flag=1 #check into the file and see if theres any missing
            else:
                #overwrite the existing .log file
                OUT1=open(repeq_dir+'/'+net_sta_key+'.log','w')
        else:
            #log file not exist, make a new file
            OUT1=open(repeq_dir+'/'+net_sta_key+'.log','w')

        sav_Date=[]
        n=0
        pre_sampr=[]      #make sure sampling rate for a same station always keeps the same. If not, don't use it
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
                continue
            sav_Date.append(Date)
            #read the data path by obspy
            data=obspy.read(n_date+'/waveforms/'+net_sta_key+'*.mseed')
            tmp_sampr=data[0].stats.sampling_rate #check sampling rate
            if pre_sampr==[]:
                pre_sampr=tmp_sampr #set the sampr for first station
            if tmp_sampr!=pre_sampr:
                print('Sampling rate inconsistent')
                continue
            #-------------get info------------
            stlon,stlat=get_staloc(net_sta_key,n_date)
            #get travel time
            tP,tS,GCARC=get_traveltime(stlon,stlat,eqlon,eqlat,eqdep,model_name='iasp91')
            sav_Parrivl.append(tP)   #P wave travel time
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
                print('detect Nan value in the data, continue to the next')
                continue
            #cut data by P-arrival
            idxt=np.where( (tP-p_wind[0]<=t) & (t<=tP+p_wind[1]) )[0] #find the index of P-waves
            Pwave_t=t[idxt]
            Pwave_y=y[idxt]/np.max(np.abs(y[idxt]))
            sav_data.append(Pwave_y) #sav_data: cutted P wave
            sav_t.append(Pwave_t)    #sav_t: save Pwave time. 0 at the KZtime.  plt.plot(sav_t[0],sav_data[0]) should make sense. t=0 at kztime
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
                #check whether needs to calculate
#                if (mkdirflag=='N') & (updateflag=='Y') & (logf_exist=='Y'):
                if check_file_flag:
                    #log file exist and you must have already read it into a large string, right?
                    if '%s-%s'%(sav_legend[idate],sav_legend[jdate]) in Chk_lines:
                        continue
                eqdist_degree=obspy.geodetics.locations2degrees(lat1=sav_evlat[idate],long1=sav_evlon[idate],lat2=sav_evlat[jdate],long2=sav_evlon[jdate])
                if (eqdist_degree>0.2):
                    continue
                if jdate<=idate:
                    continue
                CCC,lag=cal_CCF(sav_data[idate],sav_data[jdate])
                if np.isnan(CCC):
                    continue
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





