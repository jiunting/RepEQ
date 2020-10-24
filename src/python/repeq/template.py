#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 00:09:04 2020

@author: TimLin
"""
#tools for template matching method

import numpy as np
import pandas as pd


class Template():
    def __init__(self,home,project_name,catalog,download,sampling_rate,filter,tcs_length=[1,9],filt_CC=0.2,filt_nSTA=5,plot_check=False):
        '''
            home: path of the home working directory <str>
            project_name: project name <str>
            catalog: catalog name, can be absolute path or file name in home/project_name/catalog/ <str>
            download: download waveform through libcomcat or load from previously downloaded waveforms (in home/project_name/waveforms_template) <boolean>
            sampling_rate: sampling rate of data <int>
            filter: bandpass filter frequency range [min_freq,max_freq] <list with length=2>
            tcs_length: download data length center at the arrival time [t1,t2] i.e. from T-t1 to T+t2 <list with length=2>
            filt_CC: individual CC value larger than this threshold <float>
            filt_nSTA: N station larger than this threshold <int>
            plot_check: plot detailed figure for checking <boolean>
        '''
        self.home = home
        self.project_name = project_name
        if not '/' in catalog:
            catalog = home+'/'+project_name+'/catalog/'+catalog
        self.catalog = catalog #full path of catalog
        self.download = download
        self.sampling_rate = sampling_rate
        self.filter = filter
        self.tcs_length = tcs_length
        self.filt_CC = filt_CC
        self.filt_nSTA = filt_nSTA
        self.plot_check = plot_check
        self.ms = None
        #read catalog
        #cat = np.genfromtxt(catalog, delimiter=',', skip_header=0,usecols=(0,1,2,3,4,10,11), dtype=("|U19",float,float,float,float,"|U2","|U19")) #accuracy to sec
        cat = np.genfromtxt(catalog, delimiter=',', skip_header=0,usecols=(0,1,2,3,4,10,11), dtype=("|U23",float,float,float,float,"|U2","|U23")) #accuracy to ms
        df = pd.DataFrame(cat, columns=['ID','Time','Magnitude','Lat','Lon','Depth','Regional'])
        for ii in range(len(cat)):
            df = df.append({'ID': cat[ii][6][2:],'Time': cat[ii][0][11:], 'Magnitude': cat[ii][4], 'Date': cat[ii][0][:10],
                           'Lat': cat[ii][1], 'Lon': cat[ii][2], 'Depth': cat[ii][3], 'Regional': cat[ii][5]}, ignore_index=True)
        self.catalog = df

    def template_load(self):
        import os,glob
        from repeq import download_tools
        home = self.home
        project_name = self.project_name
        if self.download:
            #start download templates
            if not(os.path.exists(home+'/'+project_name+'/waveforms_template')):
                os.makedirs(home+'/'+project_name+'/waveforms_template')
            download_tools.bulk_make_template(home,project_name,self.catalog,self.sampling_rate,self.filter,self.tcs_length)
            MS = glob.glob(home+'/'+project_name+'/waveforms_template/template_*.ms')
            MS.sort()
            self.ms = MS
        else:
            #read previously downloaded templetes
            if not(os.path.exists(home+'/'+project_name+'/waveforms_template')):
                print('Make sure templates are in:',home+'/'+project_name+'/waveforms_template')
                print('If not, run .template_load() with the attribute download=True first')
                return False
            else:
                MS = glob.glob(home+'/'+project_name+'/waveforms_template/template_*.ms')
                MS.sort()
                self.ms = MS  #this is template paths

    def xcorr_cont(self,save_CCF=False,fmt=1):
        '''
            save_CCF: save individual CCFs in project_name/output/Template_match/CCF_records/
            fmt: output format if
                fmt=1:
                    #OriginTime meanCC nSTA templateIDX
                fmt=2:
                    #OriginTime meanCC nSTA templateIDX mean_maxCCC
        '''
        from obspy import UTCDateTime,read,Stream,Trace
        import glob
        from scipy import signal
        from obspy.signal.cross_correlation import correlate_template
        import matplotlib
        matplotlib.use('pdf') #instead using interactive backend
        import matplotlib.pyplot as plt
        
        home = self.home
        project_name = self.project_name
        
        def cal_CCF(data1,data2):
            #calculate normalize CCF, find max CCC, and lag idx
            tmpccf=signal.correlate(data1,data2,'full')
            auto1=signal.correlate(data1,data1,'full')
            auto2=signal.correlate(data2,data2,'full')
            tmpccf=tmpccf/np.sqrt(np.max(auto1)*np.max(auto2))
            maxCCC=np.max(tmpccf)
            lag=tmpccf.argmax()
            return(maxCCC,lag)
        
        if self.ms == None:
            print('Run .template_load() first')
            return False
        else:
            #loop the templates (paths)
            for i_tmp in self.ms:
                print('----------------------------------------------')
                print('In template: %s'%(i_tmp))
                tmp_idx = int(i_tmp.split('/')[-1].split('_')[-1].split('.')[0])
                OUT1 = open(home+'/'+project_name+'/output/Template_match/Detections/Detected_tmp_%05d.txt'%(tmp_idx),'w') #output earthquake origin time
                if fmt==1:
                    OUT1.write('#OriginTime meanCC stdCC nSTA templateIDX\n')
                elif fmt==2:
                    OUT1.write('#OriginTime meanCC stdCC nSTA templateIDX mean_maxCCC std_maxCCC\n') # mean(max CCC for each stations), so that shift in each sta is negletable
                origintime = UTCDateTime(self.catalog.iloc[tmp_idx].Date+'T'+self.catalog.iloc[tmp_idx].Time)
                st = read(i_tmp) #read template in
                
                #read all directories of daily data
                dayst_paths = glob.glob(home+'/'+project_name+'/waveforms/'+'*000000')
                dayst_paths.sort()
                
                sav_mean_sh_CCF=[] #save all the daily CCF for later plotting
                sav_daily_nSTA=[] #number of stations for each daily CCF
                sav_alldays_eq_sta = {} #detailed info for CC,CCC,shifts for every station for all searched days by the same template
                #loop the daily data
                for dayst_path in dayst_paths:
                    sav_NET=[]; sav_STA=[]; sav_CHN=[]; sav_CCF=[]; sav_travel_npts=[]; sav_continuousdata=[]; sav_template=[] #initial for saving
                    YMD = dayst_path.split('/')[-1][:8]
                    print(' --Reading daily data: %s'%(dayst_path))
                    i_dayst = read(dayst_path+'/waveforms/merged.ms') #load daily data
                    #print(i_dayst.__str__(extended=True))
                    for i in range(len(st)):
                        #-----loop individual pick/station/comp of template-----
                        NET = st[i].stats.network
                        STA = st[i].stats.station
                        CHN = st[i].stats.channel
                                        
                        #in daily data... search for same station,channel,comp,sampling rate....that matches the i_th pick in particular template
                        tmp_dayst = i_dayst.select(network=st[i].stats.network,station=STA,sampling_rate=st[i].stats.sampling_rate,
                                                   channel=st[i].stats.channel,location=st[i].stats.location)
                        if len(tmp_dayst)!=1:
                            if len(tmp_dayst)==0:
                                #print('No data found:%s, skip this station'%(STA+'.'+CHN))
                                pass
                            else:
                                #print('Multiple data found:%d, probably breaking tcs, skip this station'%(len(tmp_dayst)))
                                #print(tmp_dayst) #tmp_dayst should be only one
                                pass
                            continue
                        else:
                            #find the station travel time
                            #if len(phases[phases.Channel.str.startswith(regional.upper()+'.'+STA+'.'+CHN)])==0:
                            #    continue #cannot find station shift
                            travel_time = st[i].stats.starttime + self.tcs_length[0] - origintime
                            travel_npts = int(np.round(travel_time*self.sampling_rate)) #travel time in npts
                                        
                            #get data value for template and continuous(daily) data
                            template = np.nan_to_num(st[i].data)
                            continuousdata = np.nan_to_num(tmp_dayst[0].data)
                                        
                            #run xcorr
                            CCF = correlate_template(continuousdata,template)
                            CCF = np.nan_to_num(CCF)
                            
                            #save for later checking
                            sav_NET.append(NET)
                            sav_STA.append(STA)
                            sav_CHN.append(CHN)
                            sav_travel_npts.append(travel_npts)
                            sav_CCF.append(CCF)
                            sav_continuousdata.append(continuousdata)
                            sav_template.append(template)
                                        
                    if len(sav_CCF)<self.filt_nSTA:
                        print('   Number of CCF: %d, not enough for threshold'%(len(sav_CCF)))
                        continue #not enough data available, continue to next daily data

                    #----------dealing with shifting of each CCF----------
                    travel_npts = np.array(travel_npts)
                    sh_sav_CCF = np.array(sav_CCF) #copy the original CCF
                    #shifted CCF based on the template arrival
                    for ii in range(len(sh_sav_CCF)):
                        sh_sav_CCF[ii] = np.roll(sav_CCF[ii],-int(sav_travel_npts[ii]))
                    
                    print('   Number of CCF: %d, continue searching earthquakes'%(len(sav_CCF)))
                    mean_sh_CCF = np.mean(sh_sav_CCF,axis=0) #stack/mean all the CCFs.
                    std_sh_CCF = np.std(sh_sav_CCF,axis=0) #also calculate std
                    
                    #save the individual CCF in Stream (for only debug purpose)
                    if save_CCF:
                        ST = Stream()
                        for ii,iCCF in enumerate(sh_sav_CCF):
                            tmpCCF = Trace(iCCF)
                            tmpCCF.stats.sampling_rate = i_dayst[0].stats.sampling_rate
                            tmpCCF.stats.starttime = i_dayst[0].stats.starttime
                            tmpCCF.stats.network = sav_NET[ii]
                            tmpCCF.stats.station = sav_STA[ii]
                            tmpCCF.stats.channel = sav_CHN[ii]
                            ST += tmpCCF
                        ST.write(home+'/'+project_name+'/output/Template_match/CCF_records/'+'shftCCF_template_%05d_daily_%s.ms'%(tmp_idx,YMD),format="MSEED")
                    
                    #----------Find earthquakes by the mean CCF----------
                    time = i_dayst[0].times()
                    eq_idx = np.where(mean_sh_CCF>=self.filt_CC)[0]
                    
                    sav_eq_sta = {} #save the detailed result(lag info, CCC value) for use later
                    for neqid in eq_idx:
                        #new_dayst[0].stats.starttime+time[np.argmax(mean_sh_CCF)] #find itself
                        detected_OT = i_dayst[0].stats.starttime+time[neqid]+self.tcs_length[0] #Origin time of which detection
                        detected_OT_str = detected_OT.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-4] #accuracy to 0.01 sec
                        print('    New event found:',i_dayst[0].stats.starttime+time[neqid]+self.tcs_length[0]) #find earthquakes,this is the arrival for template.st
                        if fmt == 1:
                            OUT1.write('%s %.3f %.3f %d %s\n'%(detected_OT_str,mean_sh_CCF[neqid],std_sh_CCF[neqid],len(sav_STA),'template_%05d'%(tmp_idx)))
                        elif fmt == 2:
                            #calculate CCC for individual stations
                            sav_maxCCC = []; #sav_sh_sec=[]
                            for n in range(len(sav_template)):
                                #loop in every station
                                cut_daily = sav_continuousdata[n][neqid+sav_travel_npts[n]:neqid+sav_travel_npts[n]+len(sav_template[n])]
                                maxCCC,lag = cal_CCF(sav_template[n],cut_daily)
                                if np.isnsn(maxCCC):
                                    maxCCC = 0 #this is probalby due to cross-correlate on a zero array
                                midd = (len(cut_daily))-1  #length of b?? at this idx, refdata align with target data
                                sh_sec = (lag-midd)*(1.0/self.sampling_rate) #convert to second (dt correction of P)
                                sav_maxCCC.append(maxCCC)
                                if detected_OT_str in sav_eq_sta:
                                    sav_eq_sta[detected_OT_str]['net_sta_comp'].append(sav_NET[n]+'.'+sav_STA[n]+'.'+sav_CHN[n])
                                    sav_eq_sta[detected_OT_str]['CCC'].append(maxCCC)
                                    sav_eq_sta[detected_OT_str]['CC'].append(sh_sav_CCF[n][neqid])
                                    sav_eq_sta[detected_OT_str]['shift'].append(sh_sec)
                                
                                else:
                                    #initial dictionary
                                    sav_eq_sta[detected_OT_str] = {}
                                    sav_eq_sta[detected_OT_str]['net_sta_comp'] = [sav_NET[n]+'.'+sav_STA[n]+'.'+sav_CHN[n]]
                                    sav_eq_sta[detected_OT_str]['CCC'] = [maxCCC]
                                    sav_eq_sta[detected_OT_str]['CC'] = [sh_sav_CCF[n][neqid]] #sh_sav_CCF[n][neqid]
                                    sav_eq_sta[detected_OT_str]['shift'] = [sh_sec]
    
                                        
                                #sav_sh_sec.append(sh_sec)
                            OUT1.write('%s %.3f %.3f %d %s %.3f %.3f\n'%(detected_OT_str,mean_sh_CCF[neqid],std_sh_CCF[neqid],len(sav_STA),'template_%05d'%(tmp_idx),np.mean(sav_maxCCC),np.std(sav_maxCCC)))

                    #-----Only for checking: plot the one with largest CC value and check (find itself if the template and daily are the same day)-----
                    if self.plot_check:
                        tmp_T = st[0].times()
                        for i_eqidx,neqid in enumerate(eq_idx):
                            #loop in detection
                            detected_OT = i_dayst[0].stats.starttime+time[neqid]+self.tcs_length[0] #Origin time of which detection
                            detected_OT_str = detected_OT.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-4] #accuracy to 0.01 sec
                            plt.figure(1)
                            for n in range(len(sav_template)):
                                #loop in every station
                                #cut_daily = sav_continuousdata[n][np.argmax(mean_sh_CCF)+sav_travel_npts[n]:np.argmax(mean_sh_CCF)+sav_travel_npts[n]+len(sav_template[n])] #old version only plot maximum
                                cut_daily = sav_continuousdata[n][neqid+sav_travel_npts[n]:neqid+sav_travel_npts[n]+len(sav_template[n])]
                                cut_daily = cut_daily/np.max(np.abs(cut_daily))
                                plt.plot(tmp_T,cut_daily+n,'k',linewidth=2) #time series cutted from daily time series
                                plt.plot(tmp_T,sav_template[n]/np.max(np.abs(sav_template[n]))+n,'r',linewidth=1.2) #template data
                                plt.text(tmp_T[-1],n,sav_STA[n]+'.'+sav_CHN[n])
                                #---add individual CC value and max_CCC value---
                                if fmt==1:
                                    #maxCCC,lag = cal_CCF(sav_template[n],cut_daily)
                                    #midd = (len(cut_daily))-1  #length of b?? at this idx, refdata align with target data
                                    #sh_sec = (lag-midd)*(1.0/self.sampling_rate) #convert to second (dt correction of P)
                                    plt.text(tmp_T[-1]*0.05,n,'CC=%.2f'%(sh_sav_CCF[n][neqid]))
                                elif fmt==2:
                                    maxCCC = sav_eq_sta[detected_OT_str]['CCC'][n]
                                    sh_sec = sav_eq_sta[detected_OT_str]['shift'][n]
                                    plt.text(tmp_T[-1]*0.05,n,'CC=%.2f,max_CCC=%.2f,dt=%.3f'%(sh_sav_CCF[n][neqid],maxCCC,sh_sec))
                                #Future improvement: if fmt==2, the value have been calculated, just get the value
                                #if fmt == 1:
                                #elif fmt ==2:
                                
                            #plt.title('Time:%s  CC=%5.2f'%((i_dayst[0].stats.starttime+time[neqid]+self.tcs_length[0]).strftime('%H:%M:%S'),np.max(mean_sh_CCF)))
                            plt.title('Time:%s  CC=%5.2f'%((i_dayst[0].stats.starttime+time[neqid]+self.tcs_length[0]).strftime('%H:%M:%S.%f'),mean_sh_CCF[neqid]))
                            plt.savefig(home+'/'+project_name+'/output/Template_match/Figs/'+'template_%05d_daily_%s_%03d.png'%(tmp_idx,YMD,i_eqidx))
                            plt.close()
                            if i_eqidx>99:
                                break #don't plot if more than 99 plots in the same day
                    
                    sav_mean_sh_CCF.append(mean_sh_CCF)
                    sav_daily_nSTA.append(len(sav_CCF))
                    
                    sav_alldays_eq_sta.update(sav_eq_sta) #not support for fmt=1
                        
                ##------output detailed data(lag information for each station) in .npy ---------
                #only if fmt=2, fmt=1 didnt calculate the CCC
                if fmt==2:
                    np.save(home+'/'+project_name+'/output/Template_match/Detections/'+'Detected_tmp_%05d.npy'%(tmp_idx),sav_alldays_eq_sta)
                    
                
                
                #----plot the mean_shifted_CCF for all days----
                plt.figure(1)
                for n in range(len(sav_mean_sh_CCF)):
                    plt.plot(sav_mean_sh_CCF[n]+n,linewidth=1)
                    if n==0:
                        plt.text(len(sav_mean_sh_CCF[n]),n,'N=%d'%(sav_daily_nSTA[n])) #number of stations
                    else:
                        plt.text(len(sav_mean_sh_CCF[n]),n,'%d'%(sav_daily_nSTA[n])) #number of stations
                plt.title('Mean CCF (template_%05d)'%(tmp_idx),fontsize=16)
                plt.ylabel('Days after %s'%(dayst_paths[0].split('/')[-1][:8]),fontsize=16)
                plt.savefig(home+'/'+project_name+'/output/Template_match/Figs/'+'MeanCCF_%05d.png'%(tmp_idx))
                plt.close()
                OUT1.close()






