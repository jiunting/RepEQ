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
    def __init__(self,home,project_name,catalog,download,sampling_rate,filter,tcs_length=[1,9],filt_CC=0.2,filt_nSTA=5):
        '''
            home: path of the home working directory <str>
            project_name: project name <str>
            catalog: catalog name, can be absolute path or file name in home/project_name/catalog/ <str>
            download: download waveform through libcomcat or load from previously downloaded waveforms (in home/project_name/waveforms_template) <boolean>
            sampling_rate: sampling rate of data <int>
            filter: bandpass filter frequency range [min_freq,max_freq] <list with length=2>
            tcs_length: download data length center at the arrival time [t1,t2] i.e. from T-t1 to T+t2 <list with length=2>
            filt_CC: 
            filt_nSTA:
            
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
        self.ms = None
        #read catalog
        cat = np.genfromtxt(catalog, delimiter=',', skip_header=0,usecols=(0,1,2,3,4,10,11), dtype=("|U19",float,float,float,float,"|U2","|U19"))
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

    def xcorr_cont(self):
        from obspy import UTCDateTime,read
        import glob
        from obspy.signal.cross_correlation import correlate_template
        import matplotlib
        matplotlib.use('pdf') #instead using interactive backend
        import matplotlib.pyplot as plt
        
        home = self.home
        project_name = self.project_name
        
        if self.ms == None:
            print('Run .template_load() first')
            return False
        else:
            #loop the templates (paths)
            for i_tmp in self.ms:
                print('In template: %s'%(i_tmp))
                tmp_idx = int(i_tmp.split('/')[-1].split('_')[-1].split('.')[0])
                OUT1 = open(home+'/'+project_name+'/output/Template_match/Detections/Detected_tmp_%05d.txt'%(tmp_idx),'w') #output earthquake origin time
                OUT1.write('#OriginTime meanCC nSTA templateIDX\n')
                origintime = UTCDateTime(self.catalog.iloc[tmp_idx].Date+'T'+self.catalog.iloc[tmp_idx].Time)
                st = read(i_tmp) #read template in
                
                #read all directories of daily data
                dayst_paths = glob.glob(home+'/'+project_name+'/waveforms/'+'*000000')
                dayst_paths.sort()
                
                sav_mean_sh_CCF=[] #save all the daily CCF for later plotting
                #loop the daily data
                for dayst_path in dayst_paths:
                    sav_STA=[]; sav_CHN=[]; sav_CCF=[]; sav_travel_npts=[]; sav_continuousdata=[]; sav_template=[] #initial for saving
                    YMD = dayst_path.split('/')[-1][:8]
                    print('--Reading daily data: %s'%(dayst_path))
                    i_dayst = read(dayst_path+'/waveforms/merged.ms') #load daily data
                    #print(i_dayst.__str__(extended=True))
                    for i in range(len(st)):
                        #-----loop individual pick/station/comp of template-----
                        STA=st[i].stats.station
                        CHN=st[i].stats.channel
                                        
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
                            sav_STA.append(STA)
                            sav_CHN.append(CHN)
                            sav_travel_npts.append(travel_npts)
                            sav_CCF.append(CCF)
                            sav_continuousdata.append(continuousdata)
                            sav_template.append(template)
                                        
                    if len(sav_CCF)<self.filt_nSTA:
                        print('  Number of CCF: %d, not enough for threshold'%(len(sav_CCF)))
                        continue #not enough data available, continue to next daily data

                    #----------dealing with shifting of each CCF----------
                    travel_npts = np.array(travel_npts)
                    sh_sav_CCF = np.array(sav_CCF) #copy the original CCF
                    #shifted CCF based on the template arrival
                    for ii in range(len(sh_sav_CCF)):
                        sh_sav_CCF[ii] = np.roll(sav_CCF[ii],-int(sav_travel_npts[ii]))
                    
                    print('  Number of CCF: %d, continue searching earthquakes'%(len(sav_CCF)))
                    mean_sh_CCF = np.mean(sh_sav_CCF,axis=0) #stack/mean all the CCFs.
                                        
                    #----------Find earthquakes by the mean CCF----------
                    time = i_dayst[0].times()
                    eq_idx = np.where(mean_sh_CCF>=self.filt_CC)[0]
                                        
                    for neqid in eq_idx:
                        #new_dayst[0].stats.starttime+time[np.argmax(mean_sh_CCF)] #find itself
                        print('    New event found:',i_dayst[0].stats.starttime+time[neqid]) #find earthquakes,
                        OUT1.write('%s %.3f %d %s\n'%((i_dayst[0].stats.starttime+time[neqid]).strftime('%Y-%m-%dT%H:%M:%S.%f'),
                                                      mean_sh_CCF[neqid],len(sav_STA),'template_%05d'%(tmp_idx)))

                    sav_mean_sh_CCF.append(mean_sh_CCF)
                                        
                    #-----Only for checking: plot the one with largest CC value and check (find itself if the template and daily are the same day)-----
                    plt.figure(1)
                    for n in range(len(sav_template)):
                        #print('continuous data=',sav_continuousdata[n])
                        #print('slice from:',np.argmax(mean_sh_CCF))
                        #print('plus:',sav_travel_npts[n])
                        #print('len temp:',len(sav_template[n]))
                        cut_daily = sav_continuousdata[n][np.argmax(mean_sh_CCF)+sav_travel_npts[n]:np.argmax(mean_sh_CCF)+sav_travel_npts[n]+len(sav_template[n])]
                        cut_daily = cut_daily/np.max(np.abs(cut_daily))
                        plt.plot(cut_daily+n,'k',linewidth=2) #time series cutted from daily time series
                        plt.plot(sav_template[n]/np.max(np.abs(sav_template[n]))+n,'r',linewidth=1.2) #template data
                        plt.text(400,n,sav_STA[n])
                        plt.title('CC=%5.2f'%(np.max(mean_sh_CCF)))
                    plt.title('CC=%5.2f'%(np.max(mean_sh_CCF)))
                    plt.savefig(home+'/'+project_name+'/output/Template_match/Figs/'+'tmp_%05d_daily_%s.png'%(tmp_idx,YMD))
                    plt.close()
                                        
                #----plot the mean_shifted_CCF for all days----
                plt.figure(1)
                for n in range(len(sav_mean_sh_CCF)):
                    plt.plot(sav_mean_sh_CCF[n]+n)
                plt.title('Mean CCF (template_%05d)'%(tmp_idx),fontsize=16)
                plt.ylabel('Days after %s'%(dayst_paths[0].split('/')[-1][:8]),fontsize=16)
                plt.savefig('./Figs/MeanCCF_%05d.png'%(tmp_idx))
                plt.close()
                OUT1.close()

                                        
                                        
                                        

                                        
                                        
                                        
                                        
                                        
                                        
                                        
                                        
                                        
                                        
                                        


