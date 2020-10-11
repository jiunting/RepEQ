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
        #read catalog
        cat = np.genfromtxt(catalog, delimiter=',', skip_header=0,usecols=(0,1,2,3,4,10,11), dtype=("|U19",float,float,float,float,"|U2","|U19"))
        df = pd.DataFrame(cat, columns=['ID','Time','Magnitude','Lat','Lon','Depth','Regional'])
        for ii in range(len(cat)):
            df = df.append({'ID': cat[ii][6][2:],'Time': cat[ii][0][11:], 'Magnitude': cat[ii][4], 'Date': cat[ii][0][:10],
                           'Lat': cat[ii][1], 'Lon': cat[ii][2], 'Depth': cat[ii][3], 'Regional': cat[ii][5]}, ignore_index=True)
        self.catalog = df

    def template_load(self):
        import os
        from repeq import download_tools
        home = self.home
        project_name = self.project_name
        if self.download:
            #start download templates
            if not(os.path.exists(home+'/'+project_name+'/waveforms_template')):
                os.makedirs(home+'/'+project_name+'/waveforms_template')
            download_tools.bulk_make_template(home,project_name,self.catalog,self.sampling_rate,self.filter,self.tcs_length)
        else:
            #read previously downloaded templetes
            if not(os.path.exists(home+'/'+project_name+'/waveforms_template')):
                print('Make sure templates are in:',home+'/'+project_name+'/waveforms_template')
                print('If not, run .template_load() with the download=True first')
                return False
            else:
                import glob
                MS = glob.glob(home+'/'+project_name+'/waveforms_template/template_*.ms')
                MS.sort()
                self.ms = MS









