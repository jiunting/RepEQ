#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 14:46:02 2020

@author: TimLin
"""
import numpy as np
import obspy
import glob
from obspy import UTCDateTime
import pandas as pd
import os

#--------repeating earthquake relocation tools------------

def GMD_solve(G,D):
    #calculate inversion here
    GT=np.transpose(G)
    try:
        M=np.dot(np.dot(np.linalg.inv(np.dot(GT,G)),np.transpose(G)),D)
    except:
        print('encounter singular matrix, use pinv instead')
        M=np.dot(np.dot(np.linalg.pinv(np.dot(GT,G)),np.transpose(G)),D)
    return M



def travel_time(model,stlon,stlat,phase,eqlon,eqlat,eqdep):
    '''
        Input:
            model: TauPyModel object (i.e. model = TauPyModel(model=model_path) )
            stlon/stlat: stlon/stlat in array or list
            phase: phase type at the stlon/stlat ('P' or 'S'). same length as stlon/stlat
            eqlon/eqlat/eqdep: EQ location. Note that eqdep should be positive [km]
        Output:
            predicted travel time at stlon/stlat [sec]
    '''
    sav_time = []
    for ist in range(len(stlon)):
        #calculate distance in degree
        dist = obspy.geodetics.locations2degrees(lat1=eqlat,long1=eqlon,lat2=stlat[ist],long2=stlon[ist])
        if phase[ist]=='P':
            PS = model.get_travel_times(source_depth_in_km=eqdep, distance_in_degree=dist, phase_list=('P','p'), receiver_depth_in_km=0)
        else:
            PS = model.get_travel_times(source_depth_in_km=eqdep, distance_in_degree=dist, phase_list=('S','s'), receiver_depth_in_km=0)
        try:
            sav_time.append(PS[0].time)
        except:
            #ray unavailable
            sav_time.append(np.nan)
    sav_time = np.array(sav_time)
    return sav_time


def cal_derivative(model_path,eqlon,eqlat,eqdep,stlon,stlat,phase,dx=0.1,dy=0.1,dz=0.1,dt=0.04):
    from obspy.taup import TauPyModel
    from obspy.geodetics import kilometers2degrees
    #model_path: path of velocity npz file
    #eqlon,eqlat,eqdep:longitude/lat/depth(km) of earthquake
    #stlon,stlat: array of stations lon/lat
    #phase: calculate in P or S phase (array same length as stlon/stlat)
    #dx,dy,dz: perturtabation used for calculating gradient. [input unit in km]
    model = TauPyModel(model=model_path)
    #convert km to degree
    dx = kilometers2degrees(dx)
    dy = kilometers2degrees(dy)
    G = []
    #avail_idx = [] #all the available index. sometime the P or S arrival is not available,
    T0 = travel_time(model,stlon,stlat,phase,eqlon,eqlat,eqdep)
    T_dx = travel_time(model,stlon,stlat,phase,eqlon+dx,eqlat,eqdep)
    T_dy = travel_time(model,stlon,stlat,phase,eqlon,eqlat+dy,eqdep)
    T_dz = travel_time(model,stlon,stlat,phase,eqlon,eqlat,eqdep+dz)
    #time changes
    dT_dx = T_dx-T0
    dT_dy = T_dy-T0
    dT_dz = T_dz-T0
    #combine all into a G array
    G = np.hstack([dT_dx.reshape(-1,1),dT_dy.reshape(-1,1),dT_dz.reshape(-1,1)])
    return G
    



def get_lonlat(sta_table,sta):
    #sta_table is all station lon, lat, name in pandas format
    #sta: list of net.staname or net.stname.comp (the comp will be ignored)
    #find sta lon,lat
    sta = [i.split('.')[0]+'.'+i.split('.')[1] for i in sta]
    sav_lon = []
    sav_lat = []
    for ista in sta:
        idx = np.where(sta_table.stname==ista)[0][0]
        sav_lon.append(sta_table.stlon[idx])
        sav_lat.append(sta_table.stlat[idx])
    return np.array(sav_lon),np.array(sav_lat)

def get_all_sta(detc):
    #input detc, output all the net.staname in array
    all_sta = []
    for k in detc.keys():
        all_sta += detc[k]['net_sta_comp']
    return list(set(all_sta))

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


def cal_VR(y,yhat):
    return 1-(np.sum((y-yhat)**2)/np.sum(y**2))

def EQreloc(home,project_name,catalog,vel_model,filter_detc,fiter_inv,T0):
    '''
        Main relocation
        fiter_inv={
        'CCC_threshold':0.3,  #use observed shift with CCC greater than this threshold
        'min_stan':5,         #minumim observations
        'max_shift':2,        #maximum shift seconds(observation)
        'VR':0.5,             #after inversion, check the VR
        }
        T0: reference time in UTCDateTime format. Used in output data
    '''
    from repeq import analysis
    from obspy.geodetics import kilometers2degrees
    
    #load catalog in pd
    catalog=home+'/'+project_name+'/catalog/'+catalog.split('/')[-1]
    cat = np.genfromtxt(catalog, delimiter=',', skip_header=0,usecols=(0,1,2,3,4,10,11), dtype=("|U23",float,float,float,float,"|U2","|U23")) #accur
    df = pd.DataFrame(cat, columns=['ID','Time','Magnitude','Lat','Lon','Depth','Regional'])
    for ii in range(len(cat)):
        df = df.append({'ID': cat[ii][6][2:],'Time': cat[ii][0][11:], 'Magnitude': cat[ii][4], 'Date': cat[ii][0][:10],'Lat': cat[ii][1], 'Lon': cat[ii][2], 'Depth': cat[ii][3], 'Regional': cat[ii][5]}, ignore_index=True)

    #make velocity model (.mod to .npz) if not exist
    if not(os.path.exists(home+'/'+project_name+'/structure/'+vel_model.replace('mod','npz'))):
        model_path = analysis.build_TauPyModel(home,project_name,vel_model) #make .npz file
    else:
        #npz file already exist
        model_path = home+'/'+project_name+'/structure/'+vel_model.replace('mod','npz')

    ##load station table in home/project_nam/stations/stations.txt
    sta_table = pd.read_table(home+'/'+project_name+'/stations/'+'stations.txt',header=None,names=['stlon','stlat','stelev','stname'],sep=' ')

    #=====small step for calculating derivative=====
    #change these only if necessarily
    #dx = 0.1 #[unit:km]
    #dy = 0.1
    #dz = 0.1
    #dt = 0.04
    dx = 1 #[unit:km]
    dy = 1
    dz = 1
    dt = 0.4 #don't use here
    #========================================

    #load all detailed detection .npy files
    detect_details = glob.glob(home+'/'+project_name+'/output/Template_match/Detections/'+'Detected_tmp_*.npy')
    detect_details.sort()
    OUT1 = open(home+'/'+project_name+'/output/Template_match/Detections/'+'EQreloc_info.txt','w')
    OUT1.write('#dx dy dz dt VR nstan refT ref_tempT eqlon eqlat eqdep templateID\n')
    OUT1.close()
    for ndet,detect_detail in enumerate(detect_details):
        #detect_detail=detect_details[3]
        #detect_detail=detect_details[18]
        print('=========Now in:%s==========='%(detect_detail))
        #detailed detection file(CCC value and time shift)
        detc = np.load(detect_detail,allow_pickle=True)
        detc = detc.item()
        detc = clean_detc(detc,filter_detc)
        if detc=={}:
            continue
        neqid = int(detect_detail.split('_')[-1].split('.')[0])
        #get eqinfo
        eqlon = df.iloc[neqid].Lon
        eqlat = df.iloc[neqid].Lat
        eqdep = df.iloc[neqid].Depth
        OT = UTCDateTime(df.iloc[neqid].Date+'T'+df.iloc[neqid].Time)
        if eqdep<0:
            eqdep=0.01
        #all station in the dictionary
        #all_sta = get_all_sta(detc)
        sav_M_km = [] #save inversion parameters
        sav_reloc = [] #final result of relocation
        sav_VR = []
        sav_nstan = []
        sav_date = [] #tims relative to the T0
        for k in detc.keys():
            stlon,stlat = get_lonlat(sta_table,detc[k]['net_sta_comp'])
            #pick information (P or S?)
            phase = np.array(detc[k]['phase'])
            #time shift. positive means arrives earlier, means timeline moves negative
            #i.e. arrival of template=1.3sec, observation=1.1sec, then shift=0.2, while observation-template=-0.2
            shifts = np.array(detc[k]['shift'])
            CCC = np.array(detc[k]['CCC'])
            #QC1. with min_stations and CCC threshold
            use_idx = np.where( (CCC>=fiter_inv['CCC_threshold']) & (np.abs(shifts)<=fiter_inv['max_shift']))[0]
            if len(use_idx) < fiter_inv['min_stan']:
                continue
            #only keep whe measurement with high CCC value
            #print('k=',k)
            stlon = stlon[use_idx]
            stlat = stlat[use_idx]
            phase = phase[use_idx]
            shifts = shifts[use_idx]
            CCC = CCC[use_idx]
            #build G and do inversion
            #note that this only invert once! while most of the methods use iterative inversion
            G,avail_idx = cal_derivative(model_path,eqlon,eqlat,eqdep,stlon,stlat,phase,dx=dx,dy=dy,dz=dz,dt=dt)
            M = GMD_solve(G,shifts[avail_idx]*-1) #note that positive time shift means arrives earlier
            #QC2. filter with VR
            VR = cal_VR(shifts[avail_idx]*-1,np.dot(G,M))
            if not np.isnan(VR):
                if VR<fiter_inv['VR']:
                    continue
            else:
                VR = 1 # ()/sum(y**2) is nan means all y are zeros, and M=0
            #convert to original scale
            #M_km = M*np.array([dx,dy,dz,dt]) #M to the original scale(km)
            M_km = M*np.array([dx,dy,dz]) #M to the original scale(km)
            
            #convert shifted_km to shifted_deg
            sh_deg = kilometers2degrees(M_km[:2])
            #QC3. filter with unreasonable large M (almost impossible)
            if np.abs(sh_deg[0])>2 or np.abs(sh_deg[1])>2 or np.abs(M_km[2])>200:
                #import sys
                #sys.exit() # exit and see what causes this!
                continue
            new_lon = eqlon+sh_deg[0]
            new_lat = eqlat+sh_deg[1]
            sav_M_km.append(M_km)
            sav_reloc.append(np.array([new_lon,new_lat,eqdep+M_km[2]]))
            sav_date.append( (UTCDateTime(k)-T0)/86400.0)
            sav_VR.append(VR)
            sav_nstan.append(len(shifts))
        #finally make all list to array
        sav_M_km = np.array(sav_M_km)
        sav_reloc = np.array(sav_reloc)
        sav_date = np.array(sav_date)
        if len(sav_reloc)>0:
            #save dX,dY,dZ,dT to a txt file
            OUT1 = open(home+'/'+project_name+'/output/Template_match/Detections/'+'EQreloc_info.txt','a')
            if sav_M_km.ndim==1:
                OUT1.write('%f %f %f  %f %d %f %f  %f %f %f %05d\n'%(sav_M_km[0],sav_M_km[1],sav_M_km[2],sav_VR[0],sav_nstan[0],sav_date[0],(OT-T0)/86400.0,eqlon,eqlat,eqdep,neqid))
            else:
                for ii in range(len(sav_M_km)):
                    OUT1.write('%f %f %f  %f %d %f %f  %f %f %f %05d\n'%(sav_M_km[ii][0],sav_M_km[ii][1],sav_M_km[ii][2],sav_VR[ii],sav_nstan[ii],sav_date[ii],(OT-T0)/86400.0,eqlon,eqlat,eqdep,neqid))
            OUT1.close()















