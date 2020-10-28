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
    avail_idx = [] #all the available index. sometime the P or S arrival is not available,
    #loop in every station
    for ist in range(len(stlon)):
        #calculate distance in degree
        dist = obspy.geodetics.locations2degrees(lat1=eqlat,long1=eqlon,lat2=stlat[ist],long2=stlon[ist])
        #distance perturbation
        dist_add_dx = obspy.geodetics.locations2degrees(lat1=eqlat,long1=eqlon+dx,lat2=stlat[ist],long2=stlon[ist])
        dist_add_dy = obspy.geodetics.locations2degrees(lat1=eqlat+dy,long1=eqlon,lat2=stlat[ist],long2=stlon[ist])
        dist_add_dz = eqdep+dz
        #phase type P or S
        #print(stlon[ist],stlat[ist],dist)
        if phase[ist]=='P':
            PS = model.get_travel_times(source_depth_in_km=eqdep, distance_in_degree=dist, phase_list=('P','p'), receiver_depth_in_km=0)
            PS_dx = model.get_travel_times(source_depth_in_km=eqdep, distance_in_degree=dist_add_dx, phase_list=('P','p'), receiver_depth_in_km=0)
            PS_dy = model.get_travel_times(source_depth_in_km=eqdep, distance_in_degree=dist_add_dy, phase_list=('P','p'), receiver_depth_in_km=0)
            PS_dz = model.get_travel_times(source_depth_in_km=dist_add_dz, distance_in_degree=dist_add_dy, phase_list=('P','p'), receiver_depth_in_km=0)
        else:
            PS = model.get_travel_times(source_depth_in_km=eqdep, distance_in_degree=dist, phase_list=('S','s'), receiver_depth_in_km=0)
            PS_dx = model.get_travel_times(source_depth_in_km=eqdep, distance_in_degree=dist_add_dx, phase_list=('S','s'), receiver_depth_in_km=0)
            PS_dy = model.get_travel_times(source_depth_in_km=eqdep, distance_in_degree=dist_add_dy, phase_list=('S','s'), receiver_depth_in_km=0)
            PS_dz = model.get_travel_times(source_depth_in_km=dist_add_dz, distance_in_degree=dist_add_dy, phase_list=('S','s'), receiver_depth_in_km=0)
        #T derivates
        try:
            dTdx = PS_dx[0].time-PS[0].time
            dTdy = PS_dy[0].time-PS[0].time
            dTdz = PS_dz[0].time-PS[0].time
            dTdt = dt #the linear part in G
        except:
            print('fail calculating traveltime: depth=%f km, dist=%f deg'%(eqdep,dist))
            continue
        avail_idx.append(ist)
        G_row = np.array([dTdx,dTdy,dTdz,dTdt])
        try:
            G = np.vstack([G,G_row])
        except:
            G = G_row.copy()
    avail_idx = np.array(avail_idx)
    return G,avail_idx



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




def EQreloc(home,project_name,catalog,station_table,vel_model,fiter_inv,T0):
    '''
        Main relocation
        fiter_inv={
        'CCC_threshold':0.3,  #use observed shift with CCC greater than this threshold
        'min_stan':5,         #minumim observations
        'max_shift':2,        #maximum shift seconds(observation)
        'VR':0.5,             #after inversion, check the VR
        }
        T0: reference time in UTCDateTime format
    '''
    #load catalog in pd
    catalog=home+'/'+project_name+'/catalog/'+catalog.split('/')[-1]
    cat = np.genfromtxt(catalog, delimiter=',', skip_header=0,usecols=(0,1,2,3,4,10,11), dtype=("|U23",float,float,float,float,"|U2","|U23")) #accur
    df = pd.DataFrame(cat, columns=['ID','Time','Magnitude','Lat','Lon','Depth','Regional'])
    for ii in range(len(cat)):
        df = df.append({'ID': cat[ii][6][2:],'Time': cat[ii][0][11:], 'Magnitude': cat[ii][4], 'Date': cat[ii][0][:10],'Lat': cat[ii][1], 'Lon': cat[ii][2], 'Depth': cat[ii][3], 'Regional': cat[ii][5]}, ignore_index=True)

    #make velocity model (.mod to .npz)
    TauPy_name = analysis.build_TauPyModel(home,project_name,vel_model) #make .npz file
    #model = TauPyModel(model=TauPy_name)
    model_path = home+'/'+project_name+'/structure/'+vel_model.replace('mod','npz')


    detect_details = glob.glob(home+'/'+project_name+'/output/Template_match/Detections/'+'Detected_tmp_*.npy')
    detect_details.sort()
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
















