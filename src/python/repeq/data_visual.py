#functions that show the data

import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime

#def plot_map(home,project_name,cata_file,create_plot=True,show=False):
#    if create_plot:
#        fig = plt.figure()
#    #load catalog in pandas
#    pnsn_cat = np.genfromtxt('/Users/timlin/Documents/Project/Hawaii_cont/foreshock2018/catalog/area2_1.0.cat', delimiter=',', skip_header=0,usecols=(0,1,2,3,4,10,11), dtype=("|U22",float,float,float,float,"|U2","|U22"))
##cat = np.genfromtxt('/Users/timlin/Documents/Project/Hawaii_cont/foreshock2018/catalog/area1_1.0.cat', delimiter=',', skip_header=0,usecols=(0,1,2,3,4,10,11), dtype=("|U23",float,float,float,float,"|U2","|U23")) #accuracy to ms
#df = pd.DataFrame(pnsn_cat, columns=['ID','Time','Magnitude','Lat','Lon','Depth','Regional'])
#for ii in range(len(pnsn_cat)):
#    df = df.append({'ID': pnsn_cat[ii][6][2:],'Time': pnsn_cat[ii][0][11:], 'Magnitude': pnsn_cat[ii][4], 'Date': pnsn_cat[ii][0][:10],'Lat': pnsn_cat[ii][1], 'Lon': pnsn_cat[ii][2], 'Depth': pnsn_cat[ii][3], 'Regional': pnsn_cat[ii][5]}, ignore_index=True)


def plot_detc_tcs(daily_cut,template,filter_detc,outname):
    '''
        daily_cut: cutted daily data from the data_proc.cut_dailydata
        template: template .ms data in waveforms_template
        filter_detc: filter dictionary before plot
        outname: output name
    '''
    import obspy
    import numpy as np
    import matplotlib
    matplotlib.use('pdf') #instead using interactive backend
    import matplotlib.pyplot as plt
    from obspy import UTCDateTime
    from repeq import data_proc
    if type(template)==str:
        temp = obspy.read(template)
    else:
        temp = template
    if type(daily_cut)==str:
        daily_cut = np.load(daily_cut,allow_pickle=True)
        daily_cut = daily_cut.item()
    #apply filter
    daily_cut = data_proc.clean_data_cut(daily_cut,filter_detc)
    if len(daily_cut['detc_tcs'].keys())==0:
        return 1 #nothing left, just return
    OT_temp = UTCDateTime(daily_cut['OT_template']) #origin time for template
    for ik in daily_cut['detc_tcs'].keys():
        D = daily_cut['detc_tcs'][ik]
        phase = daily_cut['phase'][ik] # assume D and phase have the same order
        OT_D = UTCDateTime(ik) #origin time of Detection (cut from daily data)
        XLIM=[]
        #create figure based on how many traces
        print('Ntraces=',len(D))
        #if 0<len(D)<=20:
        #    fig = plt.figure(figsize=(8.5,5.5))
        #elif 20<len(D)<=30:
        #    fig = plt.figure(figsize=(8.5,6.5))
        #elif 30<len(D):
        fig = plt.figure(figsize=(8.5,8.5)) #all with the same size
        for ista in range(len(D)):
            net = D[ista].stats.network
            sta = D[ista].stats.station
            channel = D[ista].stats.channel
            location = D[ista].stats.location
            PS = phase[ista] #'P' or 'S'
            selected_temp = temp.select(network=net,station=sta,channel=channel,location=location)
            selected_temp = selected_temp.copy()
            #most of the case should return only 1 data, but if there's P and S in 1 station...
            if len(selected_temp)!=1:
                t1 = selected_temp[0].stats.starttime
                t2 = selected_temp[1].stats.starttime
                print('phase=',PS)
                if t2-t1>0:
                    if PS=='P':
                        selected_temp = obspy.Stream(selected_temp[0])
                        print('return first one')
                    elif PS=='S':
                        selected_temp = obspy.Stream(selected_temp[1])
                        print('return second one')
                else:
                    if PS=='P':
                        selected_temp = obspy.Stream(selected_temp[1])
                    elif PS=='S':
                        selected_temp = obspy.Stream(selected_temp[0])
                print('multiple data selected, return data based on basic PS wave assumption') #have to check this!
                #continue #!!!!!!!! deal with this later!!!!!!!!!!
            #dealing with time
            T_D = D[ista].times()
            T_temp = selected_temp[0].times() #length should only be 1, unless P/S in same data
            #Time relative to origin, so that at origin is zero
            dt_D = D[ista].stats.starttime-OT_D
            T_D = T_D+dt_D
            dt_temp = selected_temp[0].stats.starttime-OT_temp
            T_temp = T_temp+dt_temp
            #normalize data
            data_D = D[ista].data/np.max(D[ista].data)
            #data_D = D[ista].data/np.max(selected_temp[0].data) #normalize the data based on template amplitude, not daily data amplitude
            data_temp = selected_temp[0].data/np.max(selected_temp[0].data)
            #data_temp = selected_temp[0].data/np.max(D[ista].data)
            #plot both tcs
            plt.plot(T_D,data_D+ista*1.5,'k')
            if PS=='P':
                plt.plot(T_temp,data_temp+ista*1.5,'r')
            else:
                plt.plot(T_temp,data_temp+ista*1.5,'b')
            #get xlim bound
            if ista==0:
                XLIM.append(T_temp[0]-1)
        XLIM.append(T_temp[-1]+1)
        YLIM = plt.ylim()
        YLIM = [-1,ista*1.5+1]
        YLIM = [YLIM[0],YLIM[1]+0.08*(YLIM[1]-YLIM[0]) ]
        #add text
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        text_xloc = (XLIM[1]-XLIM[0])*0.04+XLIM[0]
        text_yloc = (YLIM[1]-YLIM[0])*0.86+YLIM[0]
        text_yloc_temp = (YLIM[1]-YLIM[0])*0.94+YLIM[0]
        plt.text(text_xloc,text_yloc,ik,fontsize=12,bbox=props)
        plt.text(text_xloc,text_yloc_temp,OT_temp.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-4],fontsize=12,color=[1,0,0],bbox=props)
        plt.xlabel('Origin time (s)',fontsize=15,labelpad=0)
        plt.xticks(fontsize=12)
        plt.yticks([],[])
        #add title
        plt.title('CC=%.2f'%(daily_cut['meanCC'][ik]))
        ax1 = plt.gca()
        ax1.tick_params(pad=1) #make axis closer
        plt.xlim(XLIM)
        plt.ylim(YLIM)
        #savName = template.split('_')[-1].split('.')[0] #this is the template ID
        if outname:
            print('save fig:',outname+ik+'.png')
            plt.savefig(outname+ik.replace(':','')+'.png',dpi=300)
        #plt.show()
        plt.close()




def bulk_plot_detc_tcs(home,project_name,filter_detc):
    #make figure and save in home/project_name/output/Template_match/Fig
    import glob
    daily_cuts = glob.glob(home+'/'+project_name+'/output/Template_match/Data_detection_cut/Detected_data_*.npy')
    daily_cuts.sort()
    #loop in daily cut data
    for daily_cut in daily_cuts:
        #find the ID first
        tempID = int(daily_cut.split('_')[-1].split('.')[0])
        #find its corresponding template
        template = home+'/'+project_name+'/waveforms_template/template_%05d.ms'%(tempID)
        outName = home+'/'+project_name+'/output/Template_match/Figs/template_%05d_'%(tempID)
        plot_detc_tcs(daily_cut,template,filter_detc,outName)




def plot_accNumber(home,project_name,cata_name,filter_detc,min_inter,time1,time2):
    #plot accumulated number of EQ in catalog v.s. detections
    '''
        min_inter: minimum inter event time (s)
        time1,time2: plot data between the range
    '''
    import glob
    from repeq import data_proc
    from obspy import UTCDateTime
    import datetime
    from repeq import data_proc
    import matplotlib
    matplotlib.use('pdf') #instead using interactive backend
    import matplotlib.pyplot as plt
    '''
    filter_detc = {
        'min_stan':9, #number of non-zero CC measurements
        'min_CC':0.5, #min mean(CC) value
        'diff_t':60, #time difference between events should larger than this
    }
    '''
    #load catalog and get their time
    df = data_proc.cat2pd(home+'/'+project_name+'/catalog/'+cata_name)
    template_time = [UTCDateTime(df.Date[i]+'T'+df.Time[i]) for i in range(len(df))]
    template_time = np.array(template_time)

    #load detections and get their time
    detcs = glob.glob(home+'/'+project_name+'/'+'output/Template_match/Detections/'+'Detected_tmp_*.npy')
    detcs.sort()
    detc_time = []
    for detc_path in detcs:
        detc = np.load(detc_path,allow_pickle=True)
        detc = detc.item()
        detc = data_proc.clean_detc(detc,filter_detc)
        detc_time += detc.keys()

    detc_time.sort()
    detc_time = np.array(detc_time)
    detc_time = [UTCDateTime(i) for i in detc_time]

    #set min-interevent time to remove redundant data
    clean_template_time = data_proc.clean_events_time(template_time,min_time=min_inter)
    clean_detc_time = data_proc.clean_events_time(detc_time,min_time=min_inter)

    t_temp, accnum_temp = data_proc.cal_accum(clean_template_time,time1,time2,dt=3600)
    t_detc, accnum_detc = data_proc.cal_accum(clean_detc_time,time1,time2,dt=3600)

    main_OT = UTCDateTime("2018-05-04T22:32:54.650Z").datetime #mainshock OT
    #convert UTCDateTime to datetime for plotting
    t_temp = [i.datetime for i in t_temp]
    t_detc = [i.datetime for i in t_detc]
    plt.figure(figsize=(10,4.5))
    plt.plot(t_temp,accnum_temp,'k')
    plt.plot(t_detc,accnum_detc,'r')
    print('***manually add something in plot function***')
    plt.plot([main_OT,main_OT],[0,np.max(accnum_detc)],'r--')
    plt.ylim([0,np.max(accnum_detc)])
    plt.xlim([UTCDateTime(time1).datetime,UTCDateTime(time2).datetime])
    plt.xlabel('Date',fontsize=14)
    plt.ylabel('Accumulated number',fontsize=14)
    plt.savefig(home+'/'+project_name+'/'+'output/Template_match/Detections/'+'detections.png')
    plt.close()
    #plt.show()



def plot_reptcs(home,project_name,tempID,staChn,phs,cut_window,ref_OT="2018-05-04T22:32:54.650"):
    '''
        #plot detected tcs by template ID
        tempID:template ID(e.g. '00836')
        staChn: station_name.channel (e.g. JOKA.HHZ)
        phs: phase name in case both P/S in same staChn
        cut_window: window same as when using data_proc.cut_dailydata or data_proc.bulk_cut_dailydata
        ref_OT: set y at ref_OT=0
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    import numpy as np
    from obspy import UTCDateTime
    import os
    ref_OT = UTCDateTime(ref_OT)
    #find the file to be plotted
    tcs_cut = home+'/'+project_name+'/'+'output/Template_match/Data_detection_cut/Detected_data_'+tempID+'.npy'
    file_lag = home+'/'+project_name+'/'+'output/Template_match/Measure_lag/measure_lag_temp_'+tempID+'.npy'
    flag_plotall = True
    #to see if the file exist
    if not(os.path.exist(tcs_cut)):
        print('File:%s do not exist! break'%(tcs_cut))
        return
    if not(os.path.exist(file_lag)):
        print(' Lag measurement file:%s do not exist, only show time series'%(file_lag))
        flag_plotall = False #only one subplot

    #load tcs_cut file
    D = np.load(tcs_cut,allow_pickle=True)
    D = D.item()

    #load lag measurement if file exist
    if os.path.exist(file_lag):
        MeasLag = np.load(file_lag,allow_pickle=True)
        MeasLag = MeasLag.item()

    #some plot setting(normalize, get range of y etc.)
    y_range = [i for i in D['detc_tcs'].keys()]
    y_range.sort()
    num_y = len(y_range)
    y_range = [UTCDateTime(y_range[0]),UTCDateTime(y_range[-1])] #the data spanning from y_range[0] to y_range[1]
    dy_range = (y_range[1]-y_range[0])/86400.0
    data_mul = dy_range/num_y #tcs with this amplitude should be good

    fig = plt.figure(constrained_layout=True)
    props = dict(boxstyle='round', facecolor='white', alpha=0.8)
    gs = fig.add_gridspec(3,1)
    if os.path.exist(file_lag):
        f3_ax1 = fig.add_subplot(gs[:2, 0]) #merge the first two row
    else:
        f3_ax1 = fig.add_subplot(gs[:, 0]) #merge all the three row

    #####start plotting#####
    for ik in D['detc_tcs'].keys():
        DD = D['detc_tcs'][ik].select(station=staChn.split('.')[0],channel=staChn.split('.')[1])
        if len(DD)!=1:
            #selected two phases, add phs condition and select data again
            phases = D['phase'][ik]
            #***make sure the order if D and phases are the same
            for i in range(len(DD)):
                if D[i].stats.station==staChn.split('.')[0] & D[i].stats.channel==staChn.split('.')[1]:
                    if phases[i]==phs:
                        #also check the phase
                        DD = D[i].copy()
                        break
        #selected the data, start plotting data
        time = DD[0].times()
        data = DD[0].data
        data_norm = data/np.max(data)*data_mul
        shft_plt = (UTCDateTime(ik)-ref_OT)/86400.0 # reference time in days
        f3_ax1.plot(time-cut_window[0],data_norm+shft_plt,'k',linewidth=1.0)
        f3_ax1.fill_between(time-cut_window[0],np.zeros_like(time)+shft_plt,data_norm+shft_plt,where=np.zeros_like(time)+shft_plt>data_norm+shft_plt,color=[0,0,0],interpolate=True)
    #after plotting tcs, add text and adjust axis, plot label
    x_pos = ((cut_window[1]-cut_window[0]*-1))*0.06 + +cut_window[0]*-1
    y_pos = f3_ax1.get_ylim()
    y_pos = (y_pos[1]-y_pos[0])*0.9+y_pos[0]
    f3_ax1.text(x_pos,y_pos,staChn,fontsize=12,bbox=props)
    f3_ax1.set_xlim([cut_window[0]*-1,cut_window[1]])
    f3_ax1.set_ylabel('Day relative to mainshock',fontsize=14)
    if (os.path.exist(file_lag)):
        #two subplots
        f3_ax1.set_xticklabels([]) #remove xlabel in the first subplot
    else:
        #if only one subplot
        f3_ax1.set_xlabel('Origin time (s)',fontsize=14)
        plt.savefig(home+'/'+project_name+'/'+'output/Template_match/Figs/template_'+tempID+'_'+staChn+'.'+phs+'.png')
        plt.close()
        return

    #two subplots case
    if os.path.exist(file_lag):
        f3_ax2 = fig.add_subplot(gs[-1, 0])
        f3_ax2.set_xlim([cut_window[0]*-1,cut_window[1]])
        f3_ax2.set_xlabel('Origin time (s)',fontsize=14)

    plt.savefig(home+'/'+project_name+'/'+'output/Template_match/Figs/template_'+tempID+'_'+staChn+'.'+phs+'.png')















