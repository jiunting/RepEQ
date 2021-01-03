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



def my_colormap():
    #define my_colormap
    '''
        regular seismic change from blue-white-white-red
        make new seismic change from white-blue-red-white
    '''
    import matplotlib
    N = 1000
    tmp_seis = plt.cm.get_cmap('bwr',N) #from blue-white-white-red
    #make the new seis to be white-blue-red-white
    tmpRGB = tmp_seis(range(N))
    R1 = np.flipud(tmpRGB[:N//2,0])
    G1 = np.flipud(tmpRGB[:N//2,1])
    B1 = np.flipud(tmpRGB[:N//2,2])

    #tmp_seis = plt.cm.get_cmap('spring',N//2) #concate whole spring
    #tmp_seis = plt.cm.get_cmap('hsv',N) #concate hsv begin half
    tmp_seis = plt.cm.get_cmap('hot',N) #concate hot last half
    tmpRGB = tmp_seis(range(N))
    R2 = tmpRGB[N//2:,0]
    G2 = tmpRGB[N//2:,1]
    B2 = tmpRGB[N//2:,2]
#    R2 = tmpRGB[N//2:,0]
#    G2 = tmpRGB[N//2:,1]
#    B2 = tmpRGB[N//2:,2]
#    R2 = np.flipud(tmpRGB[N//2:,0])
#    G2 = np.flipud(tmpRGB[N//2:,1])
#    B2 = np.flipud(tmpRGB[N//2:,2])
    print('Length R1-R2=',len(R1),len(R2))
    R = np.hstack([R1,R2])
    G = np.hstack([G1,G2])
    B = np.hstack([B1,B2])

    cdict = {'red':[],'green':[],'blue':[]}
    #fill in the RGB value into cdict
    pre_R = pre_G = pre_B = 0.0
    for i in range(len(R)):
        cdict['red'].append(((i)/len(R), pre_R, R[i]))
        cdict['green'].append(((i)/len(G), pre_G, G[i]))
        cdict['blue'].append(((i)/len(B), pre_B, B[i]))
        pre_R = R[i]
        pre_G = G[i]
        pre_B = B[i]

    cdict['red'].append((1,R[-1],R[-1]))
    cdict['green'].append((1,G[-1],G[-1]))
    cdict['blue'].append((1,B[-1],B[-1]))

    my_colormap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,N)
    return my_colormap


def plot_reptcs(home,project_name,tempID,NetStaChnLoc,phs,cut_window,v_minmax=[-3,3],ref_OT="2018-05-04T22:32:54.650"):
    '''
        #plot detected tcs by template ID
        tempID:template ID(e.g. '00836')
        NetStaChnLoc: net.station_name.channel.loc (e.g. HV.JOKA.HHZ. )
        phs: phase name in case both P/S in same NetStaChnLoc
        cut_window: window same as when using data_proc.cut_dailydata or data_proc.bulk_cut_dailydata
                    Note that time information of lag measurements is provided in measure_lag_temp*.npy file
        v_minmax: range of colormap plotting for lag measurement
        ref_OT: set y at ref_OT=0
    '''
    import numpy as np
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    import numpy as np
    from obspy import UTCDateTime,Stream
    import os
    
    ref_OT = UTCDateTime(ref_OT)
    fullName = NetStaChnLoc+'.'+phs
    #find the file to be plotted
    tcs_cut = home+'/'+project_name+'/'+'output/Template_match/Data_detection_cut/Detected_data_'+tempID+'.npy'
    file_lag = home+'/'+project_name+'/'+'output/Template_match/Measure_lag/measure_lag_temp_'+tempID+'.npy'
    flag_plotall = True
    #to see if the file exist
    if not(os.path.exists(tcs_cut)):
        print('File:%s do not exist! break'%(tcs_cut))
        return
    if not(os.path.exists(file_lag)):
        print(' Lag measurement file:%s do not exist, only show time series'%(file_lag))
        flag_plotall = False #only one subplot

    #load tcs_cut file
    D = np.load(tcs_cut,allow_pickle=True)
    D = D.item()

    #load lag measurement if file exist
    if os.path.exists(file_lag):
        MeasLag = np.load(file_lag,allow_pickle=True)
        MeasLag = MeasLag.item()

    #some plot setting(normalize, get range of y etc.)
    y_range = [i for i in D['detc_tcs'].keys()]
    y_range.sort()
    num_y = len(y_range)
    if len(y_range)==1:
        y_range = [UTCDateTime(y_range[0])-86400,UTCDateTime(y_range[0])+86400]
    else:
        y_range = [UTCDateTime(y_range[0]),UTCDateTime(y_range[-1])] #the data spanning from y_range[0] to y_range[1]
    dy_range = (y_range[1]-y_range[0])/86400.0
    data_mul = dy_range/num_y #tcs with this amplitude should be good

    fig = plt.figure(constrained_layout=True)
    props = dict(boxstyle='round', facecolor='white', alpha=0.8)
    gs = fig.add_gridspec(3,1)
    if os.path.exists(file_lag):
        f3_ax1 = fig.add_subplot(gs[:2, 0]) #merge the first two row
    else:
        f3_ax1 = fig.add_subplot(gs[:, 0]) #merge all the three row

    #####start plotting#####
    for ik in D['detc_tcs'].keys():
        #print('--select:',NetStaChnLoc,'from D["detc_tcs"][%s]'%(ik))
        DD = D['detc_tcs'][ik].select(network=NetStaChnLoc.split('.')[0],station=NetStaChnLoc.split('.')[1],channel=NetStaChnLoc.split('.')[2],location=NetStaChnLoc.split('.')[3])
        if len(DD)==0:
            continue #note that NetStaChnLoc doesnt always in every D['detc_tcs'][ik]
        if len(DD)!=1:
            #selected two or more phases, add phs condition and select data again
            phases = D['phase'][ik]
            #***make sure the order if D and phases are the same
            for i in range(len(DD)):
                if (DD[i].stats.network==NetStaChnLoc.split('.')[0]) & (DD[i].stats.station==NetStaChnLoc.split('.')[1]) & (DD[i].stats.channel==NetStaChnLoc.split('.')[2]) & (DD[i].stats.location==NetStaChnLoc.split('.')[3]):
                    if phases[i]==phs:
                        #also check the phase
                        DD = Stream(DD[i].copy())
                        break
        #selected the data, start plotting data
        #print('----data selected:',DD)
        time = DD[0].times()
        data = DD[0].data
        data_norm = data/np.max(data)*data_mul
        shft_plt = (UTCDateTime(ik)-ref_OT)/86400.0 # reference time in days
        f3_ax1.plot(time-cut_window[0],data_norm+shft_plt,'k',linewidth=1.0)
        f3_ax1.fill_between(time-cut_window[0],np.zeros_like(time)+shft_plt,data_norm+shft_plt,where=np.zeros_like(time)+shft_plt>data_norm+shft_plt,color=[0,0,0],interpolate=True)
    #after plotting tcs, add text and adjust axis, plot label
    x_pos = ((cut_window[1]-cut_window[0]*-1))*0.03 + +cut_window[0]*-1
    y_pos = f3_ax1.get_ylim()
    y_pos = (y_pos[1]-y_pos[0])*0.9+y_pos[0]
    f3_ax1.text(x_pos,y_pos,fullName,fontsize=12,bbox=props)
    f3_ax1.set_xlim([cut_window[0]*-1,cut_window[1]])
    f3_ax1.set_ylabel('Day relative to mainshock',fontsize=14)
    if (os.path.exists(file_lag)):
        #two subplots
        f3_ax1.set_xticklabels([]) #remove xlabel in the first subplot
    else:
        #if only one subplot
        f3_ax1.set_xlabel('Arrival time (s)',fontsize=14)
        plt.savefig(home+'/'+project_name+'/'+'output/Template_match/Figs/reptcs_'+tempID+'_'+fullName+'.png')
        print('figure saved:',home+'/'+project_name+'/'+'output/Template_match/Figs/reptcs_'+tempID+'_'+fullName+'.png')
        plt.close()
        return

    #use vmin and vmax from input
    vmin = v_minmax[0]
    vmax = v_minmax[1]

    #if not returned, continue to two subplots case
    f3_ax2 = fig.add_subplot(gs[-1, 0])
    f3_ax2.set_xlim([cut_window[0]*-1,cut_window[1]])
    f3_ax2.set_ylim([-0.2,0.2])
    f3_ax2.set_xlabel('Arrival time (s)',fontsize=14)
    f3_ax2.set_ylabel(r'$\tau$ (s)',fontsize=14)

    #get range of day relative to reftime (for different color)
    iks = [ik for ik in MeasLag['detc_OT'].keys()]
    iks.sort()
    iks_ref = np.array([(UTCDateTime(ik)-ref_OT)/86400.0 for ik in iks])
    print('iks_ref=',iks_ref)
    print('***Fix the vmin,vmax to %f,%f***'%(vmin,vmax))
    #cmap_ref = plt.cm.seismic(plt.Normalize(iks_ref[0],iks_ref[-1])(iks_ref)) #use the vminmax from data
    
    #cmap_ref = plt.cm.seismic(plt.Normalize(-10,10)(iks_ref)) #use the define seismic
    my_color = my_colormap()
    print('use manual seismic colormap')
    cmap_ref = my_color(plt.Normalize(vmin,vmax)(iks_ref))

    ik_color = {} #make color table
    for i in range(len(iks_ref)):
        ik_color[iks[i]] = cmap_ref[i]


    plt.plot([cut_window[0]*-1,cut_window[1]],[0,0],'k--',linewidth=0.3) #plot horizontal line
    #loop all the available measurements
    for ik in MeasLag['detc_OT'].keys():
        if fullName in MeasLag['detc_OT'][ik]:
            lag_time = MeasLag['detc_OT'][ik][fullName]['time']
            lag_shift = MeasLag['detc_OT'][ik][fullName]['shift']
            lag_CCC = MeasLag['detc_OT'][ik][fullName]['CCC']
            #if np.mean(lag_CCC)<0.5:
            #    continue
            #plt_idx = np.where(lag_CCC>=0.5)[0]
            if (UTCDateTime(ik)-ref_OT)/86400.0 < vmin or (UTCDateTime(ik)-ref_OT)/86400.0>vmax:
                continue
            plt.plot(lag_time,lag_shift,color=ik_color[ik],linewidth=0.5)
            #plt.plot(lag_time[plt_idx],lag_shift[plt_idx],'.-',color=ik_color[ik],linewidth=0.3,markersize=0.1)

    #add colormap
    #norm = matplotlib.colors.Normalize(vmin=iks_ref[0], vmax=iks_ref[-1])
    print('***Fix the vmin,vmax to %f,%f***'%(vmin,vmax))
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    #cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap='seismic')
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=my_color)
    cmap.set_array([])

    #These two lines mean put the bar inside the plot
    cbaxes = fig.add_axes([0.73, 0.3, 0.12, 0.022])
    clb=plt.colorbar(cmap,cax=cbaxes, orientation='horizontal',label='Day')
    clb.set_label('Day', rotation=0,labelpad=0)

    plt.savefig(home+'/'+project_name+'/'+'output/Template_match/Figs/reptcs_'+tempID+'_'+fullName+'.png')
    print('figure-2subplots saved:',home+'/'+project_name+'/'+'output/Template_match/Figs/reptcs_'+tempID+'_'+fullName+'.png')



def plot_lag_all(home,project_name,cata_name,sta_name,filter_slope,ref_OT="2018-05-04T22:32:54.650",coast_path=''):
    import matplotlib
    matplotlib.use('pdf') #instead using interactive backend
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    from repeq import data_proc
    from repeq.EQreloc import get_lonlat
    '''
        read home/project_name/output/Template_match/Measure_lag/measure_lag_all.npy and plot slope of shift measurements
        cata_name: catalog name
        sta_name: station table created from data_proc.make_sta_table
        filter_slope: parameter for measuring slope
        filter_slope = {
        'diff_t':60,    #minimum dt between template and detected_OT (note the definition is different from filter_detc used by data_proc.bulk_cut_dailydata)
        'aligned_CC':0.7,
        'measured_CC':0.5,
        'max_shift':0.5, #drop the shift larger than this number (very large shift due to cycle slip)
        'min_length':0.7, #length of the available time series pass the above criteria (0~1)
        }
    '''
    sns.set()
    sns.set_palette('husl',n_colors=10)

    # load coast data if given
    if coast_path:
        coast=np.genfromtxt(coast_path)

    # load catalog
    df = data_proc.cat2pd(home+'/'+project_name+'/catalog/'+cata_name)

    # load station table in home/project_nam/stations/stations.txt
    sta_table = pd.read_table(home+'/'+project_name+'/catalog/'+sta_name,header=None,names=['stlon','stlat','stelev','stname'],sep=' ')

    #load all lag measurements (this is a huge file, make sure memory fit)
    lag_all = np.load(home+'/'+project_name+'/output/Template_match/Measure_lag/'+'measure_lag_all.npy',allow_pickle=True)
    lag_all = lag_all.item()


    sav_slope = {} #with sta as key
    for ik in lag_all.keys():
        temp_OT = lag_all[ik]['template_OT']
        #loop all the detections
        for detc_OT in lag_all[ik]['detc_OT'].keys():
            if np.abs(UTCDateTime(temp_OT)-UTCDateTime(detc_OT))<filter_slope['diff_t']:
                continue #dont want the detection too close to template(basically itself)

            #loop every stations
            for sta in lag_all[ik]['detc_OT'][detc_OT].keys():
                #in each lag_all[ik]['detc_OT'][detc_OT][sta] there are "time","shift","CCC" keys
                time = lag_all[ik]['detc_OT'][detc_OT][sta]['time']
                shift = lag_all[ik]['detc_OT'][detc_OT][sta]['shift']
                CCC = lag_all[ik]['detc_OT'][detc_OT][sta]['CCC']

                #find the alignedCC (time closest to zero)
                zeroidx = np.where(np.min(np.abs(time)))[0][0]
                if CCC[zeroidx]<filter_slope['aligned_CC']:
                    continue #alignment is not robust
                        
                #take the time,shift and fit by a slope
                idx = np.where((CCC>=filter_slope['measured_CC']) & (shift<filter_slope['max_shift']) )[0] #shift cannt be too large otherwise is cycle slip
                if len(idx) >= len(shift)*filter_slope['min_length']:
                    #plt.plot(time[idx],shift[idx])
                    #plt.show()
                    #80% data pass threshold, then calculate slope
                    M = data_proc.cal_slope(time[idx],shift[idx])
                    #sav_slope.append(M[1]) #M[0] is intercept, M[1] is slope
                    G = np.hstack([np.ones([len(idx),1]),time[idx].reshape(-1,1)])
                    yhat = np.dot(G,M.reshape(-1,1))
                    fit_std = np.std(yhat-shift[idx]) #standard deviation of misfit
                    #sav_std.append(fit_std)
                    #get the reference time WRS to template_OT
                    ref_tempT = (UTCDateTime(detc_OT)-UTCDateTime(temp_OT))/86400.0 #relative days from template
                    #sav_reftime.append(ref_tempT)
                    #create new sta key if its not there
                    if not (sta in sav_slope):
                        sav_slope[sta] = {temp_OT:{'slope':[],'std':[],'ref_time':[],'ID':ik}}
                    #template_OT as new key
                    if not (lag_all[ik]['template_OT'] in sav_slope[sta]):
                        sav_slope[sta][temp_OT] = {'slope':[],'std':[],'ref_time':[],'ID':ik}
                    #appending data
                    sav_slope[sta][temp_OT]['slope'] = np.hstack([sav_slope[sta][temp_OT]['slope'],M[1]]) #appending data as array
                    sav_slope[sta][temp_OT]['std'] = np.hstack([sav_slope[sta][temp_OT]['std'],fit_std]) #
                    sav_slope[sta][temp_OT]['ref_time'] = np.hstack([sav_slope[sta][temp_OT]['ref_time'],ref_tempT])


    #=========plot station result============
    for sta in sav_slope.keys():
        plt.figure(figsize=(8.5,4.5))
        plt.subplot(1,2,1)
        #loop templates
        n_meas = 0 #n-measurements
        sav_tmplon = [] #save template lon
        sav_tmplat = [] #template lat
        sav_h = []
        for temp in sav_slope[sta].keys():
            #for each template, all the measurements at this station
            time = sav_slope[sta][temp]['ref_time'] #0 is the template time
            slope = sav_slope[sta][temp]['slope']
            stdn = sav_slope[sta][temp]['std']
            #add templat itself in the data
            time = np.hstack([time,0])
            slope = np.hstack([slope,0])
            stdn = np.hstack([stdn,0])
            #sort
            sor_idx = np.argsort(time)
            time = time[sor_idx]
            slope = slope[sor_idx]
            stdn = stdn[sor_idx]
            dt_main = (UTCDateTime(temp)-main_OT)/86400.0
            time += dt_main
            #check if slope measurement across the mainshock
            #if (time.min()<0) & (time.max()>0):
            if (time.max()>0):
                #plt.errorbar(time,slope+n_meas,stdn)
                h = plt.plot(time,slope+n_meas*0.01,'.-')
                sav_h.append(h[0])
                plt.plot(dt_main,0+n_meas*0.01,'rv')
                plt.plot([-10,10],[n_meas*0.01,n_meas*0.01],'k--',linewidth=0.5)
                #save template information (loc)
                sav_tmplon.append(df.iloc[ int(sav_slope[sta][temp]['ID']) ].Lon)
                sav_tmplat.append(df.iloc[ int(sav_slope[sta][temp]['ID']) ].Lat)
                n_meas += 1
        if n_meas == 0:
            plt.close()
            continue
        plt.plot([0,0],[-0.01,n_meas*0.01],'r',linewidth=0.5)
        plt.xlim([-5,5])
        plt.ylim([-0.01,n_meas*0.01])
        plt.yticks([],[])
        plt.xlabel('days from mainshock',fontsize=15,labelpad=0)
        plt.title(sta,fontsize=15)
        plt.grid(False)
        #another subplot plot map
        plt.subplot(1,2,2)
        for itmp in range(len(sav_tmplon)):
            plt.plot(sav_tmplon[itmp],sav_tmplat[itmp],'o',color=sav_h[itmp].get_color(),markeredgecolor=[0,0,0],mew=0.8,alpha=0.9)
        if coast_path:
            plt.plot(coast[:,0],coast[:,1],'k-')
        #get station lon,lat
        stlon,stlat = get_lonlat(sta_table,[sta])
        plt.plot(stlon,stlat,'^',markersize=10,color=[0,1,0],markeredgecolor=[0,0,1],mew=1)
        print('***manually plot mainshock loc, set xlim and ylim')
        plt.plot(-154.9996667,19.3181667,'r*',markersize=14,markeredgecolor=[0,0,0],mew=1,alpha=0.9)
        plt.xlim([-155.85,-154.74])
        plt.ylim([18.86,19.88])
        plt.xticks(rotation=30,fontsize=10)
        plt.savefig(home+'/'+project_name+'/output/Template_match/Figs/'+'slopeSummary_%s.png'%(sta))
        plt.close()
        #plt.show()


