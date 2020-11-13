#functions that show the data

import matplotlib.pyplot as plt
import numpy as np


#def plot_map(home,project_name,cata_file,create_plot=True,show=False):
#    if create_plot:
#        fig = plt.figure()
#    #load catalog in pandas
#    pnsn_cat = np.genfromtxt('/Users/timlin/Documents/Project/Hawaii_cont/foreshock2018/catalog/area2_1.0.cat', delimiter=',', skip_header=0,usecols=(0,1,2,3,4,10,11), dtype=("|U22",float,float,float,float,"|U2","|U22"))
##cat = np.genfromtxt('/Users/timlin/Documents/Project/Hawaii_cont/foreshock2018/catalog/area1_1.0.cat', delimiter=',', skip_header=0,usecols=(0,1,2,3,4,10,11), dtype=("|U23",float,float,float,float,"|U2","|U23")) #accuracy to ms
#df = pd.DataFrame(pnsn_cat, columns=['ID','Time','Magnitude','Lat','Lon','Depth','Regional'])
#for ii in range(len(pnsn_cat)):
#    df = df.append({'ID': pnsn_cat[ii][6][2:],'Time': pnsn_cat[ii][0][11:], 'Magnitude': pnsn_cat[ii][4], 'Date': pnsn_cat[ii][0][:10],'Lat': pnsn_cat[ii][1], 'Lon': pnsn_cat[ii][2], 'Depth': pnsn_cat[ii][3], 'Regional': pnsn_cat[ii][5]}, ignore_index=True)


def plot_detc_tcs(daily_cut,template,outname):
    '''
        daily_cut: cutted daily data from the data_proc.cut_dailydata
        template: template .ms data in waveforms_template
        outname: output name
    '''
    import obspy
    import numpy as np
    import matplotlib
    matplotlib.use('pdf') #instead using interactive backend
    import matplotlib.pyplot as plt
    from obspy import UTCDateTime
    
    temp = obspy.read(template)
    if type(daily_cut)==str:
        daily_cut = np.load(daily_cut,allow_pickle=True)
        daily_cut = daily_cut.item()
    OT_temp = UTCDateTime(daily_cut['OT_template']) #origin time for template
    for ik in daily_cut['detc_tcs'].keys():
        D = daily_cut['detc_tcs'][ik]
        phase = daily_cut['phase'][ik]
        OT_D = UTCDateTime(ik) #origin time of Detection (cut from daily data)
        XLIM=[]
        #create figure based on how many traces
        print('Ntraces=',len(D))
        if 0<len(D)<=20:
            fig = plt.figure(figsize=(8.5,5.5))
        elif 20<len(D)<=30:
            fig = plt.figure(figsize=(8.5,6.5))
        elif 30<len(D)<=40:
            fig = plt.figure(figsize=(8.5,8.5))
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
            data_temp = selected_temp[0].data/np.max(selected_temp[0].data)
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
        ax1 = plt.gca()
        ax1.tick_params(pad=1) #make axis closer
        plt.xlim(XLIM)
        plt.ylim(YLIM)
        #savName = template.split('_')[-1].split('.')[0] #this is the template ID
        if outname:
            print('save fig:',outname+ik+'.png')
            plt.savefig(outname+ik+'.png',dpi=300)
        #plt.show()
        plt.close()




def bulk_plot_detc_tcs(home,project_name):
    #make figure and save in home/project_name/output/Template_match/Fig
    import glob
    daily_cuts = glob.glob(home+'/'+project_name+'/output/Template_match/Data_detection_cut/Detected_data_*.npy')
    daily_cuts.sort()
    #loop in daily cut data
    for daily_cut in daily_cuts:
        #find the ID first
        tempID = int(daily_cut.split('_')[-1].split('.')[0])
        #find its corresponding template
        template = home+'/'+project_name+'/waveforms_template/template_%05d.npy'%(tempID)
        plot_detc_tcs(daily_cut,template)
        outName = home+'/'+project_name+'/output/Template_match/Figs/template_%05d_'%(tempID)
        plot_detc_tcs(daily_cut,template,outName)
















