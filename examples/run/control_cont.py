#control file for continuous dataset

from repeq import Repeq_starter 
from repeq import download_tools
from repeq import analysis
from obspy import UTCDateTime

init = 0             #initial project
get_catalog = 1     #download waveforms
get_waveform = 1     #download waveforms
get_template = 0

#================Parameter settings================
home = '/Users/timlin/Documents/Project/TestREPEQ'
project_name = 'QQQ'
cata_name = 'continuous.cat'
cata_name2 = 'area1_1.0.cat'  #template catalog. download this catalog by download_tools.catalog_USGS (see control.py for example)

##dealing with catalog: download or make fakecatalog
cata_times = ['20180101','20180505'] #this can be large
#--params for download only
#cata_area = [-156.357,-154.061,18.407,20.437] 
#cata_magnitude = [3.0,9.5]
#--params for making fakecatalog only
lon_lat = [-154.9996667,19.318]    #fake hypocenter
dt = 86400


##filter the catalog
cata_filters={
#'filt_times':['20170101','20190105'],
'filt_times':['19000101','20500101'],
'filt_lon':[-180, 180],
'filt_lat':[-90, 90],
'filt_dep':[0,800],
'filt_m':[0.0,10.0],
}


##waveform params for downloading based on the filtered catalog
sec_bef_aft = [0,86400] #second before and after origin (based on catalog)
Ftype = 'circ' #circ or rect
lon_lat = [120,24] #if Ftype=="circ", center of the circle. if Ftype=='rect', [lon1,lon2,lat1,lat2] 4 points boundary
range_rad = [0,6] #range of the searching circle. Must be provided if Ftype=="circ"
channel = ['BHZ','HHZ']
provider = ["IRIS"]
sampling_rate = 25
waveforms_outdir = home+'/'+project_name+'/'+'waveforms'

filter=[1,8]

#================Parameter settings END================



cata_out = home+'/'+project_name+'/catalog/'+cata_name

if init:
    Repeq_starter.create_dirs(home,project_name)

if get_catalog:
    #make fakecatalog
    download_tools.make_catalog(times=[cata_times[0],cata_times[1]],dt=dt,lon_lat=lon_lat,outname=cata_out)

if get_waveform:
    #get continuous waveform based on the fakecatalog
    download_tools.download_waves_catalog(cata_out,cata_filters,sec_bef_aft,range_rad,channel,provider,waveforms_outdir)
    #merge all the mseed data into a large merged.ms file
    from repeq import data_proc
    data_proc.merge_daily(home,project_name,sampling_rate,filter=filter,pattern='20*')

if get_template:
    #download template based on (real)catalog with manual pick
    from repeq import template
    #if download data
    T = template.Template(home,project_name,cata_name2,True,sampling_rate,filter=filter,tcs_length=[1,9],filt_CC=0.3,filt_nSTA=6,plot_check=True)
    #if data already exist
    T.download=False
    T.template_load() 

#run xcorr calculation
#T.xcorr_cont(save_CCF=False,fmt=1) #True means save all the CCF function in Template_match/CCF_records/
T.xcorr_cont(save_CCF=False,fmt=2) #fmt=2 output detailed calculation

#Or run xcorr calculation by multiprocessing
#Adjust the n-part to smaller if out-of-memory
n_part = 8
T_part = template.T_partition(T,n_part=n_part) #partitioning the T
template.T_parallel(T_part,n_part=n_part,save_CCF=False,fmt=2) #parallel for all T_part


# reading fmt=1 for T.xcorr_cont() (won't work for fmt=2)
#read all measurements in project_name/output/Template_match/Detections and make summary file based on the given filter criteria
from repeq import data_proc
filter_params={
'diff_t':60,
'min_sta':6,
'min_CC':0.3
}
data_proc.read_detections(home,project_name,filter_params) #this doesn't work for the newest format, read the data manually

#filter the dictionary by clean_detc
filter_detc = {
        'min_stan':5, #number of non-zero CC measurements
        'min_CC':0.2, #min mean(CC) value
        'diff_t':60, #time difference between events should larger than this
        }
data_proc.clean_detc(detc,filter_detc)



'''
#plot number of new detections v.s. original catalog 
from repeq import data_visual
filter_detc = {
        'min_stan':9, #number of stations/phases
        'min_CC':0.5, #min CC value
        'diff_t':60, #time difference between events should larger than this
        }
min_inter = 2.0 #sec
plot_time1 = "20180501"
plot_time2 = "20180515"
data_visual.plot_accNumber(home,project_name,cata_name,filter_detc,min_inter,plot_time1,plot_time2)
'''


#---write/cut the detection timeseries from daily data
#write detection timeseries (cut from continuous data based on the detected time)
from repeq import data_proc
filter_params={
'diff_t':60,
'min_sta':6,
'min_CC':0.3
}
cut_window=[5,20] #cut_window[t1,t2] means t1 sec "before" the pick time and t2 sec "after" the pick time
data_proc.bulk_cut_dailydata(home,project_name,filter_detc,cut_window) #the results will be saved in home+project_name/output/Template_match/Data_detection_cut

'''
#make figure from the above (cut) timeseries, result will be saved in home+project_name/output/Template_match/Figs
#detected tcs at all the stations
from repeq import data_visual
data_visual.bulk_plot_detc_tcs(home,project_name,filter_detc) #read the files in home+project_name/output/Template_match/Data_detection_cut and find template file
'''


'''
#plot detected tcs at same station
from repeq import data_visual
tempID = '00295'
tempID = '00628'
tempID = '00005'
#staChn = 'JOKA.HHZ'
#phs = 'P'
#or get all available staChn and phs if you dont know
from repeq import data_proc
staChns,phss = data_proc.get_cut_info(home,project_name,tempID)
print('Find:',staChns,phss)
staChn = staChns[0]
phs = phss[0]
cut_window = [5,20]
v_minmax = [-3,3]
data_visual.plot_reptcs(home,project_name,tempID,staChn,phs,cut_window,v_minmax=v_minmax,ref_OT="2018-05-04T22:32:54.650")
#or plot all available phase
for i,j in zip(staChns,phss):
    data_visual.plot_reptcs(home,project_name,tempID,i,j,cut_window,ref_OT="2018-05-04T22:32:54.650")
'''




#calculate lag from the detection cut
from repeq import data_proc
measure_params={
    'wind':[0.5,1],
    'mov':0.05,
    'interp':0.01,
    'taper':0.05,     #taper percentage
    }
tcs_length_temp = [1,9] #same as when download template
tcs_length_daily = [5,20] #same as when cut the data
align_wind = [[1,9],[0.5,1.5]] #(multiple layers )alignment e.g. [[1,9],[0.5,2],[0.1,1]]
#data_proc.bulk_cal_lag(home,project_name,tcs_length_temp,tcs_length_daily,align_wind,measure_params,overwrite=False,n_jobs=1,i_par=0) #original function without multiprocessing
#or use multiprocessing
data_proc.bulk_cal_lag_parallel(home,project_name,tcs_length_temp,tcs_length_daily,align_wind,measure_params,overwrite=False,n_jobs=4)






#------repEQ relocation---------
#Repeating EQ relocation
from repeq import EQreloc
#set some filtering criteria
#Filter 1: focus on each detection .npy file
filter_detc={
        'diff_t':120,  #[sec] event interval. detected events should be at least longer than this time. Only apply to the events with same templateID
        'min_stan':5,  #at least number of stations. Delete the detection if not enough station constraint
        'min_CC':0.3   #minimum mean(CC) value. Delete the detection if too low mean CC value
}

#Filter 2: focus on inversion part
invsion_params={
        'CCC_threshold':0.3,  #use observed shift with CCC greater than this threshold
        'min_stan':5,         #minumim observations
        'max_shift':2.0,        #maximum shift seconds(observation)
        'update_thres':1e-9,  # (optional) if update smaller than this, but inversion havent converged, add a perturbation
        'misfit_thres':1e-3,  # (optional) modeled arrival time close enough to observation, inversion converged
        'max_iter':20,        # (optional) maximum iterations
        'dx':0.05,            # (optional) small step in x-dir to calculate time changes
        'dy':0.05,            # (optional)
        'dz':0.05,            # (optional)
        'dt':0.04,            # (optional)
}

vel_model = 'Hawaii.litho.mod'  #1D velocity model same as fk format
T0 = UTCDateTime("2018-05-04T22:32:54.65")  #output data time will be related to this datetime
EQreloc.EQreloc(home,project_name,cata_name,vel_model,filter_detc,invsion_params,T0)  #final result will be in output/Template_match/Detections/EQreloc_info.txt
#------repEQ relocation END---------


#-----make station table-----
from repeq import data_proc
data_proc.make_sta_table(home,project_name,pattern='*000000')


#-----continuous data QC-------
from repeq import data_proc
data_proc.cal_moving_all(home,project_name,pattern='*0000',type='samp',window_pts=50000,mov_pts=25000) #the 3 files will be saved under waveforms/


#-----get CC of template-template-----
from repeq import template
import numpy as np
T = template.Template(home,project_name,cata_name,download=False,sampling_rate=sampling_rate,filter=filter,tcs_length=[1,9],filt_CC=0.2,filt_nSTA=6,plot_check=False)
T.template_load()
CC_temp = T.xcorr_temp()
np.save('CC_temp.npy',CC_temp)




'''
    #This method is too slow, please use mass_downloader instead
download_params = {
    'net':'*',
    'sta':'*',
    'comp':'BHZ,HHZ',
    'chn':'*',
    'sampl':40,
    'cent_lon':-154.9996667,
    'cent_lat':19.318,
    'min_radius':0,
    'max_radius':1,
    't1':UTCDateTime("2018-01-01"),
    't2':UTCDateTime("2018-05-05"),
    }
n_cores = 1
waveforms_outdir=home+'/'+project_name+'/'+'waveforms'
download_tools.batch_download_continuous_cent(home,project_name,download_params,n_cores=n_cores,waveforms_outdir=waveforms_outdir)

'''
