#control file for continuous dataset

from repeq import Repeq_starter 
from repeq import download_tools
from repeq import analysis
from obspy import UTCDateTime

init=0             #initial project
get_waveform=0     #download waveforms

#================Parameters settings================
home='/Users/timlin/Documents/Project/TestREPEQ'
project_name='QQQ'


#download continuous waveform
download_params = {
'sta':'*',
'comp':'*H*',
'chn':'*',
'sampl':40,
'cent_lon':
'cent_lat':
'min_radius':
'max_radius':
't1':UTCDateTime("2018-01-01"),
't2':UTCDateTime("2018-05-05"),
}

waveforms_outdir=home+'/'+project_name+'/'+'waveforms'



if init:
    Repeq_starter.create_dirs(home,project_name)

if get_waveform:
    download_tools.download_continuous_cent




##catalog params
cata_times=['20000101','20200101'] #this can be large
cata_area=[-156.357,-154.061,18.407,20.437]
cata_magnitude=[3.0,9.5]
cata_name='area1.cat'

##filter the catalog and download the data by these criteria
cata_filters={
#'filt_times':['20170101','20190105'],
'filt_times':['20000101','20180515'],
'filt_lon':[-156.357, -154.061],
'filt_lat':[18.407, 20.437],
'filt_dep':[0,30],
'filt_m':[4.0,9.0],
}

##waveform params for downloading
sec_bef_aft=[120,600] #second before and after origin (based on catalog)
Ftype='circ' #circ or rect
lon_lat=[120,24] #if Ftype=="circ", center of the circle. if Ftype=='rect', [lon1,lon2,lat1,lat2] 4 points boundary
range_rad=[0,6] #range of the searching circle. Must be provided if Ftype=="circ"
channel=['BHZ','HHZ']
provider=["IRIS"]
waveforms_outdir=home+'/'+project_name+'/'+'waveforms'

##waveform params for cross-correlation measurement
data_filters={
'freq':(0.8,2),  ## bandpass filter range
'window':(0,30), ##window(secs) before and after P arrival from vel_model prediction
'max_sepr':0.5,  ##max separation for the two events
}
startover=False   #True: re-run everything, or False: check the files that already exist to see if they need to be updated
make_fig_CC=2     #plot figure when measured CC larger than this number

##params for sequence defination
# This example means at least 2 stations with CC>=0.9
seq_filters={
'min_nsta':2,
'min_CC':0.9,
'time_sep':86400*2 #p1-p2 needs to be separated by at least n sec
}

##params for lag measurements
#-----parameters for correcting P arrival-----#
lag_params={
'filt_freq_HR':[(0.5,2),(0.5,2)],  #set n-step Pwave corrections
'p_wind':[(5,15),(3,3)], #window for Pwave correction. Seconds before(positive!) and after theoritical P arrival
'CCC_thres':0.9, #threshold for repEQ from log file
'CCsta_thres':0.9, #threshold for individual station
'min_num':1, #at least n stations got this threshold
'L_wind':(20,150), #Total(large window) data to be measured. Seconds before, after corrected P arrival
'filt_L_wind':(0.5,2), #filter for the large window
'S_wind':6, # n seconds for S(small window) of measurement each time
'mov':0.2, # moving seconds
'sampt':0.005, #interpolate to this interval
'Write_out':True, #write measured lag?
}
sequence_file='QQQ.sequence' #sequence file used for measure_lag.  Note that default name of sequence file is project_name.sequence



#================Parameters settings END================




cata_out=home+'/'+project_name+'/catalog/'+cata_name

if init:
    Repeq_starter.create_dirs(home,project_name)

if get_catalog:
    download_tools.catalog_USGS(cata_times,cata_area,cata_magnitude,cata_out)

if get_waveform:
    print('Load catalog:',cata_out)
    download_tools.download_waves_catalog(cata_out,cata_filters,sec_bef_aft,range_rad,channel,provider,waveforms_outdir)


if search:
#    from obspy.taup import TauPyModel
#    TauPy_name=analysis.build_TauPyModel(home,project_name,vel_model) #make .npz file
#    analysis.searchRepEQ(home,project_name,vel_model,cata_name,data_filters,startover=startover,make_fig_CC=make_fig_CC,QC=True,save_note=True)
#    analysis.read_logs(home,project_name) #merge all the .log file into a large summary file: project_name.summary
#    analysis.sequence(home,project_name,seq_filters) #make sequence file: project_name.summary
    analysis.measure_lag(home,project_name,lag_params,sequence_file,cata_name)


#download continuous data
#create fake catalog
#download_tools.make_catalog(times=['2018','20180506'],dt=86400,lon_lat=[-154.9996667,19.3181667],outname=cata_out)
#then download data by
#download_tools.download_waves_catalog(cata_out,cata_filters,sec_bef_aft,range_rad,channel,provider,waveforms_outdir)
