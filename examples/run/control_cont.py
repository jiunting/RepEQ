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


# reading fmt=1 for T.xcorr_cont() (won't work for fmt=2)
#read all measurements in project_name/output/Template_match/Detections and make summary file based on the given filter criteria
from repeq import data_proc
filter_params={
'diff_t':60,
'min_sta':6,
'min_CC':0.3
}
data_proc.read_detections(home,project_name,filter_params) #this doesn't work for the newest format, read the data manually


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



#-----continuous data QC-------
from repeq import data_proc
data_proc.cal_moving_all(home,project_name,pattern='*0000',type='samp',window_pts=50000,mov_pts=25000) #the 3 files will be saved under waveforms/





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
