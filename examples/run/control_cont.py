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
    data_proc.merge_daily(home,project_name,sampling_rate,filter=None)

if get_template:
    #download template based on (real)catalog with manual pick
    from repeq import template
    #if download data
    T = template.Template(home,project_name,cata_name2,True,sampling_rate,filter=[0.2,8],tcs_length=[1,9],filt_CC=0.2,filt_nSTA=5)
    #if data already exist
    T.download=False
    T.template_load() 

#run xcorr
T.xcorr_cont(save_CCF=False) #True means save all the CCF function in Template_match/CCF_records/


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
