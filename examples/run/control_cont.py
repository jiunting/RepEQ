#control file for continuous dataset

from repeq import Repeq_starter 
from repeq import download_tools
from repeq import analysis
from obspy import UTCDateTime

init = 0             #initial project
get_catalog = 1     #download waveforms
get_waveform = 1     #download waveforms

#================Parameter settings================
home = '/Users/timlin/Documents/Project/TestREPEQ'
project_name = 'QQQ'
cata_name = 'continuous.cat'


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
'filt_times':['20000101','20300101'],
'filt_lon':[-180, 180],
'filt_lat':[-90, 90],
'filt_dep':[0,300],
'filt_m':[0.0,10.0],
}


##waveform params for downloading based on the filtered catalog
sec_bef_aft = [0,86400] #second before and after origin (based on catalog)
Ftype = 'circ' #circ or rect
lon_lat = [120,24] #if Ftype=="circ", center of the circle. if Ftype=='rect', [lon1,lon2,lat1,lat2] 4 points boundary
range_rad = [0,6] #range of the searching circle. Must be provided if Ftype=="circ"
channel = ['BHZ','HHZ']
#channel = 'BHZ,HHZ' #try this on
provider = ["IRIS"]
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