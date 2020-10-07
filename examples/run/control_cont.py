#control file for continuous dataset

from repeq import Repeq_starter 
from repeq import download_tools
from repeq import analysis
from obspy import UTCDateTime

init=0             #initial project
get_waveform=1     #download waveforms

#================Parameters settings================
home='/Users/timlin/Documents/Project/TestREPEQ'
project_name='QQQ'


if init:
    Repeq_starter.create_dirs(home,project_name)

if get_waveform:
    #make fakecatalog
    download_tools.make_catalog(times=['20180101','20180505'],dt=86400,lon_lat=[-154.9996667,19.318],outname='continuous.cat')




'''
    #This method is no longer needed, please use mass_downloader instead
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
