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
'cent_lon':-154.9996667,
'cent_lat':19.318,
'min_radius':0,
'max_radius':1,
't1':UTCDateTime("2018-01-01"),
't2':UTCDateTime("2018-05-05"),
}

waveforms_outdir=home+'/'+project_name+'/'+'waveforms'



if init:
    Repeq_starter.create_dirs(home,project_name)

if get_waveform:
    download_tools.bulk_download_continuous_cent(home,project_name,download_params,waveforms_outdir=waveforms_outdir)

