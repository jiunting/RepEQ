
from repeq import Repeq_starter 
from repeq import download_tools
from repeq import analysis


init=0             #initial project
get_catalog=0      #download EQ catalog from USGS
get_waveform=0     #download waveforms
search=1           #search repEQ

#-------Parameters settings-----------
home='/Users/timlin/Documents/Project/TestREPEQ'
project_name='QQQ'
vel_model='/Users/timlin/Documents/Project/Inversion/HawaiiEQ_SSE/All_InvEQ/structure/Hawaii.litho.mod' #path of velocity model in fk format



##catalog params
cata_times=['20100102','20200101'] #this can be large
cata_area=[-156.357,-154.061,18.407,20.437]
cata_magnitude=[3.0,9.5]
cata_name='area1.cat'

##filter the catalog and download the data by these criteria

cata_filters={
'filt_times':['20170101','20190105'],
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
#'freq': bandpass filter range
#window(secs) before and after P arrival from vel_model prediction
data_filters={
'freq':(0.8,2),
'window':(0,30),
}
startover=False   #True: re-run everything, or False: check the files that already exist to see if they need to be updated
make_fig_CC=2     #plot figure when measured CC larger than this number


#-------Parameters settings END-----------

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
    analysis.searchRepEQ(home,project_name,vel_model,cata_name,data_filters,startover=startover,make_fig_CC=make_fig_CC)
    analysis.read_logs(home,project_name,data_filters,outdir='')





