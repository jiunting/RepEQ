
from repeq import Repeq_starter 
from repeq import download_tools

init=0             #initial project
get_catalog=0      #download EQ catalog from USGS
get_waveform=1     #

#-------Parameters settings-----------
home='/Users/timlin/Documents/Project/TestREPEQ'
project_name='QQQ'

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

##waveforms params
sec_bef_aft=[120,600] #second before and after origin (based on catalog)
Ftype='circ' #circ or rect
lon_lat=[120,24] #if Ftype=="circ", center of the circle. if Ftype=='rect', [lon1,lon2,lat1,lat2] 4 points boundary
range_rad=[0,6] #range of the searching circle. Must be provided if Ftype=="circ"
channel=['BHZ','HHZ']
provider=["IRIS"]
waveforms_outdir=home+'/'+project_name+'/'+'waveforms'

#-------Parameters settings END-----------

cata_out=home+'/'+project_name+'/catalog/'+cata_name



if init:
    Repeq_starter.create_dirs(home,project_name)

if get_catalog:
    download_tools.catalog_USGS(cata_times,cata_area,cata_magnitude,cata_out)

if get_waveform:
    print('Load catalog:',cata_out)
    download_tools.download_waves_catalog(cata_out,cata_filters,sec_bef_aft,range_rad,channel,provider,waveforms_outdir)


