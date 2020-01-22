#This example shows how to download waveforms from single event

from repeq import download_tools

#-------Parameters settings-----------
evtime='2000-01-01T06:58:39.78'
sec_bef_aft=[120,600]
Ftype='circ'
lon_lat=[120,24]
range_rad=[0,6]
provider=["IRIS"]
OUT='./QQQ'

#-------Run the code------------------
download_tools.download_waves(evtime,sec_bef_aft,Ftype,lon_lat,range_rad,provider,OUT)
