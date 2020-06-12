#This example shows how to get the event catalog from USGS API

from repeq import download_tools


#-------Parameters settings-----------
times=['20100102','2020']
area=[-156.357,-154.061,18.407,20.437]
magnitude=[3.0,9.5]
outname='testout.cat'



download_tools.catalog_USGS(times,area,magnitude,outname)

