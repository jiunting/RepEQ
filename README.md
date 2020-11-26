# RepEQ
### Event based or template-matching repeating earthquake searching and analyzing tool

****
## 1. Installation
> RepEQ uses phase picks from the ANSS(USGS) catalog, so make sure install [libcomcat][] first 

#### cd to the place where you want to put the source code  
```console
cd Your_Local_Path  
git clone https://github.com/jiunting/RepEQ.git
```

#### Add RepEQ to PYTHONPATH

> Go to your environval variable file (.base_profile or .bashrc)  
```console
vi ~/.bashrc  
```
> or  
```console
vi ~/.bash_profile      
```
> and add the following line in the file

```bash
#set MLARGE
export PYTHONPATH=$PYTHONPATH:YOUR_PATH_MARGE/RepEQ/src/python
```

## 2. Download catalog  
#### 2-1 RepEQ uses USGS's API independently from the libcomcat to download events  
> Simply copy example file control.py and modify the parameters for event based catalog.  
```python
#in control file
download_tools.catalog_USGS(cata_times,cata_area,cata_magnitude,cata_out)
```
>The function takes 4 inputs  

|Variable Name  |Meaning |
| :---------- | :-----------|
| cata_times   |<array or list; len=2; dtype=str or datetime> i.e. [t1,t2] the begining and ending of catalog. |
| cata_area   |<array or list; len=4; dtype=float> area defined by 4-points [lon_min, lon_max, lat_min, lat_max]   |
| cata_magnitude   |<array or list; len=2; dtype=float> magnitude range [mag_min, mag_max]   |
| cata_name   |<str> output name   |

#### 2-2 RepEQ can also generate fake catalog for downloading continuous data later
> Copy example file control_cont.py and modify the parameters. 
```python
#in control file
download_tools.make_catalog(times=[cata_times[0],cata_times[1]],dt=dt,lon_lat=lon_lat,outname=cata_out)
```
> The function is similar to example 2-1 except the dt, which controls the sampling interval of the generated time.  
> For daily data, set dt=86400.

## 3. Download waveforms
#### RepEQ download waveforms based on the catalog generated from the above
> Copy example file control.py or control.py then set the time (i.e. how long the timeseries to be downloaded) and filter (i.e. which event should be downloaded)
```python
#in control file
download_tools.download_waves_catalog(cata_out,cata_filters,sec_bef_aft,range_rad,channel,provider,waveforms_outdir)
```



[libcomcat]:https://github.com/usgs/libcomcat "libcomcat is a project designed to provide a Python equivalent to the ANSS ComCat search API"
