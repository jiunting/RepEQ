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
#### The RepEQ calls obspy's MassDownloader to download event-based or continuous data 
> Simply copy example file control.py for event based waveform downloader or control_cont.py continuous data downloader

[libcomcat]:https://github.com/usgs/libcomcat "libcomcat is a project designed to provide a Python equivalent to the ANSS ComCat search API"
