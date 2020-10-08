#some data pre-processing tools

import glob
import obspy



def merge_daily(home,project_name):
    Ds = glob.glob(home+'/'+project_name+'/waveforms/'+'*000000')
    Ds.sort()
    print('Number of dirs=%d'%(len(Ds)))
    for D in Ds:
        print('In dir:',D)


