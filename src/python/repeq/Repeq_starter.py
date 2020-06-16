#main code that initial directories 

def create_dirs(home,project_name,overwrite=False):
    import os
    import shutil

    #check file/dir exist,otherwise mkdir
    if not(os.path.exists(home+'/'+project_name)):
        os.makedirs(home+'/'+project_name) #main dir
        os.makedirs(home+'/'+project_name+'/'+'waveforms') #put all the waveforms in here
        os.makedirs(home+'/'+project_name+'/'+'structure') #put velocity model .mod file here
        main_REPEQ=os.environ['REPEQ'] #the uppest dir of the REPEQ
        shutil.copy(main_REPEQ+'/examples/model/example.mod',home+'/'+project_name+'/'+'structure/example.mod')
        os.makedirs(home+'/'+project_name+'/'+'catalog') #catalog file in USGS format
        os.makedirs(home+'/'+project_name+'/'+'output')
        os.makedirs(home+'/'+project_name+'/'+'output'+'/'+'logs')
        shutil.copy(main_REPEQ+'/examples/notes/Format_log.txt',home+'/'+project_name+'/'+'output/logs/Format_log.txt')

    else:
        if overwrite:
            print('Careful!!!!! do you really want to delete all the project?:%s'%(home+'/'+project_name))
            print('Too risky...Please delete them manually')
            return 1
        else:
            print('Directory already exist:%s'%(home+'/'+project_name))
            return 1

