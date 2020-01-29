#Make synthetic data
import obspy
import datetime
import time
import shutil
import os
import numpy as np
import multiprocessing as mp
import subprocess
import random


#########################default parameters#############################
#home='/Users/timlin/Documents/Project/GPSInv/code/make_GFs/Demo' #working dir
home=os.getcwd() #working dir
project_name='TestWaves1' #create a folder under working dir
fk_home='/Users/timlin/Documents/Project/GPSInv/code/make_GFs/fk' #where is your fk dir, careful the permission

grdfile='test.grid' #source location
srcfile='test.src' #source file
stafile='test.sta' #recording stations [lon, lat, xxxx,xxxx,xxx] only the first 1,2 columns are taken in calculation

modelpath='/Users/timlin/Documents/Project/GPSInv/GFs/TW_dense/vel1D_CWB' #full path (or name) of the velocity model
nprocess=8
SRCtype=1 #type=1 (Mw/strike/dip/rake) or type=2 (m0/M1~M6) or type3 (generate Green's function)
#########################default parameters END#############################


def distaz(evlo,evla,STAs):
    '''
    calculate dist(km) and az
    STAs must to be a matrix with [n stations by 2(stlon,stlat)]
    '''
    #If STAs is only one station
    if STAs.ndim==1:
        dist,az,baz=obspy.geodetics.base.gps2dist_azimuth(lat1=evla,lon1=evlo,lat2=STAs[1],lon2=STAs[0])
        return([dist],[az])
    else:
        #More than 1 station
        all_dist=[]
        all_az=[]
        for sta in STAs:
            stlo=sta[0]
            stla=sta[1]
            dist,az,baz=obspy.geodetics.base.gps2dist_azimuth(lat1=evla,lon1=evlo,lat2=stla,lon2=stlo)
            all_dist.append(dist*0.001) #km
            all_az.append(az)
        return(all_dist,all_az)


def run_FK_static(ng,evlo,evla,evdp,SRC,STA,home,project_name,fk_home,modelpath,SRCtype):
    '''
    run FK static displacement
    ng:#grid nodes
    SRC:source file Mw/strike/dip/rake, is a matrix with same length of GRD. or it can be M1~M6
    STA:matrix of station coordinate
    '''
    #time.sleep(10)
    print(evlo,evla,STA)
    dists,azs=distaz(evlo,evla,STA)
    rndname=datetime.datetime.utcnow().strftime('%s%f')+str(random.random()) #make dir name based on time and random number
    tmpworkdir=home+'/fk_'+rndname
    shutil.copytree(fk_home,tmpworkdir) #copy the fk and working in that new folder
    shutil.copy(modelpath,tmpworkdir) #copy the velocity model
    os.chdir(tmpworkdir)
    #loop all the dist and az and run 1.fk.pl and 2.syn
    print('Length of dist')
    print('Length of dist=',len(dists))

    for nsta,dist in enumerate(dists):
        az=azs[nsta]
        dist_str=str(dist)
        az_str=str(az)
        run_command='./fk.pl '+'-M'+modelpath.split('/')[-1]+'/'+str(evdp)+'/f '+'-N1 '+dist_str+' >tmp.ou'
        #print('command:',run_command)
        S=subprocess.Popen(run_command,stdout = subprocess.PIPE,stderr = subprocess.STDOUT,shell=True)
        S.communicate()
        if SRCtype==1:
            run_command2='cat tmp.ou | ./syn' +' -M'+str(SRC[ng][0])+'/'+str(SRC[ng][1])+'/'+str(SRC[ng][2])+'/'+str(SRC[ng][3])+' -A'+az_str+' -P'+' >>tmp.zrt'
        elif SRCtype==2:
            run_command2='cat tmp.ou | ./syn' +' -M'+str(SRC[ng][0])+'/'+str(SRC[ng][1])+'/'+str(SRC[ng][2])+'/'+str(SRC[ng][3])+'/'+str(SRC[ng][4])+'/'+str(SRC[ng][5])+'/'+str(SRC[ng][6])+' -A'+az_str+' -P'+' >>tmp.zrt'
        elif SRCtype==3:
            run_command2='cat tmp.ou | ./syn' +' -M'+'1e20'+'/1/0/0/0/0/0'+' -A'+az_str+' -P'+' >tmp1.zrt'
            S2=subprocess.Popen(run_command2,stdout = subprocess.PIPE,stderr = subprocess.STDOUT,shell=True)
            S2.communicate()
            run_command2='cat tmp.ou | ./syn' +' -M'+'1e20'+'/0/1/0/0/0/0'+' -A'+az_str+' -P'+' | awk \'{print($3,$4,$5)}\' >tmp2.zrt'
            S2=subprocess.Popen(run_command2,stdout = subprocess.PIPE,stderr = subprocess.STDOUT,shell=True)
            S2.communicate()
            run_command2='cat tmp.ou | ./syn' +' -M'+'1e20'+'/0/0/1/0/0/0'+' -A'+az_str+' -P'+' | awk \'{print($3,$4,$5)}\' >tmp3.zrt'
            S2=subprocess.Popen(run_command2,stdout = subprocess.PIPE,stderr = subprocess.STDOUT,shell=True)
            S2.communicate()
            run_command2='cat tmp.ou | ./syn' +' -M'+'1e20'+'/0/0/0/1/0/0'+' -A'+az_str+' -P'+' | awk \'{print($3,$4,$5)}\' >tmp4.zrt'
            S2=subprocess.Popen(run_command2,stdout = subprocess.PIPE,stderr = subprocess.STDOUT,shell=True)
            S2.communicate()
            run_command2='cat tmp.ou | ./syn' +' -M'+'1e20'+'/0/0/0/0/1/0'+' -A'+az_str+' -P'+' | awk \'{print($3,$4,$5)}\' >tmp5.zrt'
            S2=subprocess.Popen(run_command2,stdout = subprocess.PIPE,stderr = subprocess.STDOUT,shell=True)
            S2.communicate()
            run_command2='cat tmp.ou | ./syn' +' -M'+'1e20'+'/0/0/0/0/0/1'+' -A'+az_str+' -P'+' | awk \'{print($3,$4,$5)}\' >tmp6.zrt'
            S2=subprocess.Popen(run_command2,stdout = subprocess.PIPE,stderr = subprocess.STDOUT,shell=True)
            S2.communicate()
        #print('command:',run_command2)
        if SRCtype!=3:
            S2=subprocess.Popen(run_command2,stdout = subprocess.PIPE,stderr = subprocess.STDOUT,shell=True)
            S2.communicate() #stdout, stderr=S2.communicate()
        elif SRCtype==3:
            run_command2='paste tmp1.zrt tmp2.zrt tmp3.zrt tmp4.zrt tmp5.zrt tmp6.zrt >>tmp.zrt'
            S2=subprocess.Popen(run_command2,stdout = subprocess.PIPE,stderr = subprocess.STDOUT,shell=True)
            S2.communicate()
    shutil.move('tmp.zrt',home+'/'+project_name+'/'+'SRC%05d.zrt'%(ng+1))
    os.chdir(home)
    shutil.rmtree(tmpworkdir)



#Main control
def run_forward():
    if not(os.path.isdir(project_name)):
        os.mkdir(project_name)

    GRD=open(grdfile,'r')
    if SRCtype==3:
        SRC='' #generate GFs only
    else:
        SRC=np.genfromtxt(srcfile)
        if SRC.ndim==1:
            SRC=[SRC]

    STA=np.genfromtxt(stafile)
    ng=0 # number of grid nodes
    ps=[]
    for grdline in GRD.readlines():
        grds=grdline.strip().split()
        evlo=float(grds[0]); evla=float(grds[1]); evdp=float(grds[2])
        #run FK main program one by one
        #run_FK_static(ng,evlo,evla,SRC,STA,home,project_name,fk_home,modelpath)
        #####parallel computation carefully use! #####
        #print('running#',ng,evlo,evla)
        p = mp.Process(target=run_FK_static, args=([ng,evlo,evla,evdp,SRC,STA,home,project_name,fk_home,modelpath,SRCtype]) )
        ps.append(p)
        ng+=1
    GRD.close()
    #start running
    queue=np.zeros(len(ps))
    for i in range(len(ps)):
        if (i%10==0):
            print('now at',i,'out of',len(ps))
        ps[i].start()
        queue[i]=1
        while True:
        #check running status
            running_idx=np.where(queue==1)[0]
            for ri in running_idx:
                if not(ps[ri].is_alive()):
                    #finnish running, close
                    ps[ri].join()
                    #ps[ri]=0 #.join() still has issue...
                    queue[ri]=0
                else:
                    continue
            if len(np.where(queue==1)[0])<=nprocess:
                #add a process
                #print('number of processer=',len(np.where(queue==1)[0]),'add a process,now at',i)
                break
            else:
                #print('number of queue reaches max:',nprocess,'try again later,now at',i)
                time.sleep(0.5) #wait and try again later
    #Final check if all the processes are done
    while True:
        if np.sum([nps.is_alive() for nps in ps ])==0:
            print('Is alive?',[nps.is_alive() for nps in ps ])
            break
        else:
            time.sleep(1)


# Test code
if __name__ == "__main__":
    print('Usage: set','make_forward_wave.[home/project_name/fk_home/grdfile/srcfile/stafile/modelpath/nprocess/SRCtype]')
    print('Run model:','make_forward_wave.run_forward()')











