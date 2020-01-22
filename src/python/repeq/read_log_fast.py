# -*- coding: utf-8 -*-
"""
Created on Wed May 08 18:19:05 2019
    
@author: TimLin
"""

import glob
import obspy
import subprocess
import os,sys


#eqpath='/Users/timlin/Documents/Project/EQrequest/Hawaii/Hawaii_ALL/' #where you put your waveform data (EQ folders)
#log_dir='/home/jtlin/RepEQ/Hawaii_repEQ'
log_dir='/Users/timlin/Documents/Project/EQrequest/Hawaii/reqEQ_log_reduced'
logs=glob.glob(log_dir+'/'+'*.log')
#logs=glob.glob(log_dir+'/'+'IU.POHA.log')
#eqids=glob.glob(eqpath+'*')
#eqids.sort()
#Cn/2 for all the possible combination
#OUT1=open('pairs.out','w')

#OUT1=open('pairs_BP0.8-2_wind30s_fast.out','w')
#OUT1=open('pairs_BP0.8-2_wind30s_TEST.out','w')
OUT1=open('pairs_BP0.8-2_wind30s_new_fast.out','w')


#####save logs into a very large dictionary#####
Large_D={}
#nsta=0
for logname in logs:
    net_sta=logname.split('/')[-1].split('.')[0]+'.'+logname.split('/')[-1].split('.')[1]
    print(net_sta)
    Large_D[net_sta]=[] #create new key
    IN1=open(logname,'r')
    for line in IN1.readlines():
        Large_D[net_sta].append(line.strip())
    IN1.close()
    #nsta=nsta+1
    #if nsta==5:
    #    break
    #break

#sys.exit(2)

#This is all the keys for net_staname
ref_keys=list(Large_D.keys())

sav_D={}

for nkey in range(len(ref_keys)):
    print('-Dealing with',ref_keys[nkey])
    total_lines=len(Large_D[ref_keys[nkey]])
    for nline,line in enumerate(Large_D[ref_keys[nkey]]):
        if nline%100==0:
            print('n_line=',nline,'out of',total_lines)
        p1p2=line.split()[0] # p1-p2
        CCC=float(line.split()[-1])
        if p1p2 in sav_D:
            sav_D[p1p2].append(ref_keys[nkey])
            sav_D[p1p2].append(CCC)
        else:
            sav_D[p1p2]=[]
            sav_D[p1p2].append(ref_keys[nkey])
            sav_D[p1p2].append(CCC)

#sort the order
all_keys=list(sav_D.keys())
all_keys.sort()


#Write sav_D to output
for all_key in all_keys:
    OUT1.write('%s '%(all_key))
    for i_elem in range( int(len(sav_D[all_key])/2) ):
        OUT1.write('%s %4.2f '%(sav_D[all_key][2*i_elem],sav_D[all_key][2*i_elem+1] ))
    OUT1.write('\n')

OUT1.close()


#for nkey in range(len(ref_keys)):
#    print('-Dealing with',ref_keys[nkey])
#    total_lines=len(Large_D[ref_keys[nkey]])
#    for nline,line in enumerate(Large_D[ref_keys[nkey]]):
#        if nline%100==0:
#            print('n_line=',nline,'out of',total_lines)
#        #print('---checking line:',line)
#        ref_p1p2=line.split()[0] #based on this p1p2 find others
#        CCC=float(line.split()[-1])
#        OUT1.write('%s '%(ref_p1p2))
#        #print('write:',ref_keys[nkey],CCC)
#        OUT1.write('%s %4.2f '%(ref_keys[nkey],CCC))
#        #Now searching the same p1p2
#        sav_net_sta=[];sav_CCC=[]
#        for tar_nkey in range(len(ref_keys)):
#            if tar_nkey<=nkey:
#                continue #no need to craw into the old or current dictionary
#            #print('   searching:',ref_keys[tar_nkey])
#            for tar_nline,tar_line in enumerate(Large_D[ref_keys[tar_nkey]]):
#                if ref_p1p2 in tar_line:
#                    #print('        Find!',tar_line)
#                    CCC=float(tar_line.split()[-1])
#                    sav_CCC.append(CCC)
#                    sav_net_sta.append(ref_keys[tar_nkey])
#                    #print('            delete:',ref_keys[tar_nkey],'line:',Large_D[ref_keys[tar_nkey]][tar_nline])
#                    del Large_D[ref_keys[tar_nkey]][tar_nline]
#                    break
#        if len(sav_CCC)!=0:
#            #OUT1.write('%s '%(ref_p1p2))
#            for n_write in range(len(sav_net_sta)):
#                #print('write:',sav_net_sta[n_write],sav_CCC[n_write])
#                OUT1.write('%s %4.2f '%(sav_net_sta[n_write],sav_CCC[n_write]))
#                #OUT1.write('\n') #change line at the final
#        OUT1.write('\n') #the ending line is necessary for every loops
#        '''
#        if '20110915150600-20111020041733' in line:
#            print(sav_net_sta)
#            print(sav_CCC)
#            OUT1.close()
#            sys.exit(2)
#        '''
#
#OUT1.close()



'''
for ref_nkey in Large_D.keys():
    for nline,line in enumerate(Large_D[ref_nkey]):
        ref_p1p2=line.split()[0] #based on this p1p2 find others
        CCC=float(line.split()[-1])
        for nkey in Large_D.keys():

'''








'''
################################################
#IN1=open('search_ID.dat','r')
#for p1p2 in IN1.readlines():
for p1p2 in reversed(open('search_ID.dat').readlines()):
    p1p2=p1p2.rstrip()
    #print('searching:',p1p2)
    sav_net_sta=[];sav_CCC=[] #initial the array for saving file later
    for nkey in Large_D.keys():
        for nline,line in enumerate(Large_D[nkey]):
            if p1p2 in line:
                CCC=float(line.split()[-1])
                #Large_D[nkey].remove(line) #already found the info you want, no need this line anymore.
                #print('ready to remove:',Large_D[nkey][nline])
                del Large_D[nkey][nline] #already found the info you want, no need this line anymore.
                sav_net_sta.append(nkey)
                sav_CCC.append(CCC)
                break
    if len(sav_net_sta)==0:
        continue
    OUT1.write('%s '%(p1p2)) #write the 1st (p1-p2) column
    for n_write in range(len(sav_net_sta)):
        OUT1.write('%s %4.2f '%(sav_net_sta[n_write],sav_CCC[n_write]))
    OUT1.write('\n') #change line at the final

OUT1.close()
'''

'''
for n1,ev1 in enumerate(eqids):
    for n2,ev2 in enumerate(eqids):
        if ( (ev1==ev2) or (n2<=n1)):
            continue
        #pairs
        p1=ev1.split('/')[-1][:] #YYYYMMDDHHmmss
        p2=ev2.split('/')[-1][:]
        print(p1,p2)
        #print('search %s-%s from Large_Dictionary'%(p1,p2))
        #loop through keys
        sav_net_sta=[];sav_CCC=[] #initial the array for saving file later
        for nkey in Large_D.keys():
            for nline,line in enumerate(Large_D[nkey]):
                if '%s-%s'%(p1,p2) in line:
                    CCC=float(line.split()[-1])
                    #Large_D[nkey].remove(line) #already found the info you want, no need this line anymore.
                    #print('ready to remove:',Large_D[nkey][nline])
                    del Large_D[nkey][nline] #already found the info you want, no need this line anymore.
                    sav_net_sta.append(nkey)
                    sav_CCC.append(CCC)
                    break
        OUT1.write('%s '%(p1+'-'+p2)) #write the 1st (p1-p2) column
        for n_write in range(len(sav_net_sta)):
            OUT1.write('%s %4.2f '%(sav_net_sta[n_write],sav_CCC[n_write]))
        OUT1.write('\n') #change line at the final

OUT1.close()
'''
