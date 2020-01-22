# -*- coding: utf-8 -*-
"""
Created on Tue Jun 05 00:09:04 2019
    
@author: TimLin
"""

import datetime
import numpy as np

#determine whether p1-p2 is a pair
##parameter setting##
min_nsta=1 #at least x station have measurement
min_nsta_HiCC=5  #at least x stations has CC greater than the below threshold
min_CC=0.9
time_sep=86400*2 #p1-p2 needs to be separated by at least x sec
#####################

INfile=open('pairs_BP0.8-2_wind30s_new_fast.out','r')
OUTfile=open('pairs_BP0.8-2_wind30s_new_rob.seq','w')



def det_repeq(line,min_nsta,min_nsta_HiCC,min_CC,time_sep)->'Boolean':
    elems=line.split()
    p1p2=elems[0]
    t1=datetime.datetime.strptime(p1p2.split('-')[0],'%Y%m%d%H%M%S')
    t2=datetime.datetime.strptime(p1p2.split('-')[1],'%Y%m%d%H%M%S')
    if (t2-t1).total_seconds()<time_sep:
        return False,None,None
    count_CC=0
    count_sta=0
    for n_measu in range(int((len(elems) - 1) / 2)):
        sta=elems[1 + n_measu * 2]
        CC=float(elems[2 + n_measu * 2])
        count_sta += 1
        if CC>=min_CC:
            count_CC += 1
    if (count_sta>=min_nsta) and (count_CC>=min_nsta_HiCC):
        return True,p1p2.split('-')[0],p1p2.split('-')[1]
    else:
        return False,None,None

EQseq=[]
nseq=0 #number of sequences
for line in INfile.readlines():
    isrepEQ,p1,p2 = det_repeq(line,min_nsta,min_nsta_HiCC,min_CC,time_sep)
    if isrepEQ:
        nseq += 1
        comp_p1p2={allEQ:i for i,subset in enumerate(EQseq) for allEQ in subset} #dict of {'eqid':num of seq}
        #p1 exist, p2 not
        if (p1 in comp_p1p2):
            #check if p2 exist, if not, append p2 also
            if not(p2 in comp_p1p2):
                EQseq[comp_p1p2[p1]].append(p2)
                continue
            else:
                #both p1,p2 exist
                continue
    
        #p2 exist, p1 not
        if (p2 in comp_p1p2):
            #check if p1 exist, if not, append p1 also
            if not(p1 in comp_p1p2):
                EQseq[comp_p1p2[p2]].append(p1)
                continue
        
        #p1, p2 does not exist
        if (not(p1 in comp_p1p2)) and (not(p2 in comp_p1p2)):
            EQseq.append([p1,p2])
            continue

INfile.close()


#write EQ sequences to file
for nseq in EQseq:
    for neq in nseq:
        OUTfile.write('%s '%(neq))
    OUTfile.write('\n')


OUTfile.close()












