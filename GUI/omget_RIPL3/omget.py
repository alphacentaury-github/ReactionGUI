# -*- coding: utf-8 -*-
"""
Created on Wed May 13 13:00:38 2020

@author: Y.-H. song
"""
import sys
import os 
import subprocess 
from subprocess import (call,Popen)
import numpy as np
#---current path where this file resides 
try:
    here = os.path.dirname(os.path.realpath(__file__))
except:
    here = '.'
         
def global_omp_choose(ap,zp,at,zt,elab):
    """
    Returns
    -------
    omp_list : [(omp_id,omp_description), ... ]

    """
    if zp == 0 and ap == 1 :
        name_Targ = ' n '
    elif zp == 1 and ap == 1 :
        name_Targ = ' p '
    elif zp == 1 and ap == 2 :
        name_Targ = ' d '
    elif zp == 1 and ap == 3 :
        name_Targ = ' t '
    elif zp == 2 and ap == 3 :
        name_Targ = '3He'
    elif zp == 2 and ap == 4 :
        name_Targ = '4He'
    else :
        name_Targ = 'NONE'
#
    ompIndexFile = here+'/omp-index-RIPL-ID.txt'
    # from RIPL site
    # [0:4] Lib.No.     [6:9] Inc.Part.
    # [12:19] Model Type     [21:24] disp.     [26:29] rel.
    # [30:33] =< Z =< [34:37]     [38:41] =< A =< [42:45]     [47:51] =< E =< [52:57]
    # [60:63] Ref.No.     [65:80] 1st Author
    ff = open(ompIndexFile,'r')
    lines = ff.readlines()
    lines = list(map(lambda s: s.strip(), lines))
    omp_list=[]
    for line in lines:
        if line[6:9] == name_Targ :
            if int(line[30:33]) <= zt <= int(line[34:37]) and int(line[38:41]) <= at <= int(line[42:45]) :
                if float(line[47:51]) <= elab <= float(line[52:57]) :
#                    if line[12:28] == 'spher.   no   no' :
                    if line[12:23] == 'spher.   no' :
    #                  until now, only spherical & no dirsper
 #                       omp_list.append((line[0:4],line[0:80]))
#                        omp_list.append((line[0:4],line[6:9],line[30:33]))
                        omp_list.append((int(line[0:4]),line[6:9],int(line[30:33]),int(line[34:37]),int(line[38:41]),int(line[42:45]),float(line[47:51]),float(line[52:57]),line[65:80]))
    return omp_list

def global_omp_get(ap,zp,at,zt,elab,omp_id,omget_exe =here+"/omget.exe"):
    omp_para = {'ap': ap, 'at': at,
                'rc': 0,
                'V': [0,0,0], 'W': [0,0,0],
                'Vso': [0,0,0], 'Wso':[0,0,0],
                'Vd': [0,0,0], 'Wd':[0,0,0] }
    try:
        os.remove(here+'/ecis.inp')
        os.remove(here+'/ominput.inp' )
    except:
        print('inp files may be already deleted.')
    
    ff = open(here+'/ominput.inp','wt')
    ff.write('1\n')
    ff.write('%f\n' % elab)
    ff.write('%d %d %d -2\n' % (int(zt), int(at), int(omp_id)))
    ff.close()
    #print(omget_exe) 
    proc= Popen(omget_exe ,shell=True,cwd=here)
    out, err = proc.communicate()
    proc.wait() 
    proc.terminate() 
    #print(out)
    #print(err)
    #----------------------------------------------------------------
    # require : omget.f & om_retrieve.f
    #           gs-mass-sp.dat & om-parameter-u.dat
    # compile : gfortran -o omget -std=legacy om_retrieve.f omget.f
    #----------------------------------------------------------------
    try:
        ffr = open(here+'/ecis.inp','r')
        ffrlines = ffr.readlines()
        ffr.close()
        ompitem = []
        for i in range (8,15):
            ompimsi = ffrlines[i].split()
            ompitem = ompitem + ompimsi
        ompitem = [float(f) for f in ompitem]
        omp_para.update(ap = 0)
        omp_para.update(at = at)
        omp_para.update(rc = ompitem[18])
        omp_para.update(V = [ompitem[0],ompitem[1],ompitem[2]])
        omp_para.update(W = [ompitem[3],ompitem[4],ompitem[5]])
        omp_para.update(Vd = [ompitem[6],ompitem[7],ompitem[8]])
        omp_para.update(Wd = [ompitem[9],ompitem[10],ompitem[11]])
        omp_para.update(Vso = [ompitem[12],ompitem[13],ompitem[14]])
        omp_para.update(Wso = [ompitem[15],ompitem[16],ompitem[17]])
        return omp_para
    except: 
        print('No ecis.inp or Error while running omget')
        return 0 

    