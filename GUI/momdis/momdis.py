# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 2021

@author: Y.-H. song
!***********************************************************************
!     
!    Copyright (c) 2020, ReactionGUI GUI for reaction calculation
!                        Produced at Rare Isotope Science Project 
!                        Written by Young-Ho Song, yhsong@ibs.re.kr 
!                               and Ik-Jae Shin,  geniean@ibs.re.kr 
!                        All rights reserved.
!          
!    This file is part of ReactionGUI.
!
!    ReactionGUI is free software: you can redistribute it and/or modify it
!    under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!     
!    ReactionGUI is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!     
!    You should have received a copy of the GNU General Public License
!    along with FRESCO. If not, see <http://www.gnu.org/licenses/>.
!     
!    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
!    LICENSE
!     
!    Our Preamble Notice
!
!      A. This work was produced at the Rare Isotope Science Project 
!         of Institute for Basic Science, 
!         funded by Ministry of Science and ICT and 
!         by National Research Foundation of Korea (2013M7A1A1075764).
!      B. Neither the Rare Isotope Science Project nor any of their employees, 
!         makes any warranty, express or implied, or assumes any liability or
!         responsibility for the accuracy, completeness, or usefulness
!         of any information, apparatus, product, or process disclosed,
!         or represents that its use would not infringe privately-owned
!         rights.
!        
!***********************************************************************
"""

import sys
import os 
import subprocess 
from subprocess import (call,Popen)
import numpy as np
import matplotlib.pyplot as plt 
import time

#---current path where this file resides 
try:
    here = os.path.dirname(os.path.realpath(__file__))
except:
    here = '.'

def ps_start(executable_file,cwd=None):
    return subprocess.Popen(
        executable_file,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=cwd 
    )
def ps_read(process):
    return process.stdout.readline().decode("utf-8").strip()
def ps_write(process, message):
    process.stdin.write(f"{message.strip()}\n".encode("utf-8"))
    process.stdin.flush()
def ps_terminate(process):
    process.stdin.close()
    process.terminate()
    process.wait(timeout=0.2)
def ps_inout(process,list_input=[]):
    """
    sends inputs in the list to the process in sequence 
    """
    for i in list_input:
        print(ps_read(process))
        ps_write(process,i)     
    print(ps_read(process))    
    process.wait() 
    ps_terminate(process)    
    return     

def prepare_base_input(zp=6,ap=15,zt=4,at=9,zc=6,ac=14,eca=103.0,s2c=1.0):
    """
    create a base input file text 
    
    zp,ap,zt,at,zc,ac = charge and mass number of projectile,target,core 
    """    
    text =''
    text += ' {} {} {} {} {}\n'.format(zp,ap,zt,at,eca)
    text += ' {} {}\n'.format(zc,ac)
    text += ' {}\n'.format(s2c) 
    return text 

def prepare_bound_input(base_input=None,
                        N0=1,L0=0,J0=0.5,
                        V0=-61.849,R0=2.67,AA=0.6,
                        VS0=0.0,RS0=2.4,AAS=0.6,RC=2.67):
    """
    create input file text for bound state calculation 
    
    N0,L0,J0 = (n,l,j) of nucleon in bound state 
    others are binding potential parameters in WS form.  
    """    
    text = base_input 
    text += ' {} {} {}\n'.format(N0,J0,L0)
    text += ' {} {} {} {} {} {} {}\n'.format(V0,R0,AA,VS0,RS0,AAS,RC)
    return text 
    
def prepare_Smat_sn_input(base_input=None,ISMAT=0,IPOT=2,IPAULI=0,
                         IDPROJ=1,DPROJ_para=[0.7,0.0,0.0],
                         IDTARG=1,DTARG_para=[1.93,0.0,0.0]):
    """
    prepare input file text for nucleon-target Scattering 
    
    ISMAT=0 for nucleon-target
    IPOT =2 for t-rho-rho method of optical potential 
    IPAULI = 0 for no Pauli correction of nn cross sections
           = 1 for Pauli correction of nn cross sections     
    Density are parametrized with parameters p0,p1,p2,... as        
    * If IDPROJ=1, projectile density is (1+p1*(r/p0)^p2)*exp(-p0*r) 
    * If IDPROJ=2, projectile density is r^p2*exp(-p0*r) 
    * If IDPROJ=3, projectile density is (1+p2*(r/p1)^p3)/(1+exp((r-p0)/p1) 
    * If IDPROJ=4, projectile density is a (1/(1+exp((r-p0)/p1))^p2, 
    *                                    p3=0,1 to use derivative 
    * If IDPROJ=5, projectile density is a calculated from liquid-drop model
    * If IDPROJ=10, projectile density is from an input file
                                         p0=filename, p1= proton size correction   
    * The same applies for the target density       
    """ 
    text = base_input 
    text += ' {}\n'.format(ISMAT) 
    text += ' {}\n'.format(IPOT)
    text += ' {}\n'.format(IPAULI) 
    text += ' {}  {}\n'.format(IDPROJ,IDTARG)
    temp = len(DPROJ_para)*' {}'+'\n'
    text += temp.format(*DPROJ_para)
    temp = len(DTARG_para)*' {}'+'\n'
    text += temp.format(*DTARG_para)    
    return text 
    
def prepare_Smat_sc_input(base_input=None,ISMAT=1,IPOT=2,IPAULI=0,
                         IDPROJ=1,DPROJ_para=[1.73,1.38,1.0],
                         IDTARG=1,DTARG_para=[1.93,0.0,0.0]):
    """
    prepare input file text for nucleon-target Scattering 
    
    ISMAT=1 for core-target
    IPOT =2 for t-rho-rho method of optical potential 
    IPAULI = 0 for no Pauli correction of nn cross sections
           = 1 for Pauli correction of nn cross sections     
    Density are parametrized with parameters p0,p1,p2,... as        
    * If IDPROJ=1, projectile density is (1+p1*(r/p0)^p2)*exp(-p0*r) 
    * If IDPROJ=2, projectile density is r^p2*exp(-p0*r) 
    * If IDPROJ=3, projectile density is (1+p2*(r/p1)^p3)/(1+exp((r-p0)/p1) 
    * If IDPROJ=4, projectile density is a (1/(1+exp((r-p0)/p1))^p2, 
    *                                    p3=0,1 to use derivative 
    * If IDPROJ=5, projectile density is a calculated from liquid-drop model
    * If IDPROJ=10, projectile density is from an input file
                                         p0=filename, p1= proton size correction   
    * The same applies for the target density     
    
    Do we need to check number of parameters? 
    """ 
    text = base_input 
    text += ' {}\n'.format(ISMAT) 
    text += ' {}\n'.format(IPOT)
    text += ' {}\n'.format(IPAULI) 
    text += ' {}  {}\n'.format(IDPROJ,IDTARG)
    temp = len(DPROJ_para)*' {}'+'\n'
    text += temp.format(*DPROJ_para)
    temp = len(DTARG_para)*' {}'+'\n'
    text += temp.format(*DTARG_para)    
    return text 
        
def run_momdis(base_input=None,bound_input=None,sn_input=None,sc_input=None,
               type_sigma=1,L0=0,M0=0,
               inf_bound='a_bound.txt',outf_bound='a_bound.out',
               inf_sn='a_sn.txt',outf_sn='a_sn.out',
               inf_sc='a_sc.txt',outf_sc='a_sc.out',
               inf_sigma='a_sigma.txt',outf_sigma='a_sigma.out',
               exe=here+'/momdis.exe' ):
    """
      type_sigma =
        1 for d sigma / d p_z
        2 for d^2 sigma / d p_t^2
        3 for d sigma / d p_y
        4 for d^2 sigma / d p_z d p_t
        5 for core + target elastic cross section
        
      L0, M0 must match with bound state calculation   
    """    
    # step1 : bound state calculation 
    ff=open(inf_bound,'w')
    ff.write(bound_input)
    ff.close() 
    
    process = ps_start(exe)
    ps_inout(process,[inf_bound,'1',outf_bound,'9'])
        
    # step2 : nucleon-target scattering
    ff=open(inf_sn,'w') 
    ff.write(sn_input)
    ff.close() 
    process = ps_start(exe)
    ps_inout(process,[inf_sn,'2',outf_sn,'9'])
    
    # step3 : nucleon-target scattering
    ff=open(inf_sc,'w') 
    ff.write(sc_input)
    ff.close() 
    process = ps_start(exe)
    ps_inout(process,[inf_sc,'2',outf_sc,'9'])
    
    # step4 : sigma 
    ff=open(inf_sigma,'w') 
    ff.write(base_input)
    ff.close() 
    process = ps_start(exe)
    ps_inout(process,[inf_sigma,'3',str(type_sigma),
             ' {} {}'.format(L0,M0),
             outf_bound,outf_sn,outf_sc,outf_sigma,'9'])
    return 
             
def read_output(cross=here+'/a_sigma.out'):
    # read output of ccfull code 
    # For the moment only cross section     
    output = np.loadtxt(cross)
    return output 
    
def cleanup():
    # remove input/output files 
    list_of_files=['a_bound.txt','a_sn.txt','a_sc.txt','a_sigma.txt',
                   'a_bound.out','a_sn.out','a_sc.out','a_sigma.out',
                   'EIGEN.txt','SMAT.TXT','SMAT_CORE_FIT.TXT']
    for item in list_of_files:
        if os.path.exists(item):
            os.remove(item)

#===================================================================================================
if __name__ == "__main__":
    base = prepare_base_input(zp=6,ap=15,zt=4,at=9,zc=6,ac=14,eca=103.0,s2c=1.0)
    bound = prepare_bound_input(base_input=base,
                        N0=1,L0=0,J0=0.5,
                        V0=-61.849,R0=2.67,AA=0.6,
                        VS0=0.0,RS0=2.4,AAS=0.6,RC=2.67)
    sn   = prepare_Smat_sn_input(base_input=base,
                         ISMAT=0,IPOT=2,IPAULI=0,
                         IDPROJ=1,DPROJ_para=[0.7,0.0,0.0],
                         IDTARG=1,DTARG_para=[1.93,0.0,0.0])
    sc   = prepare_Smat_sc_input(base_input=base,
                         ISMAT=1,IPOT=2,IPAULI=0,
                         IDPROJ=1,DPROJ_para=[1.73,1.38,1.0],
                         IDTARG=1,DTARG_para=[1.93,0.0,0.0]) 
    run_momdis(exe=here+'/momdis.exe',   
               base_input=base,bound_input=bound,sn_input=sn,sc_input=sc,
               type_sigma=1,L0=0,M0=0,
               inf_bound='a_bound.txt',outf_bound='a_bound.out',
               inf_sn='a_sn.txt',outf_sn='a_sn.out',
               inf_sc='a_sc.txt',outf_sc='a_sc.out',
               inf_sigma='a_sigma.txt',outf_sigma='a_sigma.out'
               )
     
