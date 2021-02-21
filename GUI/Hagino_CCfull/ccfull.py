# -*- coding: utf-8 -*-
"""
Created on Wed May 13 13:00:38 2020

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

#---current path where this file resides 
try:
    here = os.path.dirname(os.path.realpath(__file__))
except:
    here = '.'
    
#omget_path = 'omget_RIPL3'
#---current path where this file resides
#try:
#    here = os.path.dirname(os.path.realpath(__file__))
#except:
#    here = '.'


def run_ccfull(ccfull_exe=here+'/ccfull.exe',text=None):
    # replace ccfull.inp with new text 
    # and run the ccfull code 
    if text:
        write_input(text)
    proc= Popen(ccfull_exe, shell=True,cwd=here)
    out, err = proc.communicate()
    proc.wait() 
    proc.terminate() 
    
def write_input(txt, ccfull_inp = here+'/ccfull.inp'):
    # write txt to input file 
    ff = open(ccfull_inp,'w' )
    ff.write(txt) 
    ff.close() 
            
def read_output(cross=here+'/cross.dat'):
    # read output of ccfull code 
    # For the moment only cross section     
    output = np.loadtxt(cross)
    return output 
    
def cleanup():
    # remove input/output files 
    list_of_files=['ccfull.inp','cross.dat','OUTPUT','spin.dat']
    for item in list_of_files:
        if os.path.exists(item):
            os.remove(item)
       
#===================================================================================================
if __name__ == "__main__":
    run_ccfull(text=None) 
    output = read_output()
    plt.semilogy(output[:,0],output[:,1])
    plt.show() 
     

    
