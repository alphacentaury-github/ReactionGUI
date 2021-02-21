# -*- coding: utf-8 -*-
"""
Created on Feb 11 2021

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

from PyQt5.QtWidgets import (QApplication,QMainWindow,
                             QDialog,QFileDialog,QWidget,
                             QComboBox,QLabel,QLineEdit,QCheckBox,
                             QMenu,QMenuBar,QDialogButtonBox,
                             QHBoxLayout,QVBoxLayout,QGridLayout,
                             QStackedLayout,
                             QGroupBox,QToolBox,QTabWidget,
                             QPushButton,QTextBrowser,
                             QSpacerItem,
                             QRadioButton)
from PyQt5 import uic
#from PyQt5.QtCore import *

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
        FigureCanvasQTAgg as FigureCanvas,
        NavigationToolbar2QT as NavigationToolbar)
import matplotlib.pyplot as plt

from subprocess import (call,Popen)
import numpy as np

import qt_myutil
from qt_myutil import (combined_Widgets_horizontal,combined_Widgets_vertical,
                       combined_Widgets_grid,QLabel_aligned,
                       text_Browser,WidgetMatplot,about_Dialog)

#import momdis
import momdis.momdis as momdis

form_omp = uic.loadUiType("qt_momdis.ui")[0]

element_names = ["n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl",
		 "Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
		 "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb",
		 "Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er",
		 "Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At",
		 "Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No",
		 "Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og",
		 "119","120","121","122","123","124","125","126","127","128","129","130"]


class density_input(QWidget,):
    """
    Density are parametrized with parameters p0,p1,p2,... as        
    * If IDPROJ=1, projectile density is (1+p1*(r/p0)^p2)*exp(-r/p0) 
    * If IDPROJ=2, projectile density is r^p2*exp(-p0*r) 
    * If IDPROJ=3, projectile density is (1+p2*(r/p1)^p3)/(1+exp((r-p0)/p1) 
    * If IDPROJ=4, projectile density is a (1/(1+exp((r-p0)/p1))^p2, 
    *                                    p3=0,1 to use derivative 
    * If IDPROJ=5, projectile density is a calculated from liquid-drop model             
    """
    def __init__(self,type_index=0,default_para=None):
        super().__init__() 
        self.layout = QVBoxLayout() 
        self.setLayout(self.layout)
        self.para_values = default_para 
        self.items = ['Gaussian','Yukawa','WS','WS power']        
        self.list_label = ['(1+p1*(r/p0)^p2)*exp(-r^2/p0^2)',
                           'r^p2*exp(-p0*r)/(p0*r)',
                           '(1+p2*(r/p1)^p3)/(1+exp((r-p0)/p1)',
                           '1/(1+exp((r-p0)/p1))^p2'
                           ] 
        self.label_form = QLabel(self.list_label[type_index])
        self.lineEdit_p0 = QLineEdit('0.0')
        self.lineEdit_p1 = QLineEdit('0.0')
        self.lineEdit_p2 = QLineEdit('0.0')
        self.lineEdit_p3 = QLineEdit('0.0')
        self.layout.addWidget(self.label_form)
        
        if default_para: 
            self.put_values(default_para)
            
        self.para_widget = combined_Widgets_grid(   
             [ [QLabel('p0'),self.lineEdit_p0],
               [QLabel('p1'),self.lineEdit_p1],
               [QLabel('p2'),self.lineEdit_p2],
               [QLabel('p3'),self.lineEdit_p3] ]) 
        self.layout.addWidget(self.para_widget )
                    
    def change_type(self,type_index):   
        self.label_form.setText(self.list_label[type_index])
                
    def reset_para(self,):
        self.lineEdit_p0.setText('0.0')
        self.lineEdit_p1.setText('0.0')
        self.lineEdit_p2.setText('0.0')
        self.lineEdit_p3.setText('0.0')
            
    def get_values(self,):
        value_list=[]
        value_list.append(float(self.lineEdit_p0.text()))
        value_list.append(float(self.lineEdit_p1.text()))
        value_list.append(float(self.lineEdit_p2.text()))
        value_list.append(float(self.lineEdit_p3.text()))
        return value_list  
    
    def put_values(self,value_list):
        #--opposite of read_input
        self.lineEdit_p0.setText('{}'.format(value_list[0]))
        self.lineEdit_p1.setText('{}'.format(value_list[1]))
        self.lineEdit_p2.setText('{}'.format(value_list[2]))
        self.lineEdit_p3.setText('{}'.format(value_list[3]))
        return 

class momdis_GUI(QDialog,form_omp):
    def __init__(self,zp=6,ap=15,zt=4,at=9,zc=6,ac=14,eca=103.0,s2c=1.0,
                 momdis_directory='momdis/',
                 default_values=True):
        super().__init__()
        self.setupUi(self)
        self.momdis_directory = momdis_directory 
        
        self.widg_dens_n = density_input(0,[0.7,0.0,0.0,0.0])
        self.widg_dens_c = density_input(0,[1.73,1.38,1.0,0.0])
        self.widg_dens_t = density_input(0,[1.93,0.0,0.0,0.0])
        self.gridLayout_trhorho.addWidget(self.widg_dens_n,3,1 )
        self.gridLayout_trhorho.addWidget(self.widg_dens_c,3,3 )
        self.gridLayout_trhorho.addWidget(self.widg_dens_t,3,5 )
        self.comboBox_den_n.addItems(self.widg_dens_n.items)
        self.comboBox_den_c.addItems(self.widg_dens_c.items)
        self.comboBox_den_t.addItems(self.widg_dens_t.items)
        self.comboBox_den_n.currentIndexChanged.connect(self.dens_n_change)
        self.comboBox_den_c.currentIndexChanged.connect(self.dens_c_change)
        self.comboBox_den_t.currentIndexChanged.connect(self.dens_t_change)
        
        self.list_task = ['d sigma / d p_z', 
                          'd^2 sigma / d p_t^2',
                          'd sigma / d p_y',
                          'd^2 sigma / d p_z d p_t' ,
                          'bound w.f.',
                          'S-matrix for nucleon-target',
                          'S-matrix for core-target'
                          ]
        self.comboBox_sigma_type.addItems(self.list_task)
        
        self.widg_external = qt_myutil.Widget_external(
             momdis_directory+"momdis.exe",[],
             working_directory=momdis_directory,
             before_run=self.before_momdis,
             after_run=self.after_momdis)
        
        self.verticalLayout.addWidget(self.widg_external)  
                
        self.btn_plot = QPushButton('Plot')
        self.btn_plot.clicked.connect(self.plot_momdis)
        self.verticalLayout.addWidget(self.btn_plot)
        
        self.data = None 
        
    def dens_n_change(self,):
        self.widg_dens_n.change_type(self.comboBox_den_n.currentIndex())
    def dens_c_change(self,):
        self.widg_dens_c.change_type(self.comboBox_den_c.currentIndex())
    def dens_t_change(self,):
        self.widg_dens_t.change_type(self.comboBox_den_t.currentIndex())

    def put_values(self,input_dict):
        self.lineEdit_ZP.setText('{}'.format(input_dict['ZP']))
        self.lineEdit_AP.setText('{}'.format(input_dict['AP']))
        self.lineEdit_ZT.setText('{}'.format(input_dict['ZT']))
        self.lineEdit_AT.setText('{}'.format(input_dict['AT']))
        self.lineEdit_ECA.setText('{}'.format(input_dict['ECA']))
        if input_dict['ZC']==input_dict['ZP']:
            self.comboBox_nucleon.setCurrentIndex(0)
        elif input_dict['ZC']==input_dict['ZP']-1:
            self.comboBox_nucleon.setCurrentIndex(1)
        self.lineEdit_S2C.setText('{}'.format(input_dict['S2C']))
        self.lineEdit_N0.setText('{}'.format(input_dict['N0']))
        self.lineEdit_L0.setText('{}'.format(input_dict['L0']))
        self.lineEdit_J0.setText('{}'.format(input_dict['J0']))
        self.lineEdit_V0.setText('{}'.format(input_dict['V0']))
        self.lineEdit_R0.setText('{}'.format(input_dict['R0']))
        self.lineEdit_AA.setText('{}'.format(input_dict['AA']))
        self.lineEdit_VS0.setText('{}'.format(input_dict['VS0']))
        self.lineEdit_RS0.setText('{}'.format(input_dict['RS0']))
        self.lineEdit_AAS.setText('{}'.format(input_dict['AAS']))
        self.lineEdit_RC.setText('{}'.format(input_dict['RC']))
        self.comboBox_IPAULIn.setCurrentIndex(input_dict['IPAULIn'])
        self.comboBox_IPAULIc.setCurrentIndex(input_dict['IPAULIc'])
        self.comboBox_den_n.setCurrentIndex(input_dict['den_n'])
        self.comboBox_den_c.setCurrentIndex(input_dict['den_c'])
        self.comboBox_den_t.setCurrentIndex(input_dict['den_t'])
        self.widg_dens_n.put_values(input_dict['den_n_para'])
        self.widg_dens_c.put_values(input_dict['den_c_para'])
        self.widg_dens_t.put_values(input_dict['den_t_para'])
        self.comboBox_sigma_type.setCurrentIndex(input_dict['sigma_type'])
        self.lineEdit_sig_L0.setText('{}'.format(input_dict['sig_L0']))
        self.lineEdit_sig_M0.setText('{}'.format(input_dict['sig_M0']))
        return
    
    def get_values(self,):
        output={} 
        ZP = int(self.lineEdit_ZP.text()) 
        AP = int(self.lineEdit_AP.text())
        ZT = int(self.lineEdit_ZT.text()) 
        AT = int(self.lineEdit_AT.text())
        ECA = float(self.lineEdit_ECA.text()) 
        if self.comboBox_nucleon.currentIndex()==0: #neutron
            ZC = ZP
            AC = AP-1
        elif self.comboBox_nucleon.currentIndex()==1: #proton 
            ZC = ZP-1
            AC = AP-1 
        S2C = float(self.lineEdit_S2C.text())
        N0 = int(self.lineEdit_N0.text())
        L0 = int(self.lineEdit_L0.text())
        J0 = float(self.lineEdit_J0.text())
        V0 = float(self.lineEdit_V0.text())
        R0 = float(self.lineEdit_R0.text())
        AA = float(self.lineEdit_AA.text())
        VS0 = float(self.lineEdit_VS0.text())
        RS0 = float(self.lineEdit_RS0.text())
        AAS = float(self.lineEdit_AAS.text())
        RC = float(self.lineEdit_RC.text())
        IPAULIn  = self.comboBox_IPAULIn.currentIndex() 
        IPAULIc  = self.comboBox_IPAULIc.currentIndex()
        den_n  = self.comboBox_den_n.currentIndex()
        den_c  = self.comboBox_den_c.currentIndex()
        den_t  = self.comboBox_den_t.currentIndex()
        den_n_para = self.widg_dens_n.get_values()
        den_c_para = self.widg_dens_c.get_values()
        den_t_para = self.widg_dens_t.get_values()        
        sigma_type = self.comboBox_sigma_type.currentIndex()
        sig_L0 =  int(self.lineEdit_sig_L0.text())
        sig_M0 =  int(self.lineEdit_sig_M0.text()) 
        output = {'ZP':ZP,'AP':AP,'ZT':ZT,'AT':AT,'ZC':ZC,'AC':AC,
                  'ECA':ECA,'S2C':S2C,
                  'N0':N0,'L0':L0,'J0':J0,
                  'V0':V0,'R0':R0,'AA':AA,
                  'VS0':VS0,'RS0':RS0,'AAS':AAS,
                  'RC':RC,
                  'IPAULIn':IPAULIn,'IPAULIc':IPAULIc,
                  'den_n':den_n,'den_n_para':den_n_para,
                  'den_c':den_c,'den_c_para':den_c_para,
                  'den_t':den_t,'den_t_para':den_t_para,
                  'sigma_type':sigma_type,
                  'sig_L0':sig_L0,'sig_M0':sig_M0}
        return output
    
    def before_momdis(self,): # here one prepare inputs for momdis
        data = self.get_values() 
        self.data = data 
        # prepare input texts for momdis 
        # interpret input_dict 
        base = momdis.prepare_base_input(zp=data['ZP'],ap=data['AP'],
                                    zt=data['ZT'],at=data['AT'],
                                    zc=data['ZC'],ac=data['AC'],
                                    eca=data['ECA'],s2c=data['S2C'])
        bound = momdis.prepare_bound_input(base_input=base,
                        N0=data['N0'],L0=data['L0'],J0=data['J0'],
                        V0=data['V0'],R0=data['R0'],AA=data['AA'],
                        VS0=data['VS0'],RS0=data['RS0'],AAS=data['AAS'],
                        RC=data['RC'])
        sn   = momdis.prepare_Smat_sn_input(base_input=base,
                         ISMAT=0,IPOT=2,IPAULI=data['IPAULIn'],
                         IDPROJ=data['den_n']+1,
                         DPROJ_para=data['den_n_para'],
                         IDTARG=data['den_t']+1,
                         DTARG_para=data['den_t_para'])
        sc   = momdis.prepare_Smat_sc_input(base_input=base,
                         ISMAT=1,IPOT=2,IPAULI=data['IPAULIc'],
                         IDPROJ=data['den_c']+1,
                         DPROJ_para=data['den_c_para'],
                         IDTARG=data['den_t']+1,
                         DTARG_para=data['den_t_para'])
        # prepare input files 
        ff=open(self.momdis_directory+'a_bound.txt','w')
        ff.write(bound)
        ff.close()
        ff=open(self.momdis_directory+'a_sn.txt','w')
        ff.write(sn)
        ff.close()
        ff=open(self.momdis_directory+'a_sc.txt','w')
        ff.write(sc)
        ff.close()
        ff=open(self.momdis_directory+'a_sigma.txt','w')
        ff.write(base)
        ff.close()
        # set sequence of run                 
        if data['sigma_type'] in [0,1,2,3]:
            self.widg_external.set_sequence(
              [['a_bound.txt','1','a_bound.out','9'],
               ['a_sn.txt','2','a_sn.out','9'],
               ['a_sc.txt','2','a_sc.out','9'],
               ['a_sigma.txt','3',
                '{}'.format(data['sigma_type']+1),
                ' {} {}'.format(data['sig_L0'],data['sig_M0']),
                'a_bound.out','a_sn.out','a_sc.out','a_sigma.out','9']
              ])
        elif data['sigma_type']==4: #bound w.f. 
            self.widg_external.set_sequence(
               ['a_bound.txt','1','a_bound.out','9'])
        elif data['sigma_type']==5: #'S-matrix for nucleon-target'
            self.widg_external.set_sequence(
               ['a_sn.txt','2','a_sn.out','9'])
        elif data['sigma_type']==6: #'S-matrix for nucleon-target'
            self.widg_external.set_sequence(
               ['a_sc.txt','2','a_sc.out','9'])       
        
        return
    def after_momdis(self,):
        return 
    def plot_momdis(self,): # here treat output of momdis
        sigma_type = self.comboBox_sigma_type.currentIndex()
        if sigma_type==0: #'d sigma / d p_z'
            out = np.loadtxt(self.momdis_directory+'a_sigma.out')
            plt.figure()
            plt.plot(out[:,0],out[:,1])
            plt.xlabel('p_z[MeV]')
            plt.ylabel('mb/MeV')
            plt.show() 
        elif sigma_type==1: #'d^2 sigma / d p_t^2'
            out = np.loadtxt(self.momdis_directory+'a_sigma.out')
            plt.figure()
            plt.plot(out[:,0],out[:,1])
            plt.xlabel('p_t[MeV]')
            plt.ylabel('mb/MeV^2')
            plt.show()     
        elif sigma_type==2: #'d sigma / d p_y'
            out = np.loadtxt(self.momdis_directory+'a_sigma.out')
            plt.figure()
            plt.plot(out[:,0],out[:,1])
            plt.xlabel('p_y[MeV]')
            plt.ylabel('mb/MeV')
            plt.show()     
        elif sigma_type==3: #'d^2 sigma / d p_z d p_t'
            out = np.loadtxt(self.momdis_directory+'a_sigma.out')
            # pz pt sigma 
            # plt.figure()
            #---------------------------------------------
            # to plot contour one have to re-arrange data 
            #---------------------------------------------
            #plt.contour([out[:,0],out[:,1]],out[:,2])
            #plt.xlabel('p_z[MeV]') 
            #plt.ylabel('p_t[MeV]')
            #plt.title('|d^2\sig,a/(dp_z dp_t)| [mb/MeV^2]')
            #plt.show()
            #--------------------------------------------
            # at the moment, it is not available
            #--------------------------------------------
            self.widg_external.message('====Contour Plot is not available yet=========\n'  )
            self.widg_external.message(np.array2string(out))
        elif sigma_type==4: #bound w.f.
            out = np.loadtxt(self.momdis_directory+'a_bound.out')
            plt.figure()
            plt.plot(out[:,0],out[:,1])
            plt.xlabel('r[fm]')
            plt.ylabel('u(r)=r*R(r)')
            plt.show()            
        elif sigma_type==5: #S-matrix for nucleon-target
            out = np.loadtxt(self.momdis_directory+'a_sn.out')
            plt.figure()
            plt.plot(out[:,0],out[:,1],label='Re S')
            plt.plot(out[:,0],out[:,2],label='Im S')
            plt.xlabel('b(impact parameter)')
            plt.legend() 
            plt.show()
        elif sigma_type==6: #S-matrix for core-target
            out = np.loadtxt(self.momdis_directory+'a_sc.out')
            plt.figure()
            plt.plot(out[:,0],out[:,1],label='Re S')
            plt.plot(out[:,0],out[:,2],label='Im S')
            plt.xlabel('b(impact parameter)')
            plt.legend()
            plt.show()            
        return 

#===================================================================================================
if __name__ == "__main__":
    def run_app():
        """
        launcher of Qt in Spyder
        to avoid error
        """
        if not QApplication.instance():
            app = QApplication(sys.argv)
        else:
            app = QApplication.instance()

        myWindow = momdis_GUI()
        myWindow.show()
        app.exec_()
        return myWindow
    m = run_app()
    #sys.exit(app.exec_()) # when run in console
    #app.exec_()           # when run in cmd
