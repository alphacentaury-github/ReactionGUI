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

#import ccfull 
import Hagino_CCfull.ccfull as ccfull 

form_omp = uic.loadUiType("qt_ccfull.ui")[0]
#omget_path = 'omget_RIPL3'
#---current path where this file resides
#try:
#    here = os.path.dirname(os.path.realpath(__file__))
#except:
#    here = '.'


element_names = ["n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl",
		 "Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
		 "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb",
		 "Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er",
		 "Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At",
		 "Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No",
		 "Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og",
		 "119","120","121","122","123","124","125","126","127","128","129","130"]

def AW_heavy_ion_potential_para(zp,ap,zt,at):
    """
    Akyuz-Winther heavy ion potential 
    it has only real part in WS form 
    
    returns WS potential parameters

    typically a~ 0.63 fm 
    """
    Np = ap-zp
    Nt = at-zt 
    R_p = 1.20*ap**(1./3.)-0.09
    R_t = 1.20*at**(1./3.)-0.09
    Rbar = R_p*R_t/(R_p+R_t)
    R0 = R_p + R_t
    gamma = 0.95*(1.0-1.8*(Np-zp)/ap*(Nt-zt)/at    )
    a_inv = 1.17*(1.0+0.53*(ap**(-1./3.)+at**(-1./3.)))
    V0 = 16*np.pi*gamma*Rbar/a_inv  
    return (V0,R0,1.0/a_inv) # (V>0,R,a) of WS potential 


class excitation_input(QWidget,):
    """
    class to take input of excited states for CCFULL
    
    if type = 0: inert. no input required. 
    if type = 1: vibrational case 
              omega, beta, lambda, nphonon 
              
              if nuclei='target', additional input 
              omega_2, beta_2, lambda_2, nphonon_2 
              
    if type =2 : rotational case 
              e2, beta2, beta4, nrot             
    """
    def __init__(self,nuclei='target',type_index=0):
        super().__init__() 
        self.layout = QStackedLayout() 
        self.setLayout(self.layout)
        self.nuclei = nuclei 
        #---first page 
        self.widget_first_page = QLabel('Inert ({})'.format(nuclei))
        self.layout.addWidget(self.widget_first_page )
        #--second page 
        self.widget_second_page = QWidget() 
        self.layout.addWidget(self.widget_second_page )
        self.second_layout = QGridLayout()
        self.widget_second_page.setLayout(self.second_layout)
        
        self.second_layout.addWidget( 
              QLabel('Vibrational ({})'.format(nuclei)),0,0)
        self.vib_lambda = QLineEdit('3') 
        self.vib_omega = QLineEdit('1.81') 
        self.vib_beta = QLineEdit('0.205') 
        self.vib_nphonon = QLineEdit('1') 
        self.vib_add_mod = QCheckBox('add. mod.')
        self.vib_lambda2 = QLineEdit('0') 
        self.vib_omega2 = QLineEdit('0.0') 
        self.vib_beta2 = QLineEdit('0.0') 
        self.vib_nphonon2 = QLineEdit('0') 
        
        self.vib_main_widget = combined_Widgets_grid(   
             [ [QLabel('lambda'),self.vib_lambda,
                QLabel('hbar omega'),self.vib_omega],
               [QLabel('beta_lambda'),self.vib_beta,
                QLabel('number of phonons'),self.vib_nphonon]
              ] ) 
        self.second_layout.addWidget(self.vib_main_widget )
        
        if nuclei=='target': 
            self.vib_2nd_widget = combined_Widgets_grid(
                [[QCheckBox('additional mode')],
                 [QLabel('lambda'),self.vib_lambda2,
                  QLabel('hbar omega'),self.vib_omega2 ],
                 [QLabel('beta_lambda'),self.vib_beta2,
                QLabel('number of phonons'),self.vib_nphonon2  ]
                    ])
            self.second_layout.addWidget(self.vib_2nd_widget )
            
        #--third page 
        self.widget_third_page = QWidget() 
        self.layout.addWidget(self.widget_third_page )  
        self.third_layout = QGridLayout()
        self.widget_third_page.setLayout(self.third_layout)
        
        self.linEdit_E2 = QLineEdit() 
        self.linEdit_beta2 = QLineEdit() 
        self.linEdit_beta4 = QLineEdit() 
        self.linEdit_nrot = QLineEdit()
        self.rot_widget =  combined_Widgets_grid(
            [[QLabel('Rotational ({})'.format(nuclei))],
             [QLabel('energy(2+)'),self.linEdit_E2  ],
             [QLabel('beta_2'),self.linEdit_beta2 ],
             [QLabel('beta_4'),self.linEdit_beta4 ],
             [QLabel('number of levels'),self.linEdit_nrot ]             
             ]
            )        
        self.third_layout.addWidget(self.rot_widget )
        
        #---initial page 
        self.layout.setCurrentIndex(type_index) 
            
    def change_type(self,type_index):   
        self.layout.setCurrentIndex(type_index)  
    
    def read_input(self,):
        output = {} 
        type_index = self.layout.currentIndex()
        output['type'] = type_index 
        output['nuclei'] = self.nuclei
        if type_index==0 : # inert 
            pass
        elif type_index==1: # vibrational 
            output['vib'] = self.vib_main_widget.get_values()
            if self.nuclei=='target':
                output['vib2'] = self.vib_2nd_widget.get_values()
        elif type_index==2: #rotational 
            output['rot'] = self.rot_widget.get_values()
        return output 
    
    def set_values(self,input_dict):
        #--opposite of read_input 
        self.layout.setCurrentIndex(input_dict['type'] )
        type_index = input_dict['type']
        self.nuclei = input_dict['nuclei']
        if type_index ==1:
            self.vib_main_widget.put_values(input_dict['vib'])
            if self.nuclei=='target':
                self.vib_2nd_widget.put_values(input_dict['vib2'])
        elif type_index ==2:
            self.rot_widget.put_values(input_dict['rot'])
        return 

class CCFULL_GUI(QDialog,form_omp):
    def __init__(self,ap=16,zp=8,at=144,zt=62,enrange=(55.0,72.0,1.0),
                 ccfull_directory ='Hagino_CCfull/',
                 default_values=True):
        super().__init__()
        self.setupUi(self)
                
        self.widget_proj_ex = excitation_input('projectile')
        self.widget_proj_ex.change_type(self.comboBox_proj_ex.currentIndex()) 
        self.widget_targ_ex = excitation_input('target')
        self.widget_targ_ex.change_type(self.comboBox_targ_ex.currentIndex()) 
                                  
        self.gridLayout_coup.addWidget( self.widget_proj_ex,1,0)
        self.gridLayout_coup.addWidget( self.widget_targ_ex,1,1)
        self.pushButton_AWpot.clicked.connect(self.AWpot)
        
        self.widg_external = qt_myutil.Widget_external(
             ccfull_directory+"ccfull.exe",[],
             working_directory=ccfull_directory,
             before_run=self.before_ccfull,
             after_run=self.after_ccfull)
        
        self.verticalLayout.addWidget(self.widg_external)
        
        self.btn_plot = QPushButton('Plot')
        self.btn_plot.clicked.connect(self.plot_data)
        self.verticalLayout.addWidget(self.btn_plot)
        
        #--signal/slots 
        self.comboBox_proj_ex.currentIndexChanged.connect(self.proj_ex_change)
        self.comboBox_targ_ex.currentIndexChanged.connect(self.targ_ex_change)
        
        #---fill/empty values 
        if default_values==False:
            self.lineEdit_ap.setText('{}'.format(ap)) 
            self.lineEdit_zp.setText('{}'.format(zp))
            self.lineEdit_at.setText('{}'.format(at))
            self.lineEdit_zt.setText('{}'.format(zt))
            self.lineEdit_emin.clear() 
            self.lineEdit_emax.clear() 
            self.lineEdit_dE.clear() 
            self.lineEdit_V0.clear() 
            self.lineEdit_R0.clear() 
            self.lineEdit_a0.clear() 
            self.comboBox_proj_ex.setCurrentIndex(0)
            self.comboBox_targ_ex.setCurrentIndex(0)
            
        
       
    def AWpot(self,):
        ap = float(self.lineEdit_ap.text())
        zp = float(self.lineEdit_zp.text())
        at = float(self.lineEdit_at.text())
        zt = float(self.lineEdit_zt.text())
        (V0,R0,a0) = AW_heavy_ion_potential_para(zp,ap,zt,at)
        self.lineEdit_V0.setText('{:.4f}'.format(V0)) 
        self.lineEdit_R0.setText('{:.4f}'.format(R0)) 
        self.lineEdit_a0.setText('{:.4f}'.format(a0)) 
                
    def proj_ex_change(self,):
        type_index = self.comboBox_proj_ex.currentIndex() 
        self.widget_proj_ex.change_type(type_index) 
        return 
    
    def targ_ex_change(self,):
        type_index = self.comboBox_targ_ex.currentIndex() 
        self.widget_targ_ex.change_type(type_index) 
        return 
        
    def read_input(self,):
        #read inputs and return dictionary
        self.ap = float(self.lineEdit_ap.text())
        self.zp = float(self.lineEdit_zp.text())
        self.at = float(self.lineEdit_at.text())
        self.zt = float(self.lineEdit_zt.text())
        self.emin = float(self.lineEdit_emin.text())
        self.emax = float(self.lineEdit_emax.text())
        self.dE   = float(self.lineEdit_dE.text())
        self.R_p = float(self.lineEdit_Rp.text())
        self.R_t = float(self.lineEdit_Rt.text())
        #---read input for projectile excitation 
        self.proj_ex_type =  self.comboBox_proj_ex.currentIndex() 
        proj_ex_inputs = self.widget_proj_ex.read_input()
        #---read input for target excitation 
        self.targ_ex_type =  self.comboBox_targ_ex.currentIndex() 
        targ_ex_inputs = self.widget_targ_ex.read_input()
        #---read Nuclear Pot 
        self.V0 = float(self.lineEdit_V0.text()) 
        self.R0 = float(self.lineEdit_R0.text()) 
        self.a0 = float(self.lineEdit_a0.text()) 
        self.rmatch = float(self.lineEdit_rmatch.text()) 
        self.rstep = float(self.lineEdit_rstep.text()) 
        #---read transfer coup 
        self.ntrans = self.comboBox_ntrans.currentIndex() 
        self.qtrans = float(self.lineEdit_qtrans.text() )
        self.ftrans = float(self.lineEdit_ftrans.text() )
        output= {'AP': self.ap,'ZP':self.zp,
                 'AT': self.at,'ZT':self.zt,
                 'EMIN':self.emin,'EMAX':self.emax,'dE':self.dE,
                 'R_p': self.R_p,'R_t':self.R_t,
                 'proj.ex.': proj_ex_inputs ,
                 'targ.ex.': targ_ex_inputs ,
                 'V0':self.V0,'R0':self.R0,'a0':self.a0,
                 'RMATCH':self.rmatch,'STEP':self.rstep,
                 'NTRANS':self.ntrans,'Qtrans':self.qtrans,
                 'FTR':self.ftrans}
        return output 
    
    def set_values(self,data_dict):
        #opposite of read_input
        self.lineEdit_ap.setText('{}'.format(data_dict['AP']))
        self.lineEdit_zp.setText('{}'.format(data_dict['ZP']))
        self.lineEdit_at.setText('{}'.format(data_dict['AT']))
        self.lineEdit_zt.setText('{}'.format(data_dict['ZT']))
        self.lineEdit_emin.setText('{}'.format(data_dict['EMIN']))
        self.lineEdit_emax.setText('{}'.format(data_dict['EMAX']))
        self.lineEdit_dE.setText('{}'.format(data_dict['dE']))
        self.lineEdit_Rp.setText('{}'.format(data_dict['R_p']))
        self.lineEdit_Rt.setText('{}'.format(data_dict['R_t']))
        
        self.comboBox_proj_ex.setCurrentIndex(data_dict['proj.ex.']['type']) 
        self.widget_proj_ex.layout.setCurrentIndex(data_dict['proj.ex.']['type'])
        if data_dict['proj.ex.']['type']==1: #vibration 
            self.widget_proj_ex.set_values(data_dict['proj.ex.']['vib'])
        elif data_dict['proj.ex.']['type']==2:
            self.widget_proj_ex.set_values(data_dict['proj.ex.']['rot'])
            
        self.comboBox_targ_ex.setCurrentIndex(data_dict['targ.ex.']['type']) 
        self.widget_proj_ex.layout.setCurrentIndex(data_dict['targ.ex.']['type']) 
        if data_dict['targ.ex.']['type']==1: #vibration 
            self.widget_targ_ex.set_values(data_dict['targ.ex.']['vib'])
            self.widget_targ_ex2.set_values(data_dict['targ.ex.']['vib2'])
        elif data_dict['targ.ex.']['type']==2:
            self.widget_targ_ex.set_values(data_dict['targ.ex.']['rot'])
                    
        self.lineEdit_V0.setText('{}'.format(data_dict['V0']))
        self.lineEdit_R0.setText('{}'.format(data_dict['R0']))
        self.lineEdit_a0.setText('{}'.format(data_dict['a0']))        
        self.lineEdit_rmatch.setText('{}'.format(data_dict['RMATCH']))
        self.lineEdit_rstep.setText('{}'.format(data_dict['STEP']))
        self.comboBox_ntrans.setCurrentIndex(data_dict['NTRANS']) 
        self.lineEdit_qtrans.setText('{}'.format(data_dict['Qtrans']))
        self.lineEdit_ftrans.setText('{}'.format(data_dict['FTR']))
                  
        return 
    
    def before_ccfull(self,): 
        # cleanup previous run 
        ccfull.cleanup() 
        # prepare input file
        data_dict = self.read_input()
        out_txt = self.interpret_input_ccfull(data_dict) 
        ccfull.write_input(out_txt)
        return 
    
    def after_ccfull(self,):
        # read result and plot 
        output = ccfull.read_output() 
        return 
    
    def plot_data(self,):
        #--plot result from array 
        # using plt.semilogy seems to cause 
        # "QCoreApplication::exec: The event loop is already running"
        # maybe better to use canvas 
        array = ccfull.read_output()
        plt.figure()
        plt.semilogy(array[:,0],array[:,1])
        plt.xlabel('MeV')
        plt.ylabel('mb')
        plt.show() 
        return 
    
    def interpret_input_ccfull(self,data_dict):
        # interpret the input in forms of ccfull inputs 
        AP=data_dict['AP']
        ZP=data_dict['ZP']
        AT=data_dict['AT']
        ZT=data_dict['ZT']
        RP=data_dict['R_p']
        RT=data_dict['R_t']
        IVIBROTP = data_dict['proj.ex.']['type']-1 # -1,0,1 
        IVIBROTT = data_dict['targ.ex.']['type']-1 # -1,0,1 
        
        OMEGAT = 0.0
        BETAT  = 0.0
        LAMBDAT= 0
        NPHONONT= 0
        OMEGAT2 = 0.0
        BETAT2 = 0.0
        LAMBDAT2 = 0
        NPHONONT2 = 0
            
        if data_dict['targ.ex.']['type']==1 : # vibrational     
            OMEGAT = float(data_dict['targ.ex.']['vib'][(0,3)]) 
            BETAT = float(data_dict['targ.ex.']['vib'][(1,1)] )
            LAMBDAT = int(data_dict['targ.ex.']['vib'][(0,1)] )
            NPHONONT = int(data_dict['targ.ex.']['vib'][(1,3)]) 
            if data_dict['targ.ex.']['vib2'][(0,0)]: # additional mode 
                OMEGAT2 = float(data_dict['targ.ex.']['vib2'][(1,3)]) 
                BETAT2 = float(data_dict['targ.ex.']['vib2'][(2,1)] )
                LAMBDAT2 = int(data_dict['targ.ex.']['vib2'][(1,1)] )
                NPHONONT2 = int(data_dict['targ.ex.']['vib2'][(2,3)]) 
        elif data_dict['targ.ex.']['type']==2 : #rotational 
            E2T = float(data_dict['targ.ex.']['rot'][(1,1)])
            BETA2T = float(data_dict['targ.ex.']['rot'][(2,1)])
            BETA4T =float(data_dict['targ.ex.']['rot'][(3,1)])
            NROTT = int(data_dict['targ.ex.']['rot'][(4,1)])
            
        OMEGAP = 0.0
        BETAP  = 0.0
        LAMBDAP= 0
        NPHONONP= 0
        if data_dict['proj.ex.']['type']==1 : # vibrational     
            OMEGAP = float(data_dict['proj.ex.']['vib'][(0,3)]) 
            BETAP = float(data_dict['proj.ex.']['vib'][(1,1)] )
            LAMBDAP = int(data_dict['proj.ex.']['vib'][(0,1)] )
            NPHONONP = int(data_dict['proj.ex.']['vib'][(1,3)]) 
        elif data_dict['proj.ex.']['type']==2 : #rotational 
            E2P = float(data_dict['proj.ex.']['rot'][(1,1)])
            BETA2P = float(data_dict['proj.ex.']['rot'][(2,1)])
            BETA4P =float(data_dict['proj.ex.']['rot'][(3,1)])
            NROTP = int(data_dict['proj.ex.']['rot'][(4,1)])
                
        NTRANS = int(data_dict['NTRANS'] )
        QTRANS = float(data_dict['Qtrans'] )
        FTR = float(data_dict['FTR'])
        V0 = float(data_dict['V0'])
        R0 = float(data_dict['R0'])
        A0 = float(data_dict['a0'])
        EMIN = float(data_dict['EMIN']) 
        EMAX = float(data_dict['EMAX'])
        DE = float(data_dict['dE'])
        RMAX = float(data_dict['RMATCH'])
        DR = float(data_dict['STEP'])
        
        #------write text input form 
        txt ='{},{},{},{}\n'.format(AP,ZP,AT,ZT)
        txt = txt +'{},{},{},{}\n'.format(RP,IVIBROTP,RT,IVIBROTT)
        if IVIBROTT in [-1,0]:
            txt = txt +'{},{},{},{}\n'.format(OMEGAT,BETAT,LAMBDAT,NPHONONT)
        elif IVIBROTT==1:
            txt = txt +'{},{},{},{}\n'.format(E2T,BETA2T,BETA4T,NROTT)
        txt = txt +'{},{},{},{}\n'.format(OMEGAT2,BETAT2,LAMBDAT2,NPHONONT2)     
        if IVIBROTP in [-1,0]:
            txt = txt +'{},{},{},{}\n'.format(OMEGAP,BETAP,LAMBDAP,NPHONONP)
        elif IVIBROTP==1:
            txt = txt +'{},{},{},{}\n'.format(E2P,BETA2P,BETA4P,NROTP)
        txt = txt +'{},{},{}\n'.format(NTRANS,QTRANS,FTR)
        txt = txt +'{},{},{}\n'.format(V0,R0,A0)
        txt = txt +'{},{},{}\n'.format(EMIN,EMAX,DE)
        txt = txt +'{},{}\n'.format(RMAX,DR)
        return txt 
            
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
        
        myWindow = CCFULL_GUI()
        myWindow.show()
        app.exec_()
        return myWindow
    m = run_app()
    #sys.exit(app.exec_()) # when run in console
    #app.exec_()           # when run in cmd

