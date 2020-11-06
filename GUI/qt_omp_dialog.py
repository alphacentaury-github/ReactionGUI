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
                             QMenu,QMenuBar,
                             QHBoxLayout,QVBoxLayout,QGridLayout,
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

from subprocess import (call,Popen)
import numpy as np
import myutil
import reactions
from omget_RIPL3.omget import (global_omp_choose,global_omp_get)

form_omp = uic.loadUiType("dialog_omp.ui")[0]
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

class DialogOMP_GUI(QDialog,form_omp):
    def __init__(self,ap,zp,at,zt,elab):
        super().__init__()
        self.setupUi(self)
        #-----set path data
        self.path_data = myutil.all_global_variables()
        self.path_data.load_from_file()
        #--dictionary
        self.input= {'ap':ap,'at':at,'zp':zp,'zt':zt,'elab':elab,
                     'omp_id': 0 ,'omp_para': {'rc':0.0} }
        #--button action
        self.buttonBox.accepted.connect(self.take_data)
        self.buttonBox.rejected.connect(self.reject)
        #--title text
        text = ' {}{}+{}{} at Elab={} MeV'.format(
            ap,element_names[zp],at,element_names[zt],elab)
        self.label.setText(text)
        #----------------------
        # choose kinematics 
        #-----------------------
        self.lineEdit_ap.setText('{}'.format(ap))
        self.lineEdit_zp.setText('{}'.format(zp))
        self.lineEdit_at.setText('{}'.format(at))
        self.lineEdit_zt.setText('{}'.format(zt))
        self.lineEdit_elab.setText('{}'.format(elab))
        self.button_change.clicked.connect(self.change_kinematics) 
        #----list         
        self.listWidget.itemSelectionChanged.connect(self.update)
        #---initialize 
        self.change_kinematics() 
         
    def change_kinematics(self,):
        # changes kinemtics to search omp para 
        # Note that this change only affects omp para values 
        # and does not change ap,zp,at,zt of input 
        ap = int(self.lineEdit_ap.text())
        zp = int(self.lineEdit_zp.text())
        at = int(self.lineEdit_at.text())
        zt = int(self.lineEdit_zt.text())
        elab = float(self.lineEdit_elab.text())
        self.listWidget.clear() 
        self.omp_list = global_omp_choose(ap,zp,at,zt,elab)
        if self.omp_list != []:
            self.listWidget.addItem('proj.     Ztarg          Atarg                   Elab                  Author')
            for i in self.omp_list:
                om_info='{0:8}{1:3} - {2:3}     {3:3} - {4:3}     {5:8.3f} - {6:8.3f}     {7}'.format(i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8])
                self.listWidget.addItem(om_info)
            self.listWidget.setCurrentRow(1)
        else:
            self.listWidget.addItem("Sorry, we cannot suggest any appropriate OM parameter.")
        self.update()
        
    def take_data(self,):
        # actions to store data and update parent window
        # if rejected, not update parent window
#        (omp_id,omp_txt)= self.omp_list[self.listWidget.currentRow()]
        if self.temp_para == {} or self.temp_para==0:
            pass
        else :     
            (omp_id,proj,ztmin,ztmanx,atmin,atmax,emin,emax,author)= self.omp_list[self.listWidget.currentRow()-1]
            if proj == ' n ' and self.temp_para['rc'] == 0.0 :
                self.temp_para['rc'] = 1.0
            else :
                pass
            self.input['omp_id'] = omp_id
            self.input['omp_para'] = self.temp_para #accepted
            print(self.temp_para)
            self.accept()

    def update(self,):
        # file name of omeget executable
        omget_exe = self.path_data.class_property['omegt_filename']
        self.temp_para = {}
        if self.omp_list == []:
#            self.textBrowser.append("ERROR")
            self.textBrowser.clear()
        elif self.listWidget.currentRow() == 0:
            self.textBrowser.clear()
            self.textBrowser.append("Select the correct OMP list !")
        else:
            self.textBrowser.clear()
#        (omp_id,omp_txt)= self.omp_list[self.listWidget.currentRow()]
            (omp_id,proj,ztmin,ztmanx,atmin,atmax,emin,emax,author)= self.omp_list[self.listWidget.currentRow()-1]
            #---------call global_omp_get
            self.temp_para = global_omp_get(self.input['ap'],
                            self.input['zp'],self.input['at'],self.input['zt'],
                            self.input['elab'], omp_id,
                            omget_exe = omget_exe )
            if self.temp_para == {} or self.temp_para== 0 :
                self.textBrowser.append('Sorry, we can not generate this potential temporarily..')
            else : 
                text_op = ' # {} '.format(self.listWidget.currentRow())
                text_amass = ' ap = {}, at = {} '.format(
                        self.temp_para['ap'], self.temp_para['at'])
                text_coul = ' rc ={0:9.3f} '.format(self.temp_para['rc'])
                text_vol_real = ' V0 ={0:9.3f},  rv ={1:9.3f},  av ={2:9.3f} '.format(
                        self.temp_para['V'][0], self.temp_para['V'][1], self.temp_para['V'][2])
                text_vol_imag = ' W0 ={0:9.3f},  rw ={1:9.3f},  aw ={2:9.3f} '.format(
                        self.temp_para['W'][0], self.temp_para['W'][1], self.temp_para['W'][2])
                text_surf_real = ' Vd ={0:9.3f},  rd ={1:9.3f},  ad ={2:9.3f} '.format(
                        self.temp_para['Vd'][0], self.temp_para['Vd'][1], self.temp_para['Vd'][2])
                text_surf_imag = ' Wd ={0:9.3f},  rwd ={1:9.3f},  awd ={2:9.3f} '.format(
                        self.temp_para['Wd'][0], self.temp_para['Wd'][1], self.temp_para['Wd'][2])
                text_so_real = ' Vso ={0:9.3f},  rso ={1:9.3f},  aso ={2:9.3f} '.format(
                        self.temp_para['Vso'][0], self.temp_para['Vso'][1], self.temp_para['Vso'][2])
                text_so_imag = ' Wso ={0:9.3f},  rwso ={1:9.3f},  awso ={2:9.3f} '.format(
                        self.temp_para['Wso'][0], self.temp_para['Wso'][1], self.temp_para['Wso'][2])
#
                self.textBrowser.append(text_op)
                self.textBrowser.append(text_amass)
                self.textBrowser.append(text_coul)
                self.textBrowser.append(text_vol_real)
                self.textBrowser.append(text_vol_imag)
                self.textBrowser.append(text_surf_real)
                self.textBrowser.append(text_surf_imag)
                self.textBrowser.append(text_so_real)
                self.textBrowser.append(text_so_imag)
#
                self.textBrowser.append('-+-+-+-+-+-+-+- Reference -+-+-+-+-+-+-+-')
                for i in range(0,6):
                    self.textBrowser.append(self.temp_para['desc'][i])
#

        return     self.temp_para
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
        ap=1;zp=1;at=63;zt=29;elab=30.0; #default
        myWindow = DialogOMP_GUI(ap,zp,at,zt,elab)
        myWindow.show()
        # app.exec_()
        return myWindow
    m = run_app()
    #sys.exit(app.exec_()) # when run in console
    #app.exec_()           # when run in cmd
