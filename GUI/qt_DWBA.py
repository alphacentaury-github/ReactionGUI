# -*- coding: utf-8 -*-
"""
Created on Fri May 29 14:13:34 2020

@author: Y.-H. Song
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
import shutil
from io import StringIO
from subprocess import (call,Popen)

#from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import (QApplication,QMainWindow,
                             QDialog,QFileDialog,QWidget,
                             QComboBox,QLabel,QLineEdit,QCheckBox,
                             QMenu,QMenuBar,
                             QHBoxLayout,QVBoxLayout,QGridLayout,
                             QGroupBox,QToolBox,QTabWidget,
                             QPushButton,QTextBrowser,
                             QSpacerItem,
                             QRadioButton)
from PyQt5 import (uic,QtCore)
from PyQt5.QtCore import (QUrl,QSize)

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
        FigureCanvasQTAgg as FigureCanvas,
        NavigationToolbar2QT as NavigationToolbar)

import numpy as np
from numpy import loadtxt

import json
import pickle

import qt_dialog_windows
import qt_myutil
from qt_myutil import (combined_Widgets_horizontal,combined_Widgets_vertical,
                       combined_Widgets_grid,QLabel_aligned,
                       text_Browser,WidgetMatplot,about_Dialog)
import myutil
import reactions
import run_fresco_v2
#import fresco_gui_rc
from DoubleFolding.qt_DFOLD_dialog import DialogDFOLD_GUI
from qt_omp_dialog import DialogOMP_GUI


#html_head = """
#            <html><head>
#             <script type="text/javascript"
#             src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
#             </script></head>
#             <body>
#             """
#html_tail =  "</body></html>"

#------------------------------------------------------------------
class enter_nuclei(QWidget,):
    """
    a Widget to enter A,Z, mass, spin,pairty, Excitation energy

    for nuc_type = 'Excited'
           only take spin/parity and excited state energy
    for nuc_type = 'fixed'
           take fixed value of ap,zp as nuc_name

    show_ex = True : show QlineEdit  for excitation energy
              False : do not show it.
    """
    def __init__(self,label_txt='Proj',nuc_type = None,
                 nuc_name='1H',show_ex = True,
                 editing_action = None):
        super().__init__()
        self.layout = QHBoxLayout()
        self.setLayout(self.layout)
        self.layout.setContentsMargins(0,0,0,0)
        #----variables of Widget

        self.lineEdit = QLineEdit(nuc_name)
        self.lineEdit.setMaximumWidth(50)
        self.lineEdit_mass = QLineEdit()
        self.lineEdit_mass.setMaximumWidth(100)
        self.lineEdit_jp = QLineEdit()
        self.lineEdit_jp.setMaximumWidth(50)
        self.lineEdit_ex = QLineEdit('0.0')
        self.lineEdit_ex.setMaximumWidth(100)
        self.name_entered() #--initiallize values

        #---add to layout to show
        dummy_label = QLabel(label_txt)
        dummy_label.setMinimumWidth(100)
        #self.layout.addWidget(QLabel(label_txt))
        self.layout.addWidget(dummy_label)
        if not nuc_type =='Excited':
            self.layout.addWidget(self.lineEdit) # no need to enter nuc_name
        if not nuc_type =='Excited':
            self.layout.addWidget(QLabel_aligned('mass(amu)','c'))
            self.layout.addWidget(self.lineEdit_mass)
        self.layout.addWidget(QLabel_aligned('spin/parity','c'))
        self.layout.addWidget(self.lineEdit_jp)
        if show_ex :
            self.layout.addWidget(QLabel_aligned('Ex(MeV)','c'))
            self.layout.addWidget(self.lineEdit_ex)

        if nuc_type == 'fixed':
            self.lineEdit.setEnabled(False) # restrict changing name
        #---action when name has changed
        if editing_action:
           self.lineEdit.editingFinished.connect(editing_action)
        else:
           self.lineEdit.editingFinished.connect(self.name_entered)
        #--action when Ex has changed #update info
        self.lineEdit_ex.editingFinished.connect(self.get_values)

    def change_name(self,nuc_name):
        self.lineEdit.setText(nuc_name)
        self.name_entered()

    def name_entered(self,):
        nuc_name = self.lineEdit.text()
        try:
            (a,z) = reactions.interp_nuclei_name(nuc_name)
            out   = reactions.read_nuclei(a,z)
            #out has 'A','Z','name','mass','J_P','J','P'
            self.lineEdit_mass.setText('{:.6}'.format(out['mass']))
            self.lineEdit_jp.setText(out['J_P'])
        except: # if not found
            self.lineEdit_mass.setText('')
            self.lineEdit_jp.setText('')

    def get_values(self,):
        self.data ={}
        self.data['name'] = self.lineEdit.text()
        (self.data['A'],self.data['Z'])= reactions.interp_nuclei_name(
                                                    self.data['name'])
        try: # if nuclei not found
            self.data['mass'] = float( self.lineEdit_mass.text())
            self.data['J_P']   = self.lineEdit_jp.text().strip()  # string..
            #----interpret 'J_P' string
            (self.data['J'],self.data['P']) = reactions.interp_spin_name(
                                                 self.data['J_P'])
        except:
            pass
        #-----
        self.data['Ex']   = float(self.lineEdit_ex.text() )
        return self.data

    def put_values(self,data):
        """
        inverse of get_values
        """
        self.lineEdit.setText(data['name'])
        self.lineEdit_mass.setText(str(data['mass']))
        self.lineEdit_jp.setText(data['J_P'])
        self.lineEdit_ex.setText(str(data['Ex']))
        return

#--------------------------------------------------------------------------
class Partition_GUI(QWidget,):
    """
    Class for partiotion part of GUI

    take inputs for the basic information of reactions.

    (1) projectile/target
    (2) energy,
    (3) projectitle-like particle/target-like particle
    (4) kinematic info.
    """
    def __init__(self,label_txt ='Partition' ):
        super().__init__()
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.list_of_reactions = ['Elastic',
                                'Inelastic-target-ex',
                                'Inelastic-proj-ex',
                                'Transfer']
        #----Entrance channel
        self.Widget_proj = enter_nuclei(label_txt='Proj',nuc_name='2H')
        self.Widget_targ = enter_nuclei(label_txt='Targ',nuc_name='12C')

        #---for energy and reaction type
        self.lineEdit_energy = QLineEdit('10.0')
        self.comboBox_en_type = QComboBox()
        self.comboBox_en_type.addItems(['ELab','E/A','Ecm'])
        self.comboBox_reaction = QComboBox()
        self.comboBox_reaction.addItems(self.list_of_reactions)
        self.Widgets_reaction = combined_Widgets_horizontal(
            [ QLabel('Energy(MeV)'),
              self.lineEdit_energy,
              self.comboBox_en_type,
              QLabel('Reaction'),
              self.comboBox_reaction])
        #--for inelastic #nuc_name is irrelevant
        self.Widget_excited = enter_nuclei(label_txt='Excited state',
                                           nuc_name = '12C',
                                           nuc_type='Excited')
        #--for transfer
        self.Widget_light =  enter_nuclei(label_txt='projectile-like',
                                          nuc_name = '1H',
                                editing_action=self.lightParticle_entered)
        self.Widget_heavy =  enter_nuclei(label_txt='target-like',
                                          nuc_type='fixed',nuc_name='13C')
        #-- Qvalue
        self.lineEdit_qval = QLineEdit('0.0')
        self.button_qval = QPushButton('get Q')
        self.kinematics  = QPushButton('Get Kinematics info')
        self.Widget_qval = combined_Widgets_horizontal([QLabel('Q(MeV)'),
                                    self.lineEdit_qval ,self.button_qval,
                                    self.kinematics ])
        self.button_qval.clicked.connect(self.get_Q)
        self.kinematics.clicked.connect(self.get_kinematics)

        #--add to layout
        self.layout.addWidget(self.Widget_proj)
        self.layout.addWidget(self.Widget_targ)
        self.layout.addWidget(self.Widgets_reaction)
        # self.toolbox.addItem(self.Widgets_reaction   ,'reaction' )
        self.comboBox_reaction.currentIndexChanged.connect(
                                              self.select_reaction)
        self.status_label = QLabel('')
        self.layout.addWidget(self.status_label)
        self.layout.addWidget(self.Widget_qval)

        #----Widget_excited, Widget_light, Widget_heavy are not added to
        #    lay out yet. They will be dynamically placed.

    def select_reaction(self,):
        reaction_index = self.comboBox_reaction.currentIndex()
        self.status_label.setText(self.list_of_reactions[reaction_index]
                                     +' reaction is chosen')
        try: # remove previous Widgets
            self.layout.removeWidget(self.Widget_qval) #always
            self.Widget_qval.setParent(None)
            self.lineEdit_qval.setText('0.0') #reset qval
            self.layout.removeWidget(self.Widget_excited) # if existed
            self.Widget_excited.setParent(None)
        except:
            pass
        try: # remove previous Widgets
            self.layout.removeWidget(self.Widget_light)
            self.Widget_light.setParent(None)
            self.layout.removeWidget(self.Widget_heavy)
            self.Widget_heavy.setParent(None)
        except:
            pass
        #----add Widgets to layout according to reaction choice
        if reaction_index == 0:
            self.layout.addWidget(self.Widget_qval)
        elif reaction_index in [1,2]:
            self.layout.addWidget(self.Widget_excited)
            self.layout.addWidget(self.Widget_qval)
        elif reaction_index == 3:
            self.layout.addWidget(self.Widget_light)
            self.layout.addWidget(self.Widget_heavy)
            self.layout.addWidget(self.Widget_qval)
            # re-evaluate
            out = self.Widget_proj.get_values()
            self.Widget_light.lineEdit.setText(out['name'])
            self.Widget_light.name_entered()
            out = self.Widget_targ.get_values()
            self.Widget_heavy.lineEdit.setText(out['name'])
            self.Widget_heavy.name_entered()
        else :
            self.status_label.setText(self.list_of_reactions[reaction_index]
                                     +' reaction is not available yet')

    def get_Q(self,):
        #--note that Q-value have to be obtained with exact values
        proj_data = self.Widget_proj.get_values()
        targ_data = self.Widget_targ.get_values()
        reaction_index = self.comboBox_reaction.currentIndex()

        if reaction_index == 0: # Elastic
            self.lineEdit_qval.setText('0.0')
        elif reaction_index == 1 : #target excited
            excit_data = self.Widget_excited.get_values()
            Q = targ_data['Ex'] - excit_data['Ex']
            self.lineEdit_qval.setText('{:.6}'.format(Q ))
        elif reaction_index == 2 : #projectile excited
            excit_data = self.Widget_excited.get_values()
            Q = proj_data['Ex'] - excit_data['Ex']
            self.lineEdit_qval.setText('{:.6}'.format(Q ))
        elif reaction_index == 3 : # transfer
            light_data = self.Widget_light.get_values()
            heavy_data = self.Widget_heavy.get_values()
            try:
                Q = reactions.get_Qvalue(
                    proj_data['A'],proj_data['Z'],proj_data['Ex'],
                    targ_data['A'],targ_data['Z'],targ_data['Ex'],
                    light_data['A'],light_data['Z'],light_data['Ex'],
                    heavy_data['A'],heavy_data['Z'],heavy_data['Ex']
                    )
                self.lineEdit_qval.setText('{:.6}'.format(Q ))
            except:
                self.lineEdit_qval.setText('')

    def get_kinematics(self,):
        data = self.get_values()
        output=''

        if data['reaction']['type']== 0 : # Elastic
            output = reactions.kin2(data['proj']['A'],
                                data['proj']['Z'],
                                data['targ']['A'],
                                data['targ']['Z'],
                                data['reaction']['energy']['E'],
                                EN_type=data['reaction']['energy']['En_type'])
        elif data['reaction']['type']== 1 : #target excit
            output = reactions.kin2(data['proj']['A'],
                                data['proj']['Z'],
                                data['targ']['A'],
                                data['targ']['Z'],
                                data['reaction']['energy']['E'],
                                EN_type=data['reaction']['energy']['En_type'],
                                EXR=data['excited']['Ex']
                                )
        elif data['reaction']['type']== 2 : #project exit
            output = reactions.kin2(data['proj']['A'],
                                data['proj']['Z'],
                                data['targ']['A'],
                                data['targ']['Z'],
                                data['reaction']['energy']['E'],
                                EN_type=data['reaction']['energy']['En_type'],
                                EXX=data['excited']['Ex']
                                )
        elif data['reaction']['type']== 3 : # transfer
            output = reactions.kin2(data['proj']['A'],
                                data['proj']['Z'],
                                data['targ']['A'],
                                data['targ']['Z'],
                                data['reaction']['energy']['E'],
                                EN_type=data['reaction']['energy']['En_type'],
                                AX = data['light']['A'],
                                ZX = data['light']['Z'],
                                AR=  data['heavy']['A'],
                                ZR=  data['heavy']['Z'],
                                EXP= data['proj']['Ex'],
                                EXT= data['targ']['Ex'],
                                EXX= data['light']['Ex'],
                                EXR= data['heavy']['Ex'] )
        else:
            print('Error: not available yet ')
        # show result in Dialog
        dlg = about_Dialog(label_text=output,text_type='text')
        dlg.exec_()
        return output

    def lightParticle_entered(self,):
        # To override the action for editingFinished of light particle
        # update light particle info first
        self.Widget_light.name_entered()
        # then modify heavy-particle
        proj_data = self.Widget_proj.get_values()
        targ_data = self.Widget_targ.get_values()
        light_data = self.Widget_light.get_values()
        # get remnant nuclei name
        (ap,zp) = reactions.interp_nuclei_name(proj_data['name'])
        (at,zt) = reactions.interp_nuclei_name(targ_data['name'])
        (a1,z1) = reactions.interp_nuclei_name(light_data['name'])
        a = ap+at-a1
        z = zp+zt-z1
        nuc_name = '{}{}'.format(a,reactions.element_names[z])
        # fill values
        self.Widget_heavy.lineEdit.setText(nuc_name)
        self.Widget_heavy.name_entered()
        self.get_Q()

    def get_values(self,):
        """
        read current values of GUI
        return

        data = {'proj': {'name','A','Z','J','P','J_P','Ex'},
                'targ': {'name','A','Z','J','P','J_P','Ex'},
                'reaction': {'energy'('E','En_type'),'type','Q'}
                'excited': {'J','P','J_P','Ex'},
                'light': {'name','A','Z','J','P','J_P','Ex'},
                'heavy': {'name','A','Z','J','P','J_P','Ex'}   }
        """
        self.data = {}
        self.data['proj'] = self.Widget_proj.get_values()
        self.data['targ'] = self.Widget_targ.get_values()
        self.data['reaction'] = {}
        self.data['reaction']['energy'] = {
                      'E': float(self.lineEdit_energy.text()),
                      'En_type':  self.comboBox_en_type.currentIndex()}
        if self.comboBox_reaction.currentIndex()==0:
            self.data['reaction']['type'] =(
                self.comboBox_reaction.currentIndex() )
        elif self.comboBox_reaction.currentIndex() in [1,2]:
            excit_data = self.Widget_excited.get_values()
            #---Note that excit_data does not have proper info
            #    except 'J_P' and 'Ex'
            if self.comboBox_reaction.currentIndex()==1:
                excit_data['name'] = self.data['targ']['name']
                excit_data['mass'] = self.data['targ']['mass']
                excit_data['A'],excit_data['Z'] = reactions.interp_nuclei_name(
                                                          excit_data['name'])
                excit_data['J'],excit_data['P'] = reactions.interp_spin_name(
                                                            excit_data['J_P'])
            elif self.comboBox_reaction.currentIndex()==2:
                excit_data['name'] = self.data['proj']['name']
                excit_data['mass'] = self.data['proj']['mass']
                excit_data['A'],excit_data['Z'] = reactions.interp_nuclei_name(
                                                          excit_data['name'])
                excit_data['J'],excit_data['P'] = reactions.interp_spin_name(
                                                            excit_data['J_P'])
            self.data['reaction']['type'] = self.comboBox_reaction.currentIndex()
            self.data['excited'] = excit_data
        elif self.comboBox_reaction.currentIndex()==3:
            light_data = self.Widget_light.get_values()
            heavy_data = self.Widget_heavy.get_values()
            self.data['reaction']['type'] = self.comboBox_reaction.currentIndex()
            self.data['light'] = light_data
            self.data['heavy'] = heavy_data
        self.data['reaction']['Q'] = float(self.lineEdit_qval.text() )
        return self.data

    def put_values(self, data):
        """
        inverse of get_values.
        set values of inputs from data which is the same format as
        the output of get_values.
        To do this, each elements should have similar inverse of get_values.
        """
        self.Widget_proj.put_values(data['proj'])
        self.Widget_targ.put_values(data['targ'])
        self.lineEdit_energy.setText(str(data['reaction']['energy']['E']))
        self.comboBox_en_type.setCurrentIndex(
                                     data['reaction']['energy']['En_type'])
        self.comboBox_reaction.setCurrentIndex( data['reaction']['type'] )
        if data['reaction']['type']==1:
            self.Widget_excited.put_values( data['excited'])
        elif data['reaction']['type']==2:
            self.Widget_excited.put_values( data['excited'])
        elif data['reaction']['type']==3:
            self.Widget_light.put_values(data['light'])
            self.Widget_heavy.put_values(data['heavy'])
        self.lineEdit_qval.setText(str(data['reaction']['Q']))

#------------------------------------------------------------------------------

class Structure_Model_GUI(QWidget,):
    def __init__(self,partition_info=None):
        """
        Make a GUI to select reaction model
        and enter additional information for the reaction .

        Parameters
        ----------
        partition_info : dictionary
            Dictionay for the information from the Partition_GUI

        """
        super().__init__()
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.list_of_reactions = ['Elastic',
                                'Inelastic-target-ex',
                                'Inelastic-proj-ex',
                                'Transfer']
        self.partition_info = partition_info
        self.combo_model_index_old=0
        #----either use QTextBrowser or QWebEngineView
        #   Let us change browser_choice ='web' later ..
        self.text_Browser = text_Browser(browser_choice='text')
        #self.text_Browser = text_Browser(browser_choice='html')
        #self.text_Browser = text_Browser(browser_choice='web')
        self.layout.addWidget(self.text_Browser)
        self.model_txt = ''
        #-----To select model
        self.combo_model = QComboBox()
        self.combo_model.currentIndexChanged.connect(self.select_model)
        self.combo_model.setCurrentIndex(0)
        self.layout.addWidget(self.combo_model)
        #---additional information Widgets for reaction...
        #--For Rotor Model of inelastic scattering
        #--test pp
        unicode_Delta = '\u0394'
        unicode_beta = '\u03B2'
        self.lineEdit_L_transfer = QLineEdit('0')
        self.lineEdit_L_transfer.editingFinished.connect(self.check_L_transfer)
        self.lineEdit_Coul_def_length = QLineEdit('0.0')
        self.lineEdit_Nuc_def_length = QLineEdit('0.0')
        self.rotor_model = combined_Widgets_horizontal([
                        QLabel(unicode_Delta+'L transfer'),
                        self.lineEdit_L_transfer,
                        QLabel('Coulomb '+unicode_beta+'L'),
                        self.lineEdit_Coul_def_length,
                        QLabel('Nuclear '+unicode_beta+'L'),
                        self.lineEdit_Nuc_def_length ])
        #--For Cluster Model of inelastic scattering
        self.Widget_valence = enter_nuclei(label_txt='excited valence',
                                           show_ex=False,
                                           editing_action=self.valence_entered)
        self.Widget_core    = enter_nuclei(label_txt='Core',nuc_type='fixed',
                                           show_ex=False)
        self.lineEdit_gs_n = QLineEdit('1')
        self.lineEdit_gs_l = QLineEdit('0')
        self.lineEdit_gs_j = QLineEdit('0')
        self.lineEdit_gs_be = QLineEdit('0')
        self.lineEdit_ex_n = QLineEdit('1')
        self.lineEdit_ex_l = QLineEdit('0')
        self.lineEdit_ex_j = QLineEdit('0')
        self.lineEdit_ex_be = QLineEdit('0')
        self.gs_lj_recommend = QLabel('')
        self.ex_lj_recommend = QLabel('')
        # check lj
#        self.lineEdit_gs_l.editingFinished.connect(lambda: self.check_lj_bound('gs'))
#        self.lineEdit_gs_j.editingFinished.connect(lambda: self.check_lj_bound('gs'))
#        self.lineEdit_ex_l.editingFinished.connect(lambda: self.check_lj_bound('ex'))
#        self.lineEdit_ex_j.editingFinished.connect(lambda: self.check_lj_bound('ex'))

        self.Widget_gs = combined_Widgets_horizontal([QLabel('ground state '),
                            QLabel('n'),self.lineEdit_gs_n,
                            QLabel('L'),self.lineEdit_gs_l,
                            QLabel('j'),self.lineEdit_gs_j,
                            QLabel('BE'),self.lineEdit_gs_be])
        self.Widget_ex = combined_Widgets_horizontal([QLabel('excited state'),
                            QLabel('n'),self.lineEdit_ex_n,
                            QLabel('L'),self.lineEdit_ex_l,
                            QLabel('j'),self.lineEdit_ex_j,
                            QLabel('BE'),self.lineEdit_ex_be])
        self.inelastic_cluster_model = combined_Widgets_vertical([
                            QLabel(),
                            self.Widget_valence,
                            self.Widget_core,
                            self.Widget_gs,
                            self.gs_lj_recommend,
                            self.Widget_ex,
                            self.ex_lj_recommend])
        #---For Cluster model of transfer
        self.Widget_transferred_valence = enter_nuclei(
                                label_txt='transferred valence',
                                show_ex=False,nuc_type='fixed')
        self.Widget_proj_core  = enter_nuclei(
                                label_txt='Projectile Core',
                                show_ex=False,nuc_type='fixed')
        self.Widget_targ_core  = enter_nuclei(
                                label_txt='Target Core',
                                show_ex=False,nuc_type='fixed')
        self.Widget_proj_bound = combined_Widgets_horizontal([
                                 QLabel('Projectile Bound'),
                                 QLabel('n'),QLineEdit('1'),
                                 QLabel('L'),QLineEdit('0'),
                                 QLabel('j'),QLineEdit('0'),
                                 QLabel('BE'),QLineEdit('0'),
                                 QLabel('S.A.'),QLineEdit('1')
            ])
        self.Widget_targ_bound = combined_Widgets_horizontal([
                                 QLabel('Target Bound    '),
                                 QLabel('n'),QLineEdit('1'),
                                 QLabel('L'),QLineEdit('0'),
                                 QLabel('j'),QLineEdit('0'),
                                 QLabel('BE'),QLineEdit('0'),
                                 QLabel('S.A.'),QLineEdit('1')
            ])
        self.proj_lj_recommend = QLabel('')
        self.targ_lj_recommend = QLabel('')

#        self.Widget_proj_core.get_Widget(5).editingFinished.connect(
#                                 lambda: self.generate_lj_bound('proj'))
        # check lj
#        self.Widget_proj_bound.get_Widget(4).editingFinished.connect(
#                                 lambda: self.check_lj_bound('proj'))
#        self.Widget_proj_bound.get_Widget(6).editingFinished.connect(
#                                 lambda: self.check_lj_bound('proj'))
#        self.Widget_targ_bound.get_Widget(4).editingFinished.connect(
#                                 lambda: self.check_lj_bound('targ'))
#        self.Widget_targ_bound.get_Widget(6).editingFinished.connect(
#                                 lambda: self.check_lj_bound('targ'))
        self.radio_postform = QRadioButton('post form')
        self.radio_priorform = QRadioButton('prior form')
        self.radio_postform.setChecked(True)
        self.Widget_radios = combined_Widgets_horizontal(
            [QLabel('Transfer Coupling'),
             self.radio_postform,
             self.radio_priorform ]
            )

        self.transfer_cluster_model = combined_Widgets_vertical([
                                QLabel('for transfer'  ),
                                self.Widget_transferred_valence,
                                self.Widget_proj_core,
                                self.Widget_targ_core,
                                self.Widget_proj_bound,
                                self.proj_lj_recommend,
                                self.Widget_targ_bound,
                                self.targ_lj_recommend,
                                self.Widget_radios
            ])

        self.status_label = QLabel('')

    def set_partition_info(self,partition_info):
        self.partition_info = partition_info
        self.reaction_type = self.partition_info['reaction']['type']

        if self.reaction_type == 0 : #Elastic
            self.reaction_txt = '{}({},{}){}'.format(
                  self.partition_info['targ']['name'],
                  self.partition_info['proj']['name'],
                  self.partition_info['proj']['name'],
                  self.partition_info['targ']['name'])
            #----explanation texts--------
            # self.reaction_txt = '<p> Choose the structure model for {}{} reaction.</p>'.format(
            #     self.list_of_reactions[self.reaction_type], self.reaction_txt)

            # self.reaction_txt += """
            #                      <p> In general, the elastic scattering would
            #                      be described with the optical potential. </p>
            #                      <p> In the next step you will provide an optical potential.</p>
            #                      <p> If the appropriate optical potential parameters are known,
            #                      you can enter those values directly.
            #                      Here several global optical potentials(OMP) are proposed
            #                      for the reaction of interest.
            #                      These global OMPs are originated from RIPL-3 library.
            #                      <a href="https://www-nds.iaea.org/RIPL-3/">
            #                      [https://www-nds.iaea.org/RIPL-3/]
            #                      </a> </p>
            #                      <p> Alternatively the optical potential may be calculated
            #                      using a folding model.</p>
            #                      <A href=https://www-nds.iaea.org/RIPL-3/>RIPL-3</A> library.
            #                      \[U_{opt}(r)\]

            #                      """
            with open("html/description_elastic.html", "r", encoding='utf-8') as f:
                self.reaction_txt += f.read()
                
            #pageSource = html_head+ self.reaction_txt + html_tail
            pageSource = self.reaction_txt
            self.text_Browser.set_page_txt(pageSource,txt_type='html' )    
            #self.text_Browser.set_page_txt(pageSource,txt_type='web' )    
            #----available model list
            self.model_list = ['Optical Potential']

            self.combo_model.currentIndexChanged.disconnect()
            self.combo_model.clear()
            self.combo_model.addItems(self.model_list)
            if self.combo_model_index_old < len(self.model_list ):
                self.combo_model.setCurrentIndex(self.combo_model_index_old)
            self.select_model()
            self.combo_model.currentIndexChanged.connect(self.select_model)

        elif self.reaction_type == 1 : #Target exc
            self.reaction_txt = '{}({},{}){}*'.format(
                  self.partition_info['targ']['name'],
                  self.partition_info['proj']['name'],
                  self.partition_info['proj']['name'],
                  self.partition_info['targ']['name'])
            #----explanation texts--------
            # self.reaction_txt = '<p>Choose structure model for {}{} reaction.</p>'.format(
            #     self.list_of_reactions[self.reaction_type], self.reaction_txt)
            # self.reaction_txt += """<p> In Rotor Model, treat the excitation as a rotation of a deformed nuclei.
            #                         One have to specify the transferred orbital angular momentum and deformation.
            #                         In general, Coulomb deformation and Nuclear deformation can be different.
            #                         deformation length <mathjax>\(\delta_L\)</mathjax>
            #                         is related with <mathjax>\(\beta_L\)</mathjax>
            #                         and Radius,
            #                         <mathjax>$$ \delta_{L} \simeq R\\times \\beta_L.$$ </mathjax>
            #                         </p>  """
            # self.reaction_txt += """ <p> In Cluster model, excitation is modeled
            #                     as a change of s.p. energy level of a nuclei.
            #                     One have to specify the core and valence particle
            #                     and change in valence particle level.<p>"""
            with open("html/description_inelastic.html", "r", encoding='utf-8') as f:
                self.reaction_txt += f.read()
                
            #pageSource = html_head+ self.reaction_txt + html_tail
            pageSource = self.reaction_txt
            self.text_Browser.set_page_txt(pageSource,txt_type='html' )  
            #self.text_Browser.set_page_txt(pageSource,txt_type='web' )  

            #-------------------------------------------------
            # self.model_list = ['Rotor model','Cluster model']
            self.model_list = ['Rotor model'] #temporary

            self.combo_model.currentIndexChanged.disconnect()
            self.combo_model.clear()
            self.combo_model.addItems(self.model_list)
            if self.combo_model_index_old < len(self.model_list ):
                self.combo_model.setCurrentIndex(self.combo_model_index_old)
            self.select_model()
            self.combo_model.currentIndexChanged.connect(self.select_model)

        elif self.reaction_type == 2 : # projectile exci
            self.reaction_txt = '{}({},{}*){}'.format(
                  self.partition_info['targ']['name'],
                  self.partition_info['proj']['name'],
                  self.partition_info['proj']['name'],
                  self.partition_info['targ']['name'])

            # self.reaction_txt = '<p>Choose structure model for {}{} reaction.</p>'.format(
            #     self.list_of_reactions[self.reaction_type], self.reaction_txt)
            # self.reaction_txt += """<p> In Rotor Model, treat the excitation as a rotation of a deformed nuclei.
            #                         One have to specify the transferred orbital angular momentum and deformation.
            #                         In general, Coulomb deformation and Nuclear deformation can be different.
            #                         deformation length <mathjax>\(\delta_L\)</mathjax>
            #                         is related with <mathjax>\(\\beta_L\)</mathjax>
            #                         and Raidus,
            #                         <mathjax>$$\int_0^\infty \delta_{L} \simeq R\\times \\beta_L.$$ </mathjax>
            #                         </p>  """
            # self.reaction_txt += """ <p> In Cluster model, excitation is modeled
            #                     as a change of s.p. energy level of a nuclei.
            #                     One have to specify the core and valence particle
            #                     and change in valence particle level.<p>"""
            with open("html/description_inelastic.html", "r", encoding='utf-8') as f:
                self.reaction_txt += f.read()
                
            #pageSource = html_head+ self.reaction_txt + html_tail
            pageSource = self.reaction_txt
            self.text_Browser.set_page_txt(pageSource,txt_type='html' )                     
            #self.text_Browser.set_page_txt(pageSource,txt_type='web' )                     

            self.model_list = ['Rotor model','Cluster model']
            self.combo_model.currentIndexChanged.disconnect()
            self.combo_model.clear()
            self.combo_model.addItems(self.model_list)
            if self.combo_model_index_old < len(self.model_list ):
                self.combo_model.setCurrentIndex(self.combo_model_index_old)
            self.select_model()
            self.combo_model.currentIndexChanged.connect(self.select_model)

        elif self.reaction_type == 3 : #transfer
            self.reaction_txt = '{}({},{}){}'.format(
                  self.partition_info['targ']['name'],
                  self.partition_info['proj']['name'],
                  self.partition_info['light']['name'],
                  self.partition_info['heavy']['name'])
            # self.reaction_txt = '<p> Choose structure model for {}{} reaction.</p>'.format(
            #     self.list_of_reactions[self.reaction_type], self.reaction_txt)
            # self.reaction_txt += """<p> For cluster model of transfer reaction,
            #               transfer is considered as moving a valence cluster to/from target.
            #               One needs Optical potential in entrance,exit channel,
            #               potential between core1-valence,
            #               potential between core2-valence,
            #               potential between core1-core2 cluster,
            #               and bound state informations in entrance,exit channel.</p>"""
            with open("html/description_transfer.html", "r", encoding='utf-8') as f:
                self.reaction_txt += f.read()
                
            #pageSource = html_head+ self.reaction_txt + html_tail
            pageSource = self.reaction_txt
            self.text_Browser.set_page_txt(pageSource,txt_type='html' )                  
            #self.text_Browser.set_page_txt(pageSource,txt_type='web' )
                          
                          
                          
            self.model_list = ['Finite range transfer without remnant correction',
                               'Finite range transfer with remnant correction']

            self.combo_model.currentIndexChanged.disconnect()
            self.combo_model.clear()
            self.combo_model.addItems(self.model_list)
            if self.combo_model_index_old < len(self.model_list ):
                self.combo_model.setCurrentIndex(self.combo_model_index_old)
            self.select_model()
            self.combo_model.currentIndexChanged.connect(self.select_model)
        else :
            self.reaction_txt = ' (reaction is not available yet) '

        #pageSource = html_head+ self.reaction_txt + html_tail
        #pageSource = self.reaction_txt
        #self.text_Browser.set_page_txt(pageSource,txt_type='html' )

    def select_model(self,):
        current_model = self.combo_model.currentIndex()
        self.combo_model_index_old = current_model
        try:
            reaction_type = self.partition_info['reaction']['type']
            #---remove previous Widgets---
            list_of_widgets = [self.status_label,self.rotor_model,
                               self.inelastic_cluster_model,
                               self.transfer_cluster_model]
            for i in list_of_widgets:
                try:
                    self.layout.removeWidget(i)
                    i.setParent(None)
                except:
                    pass
            #--- add Widgets according to the selection
            if reaction_type == 0 : # Elastic only
                pass
            elif reaction_type==1 : # target_excitation
                if current_model == 0:
                    self.layout.addWidget(self.rotor_model)
                    self.layout.addWidget(self.status_label )
                elif current_model ==1 :
                    self.layout.addWidget(self.inelastic_cluster_model)
                    self.inelastic_cluster_model.get_Widget(0).setText(
                                                     'for Target excitation')
                    self.layout.addWidget(self.status_label )
            elif reaction_type==2 : # proj_excitation
                if current_model == 0:
                    self.layout.addWidget(self.rotor_model)
                    self.layout.addWidget(self.status_label )
                elif current_model ==1 :
                    self.layout.addWidget(self.inelastic_cluster_model)
                    self.inelastic_cluster_model.get_Widget(0).setText(
                                                 'for Projectile excitation')
                    self.layout.addWidget(self.status_label )
            elif reaction_type==3 : # transfer
                if current_model == 0:
                    self.layout.addWidget(self.transfer_cluster_model)
                    self.fix_transfer_values()
                    self.layout.addWidget(self.status_label )
                    self.status_label.setText('')
                elif current_model == 1:
                    self.layout.addWidget(self.transfer_cluster_model)
                    self.fix_transfer_values()
                    self.layout.addWidget(self.status_label )
                    self.status_label.setText('Not available yet.')
                elif current_model == 2 : # ZR approximation
                    self.layout.addWidget(self.status_label )
                    self.status_label.setText('Not available yet.')
            else:
                print('Reaction is not available yet.')
        except:
            print('Error in select_model')

    def check_L_transfer(self,):
        self.status_label.setText('')
        available_list = self.get_deform_var_list()

        if int(self.lineEdit_L_transfer.text()) == 0:
            self.status_label.setText(
                "The monopole type inelastic transition is not available.")
        elif int(self.lineEdit_L_transfer.text()) in available_list:
            self.status_label.setText("")
        elif len(available_list) == 0:
            self.status_label.setText(
                "It's wrong. There is no recommended value.")
        else:
            self.status_label.setText(
                "It's wrong. The recommended value list is {}".format(available_list))

    def get_deform_var_list(self,):
        """
        Get the appropriate list of deformation variables
        """
        deform_var_list = []
        temp_list = [0,1,2,3,4,5,6,7,8]  # free possible integer list
        if self.partition_info['reaction']['type']==1: #target excitation
            j_before = self.partition_info['targ']['J']
            band_before = self.partition_info['targ']['P']
        elif self.partition_info['reaction']['type']==2: #proj excitation
            j_before = self.partition_info['proj']['J']
            band_before = self.partition_info['proj']['P']
        j_after = self.partition_info['excited']['J']
        band_after = self.partition_info['excited']['P']

        for i in temp_list:
            if int(2*(abs(j_after-j_before)))%2 == 0:   # only inelastic
                if abs(abs(j_before-i)-1./2.) <= j_after <= j_before+i+1./2.:
                    if (band_before*band_after) == 1-2*(i%2):   # parity
                        deform_var_list.append(i)
        return deform_var_list

    def generate_lj_bound(self,opt_bound):
        """
        opt_bound -  gs    : inelastic-target(or projectile) ground state
                     ex    : inelastic-target(or projectile) excited state
                     proj  : transfer projectile bound
                     targ  : transfer target bound
        """
        reaction_type = self.partition_info['reaction']['type']
        # to generate dictionary for each reaction
        if reaction_type == 1 or reaction_type == 2:
            core_JP = self.Widget_core.get_values()['J_P']
            valence_JP = self.Widget_valence.get_values()['J_P']
            if opt_bound == 'gs':
                if reaction_type == 1:
                    bound_JP = self.partition_info['targ']['J_P']
                elif reaction_type == 2:
                    bound_JP = self.partition_info['proj']['J_P']
                given_l = int(self.lineEdit_gs_l.text())
                given_j = float(self.lineEdit_gs_j.text())
                message_type = 'ground state'
            elif opt_bound == 'ex':
                bound_JP = self.partition_info['excited']['J_P']
                given_l = int(self.lineEdit_ex_l.text())
                given_j = float(self.lineEdit_ex_j.text())
                message_type = 'excited state'
            else:
                pass
        elif reaction_type == 3:
            valence_JP = self.Widget_transferred_valence.get_values()['J_P']
            if opt_bound == 'proj':
                core_JP = self.Widget_proj_core.get_values()['J_P']
                if self.partition_info['proj']['A']-self.partition_info['light']['A'] > 0: # stripping
                    bound_JP = self.partition_info['proj']['J_P']
                else: # pick-up
                    bound_JP = self.partition_info['light']['J_P']
                given_l = int(self.Widget_proj_bound.get_Widget(4).text())
                given_j = float(self.Widget_proj_bound.get_Widget(6).text())
                message_type = (self.Widget_transferred_valence.get_values()['name']
                                +'+'+self.Widget_proj_core.get_values()['name'])
            elif opt_bound == 'targ':
                core_JP = self.Widget_targ_core.get_values()['J_P']
                if self.partition_info['proj']['A']-self.partition_info['light']['A'] > 0: # stripping
                    bound_JP = self.partition_info['heavy']['J_P']
                else: # pick-up
                    bound_JP = self.partition_info['targ']['J_P']
                given_l = int(self.Widget_targ_bound.get_Widget(4).text())
                given_j = float(self.Widget_targ_bound.get_Widget(6).text())
                message_type = (self.Widget_transferred_valence.get_values()['name']
                                +'+'+self.Widget_targ_core.get_values()['name'])
            else:
                pass
        else:
            pass

        lj_dic = self.get_quantum_number_list(bound_JP,core_JP,valence_JP)
        available_l_list = list(lj_dic.keys())
        # to make a message
        if len(available_l_list) == 0:
                recommend_txt = "     Somethig wrong. There is no recommended value."
        else:
            recommend_txt = "     The recommended value list for {} :".format(message_type)
            for ll in available_l_list :
                recommend_txt += '\n          L = {}    and   its corresponding j = {}'.format(ll,lj_dic[ll])

        return recommend_txt

    def get_quantum_number_list(self,JPbound,JPcore,JPval):
        """
        Get the appropriate list of quantum number of vaelence in the frame of cluster model
        output : dictionary {'l1':{l_j1,_j2}, 'l2':{}, ... }
        """
        possible_lj_dic = {}
        temp_list = [0,1,2,3,4,5,6,7,8,9,10]  # free possible integer list

        (Jbound,Pbound) = reactions.interp_spin_name(JPbound)
        (Jcore,Pcore) = reactions.interp_spin_name(JPcore)
        (Jval,Pval) = reactions.interp_spin_name(JPval)

        possible_2J_list = range(abs(int(2*(float(Jbound)-float(Jcore)))),int(2*(float(Jbound)+float(Jcore)))+1,2)
        for i in temp_list :
            if int(Pbound)*int(Pcore)*int(Pval) == 1-2*(i%2) : # parity determines whether l is even or odd.
                jj_list = []
                for jj in possible_2J_list :
                    if abs(2*(i-float(Jval))) <= jj <= 2*(i+float(Jval)) :
                        jj_list.append(jj/2)
                if jj_list != [] :
                    possible_lj_dic[i] = jj_list
                else :
                    pass

        return possible_lj_dic

    def valence_entered(self,):
        # action of inelastic valence_partilce has entered
        self.Widget_valence.name_entered()
        valence_data = self.Widget_valence.get_values()
        if self.partition_info['reaction']['type']==1: #target excitation
            core_A = self.partition_info['targ']['A']-valence_data['A']
            core_Z = self.partition_info['targ']['Z']-valence_data['Z']
            core_name = str(core_A)+reactions.element_names[core_Z]

            core_JP = reactions.read_nuclei(core_A,core_Z)['J_P']
            gs_lj_dic = self.get_quantum_number_list(
                self.partition_info['targ']['J_P'],core_JP,valence_data['J_P'])
            gs_l_list = list(gs_lj_dic.keys())
            gs_l = gs_l_list[0]
            gs_j = gs_lj_dic[gs_l][0]

            gs_BE = self.get_BE(self.partition_info['targ']['name']
                                ,core_name,valence_data['name'])
            self.Widget_core.change_name(core_name)
        elif self.partition_info['reaction']['type']==2: #proj excitation
            core_A = self.partition_info['proj']['A']-valence_data['A']
            core_Z = self.partition_info['proj']['Z']-valence_data['Z']
            core_name = str(core_A)+reactions.element_names[core_Z]

            core_JP = reactions.read_nuclei(core_A,core_Z)['J_P']
            gs_lj_dic = self.get_quantum_number_list(
                self.partition_info['proj']['J_P'],core_JP,valence_data['J_P'])
            gs_l_list = list(gs_lj_dic.keys())
            gs_l = gs_l_list[0]
            gs_j = gs_lj_dic[gs_l][0]

            gs_BE = self.get_BE(self.partition_info['proj']['name']
                                ,core_name,valence_data['name'])
            self.Widget_core.change_name(core_name)

        ex_BE = gs_BE - self.partition_info['excited']['Ex']
        ex_lj_dic = self.get_quantum_number_list(
            self.partition_info['excited']['J_P'],core_JP,valence_data['J_P'])
        ex_l_list = list(ex_lj_dic.keys())
        ex_l = ex_l_list[0]
        ex_j = ex_lj_dic[ex_l][0]
#        self.lineEdit_gs_l.setText('{}'.format(gs_l))
#        self.lineEdit_gs_j.setText('{}'.format(gs_j))
#        self.lineEdit_ex_l.setText('{}'.format(ex_l))
#        self.lineEdit_ex_j.setText('{}'.format(ex_j))
        self.lineEdit_gs_be.setText('{:.6}'.format(gs_BE))
        gs_rec = self.generate_lj_bound('gs')
        self.gs_lj_recommend.setText(gs_rec)
        self.lineEdit_ex_be.setText('{:.6}'.format(ex_BE))
        ex_rec = self.generate_lj_bound('ex')
        self.ex_lj_recommend.setText(ex_rec)

    def fix_transfer_values(self,):
        val_A = abs(self.partition_info['proj']['A']-self.partition_info['light']['A'])
        val_Z = abs(self.partition_info['proj']['Z']-self.partition_info['light']['Z'])
        transferred_val_name = str(val_A)+reactions.element_names[val_Z]
        transferred_val_JP = reactions.read_nuclei(val_A,val_Z)['J_P']
        if self.partition_info['proj']['A']-self.partition_info['light']['A'] > 0 : # stripping
            coreP_name = self.partition_info['light']['name']
            coreT_name = self.partition_info['targ']['name']
            proj_BE = (self.get_BE(self.partition_info['proj']['name'],
                                  coreP_name,transferred_val_name)
                     + self.partition_info['light']['Ex']
                     - self.partition_info['proj']['Ex'])
            targ_BE = (self.get_BE(self.partition_info['heavy']['name'],
                                  self.partition_info['targ']['name'],
                                  transferred_val_name)
                     + self.partition_info['targ']['Ex']
                     - self.partition_info['heavy']['Ex'])

            proj_lj_dic = self.get_quantum_number_list(
                self.partition_info['proj']['J_P'],
                self.partition_info['light']['J_P'],transferred_val_JP)
            proj_l_list = list(proj_lj_dic.keys())
#            proj_l = proj_l_list[0]
#            proj_j = proj_lj_dic[proj_l][0]
            targ_lj_dic = self.get_quantum_number_list(
                self.partition_info['heavy']['J_P'],
                self.partition_info['targ']['J_P'],transferred_val_JP)
            targ_l_list = list(targ_lj_dic.keys())
#            targ_l = targ_l_list[0]
#            targ_j = targ_lj_dic[targ_l][0]

        else : # pick-up
            coreP_name = self.partition_info['proj']['name']
            coreT_name = self.partition_info['heavy']['name']
            proj_BE = (self.get_BE(self.partition_info['light']['name'],
                                  self.partition_info['proj']['name'],
                                  transferred_val_name)
                     + self.partition_info['proj']['Ex']
                     - self.partition_info['light']['Ex'])
            targ_BE = (self.get_BE(self.partition_info['targ']['name'],
                                  coreT_name,transferred_val_name)
                     + self.partition_info['heavy']['Ex']
                     - self.partition_info['targ']['Ex'])

            proj_lj_dic = self.get_quantum_number_list(
                self.partition_info['light']['J_P'],
                self.partition_info['proj']['J_P'],
                transferred_val_JP)
            proj_l_list = list(proj_lj_dic.keys())
#            proj_l = proj_l_list[0]
#            proj_j = proj_lj_dic[proj_l][0]
            targ_lj_dic = self.get_quantum_number_list(
                self.partition_info['targ']['J_P'],
                self.partition_info['heavy']['J_P'],
                transferred_val_JP)
            targ_l_list = list(targ_lj_dic.keys())
#            targ_l = targ_l_list[0]
#            targ_j = targ_lj_dic[targ_l][0]

#        # to make a message
#        if len(available_l_list) == 0:
#                warning_txt = "It's wrong. There is no recommended value."
#        else:
#            warning_txt = "It's wrong. The recommended value list for {} is".format(message_type)
#            for ll in available_l_list :
#                warning_txt += '\n     L = {}    and   its corresponding j = {}'.format(ll,lj_dic[ll])

        self.Widget_transferred_valence.change_name(transferred_val_name)
        self.Widget_proj_core.change_name(coreP_name)
        self.Widget_targ_core.change_name(coreT_name)

        self.Widget_proj_bound.get_Widget(0).setText('{}+{} bound'.format(
            transferred_val_name,coreP_name) )
        self.Widget_targ_bound.get_Widget(0).setText('{}+{} bound'.format(
            transferred_val_name,coreT_name) )
#        self.Widget_proj_bound.get_Widget(4).setText('{}'.format(proj_l))
#        self.Widget_proj_bound.get_Widget(6).setText('{}'.format(proj_j))
#        self.Widget_targ_bound.get_Widget(4).setText('{}'.format(targ_l))
#        self.Widget_targ_bound.get_Widget(6).setText('{}'.format(targ_j))
        self.Widget_proj_bound.get_Widget(8).setText('{:.6}'.format(proj_BE))
        Proj_rec = self.generate_lj_bound('proj')
        self.proj_lj_recommend.setText(Proj_rec)
        self.Widget_targ_bound.get_Widget(8).setText('{:.6}'.format(targ_BE))
        Targ_rec = self.generate_lj_bound('targ')
        self.targ_lj_recommend.setText(Targ_rec)

    def get_BE(self,sumNuc,nuc1,nuc2):
        (sumA,sumZ) = reactions.interp_nuclei_name(sumNuc)
        (a1,z1) = reactions.interp_nuclei_name(nuc1)
        (a2,z2) = reactions.interp_nuclei_name(nuc2)
        sumM = reactions.read_nuclei(sumA,sumZ)['mass_excess']
        m1 = reactions.read_nuclei(a1,z1)['mass_excess']
        m2 = reactions.read_nuclei(a2,z2)['mass_excess']
        return m1 + m2 - sumM

    def get_values(self,):
        """
        data = {'model': Optical,Rotor,Cluster
                'L_transfer',
                'Coul_deform','Nuc_deform',
                'valence',
                'core',
                'ground_state',
                'excited_state',
                'proj_core',
                'targ_core'
                'proj_bound_state'
                'targ_bound_state'
                }
        """
        self.data = {}
        reaction_type = self.partition_info['reaction']['type']
        current_model = self.combo_model.currentIndex()

        if reaction_type == 0: #Elastic
            self.data['model'] = 'Optical'
            self.data['current_model_index']=current_model
            return self.data
        elif reaction_type in [1,2] : # target/projectile excitation
            if current_model == 0 : # rotor model
                self.data['model'] = 'Rotor'
                self.data['L_transfer'] = int(self.lineEdit_L_transfer.text())
                self.data['Coul_deform'] = float(self.lineEdit_Coul_def_length.text())
                self.data['Nuc_deform'] = float(self.lineEdit_Nuc_def_length.text())
                self.data['current_model_index']=current_model
                return self.data
            elif current_model == 1 : #Cluster Model
                self.data['model'] = 'Cluster'
                self.data['valence'] = self.Widget_valence.get_values()
                self.data['core'] = self.Widget_core.get_values()
                self.data['ground_state'] = self.Widget_gs.get_values()
                self.data['ground_state'][6] = str(myutil.frac2float(self.data['ground_state'][6]))  ## fraction to string(float)
                self.data['excited_state'] = self.Widget_ex.get_values()
                self.data['excited_state'][6] = str(myutil.frac2float(self.data['excited_state'][6]))  ## fraction to string(float)
                self.data['current_model_index']=current_model
                return self.data
        elif reaction_type == 3 : # transfer
            if current_model == 0 :
                self.data['model'] = 'FNRNG_no_remnant'
                self.data['current_model_index']=current_model
                self.data['valence'] = self.Widget_transferred_valence.get_values()
                self.data['proj_core'] = self.Widget_proj_core.get_values()
                self.data['targ_core'] = self.Widget_targ_core.get_values()
                self.data['proj_bound_state'] = self.Widget_proj_bound.get_values()
                self.data['proj_bound_state'][6] = str(myutil.frac2float(self.data['proj_bound_state'][6]))  ## fraction to string(float)
                self.data['targ_bound_state'] = self.Widget_targ_bound.get_values()
                self.data['targ_bound_state'][6] = str(myutil.frac2float(self.data['targ_bound_state'][6]))  ## fraction to string(float)
                self.data['transfer_coupling'] = self.Widget_radios.get_values()
                return self.data
            elif current_model == 1 :
                self.data['model'] = 'FNRNG_with_remnant'
                self.data['current_model_index']=current_model
                self.data['valence'] = self.Widget_transferred_valence.get_values()
                self.data['proj_core'] = self.Widget_proj_core.get_values()
                self.data['targ_core'] = self.Widget_targ_core.get_values()
                self.data['proj_bound_state'] = self.Widget_proj_bound.get_values()
                self.data['proj_bound_state'][6] = str(myutil.frac2float(self.data['proj_bound_state'][6]))  ## fraction to string(float)
                self.data['targ_bound_state'] = self.Widget_targ_bound.get_values()
                self.data['targ_bound_state'][6] = str(myutil.frac2float(self.data['targ_bound_state'][6]))  ## fraction to string(float)
                self.data['transfer_coupling'] =self.Widget_radios.get_values()
                return self.data
            elif current_model == 2 :
                self.data['model'] = 'ZRNG_LEA'
                self.data['current_model_index']=current_model
                self.data['valence'] = self.Widget_transferred_valence.get_values()
                self.data['proj_core'] = self.Widget_proj_core.get_values()
                self.data['targ_core'] = self.Widget_targ_core.get_values()
                return self.data

    def put_values(self,data):
        """
        inverse of get_values

        need to be called after set_partition_info
        """
        try:
            reaction_type = self.partition_info['reaction']['type']
        except:
            print('set_partition_info before calling put_values')
            return
        current_model = data['current_model_index']
        self.combo_model.setCurrentIndex(current_model)
        if reaction_type==0:
            pass
        elif reaction_type in [1,2]:
            if current_model==0: # rotor model
                self.lineEdit_L_transfer.setText(str(data['L_transfer'] ) )
                self.lineEdit_Coul_def_length.setText(str(data['Coul_deform']))
                self.lineEdit_Nuc_def_length.setText(str(data['Nuc_deform']))
            elif current_model==1:
                self.Widget_valence.put_values(data['valence'])
                self.Widget_core.put_values(self.data['core'])
                self.Widget_gs.put_values(data['ground_state'])
                self.Widget_ex.put_values(data['excited_state']  )
        elif reaction_type==3:
            if current_model==0:
                self.Widget_transferred_valence.put_values(data['valence'] )
                self.Widget_proj_core.put_values(data['proj_core'] )
                self.Widget_targ_core.put_values(data['targ_core'] )
                self.Widget_proj_bound.put_values(data['proj_bound_state'])
                self.Widget_targ_bound.put_values(data['targ_bound_state'])
                self.Widget_radios.put_values(data['transfer_coupling'])
            elif current_model==1:
                self.Widget_transferred_valence.put_values(data['valence'] )
                self.Widget_proj_core.put_values(data['proj_core'] )
                self.Widget_targ_core.put_values(data['targ_core'] )
                self.Widget_proj_bound.put_values(data['proj_bound_state'])
                self.Widget_targ_bound.put_values(data['targ_bound_state'])
                self.Widget_radios.put_values(data['transfer_coupling'])
            elif current_model==2:
                print('not yet available')

#==============================================================================
class pot_term_GUI2(QWidget):
    """
    The data can be accessed by

    X.update_input() and then

    X.input = [type_index,shape_index,p1,p2,p3,p4,p5,p6]

    if deform is True, enable check_Box deform which enable/disable
                       input of deformation parameter.

    folding potential case, need external function to connect the button

    """
    def __init__(self,
                 type_index=0,shape_index=0,
                 type_fixed=False,shape_fixed=False,
                 deform=False,
                 complex_optical = True,
                 disable_folding = False,
                 external_folding_slot = None):
        super().__init__()
        self.layout = QHBoxLayout()
        self.layout.setContentsMargins(0,0,0,0)
        self.setLayout(self.layout)
        #----option choices
        self.opt_type_fixed = type_fixed
        self.opt_shape_fixed = shape_fixed
        self.opt_complex_optical = complex_optical
        self.opt_deform = deform
        self.opt_disable_folding = disable_folding
        self.opt_external_folding_slot = external_folding_slot
        #----store parameters
        self.data = []
        self.folding = {}

        self.type_list = ['Coulomb',
                          'Volume',
                          'Surface',
                          'LS-proj',
                          'LS-targ']
        self.type_shape_list=[ ['Coulomb'],
                               ['WS','WS_squared','Gaussian',
                                'Yukawa','Exponential',
                                'Folding',
                                'Folding+ i WS',
                                'External'],
                               ['WS','WS_squared','Gaussian',
                                'Yukawa','Exponential','External'],
                               ['WS','WS_squared','Gaussian',
                                'Yukawa','Exponential','External'],
                               ['not yet','g','h','i']  ]
        #----defined individual widgets
        self.combo_type = QComboBox()
        self.combo_type.addItems(self.type_list)
        self.combo_type.setCurrentIndex(type_index)
        self.combo_shape = QComboBox()
        self.combo_shape.addItems(self.type_shape_list[type_index])
        self.combo_shape.setCurrentIndex(shape_index)
        if type_fixed :
            self.combo_type.setEnabled(False)
        if shape_fixed :
            self.combo_shape.setEnabled(False)
        # Widgets other than type and shape will be appended
        self.list_of_Widgets = []
        self.lineEdit_a1 = QLineEdit('0')
        self.lineEdit_a2 = QLineEdit('1')
        self.lineEdit_rc = QLineEdit('1.2')

        self.lineEdit_P1 = QLineEdit('0')
        self.lineEdit_P2 = QLineEdit('1.25')
        self.lineEdit_P3 = QLineEdit('0.65')
        self.lineEdit_P4 = QLineEdit('0')
        self.lineEdit_P5 = QLineEdit('1.25')
        self.lineEdit_P6 = QLineEdit('0.65')

        # For Folding potential case
        self.button_folding = QPushButton('Get Folding V')
        self.lineEdit_NR= QLineEdit('1.0')
        self.lineEdit_NI= QLineEdit('0.0')

        # For external complex potential case
        # which will be stored in the same way as folding potential
        self.button_external_potential = QPushButton('External Complex Potential')

        # For deformation
        self.checkBox_deform = QCheckBox('deform')
        self.lineEdit_deform = QLineEdit('0.0')
        if self.opt_deform:
            self.checkBox_deform.setEnabled(True)
            self.lineEdit_deform.setEnabled(False)
        else:
            self.checkBox_deform.setEnabled(False)
            self.lineEdit_deform.setEnabled(False)

        #-----add to Layout
        self.layout.addWidget(self.combo_type)
        self.layout.addWidget(self.combo_shape)
        self.change_shape() # initialize P1-P6
        #----additional possible button (button_omp)
        # can be added later

        # connect signal-slots slots
        self.combo_type.currentIndexChanged.connect(self.change_type)
        self.combo_shape.currentIndexChanged.connect(self.change_shape)
        self.button_external_potential.clicked.connect(self.enter_external_potential)
        if external_folding_slot:
            self.button_folding.clicked.connect(external_folding_slot)

        if self.opt_deform:
            self.checkBox_deform.stateChanged.connect(self.use_deform)


    def change_option(self,type_fixed=False,shape_fixed=False,
                 complex_optical = True,
                 deform=False,
                 disable_folding = False):
        if not(type_fixed == self.opt_type_fixed):
            self.opt_type_fixed = type_fixed
            if type_fixed :
                self.combo_type.setEnabled(False)
            else:
                self.combo_type.setEnabled(True)
        if not(self.opt_shape_fixed == shape_fixed):
            self.opt_shape_fixed = shape_fixed
            if shape_fixed:
                self.combo_shape.setEnabled(False)
            else:
                self.combo_shape.setEnabled(True)
        if not(self.opt_complex_optical == complex_optical):
            self.opt_complex_optical = complex_optical
            self.change_shape()
        if not(self.opt_disable_folding == disable_folding):
            self.opt_disable_folding = disable_folding
            self.change_shape()
        if not(self.opt_deform == deform):
            self.opt_deform = deform
            if self.opt_deform:
                self.checkBox_deform.setEnabled(True)
                self.lineEdit_deform.setEnabled(False)
                self.checkBox_deform.stateChanged.connect(self.use_deform)
            else:
                self.checkBox_deform.setEnabled(False)
                self.lineEdit_deform.setEnabled(False)

    def use_deform(self,):
        if self.checkBox_deform.isChecked():
            self.lineEdit_deform.setEnabled(True)
        else:
            self.lineEdit_deform.setEnabled(False)

    def change_type(self,):
        #--chaging type-> change available shapes
        type_indx = self.combo_type.currentIndex()
        self.combo_shape.clear()
        self.combo_shape.addItems(self.type_shape_list[type_indx] )

    def change_shape(self,):
        type_index = self.combo_type.currentIndex()
        shape_index = self.combo_shape.currentIndex()
        try: # try remove previous widgets
            for i in self.list_of_Widgets:
                self.layout.removeWidget(i)
                i.setParent(None)
            self.list_of_Widgets = []
        except:
            pass

        if self.opt_disable_folding and shape_index > 4:
            #disable these options , no further changes
            return

        if type_index == 1 and shape_index == 5 : # folding
            if self.opt_complex_optical :
                self.list_of_Widgets = [self.button_folding,
                                    QLabel('NR'),self.lineEdit_NR,
                                    QLabel('NI'),self.lineEdit_NI ]
            else:
                self.list_of_Widgets = [self.button_folding,
                                    QLabel('NR'),self.lineEdit_NR]

        elif  type_index == 1 and shape_index == 6 : # folding+ i WS
            if self.opt_complex_optical:
                self.list_of_Widgets = [self.button_folding,
                                    QLabel('NR'),self.lineEdit_NR,
                                    QLabel('W'),self.lineEdit_P4,
                                    QLabel('rw'),self.lineEdit_P5,
                                    QLabel('aw'),self.lineEdit_P6
                                    ]
            else:
                self.list_of_Widgets = [self.button_folding,
                                    QLabel('NR'),self.lineEdit_NR
                                    ]
        elif type_index == 1 and shape_index ==7 : # external potential input
            if self.opt_complex_optical :
                self.list_of_Widgets = [self.button_external_potential,
                                    QLabel('NR'),self.lineEdit_NR,
                                    QLabel('NI'),self.lineEdit_NI ]
            else:
                self.list_of_Widgets = [self.button_external_potential,
                                    QLabel('NR'),self.lineEdit_NR]

        elif type_index==0 and shape_index==0 : #Coulomb
            self.list_of_Widgets = [QLabel('a1'),self.lineEdit_a1,
                                    QLabel('a2'),self.lineEdit_a2,
                                    QLabel('rC'),self.lineEdit_rc
                                    ]
        else : # normally
            if self.opt_complex_optical :
                self.list_of_Widgets = [QLabel('V'),self.lineEdit_P1,
                                    QLabel('r'),self.lineEdit_P2,
                                    QLabel('a'),self.lineEdit_P3,
                                    QLabel('W'),self.lineEdit_P4,
                                    QLabel('rw'),self.lineEdit_P5,
                                    QLabel('aw'),self.lineEdit_P6
                                    ]
            else:
                self.list_of_Widgets = [QLabel('V'),self.lineEdit_P1,
                                    QLabel('r'),self.lineEdit_P2,
                                    QLabel('a'),self.lineEdit_P3
                                    ]
        if not (type_index==1 and shape_index >=5): # if deformation
            self.list_of_Widgets.append(self.checkBox_deform)
            self.list_of_Widgets.append(self.lineEdit_deform)

        #--add to layout
        for (i,w) in enumerate(self.list_of_Widgets):
            self.layout.addWidget(w)

    def set_type(self,type_index=0):
        self.combo_type.setCurrentIndex(type_index)

    def set_shape(self,shape_index=0):
        self.combo_shape.setCurrentIndex(shape_index)

    def set_deform(self,deform=0.0):
        self.lineEdit_deform.setText(str(deform))

    def set_external_folding_slot(self,external_folding_slot ):
        # This should update self.folding info
        self.button_folding.clicked.connect(external_folding_slot)

    def enter_external_potential(self,):
        #--dialog to call external potential input
        label_txt =( 'Enter complex potential numerically\n'
                   +' r  U_Real  U_imag\n'   )
        dialog_window = qt_dialog_windows.data_Dialog(label_text=label_txt)
        dialog_window.exec_()
        try:
            out = dialog_window.get_numpy_array()
            self.folding['R'] = out[:,0]
            self.folding['US'] = out[:,1]+1j*out[:,2] #complex optical potential
        except:
            print('Error in external potential input')

    def add_Widget(self,widget):
        self.list_of_Widgets.append(widget)
        self.layout.addWidget(widget)

    def set_paras(self,list_of_paras):
        """
        assign list of parameters for the Widgets of lineEdit
        change each lineEdit values until the last of list_of_paras
        is used.
        """
        if len(list_of_paras) > 6:
            print('Error: number of parameters>6')
        m = 0
        for (i,widg) in enumerate(self.list_of_Widgets) :
            if isinstance(widg, QLineEdit):
                widg.setText(str(list_of_paras[m]))
                m += 1 # write to QlineEdit only
                if m==len(list_of_paras): #check length of inputs
                    break

    def get_values(self,):
        self.data ={}
        type_index = self.combo_type.currentIndex()
        shape_index = self.combo_shape.currentIndex()
        self.data['type'] = type_index
        self.data['shape'] = shape_index
        for (i,widg) in enumerate(self.list_of_Widgets) :
            if isinstance(widg, QLabel) :
                self.data[i] = widg.text()
            elif isinstance(widg, QLineEdit):
                self.data[i] = float(widg.text() )
            elif isinstance(widg, QCheckBox) :
                self.data[i] = widg.isChecked()
            elif isinstance(widg,QPushButton) :
                self.data[i] = 'PushButton'
                self.data['folding'] = self.folding
        return self.data

    def put_values(self,data):
        # inverse of get_values
        self.combo_type.setCurrentIndex(data['type'] )
        self.combo_shape.setCurrentIndex(data['shape'])
        for (i,widg) in enumerate(self.list_of_Widgets) :
            if isinstance(widg, QLabel) :
                widg.setText(data[i])
            elif isinstance(widg, QLineEdit):
                widg.setText(str(data[i])  )
            elif isinstance(widg, QCheckBox) :
                widg.setCheckState(data[i] )
            elif isinstance(widg,QPushButton) : # pushbutton folding or external
                self.folding = data['folding']
        return
#-----------------------------------------------------------------------------------
class potential_tab_GUI(QWidget,):
    """
    One tab of potential,
    It can have V_C+ V+ i W + V_S +i W_S + V_LS+i W_LS

    all terms are usually WS form.
    But, V can be a double folding potential.

    folding slot should get info from partition and reaction
    and then obtain DF potential and put it in the Widget
    as Widget.folding ={}

    channel_info = {'ap','zp','at','zt','elab',
                    'Coulomb_deformation','Nuclear_deformation'}
    """
    def __init__(self,potential_name='',
                 complex_optical=True,
                 deform=False,disable_folding = False,
                 omp_button=True,
                 channel_info=None):
        super().__init__()
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.potential_name = potential_name
        self.channel_info = channel_info
        self.folding_pot = {}
        #--options
        self.complex_optical = complex_optical
        self.deform = deform
        self.disable_folding = disable_folding
        self.omp_button = omp_button
        #--prepare widgets
        self.status_label = QLabel('Title...')
        #---input for Coulomb, entrance/exit channel
        self.Widget_Coul = pot_term_GUI2(type_index=0,shape_index=0,
                                    deform=deform,
                                    type_fixed=True,shape_fixed=True)
        self.button_clear = QPushButton('Reset')
        self.button_clear.clicked.connect(self.clear_pot)
        self.button_omp = QPushButton('Get Global OMP')
        self.button_omp.clicked.connect(self.get_omp)
        if omp_button :
            self.Widget_Coul.add_Widget(self.button_clear)
            self.Widget_Coul.add_Widget(self.button_omp )

        self.Widget_Vol = pot_term_GUI2(type_index=1,shape_index=0,
                                    type_fixed=True,shape_fixed=False,
                                    deform=deform,
                                    disable_folding = disable_folding,
                                    complex_optical=complex_optical,
                                    external_folding_slot=self.get_folding)
        self.Widget_Surf = pot_term_GUI2(type_index=2,shape_index=0,
                                    type_fixed=True,shape_fixed=False,
                                    deform=deform,
                                    disable_folding = disable_folding,
                                    complex_optical=complex_optical)
        self.Widget_SO_p = pot_term_GUI2(type_index=3,shape_index=0,
                                    type_fixed=True,shape_fixed=False,
                                    deform=False,
                                    disable_folding = disable_folding,
                                    complex_optical=complex_optical)
        self.layout.addWidget(self.status_label)
        self.layout.addWidget(self.Widget_Coul )
        self.layout.addWidget(self.Widget_Vol )
        self.layout.addWidget(self.Widget_Surf )
        self.layout.addWidget(self.Widget_SO_p )
        if channel_info:
            a1 = '{}{}'.format(
                channel_info['ap'],reactions.element_names[channel_info['zp']])
            a2 = '{}{}'.format(
                channel_info['at'],reactions.element_names[channel_info['zt']])
            txt = 'Potential between {} and {}'.format(a1,a2)
            if complex_optical:
                txt = txt+' at {} MeV'.format(channel_info['elab'])
            else :
                txt = txt+' binding. potential depth of vol term will be adjusted to binding energy.'
            self.status_label.setText(txt)

        if deform and channel_info :
            self.Widget_Coul.set_deform(deform= channel_info['Coulomb_deformation'] )
            self.Widget_Vol.set_deform(deform= channel_info['Nuclear_deformation'])
            self.Widget_Surf.set_deform(deform= channel_info['Nuclear_deformation'])

    def clear_pot(self,):
        """ reset potential depths to zero.
        """
        self.Widget_Vol.set_paras([0,1.25,0.65,0,1.25,0.65])
        self.Widget_Surf.set_paras([0,1.25,0.65,0,1.25,0.65])
        self.Widget_SO_p.set_paras([0,1.25,0.65,0,1.25,0.65])

    def change_option(self, complex_optical=True,
                      deform=False,
                      disable_folding=False):
        if not(self.complex_optical==complex_optical):
            self.complex_optical = complex_optical
        if not(self.deform == deform ):
            self.deform = deform
        if not(self.disable_folding==disable_folding ):
            self.disable_folding = disable_folding
        self.Widget_Coul.change_option(type_fixed=True,shape_fixed=True,
                                       deform=deform,
                                       disable_folding = disable_folding,
                                       complex_optical=complex_optical)
        self.Widget_Vol.change_option(type_fixed=True,deform=deform,
                                      disable_folding = disable_folding,
                                      complex_optical=complex_optical)
        self.Widget_Surf.change_option(type_fixed=True,deform=deform,
                                      disable_folding = disable_folding,
                                      complex_optical=complex_optical)
        self.Widget_SO_p.change_option(type_fixed=True,deform=False,
                                      disable_folding = disable_folding,
                                      complex_optical=complex_optical)

    def set_channel_info(self,ap,zp,at,zt,elab,
                          Coulomb_deformation=0.0,
                          Nuclear_deformation=0.0):
        self.channel_info ={'ap': ap, 'zp':zp,'at':at,'zt':zt,'elab':elab,
                            'Coulomb_deformation': Coulomb_deformation,
                            'Nuclear_deformation': Nuclear_deformation}
        #--change label of potential tab
        a1 = '{}{}'.format(ap,reactions.element_names[zp])
        a2 = '{}{}'.format(at,reactions.element_names[zt])
        txt = 'Potential between {} and {}'.format(a1,a2)
        if self.complex_optical:
            txt = txt+' at {} MeV'.format(self.channel_info['elab'])
        else :
            txt = txt+' binding. potential depth of vol term will be adjusted to binding energy.'
        self.status_label.setText(txt)
        #---set initial guess of ap,at
        #self.Widget_Coul.set_paras([ap,at])

        if self.deform :
             self.Widget_Coul.set_deform(Coulomb_deformation)
             self.Widget_Vol.set_deform(Nuclear_deformation)
             self.Widget_Surf.set_deform(Nuclear_deformation)

    def get_omp(self,):
        ap = self.channel_info['ap']
        zp = self.channel_info['zp']
        at = self.channel_info['at']
        zt = self.channel_info['zt']
        elab = self.channel_info['elab']
        dlg = DialogOMP_GUI(ap,zp,at,zt,elab)
        dlg.exec_()
        #-- get optical potential parameters
        omp_para = dlg.input['omp_para']
        if omp_para['rc']> 0.01 :
            #---update omp_para
            # Coulomb part
            self.Widget_Coul.set_paras([omp_para['ap'],
                                        omp_para['at'],omp_para['rc']] )
            # Volume part need to set shape WS
            self.Widget_Vol.set_shape(0)
            self.Widget_Vol.set_paras([*omp_para['V'],*omp_para['W']  ])
            # Surf part
            self.Widget_Surf.set_shape(0)
            self.Widget_Surf.set_paras([*omp_para['Vd'],*omp_para['Wd']  ])
            # SO_p part
            self.Widget_SO_p.set_shape(0)
            self.Widget_SO_p.set_paras([*omp_para['Vso'],*omp_para['Wso']  ])
        else:
            # no update parameters.
            print('rc=0.0 thus, no update of OMP parameters')

    def get_folding(self,):
        ap = self.channel_info['ap']
        zp = self.channel_info['zp']
        at = self.channel_info['at']
        zt = self.channel_info['zt']
        e_a = self.channel_info['elab']/ap
        try:
            dlg = DialogDFOLD_GUI(ap,zp,at,zt,e_a)
            dlg.exec_()
            #----take results of folding
            self.folding_pot['R'] = dlg.DFpot['R']
            self.folding_pot['UC'] = dlg.DFpot['Coulomb']
            self.folding_pot['US'] = dlg.DFpot['Isoscalar'] #real-valued
            self.Widget_Vol.folding = self.folding_pot # update within Widget
        except:
            print('Error!!')

    def get_values(self,):
        self.data ={}
        self.data['Coul'] = self.Widget_Coul.get_values()
        self.data['Vol'] = self.Widget_Vol.get_values()
        self.data['Surf'] = self.Widget_Surf.get_values()
        self.data['SO_p'] =self.Widget_SO_p.get_values()
        return self.data

    def put_values(self,data):
        self.Widget_Coul.put_values(data['Coul'])
        self.Widget_Vol.put_values(data['Vol'])
        self.Widget_Surf.put_values(data['Surf'])
        self.Widget_SO_p.put_values(data['SO_p'])
        return

#-----------------------------------------------------------------------------------
class potentials_GUI(QTabWidget,):
    """
    Collection of potential tabs
    """
    def __init__(self,partition_data=None,structure_data=None):
        super().__init__()
        self.list_of_potentials = [] #list of potentials to appear
        #--prepare all potentials
        # for elastic and inelastic-rotor model
        # change deform, channel_info later
        self.potential_entrance = potential_tab_GUI(
                potential_name='entrance',
                complex_optical=True,
                deform=False,
                channel_info=None)
        # for inelastic-cluster model
        self.potential_coreP_coreT = potential_tab_GUI(
                potential_name='coreP_coreT',
                complex_optical=True,
                deform=False,
                channel_info =None) # coreT + coreP optical
        self.potential_coreP_val =  potential_tab_GUI(
                potential_name='coreP_val',
                complex_optical=True,
                deform=False,
                channel_info =None)
        self.potential_coreT_val =  potential_tab_GUI(
                potential_name='coreT_val',
                complex_optical=True,
                deform=False,
                channel_info = None)
        self.potential_targ_binding = potential_tab_GUI(
                potential_name='targ_binding',
                complex_optical=False,
                disable_folding = True,
                deform=False, omp_button=False,
                channel_info= None)
        self.potential_targ_binding.Widget_Vol.set_paras([70.]) #set initial guess depth
        self.potential_proj_binding = potential_tab_GUI(
                potential_name='proj_binding',
                complex_optical=False,
                disable_folding = True,
                deform=False, omp_button=False,
                channel_info= None )
        self.potential_proj_binding.Widget_Vol.set_paras([70.])
        # for transfer
        # entrance use the same potential
        # exit channelpotential
        self.potential_exit = potential_tab_GUI(
                potential_name='exit',
                complex_optical=True,
                deform=False, omp_button=True,
                channel_info=None)
        self.potential_entrance_binding = potential_tab_GUI(
                potential_name='entrance_binding',
                complex_optical=False,
                disable_folding = True,
                deform=False,omp_button=False,
                channel_info=None) # coreP +x binding
        self.potential_entrance_binding.Widget_Vol.set_paras([70.])
        self.potential_exit_binding = potential_tab_GUI(
                potential_name='exit_binding',
                complex_optical=False,
                disable_folding = True,
                deform=False,omp_button=False,
                channel_info=None) # coreT +x binding
        self.potential_exit_binding.Widget_Vol.set_paras([70.])
        self.potential_core_core_for_remnant = potential_tab_GUI(
                potential_name='core_core',
                complex_optical=True,
                deform=False, omp_button=True,
                channel_info = None )
        #---initialize tabs---
        if partition_data and structure_data:
            self.partition_data = partition_data
            self.structure_data = structure_data
            self.set_tabs_from_data(self.partition_data,self.structure_data)
        else:
            self.list_of_potentials.append(self.potential_entrance) #need to add more.. !!!!!
            self.addTab(self.potential_entrance,'Entrance channel Potential') #

    def interpret_data(self,partition_info,structure_info):
        """
        instead of raw data dictionary,
        return interpreted values for each channel
        """
        summary ={}
        AP = partition_info['proj']['A']
        ZP = partition_info['proj']['Z']
        nameP = partition_info['proj']['name']
        MP = partition_info['proj']['mass']
        JP = partition_info['proj']['J']
        PP = partition_info['proj']['P']
        EXP = partition_info['proj']['Ex']
        AT = partition_info['targ']['A']
        ZT = partition_info['targ']['Z']
        nameT = partition_info['targ']['name']
        MT = partition_info['targ']['mass']
        JT = partition_info['targ']['J']
        PT = partition_info['targ']['P']
        EXT = partition_info['targ']['Ex']
        summary['Proj'] = partition_info['proj']
        summary['Targ'] = partition_info['targ']
        #---energy
        reaction_type = partition_info['reaction']['type']
        en = partition_info['reaction']['energy']['E']
        en_type = partition_info['reaction']['energy']['En_type']
        #---------
        if reaction_type ==0:
            AX = AP; ZX = ZP; MX = MP; JX = JP; PX = PP; EXX = EXP
            nameX = nameP
            AR = AT; ZR = ZT; MR = MT; JR = JT; PR = PT; EXR = EXT
            nameR = nameT
            #!--need to set elab
            (Qval,ELAB,EpA,ECM,ELABf,EpAf,ECMf) = reactions.kin2_simplified(
                                        AP,ZP,AT,ZT,en,AX,ZX,AR,ZR,
                                        EXP,EXT,EXX,EXR,en_type)
            summary['X'] = {'name':nameX,'A':AX,'Z':ZX,'mass':MX,
                            'J':JX,'P':PX,'Ex':EXX}
            summary['R'] ={'name':nameR,'A':AR,'Z':ZR,'mass':MR,
                            'J':JR,'P':PR,'Ex':EXR}
            summary['reaction'] ={
                'type':reaction_type,
                'model': structure_info['model'],
                'current_model_index': structure_info['current_model_index'],
                'Q': Qval,'ELAB':ELAB,'EpA':EpA,'ECM':ECM,
                'ELABf':ELABf,'EpAf':EpAf,'ECMf':ECMf}

        elif reaction_type==1: # target excitation
            Jex = partition_info['excited']['J']
            Pex = partition_info['excited']['P']
            Eex = partition_info['excited']['Ex']
            AX = AP; ZX = ZP; MX = MP; JX = JP; PX = PP; EXX = EXP
            nameX = nameP
            AR = AT; ZR = ZT; MR = MT; JR = Jex; PR = Pex; EXR = Eex
            nameR = nameT
            #!---need to set elab
            (Qval,ELAB,EpA,ECM,ELABf,EpAf,ECMf) = reactions.kin2_simplified(
                                        AP,ZP,AT,ZT,en,AX,ZX,AR,ZR,
                                        EXP,EXT,EXX,EXR,en_type)

            summary['X'] = {'name':nameX,'A':AX,'Z':ZX,'mass':MX,
                            'J':JX,'P':PX,'Ex':EXX}
            summary['R'] ={'name':nameR,'A':AR,'Z':ZR,'mass':MR,
                            'J':JR,'P':PR,'Ex':EXR}
            summary['reaction'] ={
                'type':reaction_type,
                'model': structure_info['model'],
                'current_model_index': structure_info['current_model_index'],
                'Q': Qval,'ELAB':ELAB,'EpA':EpA,'ECM':ECM,
                'ELABf':ELABf,'EpAf':EpAf,'ECMf':ECMf}

            if structure_info['model']=='Rotor':
                dL = structure_info['L_transfer']
                beta_C = structure_info['Coul_deform']
                beta_N = structure_info['Nuc_deform']
                summary['Targ']['Rotor']={'L': dL,'beta_C':beta_C,'beta_N':beta_N}
            elif structure_info['model']=='Cluster':
                AV = structure_info['valence']['A'] #V for valence
                ZV = structure_info['valence']['Z']
                MV = structure_info['valence']['mass']
                JV = structure_info['valence']['J']
                PV = structure_info['valence']['P']
                nameV = structure_info['valence']['name']
                AC_T = structure_info['core']['A'] # C for core
                ZC_T = structure_info['core']['Z']
                MC_T = structure_info['core']['mass']
                JC_T = structure_info['core']['J']
                PC_T = structure_info['core']['P']
                nameC_T = structure_info['core']['name']

                N_gs = int(structure_info['ground_state'][2])
                L_gs = int(structure_info['ground_state'][4])
                J_gs = float(structure_info['ground_state'][6])
                BE_gs = float(structure_info['ground_state'][8])

                N_exs = int(structure_info['excited_state'][2])
                L_exs = int(structure_info['excited_state'][4])
                J_exs = float(structure_info['excited_state'][6])
                BE_exs = float(structure_info['excited_state'][8])

                summary['Targ']['Cluster'] = {'V': {},'C':{},'st': {} }
                summary['Targ']['Cluster']['V']={'name':nameV,'A':AV,'Z':ZV,
                                                 'mass':MV,
                                                 'J':JV,'P':PV}
                summary['Targ']['Cluster']['C']={'name':nameC_T,'A':AC_T,'Z':ZC_T,
                                                 'mass':MC_T,
                                                 'J':JC_T,'P':PC_T}
                summary['Targ']['Cluster']['st']={'n':N_gs,'l':L_gs,'j':J_gs,'BE':BE_gs}

                summary['R']['Cluster'] = {'V': {},'C':{},'st':{} }
                summary['R']['Cluster']['V'] = summary['Targ']['Cluster']['V']
                summary['R']['Cluster']['C'] = summary['Targ']['Cluster']['C']
                summary['R']['Cluster']['st']={'n':N_exs,'l':L_exs,'j':J_exs,'BE':BE_exs}

        elif reaction_type==2: #projectile excitation
            Jex = partition_info['excited']['J']
            Pex = partition_info['excited']['P']
            Eex = partition_info['excited']['Ex']
            AX = AP; ZX = ZP; MX = MP; JX = Jex; PX = Pex; EXX = Eex
            nameX = nameP
            AR = AT; ZR = ZT; MR = MT; JR = JT; PR = PT; EXR = EXT
            nameR = nameT
            #!---need to set elab
            (Qval,ELAB,EpA,ECM,ELABf,EpAf,ECMf) = reactions.kin2_simplified(
                                        AP,ZP,AT,ZT,en,AX,ZX,AR,ZR,
                                        EXP,EXT,EXX,EXR,en_type)

            summary['X'] = {'name':nameX,'A':AX,'Z':ZX,'mass':MX,
                            'J':JX,'P':PX,'Ex':EXX}
            summary['R'] ={'name':nameR,'A':AR,'Z':ZR,'mass':MR,
                            'J':JR,'P':PR,'Ex':EXR}
            summary['reaction'] ={
                'type':reaction_type,
                'model': structure_info['model'],
                'current_model_index': structure_info['current_model_index'],
                'Q': Qval,'ELAB':ELAB,'EpA':EpA,'ECM':ECM,
                'ELABf':ELABf,'EpAf':EpAf,'ECMf':ECMf}
            if structure_info['model']=='Rotor':
                dL = structure_info['L_transfer']
                beta_C = structure_info['Coul_deform']
                beta_N = structure_info['Nuc_deform']
                summary['Proj']['Rotor']={'L': dL,'beta_C':beta_C,'beta_N':beta_N}
            elif structure_info['model']=='Cluster':
                AV = structure_info['valence']['A'] # V for valence
                ZV = structure_info['valence']['Z']
                MV = structure_info['valence']['mass']
                JV = structure_info['valence']['J']
                PV = structure_info['valence']['P']
                nameV = structure_info['valence']['name']
                AC_P = structure_info['core']['A'] # C for core
                ZC_P = structure_info['core']['Z']
                MC_P = structure_info['core']['mass']
                JC_P = structure_info['core']['J']
                PC_P = structure_info['core']['P']
                nameC_P = structure_info['core']['name']

                N_gs = int(structure_info['ground_state'][2])
                L_gs = int(structure_info['ground_state'][4])
                J_gs = float(structure_info['ground_state'][6])
                BE_gs = float(structure_info['ground_state'][8])

                N_exs = int(structure_info['excited_state'][2])
                L_exs = int(structure_info['excited_state'][4])
                J_exs = float(structure_info['excited_state'][6])
                BE_exs = float(structure_info['excited_state'][8])

                summary['Proj']['Cluster'] = {'V': {},'C':{},'st': {} }
                summary['Proj']['Cluster']['V']={'name':nameV,'A':AV,'Z':ZV,
                                                 'mass':MV,
                                                 'J':JV,'P':PV}
                summary['Proj']['Cluster']['C']={'name':nameC_P,'A':AC_P,'Z':ZC_P,
                                                 'mass':MC_P,
                                                 'J':JC_P,'P':PC_P}
                summary['Proj']['Cluster']['st']={'n':N_gs,'l':L_gs,'j':J_gs,'BE':BE_gs}

                summary['X']['Cluster'] = {'V': {},'C':{},'st':{} }
                summary['X']['Cluster']['V'] = summary['Proj']['Cluster']['V']
                summary['X']['Cluster']['C'] = summary['Proj']['Cluster']['C']
                summary['X']['Cluster']['st']={'n':N_exs,'l':L_exs,'j':J_exs,'BE':BE_exs}

        elif reaction_type==3: #transfer
            AX = partition_info['light']['A']
            ZX = partition_info['light']['Z']
            nameX = partition_info['light']['name']
            MX = partition_info['light']['mass']
            JX = partition_info['light']['J']
            PX = partition_info['light']['P']
            EXX = partition_info['light']['Ex']
            AR = partition_info['heavy']['A']
            ZR = partition_info['heavy']['Z']
            nameR = partition_info['heavy']['name']
            MR = partition_info['heavy']['mass']
            JR = partition_info['heavy']['J']
            PR = partition_info['heavy']['P']
            EXR = partition_info['heavy']['Ex']
            #--get kinematics
            (Qval,ELAB,EpA,ECM,ELABf,EpAf,ECMf) = reactions.kin2_simplified(
                                        AP,ZP,AT,ZT,en,AX,ZX,AR,ZR,
                                        EXP,EXT,EXX,EXR,en_type)
            summary['X'] = {'name':nameX,'A':AX,'Z':ZX,'mass':MX,
                            'J':JX,'P':PX,'Ex':EXX}
            summary['R'] ={'name':nameR,'A':AR,'Z':ZR,'mass':MR,
                            'J':JR,'P':PR,'Ex':EXR}
            summary['reaction'] ={
                'type':reaction_type,
                'model': structure_info['model'],
                'current_model_index': structure_info['current_model_index'],
                'Q': Qval,'ELAB':ELAB,'EpA':EpA,'ECM':ECM,
                'ELABf':ELABf,'EpAf':EpAf,'ECMf':ECMf}

            AV = structure_info['valence']['A']
            ZV = structure_info['valence']['Z']
            MV = structure_info['valence']['mass']
            JV = structure_info['valence']['J']
            PV = structure_info['valence']['P']
            nameV = structure_info['valence']['name']

            AC_P = structure_info['proj_core']['A']
            ZC_P = structure_info['proj_core']['Z']
            JC_P = structure_info['proj_core']['J']
            PC_P = structure_info['proj_core']['P']
            nameC_P = structure_info['proj_core']['name']

            AC_T = structure_info['targ_core']['A']
            ZC_T = structure_info['targ_core']['Z']
            JC_T = structure_info['targ_core']['J']
            PC_T = structure_info['targ_core']['P']
            nameC_T = structure_info['targ_core']['name']
            # x_projectile_core_valence
            NCV_P = int(structure_info['proj_bound_state'][2])
            LCV_P = int(structure_info['proj_bound_state'][4])
            JCV_P = float(structure_info['proj_bound_state'][6])
            BECV_P = float(structure_info['proj_bound_state'][8])
            SACV_P = float(structure_info['proj_bound_state'][10])

            NCV_T = int(structure_info['targ_bound_state'][2])
            LCV_T = int(structure_info['targ_bound_state'][4])
            JCV_T = float(structure_info['targ_bound_state'][6])
            BECV_T = float(structure_info['targ_bound_state'][8])
            SACV_T = float(structure_info['targ_bound_state'][10])

            if AP-AX > 0 : # projectile stripping
                summary['reaction']['transfer_type'] = 'stripping'
                summary['Proj']['Cluster'] = {'V': {},'C':{},'st':{} }
                summary['Proj']['Cluster']['V'] = {'name':nameV,'A':AV,'Z':ZV,
                                              'mass':MV,'J':JV,'P':PV}
                summary['Proj']['Cluster']['C'] = summary['X']
                summary['Proj']['Cluster']['st'] = {'n':NCV_P,'l':LCV_P,
                                                    'j':JCV_P,'BE':BECV_P,
                                                    'SA': SACV_P }
                summary['R']['Cluster'] = {'V': {},'C':{},'st':{} }
                summary['R']['Cluster']['V'] = summary['Proj']['Cluster']['V']
                summary['R']['Cluster']['C'] = summary['Targ']
                summary['R']['Cluster']['st'] = {'n':NCV_T,'l':LCV_T,
                                                 'j':JCV_T,'BE':BECV_T,
                                                 'SA': SACV_T }
            else : #projectile pickup
                summary['reaction']['transfer_type'] = 'pickup'
                summary['X']['Cluster'] = {'V': {},'C':{},'st':{} }
                summary['X']['Cluster']['V'] = {'name':nameV,'A':AV,'Z':ZV,
                                              'mass':MV,'J':JV,'P':PV}
                summary['X']['Cluster']['C'] = summary['Proj']
                summary['X']['Cluster']['st'] = {'n':NCV_P,'l':LCV_P,
                                                 'j':JCV_P,'BE':BECV_P,
                                                 'SA':SACV_P}

                summary['Targ']['Cluster'] = {'V': {},'C':{},'st':{} }
                summary['Targ']['Cluster']['V'] = summary['X']['Cluster']['V']
                summary['Targ']['Cluster']['C'] = summary['R']
                summary['Targ']['Cluster']['st'] = {'n':NCV_T,'l':LCV_T,
                                                    'j':JCV_T,'BE':BECV_T,
                                                    'SA':SACV_T}
            if structure_info['current_model_index']==0:
                if structure_info['transfer_coupling'][1]: # post form
                    summary['reaction']['transfer_coupling']= 'post'
                else:
                    summary['reaction']['transfer_coupling']= 'prior'
            elif structure_info['current_model_index']==1:
                if structure_info['transfer_coupling'][1]: # post form
                    summary['reaction']['transfer_coupling']= 'post'
                else:
                    summary['reaction']['transfer_coupling']= 'prior'
            elif structure_info['current_model_index']==2:
                print('Not available yet')
        return summary


    def set_tabs_from_data(self,partition_data,structure_data):
        """
        Prepare tabs of potentials using interpreted_data of inputs
        """
        interpreted_data = self.interpret_data(partition_data,structure_data)
        ap = interpreted_data['Proj']['A']
        zp = interpreted_data['Proj']['Z']
        at = interpreted_data['Targ']['A']
        zt = interpreted_data['Targ']['Z']
        elab = interpreted_data['reaction']['ELAB']
        reaction_type = interpreted_data['reaction']['type']
        reaction_model = interpreted_data['reaction']['model']
        reaction_model_index = interpreted_data['reaction']['current_model_index']
        #----hide/show tabs
        try:
            for i in range(5): #try to remove tabs 4 times
                self.removeTab(0)
        except:
            pass

        if reaction_type== 0 : #Elastic
            self.potential_entrance.change_option(
                 complex_optical=True,deform=False,
                 disable_folding=False)
            self.potential_entrance.set_channel_info(ap,zp,at,zt,elab)
            self.addTab(self.potential_entrance,'Entrance channel')
            self.list_of_potentials = [self.potential_entrance ]
        elif reaction_type in [1,2] : # excitation
            if reaction_model =='Rotor':
                if reaction_type ==1:
                    beta_C = interpreted_data['Targ']['Rotor']['beta_C']
                    beta_N = interpreted_data['Targ']['Rotor']['beta_N']
                elif reaction_type ==2:
                    beta_C = interpreted_data['Proj']['Rotor']['beta_C']
                    beta_N = interpreted_data['Proj']['Rotor']['beta_N']

                self.potential_entrance.change_option(
                    complex_optical=True,deform=True,disable_folding=False)
                self.potential_entrance.set_channel_info(ap,zp,at,zt,elab,
                      Coulomb_deformation = beta_C  ,
                      Nuclear_deformation = beta_N )

                self.addTab(self.potential_entrance,'Entrance channel')

                self.list_of_potentials = [self.potential_entrance]
            elif reaction_model =='Cluster':
                if reaction_type==1: #target excitation
                    #---entrance optical
                    self.potential_entrance.change_option(
                        complex_optical=True,deform=False,disable_folding=False
                        )
                    self.potential_entrance.set_channel_info(ap,zp,at,zt,elab)

                    self.addTab(self.potential_entrance,'Entrance channel')
                    #---core/core optical
                    ac = interpreted_data['Targ']['Cluster']['C']['A']
                    zc = interpreted_data['Targ']['Cluster']['C']['Z']
                    self.potential_coreP_coreT.set_channel_info(
                        ap, zp, ac ,zc ,elab ) # coreT + coreP optical
                    self.addTab(self.potential_coreP_coreT,'Core-core optical')
                    #---proj-val of target optical
                    self.potential_coreP_val.set_channel_info(
                        ap, zp, at-ac, zt-zc, elab )
                    self.addTab(self.potential_coreP_val,'valence optical')
                    #---val-core binding
                    self.potential_targ_binding.set_channel_info(
                        at-ac, zt-zc, ac, zc,
                        interpreted_data['Targ']['Cluster']['st']['BE'])
                        # coreT +x binding
                    self.addTab(self.potential_targ_binding,'binding potential')
                    self.list_of_potentials = [self.potential_entrance,
                                             self.potential_coreP_coreT,
                                             self.potential_coreP_val,
                                             self.potential_targ_binding ]
                elif reaction_type==2: # proj excitation
                    #---entrance optical
                    self.potential_entrance.change_option(
                        complex_optical=True,deform=False,disable_folding=False
                        )
                    self.potential_entrance.set_channel_info(ap,zp,at,zt, elab)
                    self.addTab(self.potential_entrance,'Entrance channel')
                    #---core/core optical
                    ac = interpreted_data['Proj']['Cluster']['C']['A']
                    zc = interpreted_data['Proj']['Cluster']['C']['Z']
                    self.potential_coreP_coreT.set_channel_info(
                        ac, zc, at, zt , elab*ac/ap ) # coreT + coreP optical
                    self.addTab(self.potential_coreP_coreT,'Core-core optical')
                    #--- val +target optical
                    self.potential_coreT_val.set_channel_info(
                        ap-ac, zp-zc, at , zt , elab*(ap-ac)/ap )
                    self.addTab(self.potential_coreT_val,'valence optical')
                    #---val-core binding
                    self.potential_proj_binding.set_channel_info(
                        ap-ac, zp-zc, ac, zc, interpreted_data['Proj']['Cluster']['st']['BE'] )
                        # coreT +x binding
                    self.addTab(self.potential_proj_binding,'binding potential')

                    self.list_of_potentials = [self.potential_entrance,
                                             self.potential_coreP_coreT,
                                             self.potential_coreT_val,
                                             self.potential_proj_binding ]

        elif reaction_type ==3 : # transfer
            if reaction_model_index==0: # FNRNG w/o remnant
                # entrance channel
                self.potential_entrance.change_option(
                    complex_optical=True,deform=False,disable_folding=False)
                self.potential_entrance.set_channel_info(ap, zp, at, zt, elab)
                self.addTab(self.potential_entrance,'Entrance channel')
                # exit channel
                ax = interpreted_data['X']['A']
                zx = interpreted_data['X']['Z']
                ar = interpreted_data['R']['A']
                zr = interpreted_data['R']['Z']
                elabf = interpreted_data['reaction']['ELABf']#E_X when X is a projectile
                self.potential_exit.set_channel_info(
                    ax, zx, ar, zr, elabf )
                self.addTab(self.potential_exit,'Exit channel')
                if interpreted_data['reaction']['transfer_type']=='stripping':
                    #---stripping case--------------------
                    #  P + T -> X +R
                    # P = X + (P-X)
                    # R = T + (P-X)
                    #-------------------------------------
                    # entrance binding
                    self.potential_entrance_binding.set_channel_info(
                        ap-ax, zp-zx, ax, zx,
                        interpreted_data['Proj']['Cluster']['st']['BE']) # coreP +x binding
                    # exit binding
                    self.potential_exit_binding.set_channel_info(
                        ap-ax, zp-zx, at, zt,
                        interpreted_data['R']['Cluster']['st']['BE']) # coreT +x binding
                    self.addTab(self.potential_entrance_binding,
                        'entrance binding potential')
                    self.addTab(self.potential_exit_binding,
                        'exit binding potential')
                    # collection of potentials
                    self.list_of_potentials = [self.potential_entrance,
                                             self.potential_exit,
                                             self.potential_entrance_binding,
                                             self.potential_exit_binding]
                elif  interpreted_data['reaction']['transfer_type']=='pickup':
                    #---pickup case--------------------
                    #  P + T -> X +R
                    # X = P + (X-P)
                    # T = R + (X-P)
                    #-------------------------------------
                    # entrance binding
                    self.potential_entrance_binding.set_channel_info(
                        ax-ap , zx-zp, ar, zr,
                        interpreted_data['Targ']['Cluster']['st']['BE']) # coreT +x binding
                    # exit binding
                    self.potential_exit_binding.set_channel_info(
                        ax-ap, zx-zp, ap , zp,
                        interpreted_data['X']['Cluster']['st']['BE']) # coreP +x binding
                    self.addTab(self.potential_entrance_binding,
                        'entrance binding potential')
                    self.addTab(self.potential_exit_binding,
                        'exit binding potential')
                    self.list_of_potentials = [self.potential_entrance,
                                             self.potential_exit,
                                             self.potential_entrance_binding,
                                             self.potential_exit_binding]
            elif reaction_model_index==1 : #FNRNG w remnant
                # entrance channel
                self.potential_entrance.change_option(
                    complex_optical=True,deform=False,disable_folding=False)
                self.potential_entrance.set_channel_info(ap, zp, at, zt, elab)
                self.addTab(self.potential_entrance,'Entrance channel')
                # exit channel
                ax = interpreted_data['X']['A']
                zx = interpreted_data['X']['Z']
                ar = interpreted_data['R']['A']
                zr = interpreted_data['R']['Z']
                elabf = interpreted_data['reaction']['ELABf']
                self.potential_exit.set_channel_info(
                    ax, zx, ar, zr, elabf )
                self.addTab(self.potential_exit,'Exit channel')
                if interpreted_data['reaction']['transfer_type']=='stripping':
                    #---stripping case--------------------
                    #  P + T -> X +R
                    # P = X + (P-X)
                    # R = T + (P-X)
                    #-------------------------------------
                    # entrance binding
                    self.potential_entrance_binding.set_channel_info(
                     ap-ax, zp-zx, ax, zx,
                     interpreted_data['Proj']['Cluster']['st']['BE']) # coreP +x binding
                    self.addTab(self.potential_entrance_binding,
                        'entrance binding potential')
                    # exit binding
                    self.potential_exit_binding.set_channel_info(
                     ap-ax, zp-zx, at, zt,
                     interpreted_data['R']['Cluster']['st']['BE']) # coreT +x binding
                    self.addTab(self.potential_exit_binding,
                        'exit binding potential')
                    # core-core potential
                    print('Need to change energy for core-core potential')
                    # if stripping/prior form
                    # E_cp = Elab*a_cp/a_p
                    # if stripping/post form
                    # E_cp = Elabf
                    self.potential_core_core_for_remnant.set_channel_info(
                        ax,zx,at,zt,
                        interpreted_data['reaction']['ELABf'],
                        'core-core optical'
                        )
                    self.addTab(self.potential_core_core_for_remnant,
                        'core-core potential')
                    # list of potentials
                    self.list_of_potentials = [self.potential_entrance,
                                        self.potential_exit,
                                        self.potential_entrance_binding,
                                        self.potential_exit_binding,
                                        self.potential_core_core_for_remnant]
                elif  interpreted_data['reaction']['transfer_type']=='pickup':
                    #---pickup case--------------------
                    #  P + T -> X +R
                    # X = P + (X-P)
                    # T = R + (X-P)
                    #-------------------------------------
                    # entrance binding
                    self.potential_entrance_binding.set_channel_info(
                        ax-ap , zx-zp, ar, zr,
                        interpreted_data['Targ']['Cluster']['st']['BE']) # coreT +x binding
                    self.addTab(self.potential_entrance_binding,
                        'entrance binding potential')
                    # exit binding
                    self.potential_exit_binding.set_channel_info(
                        ax-ap, zx-zp, ap , zp,
                        interpreted_data['X']['Cluster']['st']['BE']) # coreP +x binding

                    self.addTab(self.potential_exit_binding,
                        'exit binding potential')
                    # core-core potential
                    print('Need to change energy for core-core potential')
                    # if pickup/prior form
                    # E_cp = Elab
                    # if pickup/post form
                    # E_cp = Elabf *a_cp/a_x
                    self.potential_core_core_for_remnant.set_channel_info(
                        ap,zp,ar,zr,
                        interpreted_data['reaction']['ELABf'],
                        'core-core optical'
                        )
                    self.addTab(self.potential_core_core_for_remnant,
                        'core-core potential')
                    self.list_of_potentials = [self.potential_entrance,
                                        self.potential_exit,
                                        self.potential_entrance_binding,
                                        self.potential_exit_binding,
                                        self.potential_core_core_for_remnant]

            elif reaction_model_index==2 : # ZRNG
                print('reaction_model_index==2. No more work available')
                pass


        return interpreted_data

    def get_values(self,):
        self.data = {}
        for i,widg in enumerate(self.list_of_potentials):
            self.data[widg.potential_name] = widg.get_values()
        return self.data

    def put_values(self,data):
        for i,widg in enumerate(self.list_of_potentials):
            widg.put_values(data[widg.potential_name])
        return

#=============================================================================
class potentials_and_control(QWidget,):
    def __init__(self,partition_data=None,structure_data=None):
        """
        Widget combinding Potential_tabs and control/plot
        of calculations.
        """
        super().__init__()
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.fresco_input = run_fresco_v2.Fresco_input()
        self.sfresco_input = run_fresco_v2.SFresco_Input()
        self.path_data = myutil.all_global_variables()
        self.exp_data = '# empty'
        self.exp_data_dialog = qt_dialog_windows.ExpDataDialog(
               data=self.exp_data,
               label_text='',
               SFresco_Input_object=self.sfresco_input)
        self.list_of_data_to_plot = []

        #----potential tab
        self.potential_tabs = potentials_GUI(partition_data,structure_data)
        #----collection of control buttons
        self.button_calculate = QPushButton('Calculate')
        self.button_calculate.clicked.connect(self.compute)
        self.combobox_ratio = QComboBox()
        self.combobox_ratio.addItems(['Ratio to Rutherford','mb/sr'])
        self.combobox_ratio.currentIndexChanged.connect(self.re_draw)
        self.combobox_scale = QComboBox()
        self.combobox_scale.addItems(['linear scale','log scale'])
        self.combobox_scale.currentIndexChanged.connect(self.re_draw)
        self.button_expdata = QPushButton('exp data')
        self.button_expdata.clicked.connect(self.exp_data_clicked)
        self.button_fit_elastic = QPushButton('Fit Elastic OMP')
        self.button_fit_elastic.clicked.connect(self.fit_elastic)
        self.button_fresco_option = QPushButton('Fresco Options')
        self.button_fresco_option.clicked.connect(self.change_fresco_option)
        self.button_export_output = QPushButton('Export Result')
        self.button_export_output.clicked.connect(self.export_result)
        self.button_read_fresco_output = QPushButton('Read Fresco output')
        self.button_read_fresco_output.clicked.connect(self.show_fresco_out)
        # self.button_save_parameters = QPushButton('Save inputs to file')
        # self.button_load_parameters = QPushButton('Load inputs from file')

        self.Widget_controls = combined_Widgets_grid([
            [self.button_calculate, self.combobox_ratio,
             self.combobox_scale, self.button_expdata,
             self.button_fit_elastic     ],
            [self.button_fresco_option, self.button_export_output,
             self.button_read_fresco_output ]
            ] )

        self.Widget_text = QTextBrowser()
        #----graph
        self.Widget_plot = WidgetMatplot()

        if partition_data and structure_data:
            self.interpreted_reaction = self.potential_tabs.set_tabs_from_data(
                partition_data,structure_data)

        self.layout.addWidget(self.potential_tabs)
        #---scale factor         
        self.LineEdit_scale = QLineEdit('1.0')
        self.LineEdit_scale.editingFinished.connect(self.re_draw)
        self.Widget_scale = combined_Widgets_horizontal(
            [QLabel('Scale factor'),self.LineEdit_scale]    )
        self.layout.addWidget(self.Widget_scale)
        
        self.layout.addWidget(self.Widget_controls)

        self.Widget_text_plot = combined_Widgets_horizontal(
            [self.Widget_text, self.Widget_plot])
        self.layout.addWidget(self.Widget_text_plot)

    def set_partition_structure_info(self,partition_data=None,structure_data=None):
        self.interpreted_reaction = self.potential_tabs.set_tabs_from_data(
            partition_data,structure_data)

    def interpret_potential_input(self,):
        """
        interpret potentials for Fresco input form
        """
        raw_potential_inputs = self.potential_tabs.get_values()
        reaction_type= self.interpreted_reaction['reaction']['type']
        structure_model = self.interpreted_reaction['reaction']['model']

        potentials = {}
        # ['entrance']
        # ['entrance', 'coreP_coreT', 'coreP_val', 'targ_binding']
        # ['entrance', 'coreP_coreT', 'coreT_val', 'proj_binding']
        # ['entrance', 'exit', 'entrance_binding', 'exit_binding']
        for itab in raw_potential_inputs.keys():
            potentials[itab] = []
            # Coulomb term
            a1 = raw_potential_inputs[itab]['Coul'][1]
            a2 = raw_potential_inputs[itab]['Coul'][3]
            rc = raw_potential_inputs[itab]['Coul'][5]
            potentials[itab].append(
                {'TYPE':0,'SHAPE':0,
                 'P1': a1, 'P2': a2, 'P3': rc  })
            # Coulomb excitation
            if (reaction_type==1 and structure_model=='Rotor'
                and itab=='entrance'):  #target excitation 
                L = self.interpreted_reaction['Targ']['Rotor']['L']
                if raw_potential_inputs[itab]['Coul'][6]:
                    beta = raw_potential_inputs[itab]['Coul'][7]
                else:
                    beta = 0.0
                if L> 0 :
                    multipolarity = 'P{}'.format(L)
                else:
                    raise ValueError('Error: Monopole excitation is not available.')
                #--need to get Mn(Ek)
                zt = self.interpreted_reaction['Targ']['Z']
                zp = self.interpreted_reaction['Proj']['Z']
                # test R_c values 
                # radius = (a1**(1./3.)+a2**(1./3.))*rc # r_c(a1^1/3+a2^1/3)
                radius = (a2**(1./3.))*rc               # target   
                MnEk = 3.0*(zt*beta*radius**L)/(4.0*np.pi)

                potentials[itab].append({'TYPE': 11,'SHAPE': 10,
                                         multipolarity : MnEk} )
            elif (reaction_type==2 and structure_model=='Rotor'
                  and itab=='entrance'): #projectile excitation
                L = self.interpreted_reaction['Proj']['Rotor']['L']
                if raw_potential_inputs[itab]['Coul'][6]:
                    beta = raw_potential_inputs[itab]['Coul'][7]
                else:
                    beta = 0.0
                if L> 0 :
                    multipolarity = 'P{}'.format(L)
                else:
                    raise ValueError('Error: Monopole excitation is not available.')
                zp = self.interpreted_reaction['Proj']['Z']
                radius = (a1**(1./3.))*rc #projectile 
                MnEk = 3.0*(zp*beta*radius**L)/(4.0*np.pi)

                potentials[itab].append({'TYPE': 10,'SHAPE': 10,
                                         multipolarity : MnEk } )
            else :
                pass
            # Volume term
            if itab in ['entrance','exit','coreP_coreT', 'coreP_val','coreT_val']:
                if raw_potential_inputs[itab]['Vol']['shape'] < 5:
                    potentials[itab].append(
                    { 'TYPE': raw_potential_inputs[itab]['Vol']['type'],
                      'SHAPE': raw_potential_inputs[itab]['Vol']['shape'],
                      'P1': raw_potential_inputs[itab]['Vol'][1],
                      'P2': raw_potential_inputs[itab]['Vol'][3],
                      'P3': raw_potential_inputs[itab]['Vol'][5],
                      'P4': raw_potential_inputs[itab]['Vol'][7],
                      'P5': raw_potential_inputs[itab]['Vol'][9],
                      'P6': raw_potential_inputs[itab]['Vol'][11]}
                    )
                    if reaction_type==1 and structure_model=='Rotor':
                        L = self.interpreted_reaction['Targ']['Rotor']['L']
                        if raw_potential_inputs[itab]['Vol'][12]:
                            beta = raw_potential_inputs[itab]['Vol'][13]
                        else:
                            beta = 0.0
                        if L> 0 :
                            multipolarity = 'P{}'.format(L)
                        else:
                            raise ValueError('Error: Monopole excitation is not available.')
                        at = self.interpreted_reaction['Targ']['A']
                        ap = self.interpreted_reaction['Proj']['A']
                        #radius =  (ap**(1./3.)+ at**(1./3.)) *1.2 # fixed definition of radius
                        radius =  ( at**(1./3.)) *1.2 # fixed definition of radius
                        deformation_length = beta*radius
                        potentials[itab].append({'TYPE': 11, 'SHAPE':10,
                                          multipolarity: deformation_length } )
                    elif reaction_type==2 and structure_model=='Rotor':
                        L = self.interpreted_reaction['Proj']['Rotor']['L']
                        if raw_potential_inputs[itab]['Vol'][12]:
                            beta = raw_potential_inputs[itab]['Vol'][13]
                        else:
                            beta = 0.0
                        if L> 0 :
                            multipolarity = 'P{}'.format(L)
                        else:
                            raise ValueError('Error: Monopole excitation is not available.')
                        ap = self.interpreted_reaction['Proj']['A']
                        radius = ap**(1./3.)*1.2  # fixed definition of radius
                        deformation_length = beta*radius
                        potentials[itab].append({'TYPE': 10, 'SHAPE':10,
                                          multipolarity: deformation_length } )
                elif raw_potential_inputs[itab]['Vol']['shape'] == 5: # (NR+i NI)*folding
                    # shape should be 9 of Fresco
                    # to complex number
                    external = {'R': raw_potential_inputs[itab]['Vol']['folding']['R'],
                               'US': raw_potential_inputs[itab]['Vol']['folding']['US']*(1+1j)}
                    potentials[itab].append(
                        { 'TYPE': raw_potential_inputs[itab]['Vol']['type'],
                         'SHAPE': 9,
                         'P1': raw_potential_inputs[itab]['Vol'][2], #NR
                         'P2': raw_potential_inputs[itab]['Vol'][4], #NI
                         'external': external
                                       }
                        )
                elif raw_potential_inputs[itab]['Vol']['shape'] == 6: #(US+i WS)
                    # one external + imaginary WS potential
                    external = {'R': raw_potential_inputs[itab]['Vol']['folding']['R'],
                               'US': raw_potential_inputs[itab]['Vol']['folding']['US']*(1+1j*0.0)}

                    potentials[itab].append(
                        { 'TYPE': raw_potential_inputs[itab]['Vol']['type'],
                         'SHAPE': 9,
                         'P1': raw_potential_inputs[itab]['Vol'][2],
                         'P2': 0.0,
                         'external': external   # complex number
                         })
                    potentials[itab].append(
                        { 'TYPE': raw_potential_inputs[itab]['Vol']['type'],
                         'SHAPE': 1,
                         'P4': raw_potential_inputs[itab]['Vol'][4],
                         'P5': raw_potential_inputs[itab]['Vol'][6],
                         'P6': raw_potential_inputs[itab]['Vol'][8] } )
                elif raw_potential_inputs[itab]['Vol']['shape'] == 7: #(external complex)
                    potentials[itab].append(
                        { 'TYPE': raw_potential_inputs[itab]['Vol']['type'],
                         'SHAPE': 9,
                         'P1': raw_potential_inputs[itab]['Vol'][2], #NR
                         'P2': raw_potential_inputs[itab]['Vol'][4], #NI
                         'external': raw_potential_inputs[itab]['Vol']['folding']
                         }
                        )
                else :
                    print('Error: potential shape is not available yet')
            else: # binding potential cases
                if raw_potential_inputs[itab]['Vol']['shape']> 4:
                    raise ValueError('This shape should not be used.')
                potentials[itab].append({
                    'TYPE': raw_potential_inputs[itab]['Vol']['type'],
                    'SHAPE': raw_potential_inputs[itab]['Vol']['shape'],
                    'P1': raw_potential_inputs[itab]['Vol'][1],
                    'P2': raw_potential_inputs[itab]['Vol'][3],
                    'P3': raw_potential_inputs[itab]['Vol'][5] }  )
            if itab in ['entrance','exit','coreP_coreT', 'coreP_val','coreT_val']:
                # Surface term
                potentials[itab].append(
                    {'TYPE': raw_potential_inputs[itab]['Surf']['type'],
                     'SHAPE': raw_potential_inputs[itab]['Surf']['shape'],
                     'P1': raw_potential_inputs[itab]['Surf'][1],
                     'P2': raw_potential_inputs[itab]['Surf'][3],
                     'P3': raw_potential_inputs[itab]['Surf'][5],
                     'P4': raw_potential_inputs[itab]['Surf'][7],
                     'P5': raw_potential_inputs[itab]['Surf'][9],
                     'P6': raw_potential_inputs[itab]['Surf'][11],
                     } )
                if reaction_type==1 and structure_model=='Rotor':
                    L = self.interpreted_reaction['Targ']['Rotor']['L']
                    if raw_potential_inputs[itab]['Surf'][12]:
                        beta = raw_potential_inputs[itab]['Surf'][13]
                    else:
                        beta = 0.0
                    if L> 0 :
                        multipolarity = 'P{}'.format(L)
                    else:
                        raise ValueError('Error: Monopole excitation is not available.')
                    at = self.interpreted_reaction['Targ']['A']
                    radius = at**(1./3.)*1.2  # fixed definition of radius
                    deformation_length = beta*radius
                    potentials[itab].append({'TYPE': 11, 'SHAPE':10,
                                          multipolarity: deformation_length } )
                elif reaction_type==2 and structure_model=='Rotor':
                    L = self.interpreted_reaction['Proj']['Rotor']['L']
                    if raw_potential_inputs[itab]['Surf'][12]:
                        beta = raw_potential_inputs[itab]['Surf'][13]
                    else:
                        beta = 0.0
                    if L> 0 :
                        multipolarity = 'P{}'.format(L)
                    else:
                        raise ValueError('Error: Monopole excitation is not available.')
                    ap = self.interpreted_reaction['Proj']['A']
                    radius = ap**(1./3.)*1.2  # fixed definition of radius
                    deformation_length = beta*radius
                    potentials[itab].append({'TYPE': 10, 'SHAPE':10,
                                          multipolarity: deformation_length } )
                # Spin-Orbit-Projectile
                potentials[itab].append(
                    {'TYPE': raw_potential_inputs[itab]['SO_p']['type'],
                     'SHAPE': raw_potential_inputs[itab]['SO_p']['shape'],
                     'P1': raw_potential_inputs[itab]['SO_p'][1],
                     'P2': raw_potential_inputs[itab]['SO_p'][3],
                     'P3': raw_potential_inputs[itab]['SO_p'][5],
                     'P4': raw_potential_inputs[itab]['SO_p'][7],
                     'P5': raw_potential_inputs[itab]['SO_p'][9],
                     'P6': raw_potential_inputs[itab]['SO_p'][11],
                     } )
            else:
                # Surface term
                potentials[itab].append(
                    {'TYPE': raw_potential_inputs[itab]['Surf']['type'],
                     'SHAPE': raw_potential_inputs[itab]['Surf']['shape'],
                     'P1': raw_potential_inputs[itab]['Surf'][1],
                     'P2': raw_potential_inputs[itab]['Surf'][3],
                     'P3': raw_potential_inputs[itab]['Surf'][5]
                     } )
                # Spin-Orbit-Projectile
                potentials[itab].append(
                    {'TYPE': raw_potential_inputs[itab]['SO_p']['type'],
                     'SHAPE': raw_potential_inputs[itab]['SO_p']['shape'],
                     'P1': raw_potential_inputs[itab]['SO_p'][1],
                     'P2': raw_potential_inputs[itab]['SO_p'][3],
                     'P3': raw_potential_inputs[itab]['SO_p'][5]
                     } )
        self.interpreted_potentials = potentials
        return potentials

    def compute(self,):
        """
        slots for 'calculate' button

        (1) gather all inputs
        (2) prepare fresco/sfresco inputs
        (3) call fresco code
        (4) plot results
        """
        #-----for fresco/sfresco inputs
        # reset inputs
        # self.fresco_input = run_fresco_v2.Fresco_input()
        self.fresco_input.partitions = run_fresco_v2.Fresco_Partitions()
        self.fresco_input.potentials = run_fresco_v2.Fresco_Potentials()
        self.fresco_input.overlaps = run_fresco_v2.Fresco_Overlaps()
        self.fresco_input.couplings = run_fresco_v2.Fresco_Couplings()
        
        #---list of plots
        self.list_of_data_to_plot = []
        self.Widget_plot.reset_plot()
        #---read GUI inputs
        reaction_type = self.interpreted_reaction['reaction']['type']
        structure_model = self.interpreted_reaction['reaction']['model']
        structure_model_index = self.interpreted_reaction['reaction']['current_model_index']

        self.interpret_potential_input()

        (aP,zP,aT,zT) = (self.interpreted_reaction['Proj']['A'],
                         self.interpreted_reaction['Proj']['Z'],
                         self.interpreted_reaction['Targ']['A'],
                         self.interpreted_reaction['Targ']['Z'])
        #-------------------------
        # setup head
        #-------------------------
        if reaction_type==0:
            form_txt = '{}({},{}){} at {} MeV'
        elif reaction_type==1:
            form_txt = '{}({},{}){}* at {} MeV'
        elif reaction_type==2:
            form_txt = '{}({},{}*){} at {} MeV'
        elif reaction_type==3:
            form_txt = '{}({},{}){} at {} MeV'

        title_txt = form_txt.format(
            self.interpreted_reaction['Targ']['name'],
            self.interpreted_reaction['Proj']['name'],
            self.interpreted_reaction['X']['name'],
            self.interpreted_reaction['R']['name'],
            self.interpreted_reaction['reaction']['ELAB']
            )
        self.fresco_input.head.set_items(TITLE=title_txt,
                        ELAB= self.interpreted_reaction['reaction']['ELAB'])

        #--------------------------
        # setup Partition
        #-------------------------
        self.fresco_input.partitions.change_partition(0,
              NAMEP = self.interpreted_reaction['Proj']['name'],
              MASSP = self.interpreted_reaction['Proj']['mass'],
              ZP    = self.interpreted_reaction['Proj']['Z'],
              NAMET = self.interpreted_reaction['Targ']['name'],
              MASST = self.interpreted_reaction['Targ']['mass'],
              ZT    = self.interpreted_reaction['Targ']['Z']
              )
        self.fresco_input.partitions.add_state(0,
              JP = self.interpreted_reaction['Proj']['J'],
              BANDP = self.interpreted_reaction['Proj']['P'],
              EP = self.interpreted_reaction['Proj']['Ex'],
              CPOT = 1,
              JT = self.interpreted_reaction['Targ']['J'],
              BANDT = self.interpreted_reaction['Targ']['P'],
              ET = self.interpreted_reaction['Targ']['Ex'])
        # if inelastic case
        if reaction_type==0:
            pass
        elif reaction_type==1 : # target excitation
            self.fresco_input.partitions.add_state(0,
                COPYP = 1,
                CPOT = 1,
                JT = self.interpreted_reaction['R']['J'],
                BANDT = self.interpreted_reaction['R']['P'],
                ET = self.interpreted_reaction['R']['Ex'])
            if structure_model=='Cluster':
                # target= core+ valence
                aV = self.interpreted_reaction['Targ']['Cluster']['V']['A']
                zV = self.interpreted_reaction['Targ']['Cluster']['V']['Z']
                nucX_tr = reactions.read_nuclei(aP+aV,zP+zV)
                nucR_tr = reactions.read_nuclei(aT-aV,zT-zV)
                (Q,ELAB,EpA,ECM,ELABf,EpAf,ECMf) = reactions.kin2_simplified(
                              aP,zP,aT,zT,1.0,aP+aV,zP+zV)
                self.fresco_input.partitions.add_partition(
                    NAMEP = nucX_tr['name'],
                    MASSP = nucX_tr['mass'],
                    ZP    = nucX_tr['Z'],
                    NAMET = nucR_tr['name'],
                    MASST = nucR_tr['mass'],
                    ZT    = nucR_tr['Z'],
                    QVAL  = Q )
                self.fresco_input.partitions.add_state(1,
                    JP = nucX_tr['J'], BANDP=nucX_tr['P'], EP=0.0,
                    JT= nucR_tr['J'], BANDT=nucR_tr['P'], ET =0.0)
        elif reaction_type==2: # projtile excitation
            self.fresco_input.partitions.add_state(0,
                JP = self.interpreted_reaction['X']['J'],
                BANDP = self.interpreted_reaction['X']['P'],
                ET = self.interpreted_reaction['X']['Ex'],
                CPOT = 1,
                COPYT=1)
            if structure_model=='Cluster':
                # proj=core +valence
                aV = self.interpreted_reaction['Proj']['Cluster']['V']['A']
                zV = self.interpreted_reaction['Proj']['Cluster']['V']['Z']
                nucX_tr = reactions.read_nuclei(aP-aV,zP-zV)
                nucR_tr = reactions.read_nuclei(aT+aV,zT+zV)
                (Q,ELAB,EpA,ECM,ELABf,EpAf,ECMf) = reactions.kin2_simplified(
                              aP,zP,aT,zT,1.0,aP-aV,zP-zV)
                self.fresco_input.partitions.add_partition(
                    NAMEP = nucX_tr['name'],
                    MASSP = nucX_tr['mass'],
                    ZP    = nucX_tr['Z'],
                    NAMET = nucR_tr['name'],
                    MASST = nucR_tr['mass'],
                    ZT    = nucR_tr['Z'],
                    QVAL  = Q )
                self.fresco_input.partitions.add_state(1,
                    JP = nucX_tr['J'], BANDP=nucX_tr['P'], EP=0.0,
                    JT= nucR_tr['J'], BANDT=nucR_tr['P'], ET =0.0)
        elif reaction_type==3 :
            self.fresco_input.partitions.add_partition(
                NAMEP = self.interpreted_reaction['X']['name'],
                MASSP = self.interpreted_reaction['X']['mass'],
                ZP    = self.interpreted_reaction['X']['Z'],
                NAMET = self.interpreted_reaction['R']['name'],
                MASST = self.interpreted_reaction['R']['mass'],
                ZT    = self.interpreted_reaction['R']['Z'],
                QVAL  = self.interpreted_reaction['reaction']['Q'])
            self.fresco_input.partitions.add_state(1,
                JP = self.interpreted_reaction['X']['J'],
                BANDP = self.interpreted_reaction['X']['P'],
                EP = self.interpreted_reaction['X']['Ex'],
                CPOT = 2,
                JT = self.interpreted_reaction['R']['J'],
                BANDT = self.interpreted_reaction['R']['P'],
                ET = self.interpreted_reaction['R']['Ex'])
        else:
            print('Error: not available yet')
        #-------------------------
        # prepare potential
        #
        # kp=1 : entrance
        # kp=2 : exit
        # kp=3 : coreP+coreT
        # kp=4 : coreP+val
        # kp=5 : coreT+val
        #-------------------------
        #---routine to add each potentials to fresco_input
        def add_each_potential_to_fresco_input(kp,kp_indx):
            # add potential to fresco input form
            # kp is a dictionary
                self.fresco_input.potentials.add_potential(P1=kp[0]['P1'],
                                            P2=kp[0]['P2'],P3=kp[0]['P3'])
                for i in range(1,len(kp)):
                    # in case of folding potential
                    if kp[i]['TYPE']==1 and kp[i]['SHAPE'] >=  7: # Fresco potential shape 7,8,9
                        r_points = kp[i]['external']['R']
                        pot_array = kp[i]['external']['US']
                        self.fresco_input.potentials.set_external_potential(
                            kp_indx,i,r_points,pot_array)
                        self.fresco_input.potentials.add_term(
                            kp_indx,
                            TYPE=kp[i]['TYPE'],SHAPE=kp[i]['SHAPE'],
                            P1=kp[i]['P1'],P2=kp[i]['P2'])
                    else:
                        self.fresco_input.potentials.add_term(kp_indx,**kp[i])
                return
        #-------------------------------------------------
        if reaction_type ==0:
            kp1 = self.interpreted_potentials['entrance']
            add_each_potential_to_fresco_input(kp1,0)
        elif reaction_type==1 : # target excitation in rotor model
            kp1 = self.interpreted_potentials['entrance']
            add_each_potential_to_fresco_input(kp1,0)
            if structure_model=='Cluster':
                kp2 = self.interpreted_potentials['coreP_coreT']
                kp3 = self.interpreted_potentials['coreP_val']
                kp4 = self.interpreted_potentials['targ_binding']
                add_each_potential_to_fresco_input(kp2,1)
                add_each_potential_to_fresco_input(kp3,2)
                add_each_potential_to_fresco_input(kp4,3)
        elif  reaction_type==2: #projectile excitation
            kp1 = self.interpreted_potentials['entrance']
            add_each_potential_to_fresco_input(kp1,0)
            if structure_model=='Cluster':
                kp2 = self.interpreted_potentials['coreP_coreT']
                kp3 = self.interpreted_potentials['coreT_val']
                kp4 = self.interpreted_potentials['proj_binding']
                add_each_potential_to_fresco_input(kp2,1)
                add_each_potential_to_fresco_input(kp3,2)
                add_each_potential_to_fresco_input(kp4,3)
        elif reaction_type==3: #transfer
            if structure_model_index==0: #no remnant
                kp1 = self.interpreted_potentials['entrance']
                kp2 = self.interpreted_potentials['exit']
                kp3 = self.interpreted_potentials['entrance_binding']
                kp4 = self.interpreted_potentials['exit_binding']
                add_each_potential_to_fresco_input(kp1,0)
                add_each_potential_to_fresco_input(kp2,1)
                add_each_potential_to_fresco_input(kp3,2)
                add_each_potential_to_fresco_input(kp4,3)
            elif structure_model_index==1: #with remnant
                kp1 = self.interpreted_potentials['entrance']
                kp2 = self.interpreted_potentials['exit']
                kp3 = self.interpreted_potentials['entrance_binding']
                kp4 = self.interpreted_potentials['exit_binding']
                kp5 = self.interpreted_potentials['core_core']
                add_each_potential_to_fresco_input(kp1,0)
                add_each_potential_to_fresco_input(kp2,1)
                add_each_potential_to_fresco_input(kp3,2)
                add_each_potential_to_fresco_input(kp4,3)
                add_each_potential_to_fresco_input(kp5,4)
            elif structure_model_index==2: #with remnant
                print('Not available')

        #---------------------------
        # Prepare Overlaps
        #
        # cluster inelastic excitation is not done yet !!!
        #---------------------------
        if reaction_type==1 and structure_model=='Cluster': #target excitation
            # ground state
            self.fresco_input.overlaps.add_overlap(KN1=1,IC1=1,IC2=2,IN=2,KIND=0,
                    NN=self.interpreted_reaction['Targ']['Cluster']['st']['n'] ,
                    L=self.interpreted_reaction['Targ']['Cluster']['st']['l'],
                    SN=self.interpreted_reaction['Targ']['Cluster']['V']['J'],
                    J=self.interpreted_reaction['Targ']['Cluster']['st']['j'],
                    IA=1,IB=1,KBPOT= 4,
                    BE=self.interpreted_reaction['Targ']['Cluster']['st']['BE'],
                    ISC=1,IPC=0)
            # excited state
            self.fresco_input.overlaps.add_overlap(KN1=2,IC1=1,IC2=2,IN=2,KIND=0,
                    NN=self.interpreted_reaction['R']['Cluster']['st']['n'] ,
                    L=self.interpreted_reaction['R']['Cluster']['st']['l'],
                    SN=self.interpreted_reaction['R']['Cluster']['V']['J'],
                    J=self.interpreted_reaction['R']['Cluster']['st']['j'],
                    IA=1,IB=2,KBPOT= 4,
                    BE=self.interpreted_reaction['R']['Cluster']['st']['BE'],
                    ISC=1,IPC=0)
        elif reaction_type==2 and structure_model=='Cluster': #projectile excitation
             # ground state
            self.fresco_input.overlaps.add_overlap(KN1=1,IC1=1,IC2=2,IN=1,KIND=0,
                    NN=self.interpreted_reaction['Proj']['Cluster']['st']['n'] ,
                    L=self.interpreted_reaction['Proj']['Cluster']['st']['l'],
                    SN=self.interpreted_reaction['Proj']['Cluster']['V']['J'],
                    J=self.interpreted_reaction['Proj']['Cluster']['st']['j'],
                    IA=1,IB=1,KBPOT= 4,
                    BE=self.interpreted_reaction['Proj']['Cluster']['st']['BE'],
                    ISC=1,IPC=0)
            # excited state
            self.fresco_input.overlaps.add_overlap(KN1=2,IC1=1,IC2=2,IN=1,KIND=0,
                    NN=self.interpreted_reaction['X']['Cluster']['st']['n'] ,
                    L=self.interpreted_reaction['X']['Cluster']['st']['l'],
                    SN=self.interpreted_reaction['X']['Cluster']['V']['J'],
                    J=self.interpreted_reaction['X']['Cluster']['st']['j'],
                    IA=1,IB=2,KBPOT= 4,
                    BE=self.interpreted_reaction['X']['Cluster']['st']['BE'],
                    ISC=1,IPC=0)
        elif reaction_type==3 : #transfer
            # add projectile/target bound overlap
            transfer_type=self.interpreted_reaction['reaction']['transfer_type']
            transfer_coupling = self.interpreted_reaction['reaction']['transfer_coupling']
            # bound states
            if transfer_type=='stripping':
                #---entrance_binding
                self.fresco_input.overlaps.add_overlap(
                    KN1=1,IC1=1,IC2=2,IN=1,KIND=0,
                    NN=self.interpreted_reaction['Proj']['Cluster']['st']['n'] ,
                    L=self.interpreted_reaction['Proj']['Cluster']['st']['l'],
                    SN=self.interpreted_reaction['Proj']['Cluster']['V']['J'],
                    J=self.interpreted_reaction['Proj']['Cluster']['st']['j'],
                    IA=1,IB=1,KBPOT= 3,
                    BE=self.interpreted_reaction['Proj']['Cluster']['st']['BE'],
                    ISC=1,IPC=0)
                #---exit_binding
                self.fresco_input.overlaps.add_overlap(
                    KN1=2,IC1=2,IC2=1,IN=2,KIND=0,
                    NN=self.interpreted_reaction['R']['Cluster']['st']['n'] ,
                    L=self.interpreted_reaction['R']['Cluster']['st']['l'],
                    SN=self.interpreted_reaction['R']['Cluster']['V']['J'],
                    J=self.interpreted_reaction['R']['Cluster']['st']['j'],
                    IA=1,IB=1,KBPOT= 4,
                    BE=self.interpreted_reaction['R']['Cluster']['st']['BE'],
                    ISC=1,IPC=0)
            elif transfer_type=='pickup':
                #---entrance_binding
                self.fresco_input.overlaps.add_overlap(
                    KN1=1,IC1=1,IC2=2,IN=2,KIND=0,
                    NN=self.interpreted_reaction['Targ']['Cluster']['st']['n'] ,
                    L=self.interpreted_reaction['Targ']['Cluster']['st']['l'],
                    SN=self.interpreted_reaction['Targ']['Cluster']['V']['J'],
                    J=self.interpreted_reaction['Targ']['Cluster']['st']['j'],
                    IA=1,IB=1,KBPOT= 3,
                    BE=self.interpreted_reaction['Targ']['Cluster']['st']['BE'],
                    ISC=1,IPC=0)
                #---exit binding
                self.fresco_input.overlaps.add_overlap(
                    KN1=2,IC1=2,IC2=1,IN=1,KIND=0,
                    NN=self.interpreted_reaction['X']['Cluster']['st']['n'] ,
                    L=self.interpreted_reaction['X']['Cluster']['st']['l'],
                    SN=self.interpreted_reaction['X']['Cluster']['V']['J'],
                    J=self.interpreted_reaction['X']['Cluster']['st']['j'],
                    IA=1,IB=1,KBPOT= 4,
                    BE=self.interpreted_reaction['X']['Cluster']['st']['BE'],
                    ISC=1,IPC=0)
        #---------------------------
        # Prepare Couplings
        #---------------------------
        if reaction_type==1 and structure_model=='Cluster':
            self.fresco_input.couplings.add_a_coupling(
                KIND=4,IP1=10, IP2=0, IP3=0,P1=3,P2=2)
        elif reaction_type==2 and structure_model=='Cluster':
            self.fresco_input.couplings.add_a_coupling(
                KIND=3,IP1=10, IP2=0, IP3=0,P1=3,P2=2)
        elif reaction_type==3: # transfer
            transfer_model_index = self.interpreted_reaction['reaction']['current_model_index']
            transfer_type=self.interpreted_reaction['reaction']['transfer_type']
            transfer_coupling = self.interpreted_reaction['reaction']['transfer_coupling']
            if transfer_model_index==0: #finite range transfer coupling no remnant
                print('need to check ICFROM and ICTO')
                if transfer_coupling=='post':
                    self.fresco_input.couplings.add_a_coupling(
                        KIND=7,ICTO=2,ICFROM=1,IP1=0,IP2=0,IP3=0)
                elif transfer_coupling=='prior':
                    self.fresco_input.couplings.add_a_coupling(
                        KIND=7,ICTO=2,ICFROM=1,IP1=1,IP2=0,IP3=0)
            elif transfer_model_index==1: #finite range transfer with remnant
                print('not sure about IP2=1 or -1.')
                if transfer_coupling=='post':
                    self.fresco_input.couplings.add_a_coupling(
                        KIND=7,ICTO=2,ICFROM=1,IP1=0,IP2=-1,IP3=5)
                elif transfer_coupling=='prior':
                    self.fresco_input.couplings.add_a_coupling(
                        KIND=7,ICTO=2,ICFROM=1,IP1=1,IP2=-1,IP3=5)
            elif transfer_model_index==2:
                print('Need input for P1 and P2 ')
                if transfer_coupling=='post':
                    self.fresco_input.couplings.add_a_coupling(
                        KIND=5,ICTO=2,ICFROM=1,P1=1.0,P2=1.0)
                elif transfer_coupling=='prior':
                    self.fresco_input.couplings.add_a_coupling(
                        KIND=5,ICTO=2,ICFROM=1,P1=1.0,P2=1.0)
            #add CFP
            if transfer_type=='stripping':
                sa_entrance = self.interpreted_reaction['Proj']['Cluster']['st']['SA']
                sa_exit = self.interpreted_reaction['R']['Cluster']['st']['SA']

                self.fresco_input.couplings.add_cfp(0,IN=1,IB=1,IA=1,KN=1,
                                                    A=sa_entrance) #entrance
                self.fresco_input.couplings.add_cfp(0,IN=2,IB=1,IA=1,KN=2,
                                                    A=sa_exit) #exit
            elif transfer_type=='pickup':
                sa_entrance = self.interpreted_reaction['Targ']['Cluster']['st']['SA']
                sa_exit = self.interpreted_reaction['X']['Cluster']['st']['SA']
                self.fresco_input.couplings.add_cfp(0,IN=2,IB=1,IA=1,KN=1,
                                                    A=sa_entrance) #entrance
                self.fresco_input.couplings.add_cfp(0,IN=1,IB=1,IA=1,KN=2,
                                                    A=sa_exit) #exit
        self.Widget_text.clear()
        self.Widget_text.append(self.fresco_input.write())

        #---------------------------
        # Run Fresco
        #---------------------------
        fresco_input_txt = self.fresco_input.write()

        fresco_exe = self.path_data.class_property['fresco_filename']
        fresco_frin = self.path_data.class_property['frin_filename']
        fresco_frout = self.path_data.class_property['frout_filename']

        out = run_fresco_v2.run_fresco_from_input_txt(
                   fresco_input_txt,
                   fresco_path = fresco_exe,
                   fresco_input_path = fresco_frin ,
                   fresco_output_path= fresco_frout)

        # check result
        if run_fresco_v2.chck_fresco_out(fname=fresco_frout)==0:
            self.Widget_text.append('Fresco run finished.\n' )
        else:
            self.Widget_text.append('Error in Fresco run.\n' )
            ff=open(fresco_frout,'r')
            ll=ff.readlines()
            ff.close()
            # show last few lines of output
            self.Widget_text.append(ll[-3]+ll[-2]+ll[-1])
            return
        #---------------------------
        # Get results
        #---------------------------
        elastic_out = run_fresco_v2.get_elastic_result_from_fresco_out(
                fname=fresco_frout)
        angle = elastic_out[:,0]
        xs = elastic_out[:,1]
        ratio = elastic_out[:,2]
        out[0] = elastic_out # replace elastic part
        # neutral scattering case neutral_case=0
        # and np.sum(ratio) = 0.
        neutral_case=( self.interpreted_reaction['Proj']['Z']
                     *self.interpreted_reaction['Targ']['Z'] )
        #--------
        # Store results for later plotting
        #--------
        self.fresco_results = out
        #---------------------------
        # Plot Results
        #---------------------------
        scale_factor = 1.0 # by default
        try:
            scale_factor = float(self.LineEdit_scale.text())
        except:
            print('Error in scale factor')
            scale_factor=1.0 
            
        if reaction_type==0: # default plot is ratio ,linear
            self.combobox_ratio.setEnabled(True)
            self.button_fit_elastic.setEnabled(True)
            # special case of neutral particle
            if self.combobox_ratio.currentIndex()==0:
                if (neutral_case)>0:
                    self.current_plot = {'x': angle, 'y': ratio*scale_factor,
                                      'fmt':'-' ,'label':'Elastic'}
                    ylabel = 'ratio'
                else : #neutral case 
                    self.current_plot = {'x': angle, 'y': ratio*scale_factor,
                    'fmt':'-' ,'label':'ratio un-defined for neutral particle!'}
                    ylabel = 'ratio'
            else:
                self.current_plot = {'x': angle, 'y': xs*scale_factor,
                                      'fmt':'-' ,'label':'Elastic'}
                ylabel = 'mb/sr'
            # construct list of plots
            # list[0] = 'elastic' or 'dwba' cross section
            # list[1] = experimental data
            self.list_of_data_to_plot.append(self.current_plot)
            # if exp_data exists
            if self.exp_data.strip() : # what if comments exits?
                try:
                    out = myutil.txt_to_array(self.exp_data)
                    self.list_of_data_to_plot.append({'x': out[:,0] ,
                                          'y': out[:,1] ,
                                          'yerr': out[:,2] ,
                                          'fmt': 'o',
                                          'label':'exp'})
                except:
                    print('error in exp_data_plot')

            if self.combobox_scale.currentIndex() == 0 :
                yscale = 'linear'
            else:
                yscale = 'log'
            self.Widget_plot.add_plot_by_data(
                    self.list_of_data_to_plot,
                    show_grid=False,yscale=yscale,
                    xlabel='c.m. angle[deg]',
                    ylabel=ylabel)
        else: # DWBA case
            self.combobox_ratio.setEnabled(False)
            self.button_fit_elastic.setEnabled(False)
            self.current_plot = {'x': self.fresco_results[1][:,0],
                                 'y': self.fresco_results[1][:,1]*scale_factor,
                                  'fmt': '-','label':'DWBA'}
            # construct list of plots
            # construct list of plots
            # list[0] = 'elastic' or 'dwba' cross section
            # list[1] = experimental data
            self.list_of_data_to_plot.append(self.current_plot)
            # if exp_data exists
            if self.exp_data.strip() : # what if comments exits?
                try:
                    out = myutil.txt_to_array(self.exp_data)
                    self.list_of_data_to_plot.append({'x': out[:,0] ,
                                          'y': out[:,1] ,
                                          'yerr': out[:,2] ,
                                          'fmt': 'o',
                                          'label':'exp'})
                except:
                    print('error in exp_data_plot')

            if self.combobox_scale.currentIndex()==0:
                yscale='linear'
            else:
                yscale='log'
            self.Widget_plot.add_plot_by_data(self.list_of_data_to_plot,
                          show_grid=False,yscale=yscale,
                          xlabel='c.m. angle[deg]',
                          ylabel='mb/sr')

    def re_draw(self,):
        # # change of plot_option (ratio or scale)
        # only replace list_of_plot[0]
        reaction_type = self.interpreted_reaction['reaction']['type']
        scale_factor = 1.0 # by default
        try:
            scale_factor = float(self.LineEdit_scale.text())
        except:
            print('Error in scale factor')
            scale_factor=1.0 
            
        if reaction_type==0: # default plot is ratio ,linear
            neutral_case=( self.interpreted_reaction['Proj']['Z']
                     *self.interpreted_reaction['Targ']['Z'] ) 
            if self.combobox_ratio.currentIndex()==0:
                if neutral_case > 0: 
                    self.current_plot = {'x': self.fresco_results[0][:,0],
                                      'y': self.fresco_results[0][:,2]*scale_factor,
                                      'fmt':'-' ,'label':'Elastic'}
                else:
                    self.current_plot = {'x': self.fresco_results[0][:,0],
                                      'y': self.fresco_results[0][:,2]*scale_factor,
                                      'fmt':'-' ,
                                      'label':'ratio undefined for neutral particle'}
                ylabel = 'ratio'
            else:
                self.current_plot = {'x': self.fresco_results[0][:,0],
                                      'y': self.fresco_results[0][:,1]*scale_factor,
                                      'fmt':'-' ,'label':'Elastic'}
                ylabel = 'mb/sr'

            if len(self.list_of_data_to_plot) > 0:    
                self.list_of_data_to_plot[0] = self.current_plot
            else: #empty 
                self.list_of_data_to_plot.append( self.current_plot) 

            if self.combobox_scale.currentIndex() == 0 :
                yscale = 'linear'
            else:
                yscale = 'log'

            self.Widget_plot.add_plot_by_data(self.list_of_data_to_plot,
                          show_grid=False,yscale=yscale,
                          xlabel='c.m. angle[deg]',
                          ylabel=ylabel)
        else:
            self.current_plot = {'x': self.fresco_results[1][:,0],
                                  'y': self.fresco_results[1][:,1]*scale_factor,
                                  'fmt': '-','label':'DWBA'}
            self.list_of_data_to_plot[0] = self.current_plot
            if self.combobox_scale.currentIndex()==0:
                yscale='linear'
            else:
                yscale='log'
            self.Widget_plot.add_plot_by_data(self.list_of_data_to_plot,
                          show_grid=False,yscale=yscale,
                          xlabel='c.m. angle[deg]',
                          ylabel='mb/sr')
        return

    def exp_data_clicked(self,):
        """
        exp_data button clicked
        """
        label_txt = ( 'Enter data :\n'
                     +'angle mb error\n'
                     +'# for comments\n'
                     +'*** Under construction.\n'
                     +' At the moment, data conversion is not available\n')
        dlg = self.exp_data_dialog
        self.exp_data_dialog.modify(data=self.exp_data,
                           label_text=label_txt,
                           SFresco_Input_object=self.sfresco_input)
        check = dlg.exec_()
        if check==1 : # accepted
            self.exp_data = dlg.data #update
            self.Widget_text.append(self.exp_data)
            # plot with exp data...
            out = myutil.txt_to_array(self.exp_data)
            # empty or error in data
            if len(out)==0 :
                print('converting txt to array error')
                if len(self.list_of_data_to_plot)<=1:
                    self.Widget_text.append('Error: could not draw experimental data.')
                    return
                else:
                    # remove all plots except calculated one
                    self.list_of_data_to_plot = [ self.list_of_data_to_plot[0] ]
                    self.re_draw()
                    return
            # if okay
            if len(self.list_of_data_to_plot)<=1:
                self.list_of_data_to_plot.append( {'x': out[:,0] ,
                                          'y': out[:,1] ,
                                          'yerr': out[:,2] ,
                                          'fmt': 'o',
                                          'label':'exp'})
            else:
                self.list_of_data_to_plot[1]= {'x': out[:,0] ,
                                          'y': out[:,1] ,
                                          'yerr': out[:,2] ,
                                          'fmt': 'o',
                                          'label':'exp'}
            try:
                self.re_draw()
            except:
                self.Widget_text.append('Error: could not draw figure.')
                self.Widget_text.append('   Either error in the experimental data'
                                   +' or no reaction result is available.' )
        return

    def show_fresco_out(self,):
        fname=self.path_data.class_property['frout_filename']
        try:
            print(fname)
            ff = open(fname,'r')
            output = ff.read()
            print(output)
            ff.close()
            dlg = about_Dialog(label_text= output  ,text_type='text')
            dlg.exec_()
        except:
            print('Error in show_fresco_out')
        return

    def export_result(self,):
        try:
            if len(self.fresco_results)==0:
                print('no fresco results')
                return
            elif len(self.fresco_results)==1: #elastic
                elastic = self.fresco_results[0]
            elif len(self.fresco_results) > 1: #elastic+DWBA
                elastic = self.fresco_results[0]
                dwba = self.fresco_results[1]
        except:
            print('error reading fresco results')
            return
        options = QFileDialog.Options()
        fileName, _filter = QFileDialog.getSaveFileName(self,
                    "Save file",
                    "",
                    "All Files (*)",
                    options=options)
        if fileName : 
            print('Save Results to {}'.format(fileName))
            #---write to fileName
            txt = '# Elastic channel \n'
            txt += '#  theta                   mb/sr                   ratio \n'
            txt +=  myutil.array_to_txt(elastic)+'\n\n'
            if len(self.fresco_results) > 1:
                txt += '# DWBA \n'
                txt += '#  theta                   mb/sr   \n'
                txt += myutil.array_to_txt(dwba)
            print(txt)
            ff=open(fileName,'w')
            ff.write(txt)
            ff.close()
            #---copy fresco fort.16x to fileName
            #shutil.copy('fort.16x',fileName)
        else:
            return

    def change_fresco_option(self,):
        # change head of fresco into namelist forms
        old_head = self.fresco_input.head.data
        # separate TITLE part
        old_head = {key:val for key, val in old_head.items() if key != 'TITLE'}
        old_head_txt = myutil.write_namelist_form(
            old_head,upper_case_item=True,separator='\n')
        label_txt =( 'Enter Fresco options \n'   )
        dialog_window = qt_dialog_windows.data_Dialog(
            data = old_head_txt, label_text=label_txt)
        dialog_window.exec_()
        new_head_txt = dialog_window.data
        new_head = myutil.read_namelist_input(new_head_txt,
                    upper_case_item=True)
        for key in new_head.keys() :
            self.fresco_input.head.set_items(**{key : new_head[key]})
        return

    def fit_elastic(self,):
        # calls sfresco dialog
        reaction_type = self.interpreted_reaction['reaction']['type']
        if reaction_type==0:
            dlg = qt_dialog_windows.Elastic_Fit_Dialog(
                fresco_input_object=self.fresco_input,
                exp_data=self.exp_data,
                sfresco_input_object=self.sfresco_input,
                path_data=self.path_data.class_property)
            dlg_accepted =dlg.exec_()
            #--closed dialog window
            if (dlg_accepted):
                try:
                    new_para_text = dlg.fit_result['text']
                    new_para_list = dlg.fit_result['para']
                    new_result_array = dlg.fit_result['plot'][1]
                    self.Widget_text.append('new parameters: \n'+new_para_text)
                    print('To do: update parameters after fitting. ')
                    print('To do: update plot with new one')
                    self.list_of_data_to_plot.append({
                         'x': new_result_array[:,0],
                         'y': new_result_array[:,1],
                         'fmt': '-',
                         'label': 'fit'})
                    self.re_draw()

                except:
                    print('Error in ftting? ')
            # fitting is done
            return
        else:
            self.Widget_text.append(
                'Fitting is not available for non-elastic reactions yet.\n')
            return



#===============================================================================
class TestToolBox(QToolBox,):
    def __init__(self,):
        super().__init__()
        #self.layout = QVBoxLayout()
        #self.setLayout(self.layout)
        #size = QSize(720,690)
        #self.resize(size)
        self.partition = Partition_GUI()
        self.model = Structure_Model_GUI()
        self.pot_control_plot = potentials_and_control()

        self.addItem(self.partition,'Partition')
        self.addItem(self.model,'Structure Model')
        self.addItem(self.pot_control_plot,'Potentials/compute/plot')
        self.setCurrentIndex(0)
        self.currentChanged.connect(self.toolbox_changed)

        # additional signal/slots
        #self.pot_control_plot.button_save_parameters.clicked.connect(self.save_file )
        #self.pot_control_plot.button_load_parameters.clicked.connect(self.load_file )

        print('To do: better description of each model')
        print('To do: Test with various reaction examples. ')
        print('To do: check consistency of elab values in each potential tab')
        print('To do: Error handling.')
        print('To do: transfer model in zero-range approximation')
        print('To do: enable radiative capture and fusion reaction ')

    def toolbox_changed(self,):
        current_index = self.currentIndex()
        if current_index == 0 :
            #do something?
            pass
        elif current_index == 1 :
            partition_data = self.partition.get_values()
            self.model.set_partition_info(partition_data)
        elif current_index == 2:
            partition_data = self.partition.get_values()
            structure_data = self.model.get_values()
            self.pot_control_plot.set_partition_structure_info(
                partition_data,structure_data)

    def save_para(self,):
        """
        return input parameters as dictionary
        """
        partition_data = self.partition.get_values()
        structure_data = self.model.get_values()
        potential_data = self.pot_control_plot.potential_tabs.get_values()
        input_paras = {'partition': partition_data,
                       'structure': structure_data,
                       'potential': potential_data,
                       'exp_data' : self.pot_control_plot.exp_data }
        return input_paras

    def load_para(self,data):
        """
        read data and put them in the GUI

        dictionary data have the same structure as output of save_data
        """
        self.partition.put_values(data['partition'])
        self.model.set_partition_info(data['partition'])
        self.model.put_values(data['structure'])
        self.pot_control_plot.potential_tabs.set_tabs_from_data(
            data['partition'], data['structure'])
        self.pot_control_plot.potential_tabs.put_values(data['potential'])
        try: # set exp_data
            self.pot_control_plot.exp_data = data['exp_data']
        except:
            print('No exp_data saved')
        # self.pot_control_plot.potential_tabs.show()

    def save_file(self,):
        options = QFileDialog.Options()
        fileName, _filter = QFileDialog.getSaveFileName(self,
                    "Save file",
                    "",
                    "All Files (*)",
                    options=options)
        para_dict   = self.save_para()
        # save dictionary to file
        #>> using json
        #jj = json.dumps(para_dict,indent=2)
        #ff = open(fileName,'w')
        #ff.write(jj)
        #ff.close() 
        #>> or pickle
        ff= open(fileName,'wb')
        pickle.dump(para_dict,ff)
        ff.close()
        return para_dict

    def load_file(self,):
        options = QFileDialog.Options()
        fileName, _filter = QFileDialog.getOpenFileName(self,
                    "Load file",
                    "",
                    "All Files (*)",
                    options=options)
        # load file as dictionary
        #>> using json
        #ff = open(fileName,'r')
        #jj = ff.read()
        #ff.close()
        #para_dict = json.loads(jj) # convert string to dictionary
        # print('Json interpret numbered key as string. Problem!  ')
        #>> using pickle
        ff= open(fileName,'rb')
        para_dict = pickle.load(ff)
        ff.close()
        self.load_para(para_dict)
        self.setCurrentIndex(0)
        self.setCurrentIndex(1)
        self.setCurrentIndex(2)
        return  para_dict

#============================================================================
class MyWindow(QMainWindow, uic.loadUiType("main_window.ui")[0]):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        size = QSize(750,700)
        self.resize(size)
        self.setWindowTitle('GIRL: Gui for Intuitive Reaction Learning')
        # data for main window
        self.path_data = myutil.all_global_variables()
        self.load_path_info()
        #---connect Menu action
        self.actionSave.triggered.connect(self.save_file)
        self.actionOpen.triggered.connect(self.open_file)
        self.actionPath_fresco.triggered.connect(lambda: self.set_path('fresco_filename'))
        self.actionPath_sfresco.triggered.connect(lambda: self.set_path('sfresco_filename'))
        self.actionPath_omget.triggered.connect(lambda: self.set_path('omegt_filename'))
        self.actionAbout.triggered.connect(self.show_about)
        self.actionDocumentation.triggered.connect(self.show_documents)
        self.actionBug_Report.triggered.connect(self.show_bugreport)
        #---add Widgets
        self.DWBA = TestToolBox()
        self.verticalLayout.addWidget(self.DWBA)
        
    def show_documents(self,):
        external_web_browser(url ='https://github.com/alphacentaury-github/ReactionGUI')
        
    def show_bugreport(self,):
        external_web_browser(url ='https://github.com/alphacentaury-github/ReactionGUI/issues')

    def load_path_info(self,):
        # try load Fresco path
        try:
            self.path_data.load_from_file()
        except:
            # if file does not exist, create
            self.path_data.save_to_file()

    def set_path(self,option='fresco_filename'):
        """
        set pth for fresco for each option case
        possibly later change it to cover other path
        """
        options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,
                      "QFileDialog.getOpenFileName()",
                      "","All Files (*);;Python Files (*.py)",
                      options=options)
        if fileName:
            self.path_data.set_global_value(**{ option:  fileName})
            self.path_data.save_to_file()

    def save_file(self,):
        """
        save information of the window in a file
        exact action is not decided..

        """
        options = QFileDialog.Options()
        fileName, _filter = QFileDialog.getSaveFileName(self,
                    "Save file",
                    "",
                    "All Files (*)",
                    options=options)
        para_dict   = self.DWBA.save_para()
        if fileName:
            #>> pickle 
            ff= open(fileName,'wb')
            pickle.dump(para_dict,ff)
            ff.close()
            #>> json 
            #jj = json.dumps(para_dict,indent=2)
            #ff = open(fileName,'w')
            #ff.write(jj)
            #ff.close() 
            return para_dict
        else :
            print('No file is chosen')
            return

    def open_file(self,):
        """
        open previous saved file

        """
        options = QFileDialog.Options()
        fileName, _filter = QFileDialog.getOpenFileName(self,
                    "Load file",
                    "",
                    "All Files (*)",
                    options=options)

        if fileName:
            #>> using pickle
            ff= open(fileName,'rb')
            para_dict = pickle.load(ff)
            ff.close()
            #>> using json 
            #ff = open(fileName,'r')
            #jj = ff.read()
            #ff.close()
            #para_dict = json.loads(jj) # convert string to dictionary
        
            self.DWBA.load_para(para_dict)
            self.DWBA.setCurrentIndex(0)
            self.DWBA.setCurrentIndex(1)
            self.DWBA.setCurrentIndex(2)
            return  para_dict
        else:
            print('No file is chosen')
            return

    def show_about(self,):
        # testing new application
        dlg = about_Dialog()
        dlg.exec_()
        return
#-----external webbrowser
def external_web_browser(url = 'html\description_transfer.html'):
    # open html file in external web-browser. 
    import webbrowser
    webbrowser.open(url,new=2)
    return 


#============================================================================
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
        #myWindow = TestToolBox()
        myWindow = MyWindow()
        myWindow.show()
        #sys.exit(app.exec_() )
        return myWindow
    m = run_app()
