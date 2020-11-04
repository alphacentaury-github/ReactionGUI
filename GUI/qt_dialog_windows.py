# -*- coding: utf-8 -*-
"""
Created on Mon May 18 10:54:34 2020

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

from PyQt5.QtWidgets import (QApplication,QMainWindow,QDialog,QFileDialog,QWidget)
from PyQt5 import uic
#from PyQt5.QtCore import *
from subprocess import (call,Popen)

from io import StringIO
import numpy as np
import myutil 
import reactions
import run_fresco_v2 

class data_Dialog(QDialog, uic.loadUiType("dialog_data.ui")[0]):    
    def __init__(self,data=None,label_text=None):
        super().__init__()
        self.setupUi(self)
        self.data = data 
        self.accepted = False 
        self.pushButton_2.clicked.connect(self.clear_data)
        self.buttonBox.accepted.connect(self.take_data)
        self.buttonBox.rejected.connect(self.reject)
        self.plainTextEdit.setPlainText(self.data)
        self.label.setText(label_text)   
        
    def take_data(self) :
        # make the text ends with \n 
        mytext = (self.plainTextEdit.toPlainText()).strip()+'\n'
        self.data = mytext
        self.accepted= True 
        self.accept()

    def clear_data(self):
        self.plainTextEdit.clear()
        
    def get_numpy_array(self,):
        """
        convert text into numpy array 
        """
        ff = StringIO(self.data)
        numpy_array = np.loadtxt(ff)
        return numpy_array

#---------------------------------------------------------------------------------------
class ExpDataDialog(QDialog,uic.loadUiType("dialog_exp_data.ui")[0]):    
    """
    .data = input string

    data is supposed to be angle mb/sr error 
    if no error is included... 
    """
    def __init__(self,data=None,
                 label_text="Type of Data and Error",
                 SFresco_Input_object=None):
        super().__init__()
        if SFresco_Input_object:
            self.SFresco_input = SFresco_Input_object
        else: 
            self.SFresco_input = run_fresco_v2.SFresco_Input()
        self.setupUi(self)
        if data:
            self.data = data # string 
            self.plainTextEdit.setPlainText(self.data) #show data 
        else: 
            self.data = self.plainTextEdit.toPlainText() 
        # take care of the case without error     
        self.data = myutil.add_error_to_datatext(self.data)
        self.data_full_info ={} # full informations 
        
        self.pushButton_clear.clicked.connect(self.clear_data)
        self.label.setText(label_text)   
        
        self.buttonBox.accepted.connect(self.take_data)
        self.buttonBox.rejected.connect(self.reject_data)
   
    def modify(self,data='',
                  label_text="Type of Data and Error",
                  SFresco_Input_object=None):  
        self.data= data 
        self.data= myutil.add_error_to_datatext(self.data)
        self.plainTextEdit.setPlainText(self.data)
        self.data_full_info ={} # full informations 
        self.label.setText(label_text)   
        self.SFresco_input = SFresco_Input_object

    def update_input(self,): 
        """
         update inputs in the Dialog  
        """
        
        if self.comboBox_LAB.currentIndex()==0:
            lab  =  'T'
        else: 
            lab  =  'F'
        if self.comboBox_IDIR_ISCALE.currentIndex()==0: # ratio 
            idir = 1 
            iscale = -1
        else : 
            idir = 0  # absolute 
            iscale = self.comboBox_IDIR_ISCALE.currentIndex()-1
        if self.comboBox_ABSERR.currentIndex()==0:    
            abserr = 'T'
        else: 
            abserr = 'F'
            
        self.SFresco_input.set_data_info(TYPE= self.comboBox_TYPE.currentIndex(),
                                         LAB =  lab,  
                                         IDIR = idir,
                                         ISCALE = iscale, 
                                         ABSERR = abserr )  
        #---How to use this....
        self.data_full_info ={'LAB': self.comboBox_LAB.currentIndex(),
                      'IDIR': self.comboBox_IDIR_ISCALE.currentIndex(),
                      'ABSERR': self.comboBox_ABSERR.currentIndex()}
        #---
        #----change the text to end with \n 
        self.data = (self.plainTextEdit.toPlainText()).strip()+'\n' 
        data = myutil.txt_to_array(self.data)
        # take care of the case without error     
        if data.shape[0]==0: # no data 
            self.SFresco_input.set_exp_data(self.data)
            return 
        elif data.shape[1]==1: # angle only. absurd 
            self.SFresco_input.set_exp_data(self.data)
            return 
        elif data.shape[1]==2:   # data but no error bar  
            # add error to data 
            self.data = myutil.add_error_to_datatext(self.data)
            self.SFresco_input.set_exp_data(self.data)
            return 
        else : # all good ? 
            self.SFresco_input.set_exp_data(self.data)
            return 
    def convert_data(self,):
        # need to store setting of combobox and 
        # convert them into correct ones 
        return 
        
    def take_data(self) :        
        self.update_input() 
        self.accept() # return 1
        
    def reject_data(self,):
        self.data = None
        self.reject()  # return 0 

    def clear_data(self):
        self.plainTextEdit.clear()

#------------------------------------------------------------------------------------
class Elastic_Fit_Dialog(QDialog,uic.loadUiType("dialog_fit.ui")[0]):
    def __init__(self,fresco_input_object=None,exp_data=None,
                 sfresco_input_object=None,path_data=None ):
        """
        fresco_input_object : run_fresco_v2.fresco_input object 
        exp_data : exp_data in text string
        sfresco_input_object: run_fresco_v2.sfresco_input object 
        path_data:  dictionary with 
              path_data['frin_filename']
              path_data['frout_filename']
              path_data['search_filename']
              path_data['minunit_filename']
        """
        super().__init__()
        self.setupUi(self)
        if fresco_input_object:
            self.input = fresco_input_object 
        else: 
            self.input = run_fresco_v2.Fresco_input()
        if sfresco_input_object:
            self.sfresco_input = sfresco_input_object
        else: 
            self.sfresco_input = run_fresco_v2.SFresco_Input()
            
        self.exp_data = exp_data # text   
        self.path_data = path_data 
        self.fit_result ={} 
         #----check buttons      
         # I need a better way to organize these...     
        self.check_buttons = [self.checkBox_V0,self.checkBox_r0, self.checkBox_a0,
                              self.checkBox_W0,self.checkBox_rw,self.checkBox_aw,
                              self.checkBox_Wd,self.checkBox_rwd,self.checkBox_awd,
                              self.checkBox_Vso,self.checkBox_rso,self.checkBox_aso
                              ]
         # each one [kind,name,kp,pline, col]
         # Need a better way ... 
        self.check_buttons_associated = [ [1,'v0',1,2,1],[1,'r0',1,2,2],[1,'a0',1,2,3],
                                          [1,'w0',1,2,4],[1,'rw',1,2,5],[1,'aw',1,2,6],
                                          [1,'wd',1,3,4],[1,'rwd',1,3,5],[1,'awd',1,3,6],
                                          [1,'vso',1,4,1],[1,'rso',1,4,2],[1,'rso',1,4,3]
                                          ] 
        #---connect  
        self.pushButton.clicked.connect(self.update_and_run)
        self.buttonBox.accepted.connect(self.take_data)
        self.buttonBox.rejected.connect(self.reject_data)
        
    def update_and_run(self, ):
        """
        May need to separate into several steps 
        
        update input
        prepare input files 
        run Sfresco 
        read result of fitting 

        """        
        # updated check box values 
        self.sfresco_input.reset_search_variables() 
        self.textBrowser.clear() 
        for i in range(len(self.check_buttons)) : 
            if self.check_buttons[i].isChecked(): 
                self.sfresco_input.add_search_variable(
                    KIND=self.check_buttons_associated[i][0],
                    NAME=self.check_buttons_associated[i][1],
                    KP=self.check_buttons_associated[i][2],
                    PLINE=self.check_buttons_associated[i][3],
                    COL= self.check_buttons_associated[i][4] 
                    )
                 
        # Prepare Sfresco input files 
        frin_fname = self.path_data['frin_filename']
        frout_fname = self.path_data['frout_filename']
        search_fname = self.path_data['search_filename']
        min_fname = self.path_data['minunit_filename']
        fname_head = min_fname.split('.')[0]
        fit_result_filename = fname_head+'_fit.plot'
        
        #----prepare file 
        if len(self.sfresco_input.search_var)==0: # nothing is checked
            self.textBrowser.append('Fitting variables are not chosen.\n')
            return 
        #---- if checked.. 
        #--remove previous files 
        try:
            os.remove(search_fname)
            os.remove(min_fname)
            os.remove(fit_result_filename)
        except:
            print('Error deleting search/minuit files ')

        # prepare FR input file 
        ff=open(frin_fname,'w')
        ff.write(self.input.write() )
        ff.close() 
        #----prepare minuit file 
        self.sfresco_input.write_minuit(
            fname=min_fname,search_file=search_fname)
        # prepare search file 
        txt = self.sfresco_input.write_search( 
            search_fname= search_fname,
            frin_fname= frin_fname,
            frout_fname = frout_fname )
        self.textBrowser.append(txt) # for test, remove later 
        if self.sfresco_input.data=='':
            self.textBrowser.append('No data! stop\n')
            return 
        #----run sfresco 
        try: 
            # Let us use Popen instead of call 
            #call(self.path_data['sfresco_filename']+" < "+ min_fname ,shell=True)
            proc = Popen(self.path_data['sfresco_filename']+" < "+ min_fname ,shell=True)
            proc.wait()
            proc.terminate() 
        except: 
            self.textBrowser.append('Something wrong while running Sfresco \n')
        #----read results of fitting 
        # check result?? 
        ff=open(fit_result_filename,'r')
        ll=ff.readlines()
        ff.close() 
        
        # Show result.. 
        self.textBrowser.append('==Finished SFresco fitting\n')
        
        # get fitted parameter values 
        self.fit_result['para']=[]  
        self.fit_result['text']=''
        for i in range(len(self.sfresco_input.search_var) ):
            ll[i]=ll[i].replace(',',' ') # remove , 
            ww= ll[i].split() 
            if 'Var'==ww[1]:
                #print(ww[4])
                self.fit_result['para'].append(float(ww[4])) 
                self.fit_result['text'] += ' {} ={}\n'.format(
                    self.sfresco_input.search_var[i]['NAME'],
                    float(ww[4]) )
                
                self.textBrowser.append(' {} ={}'.format(
                    self.sfresco_input.search_var[i]['NAME'],
                    float(ww[4]) ) )
        # get fitted plot data           
        myutil.clean_comm(fit_result_filename)
        self.fit_result['plot'] = myutil.read_fresco_res(fit_result_filename+'x')
        # convert to array 
        self.fit_result['plot'][0] = np.array(self.fit_result['plot'][0])
        self.fit_result['plot'][1] = np.array(self.fit_result['plot'][1])
        return self.fit_result    
    
    #self.buttonBox.accepted.connect(self.take_data)
    #self.buttonBox.rejected.connect(self.reject_data)
    def take_data(self,):        
        self.accept() #return 1
    def reject_data(self,):
        self.fit_result={} 
        self.reject() #return 0 
        
        
        
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
        myWindow = data_Dialog()  
        #myWindow = ExpDataDialog() 
        myWindow.show()
        # app.exec_()
        return myWindow  
    m = run_app() 
