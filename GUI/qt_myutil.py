# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 13:23:10 2020

@author: Y.-H. Song

Collection of useful class/functions using Qt
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
import re 
#from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import (QApplication,QMainWindow,QDialog,QFileDialog,QWidget,
                             QDialogButtonBox,
                             QComboBox,QLabel,QLineEdit,QCheckBox,
                             QHBoxLayout,QVBoxLayout,QGridLayout,
                             QGroupBox,QToolBox,QTabWidget,QPushButton,QTextBrowser,
                             QPlainTextEdit,QProgressBar,
                             QSpacerItem,QRadioButton)
from PyQt5 import uic
from PyQt5 import QtCore
from PyQt5.QtCore import (QProcess,QByteArray,QUrl) 

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
        FigureCanvasQTAgg as FigureCanvas,
        NavigationToolbar2QT as NavigationToolbar)

from numpy import loadtxt

from io import StringIO
from subprocess import (call,Popen)

import numpy as np
import myutil
import reactions
import run_fresco_v2

html_head = """
            <html><head>
             <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
             </script></head>
             <body>
             """
html_tail =  "</body></html>"

#------------------------------------------------------------------------------------------------
class combined_Widgets_horizontal(QWidget,):
    """
    place list of Widgets Horizontally and
    treat the collection as a Widget.
    
    one can either use 'json' or 'pickle' form
    to save dictionary.
    (difference is that json use string index
     while pickle can use number index)
    """
    def __init__(self,list_of_Widgets,opt='pickle'):
        super().__init__()
        self.layout = QHBoxLayout()
        self.setLayout(self.layout)
        self.save_type = opt #pickle or json 

        self.list_of_Widgets = list_of_Widgets
        for i in list_of_Widgets :
            self.layout.addWidget(i)
            
    def change_opt(self,opt='pickle'):
        self.save_type = opt 

    def get_Widget(self,index):
        return self.list_of_Widgets[index]

    def set_signal_solt_Widget(self,index,signal='',slot_function=''):
        # Is it necessary ? Actually, one can directly access a Widget
        # from outside by using x.list_of_Widgets[index]
        signal_to_call = getattr(self.list_of_Widgets[index],signal)
        signal_to_call.connect(slot_function)

    def set_attribute_Widget(self,index,attribute='',value=''):
        # Is it necessary ? Actually, one can directly access a Widget
        # from outside by using x.list_of_Widgets[index]
        widget_attribute= getattr(self.list_of_Widgets[index],attribute)
        widget_attribute(value)

    def get_value_Widget(self,index,attribute=''):
        # Is it necessary ? Actually, one can directly access a Widget
        # from outside by using x.list_of_Widgets[index]
        widget_attribute= getattr(self.list_of_Widgets[index],attribute)
        return widget_attribute()

    def get_values(self,):
        self.data = {}
        for i in range(len(self.list_of_Widgets)):
            widg = self.list_of_Widgets[i]
            str_i = str(i)
            #>> using json, string keys 
            if self.save_type=='json':
                index = str_i  
            elif self.save_type=='pickle':
                index = i
            if isinstance(widg, QLabel):
                self.data[index] = widg.text()
            elif isinstance(widg, QLineEdit):
                self.data[index] = widg.text()
            elif isinstance(widg, QComboBox):
                self.data[index] = widg.currentIndex()
            elif isinstance(widg, QRadioButton):
                self.data[index] = widg.isChecked()
            elif isinstance(widg, QCheckBox):
                self.data[index] = widg.isChecked()     
            else: #other widgets are skipped
                pass 
        return self.data

    def put_values(self,data):
        # inverse of get_values
        # key is number with get_values. 0, 1
        # however, json makes the key as string, '0','1'
        for i in range(len(self.list_of_Widgets)):
            widg = self.list_of_Widgets[i]
            str_i = str(i)
            # >> json case # >> pickle case 
            if self.save_type=='json':
                index = str_i 
            elif self.save_type=='pickle':
                index = i 
            if isinstance(widg, QLabel):
                widg.setText(data[index]) # all data are strings
            elif isinstance(widg, QLineEdit):
                widg.setText(data[index])
            elif isinstance(widg, QComboBox):
                widg.setCurrentIndex(data[index])
            elif isinstance(widg,QRadioButton):
                widg.setChecked(data[index])
            elif isinstance(widg, QCheckBox):
                widg.setChecked(data[index])
            else: #other widgets are skipped
                pass

#--------------------------------------------------------------------------------------
class combined_Widgets_vertical(QWidget,):
    """
    place list of Widgets Horizontally and
    treat the collection as a Widget.
    """
    def __init__(self,list_of_Widgets,opt='pickle'):
        super().__init__()
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.list_of_Widgets = list_of_Widgets
        self.save_type= opt 
        for i in list_of_Widgets :
            self.layout.addWidget(i)
            
    def change_opt(self,opt='pickle'):
        self.save_type = opt 
        
    def get_Widget(self,index):
        return self.list_of_Widgets[index]

    def set_signal_slot_Widget(self,index,signal='',slot_function=''):
        signal_to_call = getattr(self.list_of_Widgets[index],signal)
        signal_to_call.connect(slot_function)

    def get_values(self,):
        self.data = {}
        for i in range(len(self.list_of_Widgets)):
            widg = self.list_of_Widgets[i]
            str_i = str(i)
            # pickle or json 
            if self.save_type=='json':
                index = str_i 
            elif self.save_type=='pickle':
                index = i 

            if isinstance(widg, QLabel):
                self.data[index] = widg.text()
            elif isinstance(widg, QLineEdit):
                self.data[index] = widg.text()
            elif isinstance(widg, QComboBox):
                self.data[index] = widg.currentIndex()
            elif isinstance(widg, QRadioButton):
                self.data[index] = widg.isChecked()  
            elif isinstance(widg, QCheckBox):
                self.data[index] = widg.isChecked()    
            else: #other widgets are skipped
                pass
        return self.data

    def put_values(self,data):
        # inverse of get_values
        for i in range(len(self.list_of_Widgets)):
            str_i = str(i)
            widg = self.list_of_Widgets[i]
            # pickle or json 
            if self.save_type=='json':
                index = str_i 
            elif self.save_type=='pickle':
                index = i 
            
            if isinstance(widg, QLabel):
                widg.setText(str(data[index]))
            elif isinstance(widg, QLineEdit):
                widg.setText(str(data[index]))
            elif isinstance(widg, QComboBox):
                widg.setCurrentIndex(data[index])
            elif isinstance(widg, QRadioButton):
                widg.setChecked(data[index])
            elif isinstance(widg, QCheckBox):
                widg.setChecked(data[index])
            else: #other widgets are skipped
                pass


#------------------------------------------------------------------------------
class combined_Widgets_grid(QWidget,):
    def __init__(self,list_of_Widgets,opt='pickle'):
        super().__init__()
        self.layout = QGridLayout()
        self.setLayout(self.layout)

        self.list_of_Widgets = list_of_Widgets
        self.save_type = opt 
        for i in range(len(self.list_of_Widgets)): #row
            for j in range(len(self.list_of_Widgets[i] )): #column
                self.layout.addWidget(self.list_of_Widgets[i][j],i,j)
                
    def change_opt(self,opt='pickle'):
        self.save_type = opt 
                    
    def get_Widget(self,row,column):
        return self.list_of_Widgets[row][column]

    def set_signal_slot_Widget(self,row,column,signal='',slot_function=''):
        signal_to_call = getattr(self.list_of_Widgets[row][column],signal)
        signal_to_call.connect(slot_function)

    def get_values(self,):
        """
        return dictionary

        """
        self.data ={}
        for row in range(len(self.list_of_Widgets)): #row
            for col in range(len(self.list_of_Widgets[row] )): #column
                widg = self.list_of_Widgets[row][col]
                str_row = str(row)
                str_col = str(col)
                # pickle or json 
                if self.save_type=='json':
                  row_index = str_row 
                  col_index = str_col  
                elif self.save_type=='pickle':
                  row_index = row 
                  col_index = col 
                
                if isinstance(widg, QLabel):
                    self.data[(row_index,col_index)] = widg.text()
                elif isinstance(widg, QLineEdit):
                    self.data[(row_index,col_index)] = widg.text()
                elif isinstance(widg, QComboBox):
                    self.data[(row_index,col_index)] = widg.currentIndex()
                elif isinstance(widg, QRadioButton):
                    self.data[(row_index,col_index)] = widg.isChecked() 
                elif isinstance(widg, QCheckBox):
                    self.data[(row_index,col_index)] = widg.isChecked()     
                else: #other widgets are skipped
                    pass
        return self.data

    def put_values(self,data):
        for row in range(len(self.list_of_Widgets)): #row
            for col in range(len(self.list_of_Widgets[row] )): #column
                widg = self.list_of_Widgets[row][col]
                str_row = str(row)
                str_col = str(col)
                # pickle or json 
                if self.save_type=='json':
                  row_index = str_row 
                  col_index = str_col  
                elif self.save_type=='pickle':
                  row_index = row 
                  col_index = col 
                
                if isinstance(widg, QLabel):
                    widg.setText(data[(row_index,col_index)])
                elif isinstance(widg, QLineEdit):
                    widg.setText(data[(row_index,col_index)])
                elif isinstance(widg, QComboBox):
                    widg.setCurrentIndex(data[(row_index,col_index)])
                elif isinstance(widg, QRadioButton):
                    widg.setChecked(data[(row_index,col_index)])
                elif isinstance(widg, QCheckBox):
                    widg.setChecked(data[(row_index,col_index)])
                else: #other widgets are skipped
                    pass
        return

#----------------------------------------------------------------
class QLabel_aligned(QLabel,):
    def __init__(self,label_txt,align='c'):
        super().__init__()

        self.setText(label_txt)
        if align =='r':
            self.setAlignment(QtCore.Qt.AlignRight )
        elif align=='c':
            self.setAlignment(QtCore.Qt.AlignCenter)
        elif align=='l':
            self.setAlignment(QtCore.Qt.AlignLeft)

#-------------------------------------------------------------------
class text_Browser(QWidget,):
    """
    Change between Text Browser and Web Browser
    """
    def __init__(self,browser_choice='text',initial_page=''):
        """

        browser_choice =='text' for text_browser
                               but this can also show simple Html text
                       =='web'  for web browser

        set_page_txt puts text either plain string text
        or HTML format

        """
        super().__init__()
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.browser_choice = browser_choice
        self.page_txt = initial_page
        if self.browser_choice in ['text','html']:
            self.browser = QTextBrowser()
            self.browser.setAcceptRichText(True)
            self.browser.setOpenExternalLinks(True)
        elif self.browser_choice =='web':
            #self.browser = QWebEngineView()
            pass
        #---initial view
        self.layout.addWidget(self.browser)
        self.set_page_txt(self.page_txt)

    def set_page_txt(self,page_txt,txt_type='text'):
        self.page_txt = page_txt
        if self.browser_choice =='text':
            if txt_type=='text':
                self.browser.setPlainText(page_txt)
            elif txt_type=='html':
                self.browser.setHtml(page_txt)
        elif self.browser_choice =='web':
            self.browser.setHtml(page_txt)
#------------------------------------------------------------------------
class about_Dialog(QDialog,uic.loadUiType("dialog_about.ui")[0]):
    def __init__(self,label_text=None,text_type='html'):
        super().__init__()
        self.setupUi(self)
        if label_text and text_type=='html':
            self.textBrowser.setHtml(label_text)
        elif label_text and text_type=='text':
            self.textBrowser.setPlainText(label_text)


#-------------------------------------------------------------------------
class WidgetMatplot(QWidget,):
    def __init__(self,):
        super().__init__()
        self.layout =  QVBoxLayout()
        self.setLayout(self.layout)
        # defualt figure
        self.figure_default = Figure()
        ax = self.figure_default.add_subplot(111)
        x = np.arange(0.,2*np.pi,0.1)
        ax.plot(x,np.sin(x))
        ax.text( 0.5, 0.0 ,'This is a plot area')
        # setup Widgets 
        self.add_plot(self.figure_default)
        #self.canvas = FigureCanvas(self.figure_default)
        #self.layout.addWidget(self.canvas)
        #self.toolbar = NavigationToolbar(self.canvas,
        #        self, coordinates =True)
        #self.layout.addWidget(self.toolbar)
        self.axes = None  #this is defined when figure/plot is added 
        
    def draw_axes(self,x,y,label=None,fmt=None):
        """
        Instead of replace whole figure,
        add a line in existing axes  
        
        Seems to be not working !! 
        """
        if fmt:
            self.axes.plot(x,y,fmt,label=label)
        else: 
            self.axes.plot(x,y,label=label)
        
    def reset_plot(self,):
        # show default figure 
        self.rm_plot() 
        self.add_plot(self.figure_default)
   
    def add_plot(self,fig):
        # add a figure in the canvas.  
        # figure is external
        # add canvas and toolbar to layout
        # To replace it one have to remove previous one first. 
        self.axes = fig.get_axes()[0] # suppose only one axes..
        self.canvas = FigureCanvas(fig)
        self.layout.addWidget(self.canvas) # add a canvas widget into layout
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas,
                self, coordinates =True)
        self.layout.addWidget(self.toolbar)

    def rm_plot(self,):
        # remove canvas and toolbar
        self.layout.removeWidget(self.canvas)
        self.canvas.close()
        #self.canvas.setParent(None)
        self.layout.removeWidget(self.toolbar)
        self.toolbar.close()
        #self.toolbar.setParent(None)
        
    def add_plot_by_data(self,list_of_data_to_plot,
                         show_grid=False,yscale='linear',
                         xlabel='c.m. angle[deg]',
                         ylabel='mb'):
        """
        Parameters
        ----------
        list_of_data_to_plot : list of dictionary to plot
            [plot1={'x': array ,'y':array,'yerr': array,
                    'label' : legend txt, 'fmt': matplotlib fmt},
             plot2={},..]


        show_grid : TYPE, optional
            DESCRIPTION. The default is False.
        yscale : TYPE, optional
            DESCRIPTION. The default is 'linear'.
        xlabel : TYPE, optional
            DESCRIPTION. The default is 'c.m. angle[deg]'.
        ylabel : TYPE, optional
            DESCRIPTION. The default is 'mb'.

        Returns
        -------
        None.

        """
        # remove previous plot
        self.rm_plot()
        fig = Figure()
        ax = fig.add_subplot(111)
        for p in list_of_data_to_plot:
            try: # check whether a yerr exits in data
                yerr = p['yerr']
                plot_err= True
            except:
                plot_err = False
            if plot_err: #error exist
                ax.errorbar(p['x'],p['y'],yerr=yerr,fmt=p['fmt'],label=p['label']  )
            else :
                ax.plot(p['x'],p['y'],p['fmt'],label=p['label'])
        ax.legend()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_yscale(yscale)
        if show_grid :
            ax.grid()
        self.add_plot(fig)
        return
    
#============================================================================    
class Widget_external(QWidget):
    # Widget to call external command line program 
    #  for example: Widget_external("python",["dummy_script.py"])
    """
     At the moment, setWorkingDirectory seems to be not working... 
     
     before_run is a function to be done before running program 
     after_run is a function to be done after running program
     
     one can send signal to external program 
     (1) interactive
     (2) series of inputs prepared 
         with list_sequence 
    """
    def __init__(self,program_name,argument_list=[],
                 before_run =None,
                 after_run = None,  
                 working_directory=None,
                 list_sequence=None):
        super().__init__()
        self.btn = QPushButton("Execute")
        self.btn.pressed.connect(self.start_process)
        # for series of inputs 
        self.list_sequence = None
        self.multi_run = None # if sequence requires restart
        self.set_sequence(list_sequence)
        self.ans = None 
         
        self.text = QPlainTextEdit()
        self.text.setReadOnly(True)
        self.input = QLineEdit() 
        
        self.enter_input = QPushButton('Enter') 
        self.enter_input.clicked.connect(self.handle_input) 
        answer_area=combined_Widgets_horizontal([QLabel('answer:'),
                                                 self.input,
                                               self.enter_input])

        l = QVBoxLayout()
        l.addWidget(self.btn)
        l.addWidget(self.text)
        l.addWidget(answer_area)
        
        self.setLayout(l)

        self.p = None 
        self.external_program = program_name
        self.external_argument_list = argument_list
        self.working_directory = working_directory
        self.before_run = before_run 
        self.after_run = after_run   
            
    def message(self,s):
        self.text.appendPlainText(s)
        
    def set_sequence(self,list_sequence):    
        # if list is a nested list of series of inputs [['1','2']]
        # set it as a multi_run 
        self.list_sequence = list_sequence
        if list_sequence:
          if any(isinstance(i, list) for i in list_sequence):
              self.multi_run = True  
          else:
              self.multi_run = False 
                
    def start_process(self,):
        # note that though with setWorkingDirectory 
        # external program name should contain path.
        self.text.clear()
        # prepare run 
        if self.before_run :
            try:
                self.before_run() 
            except: 
                self.text.appendPlainText(
                        'Error in preparation/input\n')
                return
        # start sequence     
        if self.list_sequence is None: # normal case
            self.run_process()  
        elif self.multi_run : # multiple run with answers
            for tasks in self.list_sequence:
                if self.p is None:
                    self.run_process()
                    if self.p.state()==2:
                        for ans in tasks:
                            self.ans = ans
                            self.handle_input() 
                            self.ans = None
                    self.p.waitForFinished()                    
        elif self.list_sequence: # a single run with answers      
            self.run_process()
            if self.p.state()==2:
                for ans in self.list_sequence: 
                    self.ans = ans 
                    self.handle_input()
                    self.ans = None # reset 
            self.p.waitForFinished()        
        # all process done  
        if self.after_run:
            try:
                self.after_run()
            except:
                self.text.appendPlainText(
                        'Error in calculation/output\n')
                return
        
    def run_process(self,):
        # note that though with setWorkingDirectory 
        # external program name should contain path.          
        if self.p is None: # No process is running 
            #----connect slots     
            self.message("Execute process.")
            self.p = QProcess()
            self.p.setProgram(self.external_program)
            self.p.setArguments(self.external_argument_list)
            if self.working_directory:
                self.p.setWorkingDirectory(self.working_directory)
            #---signals/slots     
            self.p.readyReadStandardOutput.connect(self.handle_stdout)
            self.p.readyReadStandardError.connect(self.handle_stderr)
            self.p.stateChanged.connect(self.handle_state) 
            self.p.finished.connect(self.process_finished)
            self.p.start() 
                            
    def process_finished(self):
        self.message("================Process finished======================")
        self.p = None 
               
    def handle_stderr(self):
        # get stderr message of external program 
        # Because readAllStandardError() return 
        # data as bytes,wrapped in a Qt object 
        # thus, convert it into python bytes() object
        # and then decode into a string 
        data = self.p.readAllStandardError()
        stderr = bytes(data).decode("utf8")
        self.message(stderr)            
        
    def handle_stdout(self):
        data = self.p.readAllStandardOutput()
        stdout = bytes(data).decode("utf8")
        self.message(stdout)       
        
    def handle_state(self,state):
        """
        state is signal of stateChanged
        
        QProcess.NotRunning	0	The process is not running.
        QProcess.Starting	1	The process is starting, but the program has not yet been invoked.
        QProcess.Running	2	The process is running and is ready for reading and writing.
        """
        states = {   
            QProcess.NotRunning: 'Not running',
            QProcess.Starting: 'Starting',
            QProcess.Running: 'Running',
        }
        state_name = states[state]
        self.message(f"State changed: {state_name}")
    
    def process_state(self,):    
        if self.p:
            return self.p.state() 
        else:
            return -1 # no process 
        
    def handle_input(self,):
        # send input to external program 
        # Note that Qprocess.write() takes bytes input
        if self.ans is None:    
            ans = self.input.text()
        else:
            ans = self.ans         
        a = QByteArray()
        a.append(ans+'\n') # input must end with newline
        if self.p : #process is not finished
            if self.p.state()==2: #process is ready
                #self.handle_stdout()
                self.p.write(a)
                self.message(ans+' is sent to program') 
                #--clear entered input  
                self.input.clear()
