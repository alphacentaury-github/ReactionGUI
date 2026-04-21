# -*- coding: utf-8 -*-
"""
Created on Feb 2024

@author: Y.-H. song

SIMPLE GUI using FRESCO input file and plot results.
Coupled Channel calculation is possible.  


For later use: plot update can be done
 (1) clearing and redrawing
     self.canvas.axes.cla() 
     self.canvas.axes.plot()
     self.canvas.draw() 
 (2) keep a reference to the plotted line and update the data 
        (for multiple lines , list or dictionary to store data )
        # add points to the end of old data
        self.ydata = self.ydata[1:] + [random.randint(0, 10)]
        # without clear axis 
        if self._plot_ref is None:
            # First time we have no plot reference, so do a normal plot.
            # .plot returns a list of line <reference>s, as we're
            # only getting one we can take the first element.
            plot_refs = self.canvas.axes.plot(self.xdata, self.ydata, 'r')
            self._plot_ref = plot_refs[0] # here _plot_ref is a variable to hold to the plotted line 
        else:
            # We have a reference, we can use it to update the data for that line.
            self._plot_ref.set_ydata(self.ydata)

        # Trigger the canvas to update and redraw.
        self.canvas.draw() 
    
"""
import sys
import os 
from io import StringIO
#---current path where this file resides 
try:
    here = os.path.dirname(os.path.realpath(__file__))
except:
    here = '.'
sys.path.insert(0,here) # to import DFPot3 in this directory 

from PyQt5.QtWidgets import (QApplication,QMainWindow,QDialog,QFileDialog
                             ,QWidget,QComboBox,QVBoxLayout,QLineEdit,QFormLayout)
from PyQt5.QtWidgets import QPushButton
 
from PyQt5 import uic
from PyQt5 import QtCore, QtWidgets

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
        FigureCanvasQTAgg as FigureCanvas,
        NavigationToolbar2QT as NavigationToolbar)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

from subprocess import (call,Popen)
import numpy as np
from scipy.optimize import curve_fit 

import matplotlib.pyplot as plt
import sys
import matplotlib

import myutil
import reactions
import run_fresco_v2
import qt_dialog_windows

matplotlib.use('Qt5Agg')

form_omp = uic.loadUiType(here+'/'+"qt_fresco.ui")[0]

element_names = ["n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl",
		 "Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
		 "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb",
		 "Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er",
		 "Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At",
		 "Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No",
		 "Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og",
		 "119","120","121","122","123","124","125","126","127","128","129","130"]
          
def F_Gaussian(x,A,P0,sigma):
    # A*exp(-(x-P0)^2/(2*sigma^2))
    # x, P0, sigma are all MeV units 
    return A*np.exp(-(x-P0)**2/2.0/sigma**2)

def F_AsymGaussian(x,A,P0,sig_L,sig_R):
    # asymmetric Gaussian
    # A*exp(-(x-P0)^2/(2*sig_L^2)) if x < P0
    # A*exp(-(x-P0)^2/(2*sig_R^2)) if x >= P0
    return (A*np.exp(-(x-P0)**2/(2.0*sig_L**2))*(x-P0 < 0.0) 
            +A*np.exp(-(x-P0)**2/(2.0*sig_R**2))*(x-P0 >= 0.0) )

def F_CauchyHazard(x,A,P0,a,b):
    #-----modified CauchyHazard function 
    #     for the slow tail on left. 
    x= (x-P0)/b; #shift centre 
    return A/(1+x**2)/(0.5*np.pi*a - np.arctan(-x))

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)

class addData_Dialog(QDialog, uic.loadUiType("dialog_data.ui")[0]):    
    def __init__(self,data=None,label_text=None):
        super().__init__()
        self.setupUi(self)
        self.data = data #text data 
        self.accepted = False 
        self.exp_data_dict= {} # dictionary  
        #---add 3 buttons
        self.pushButton_exp0 = QPushButton('set data as 0')
        self.pushButton_exp1 = QPushButton('set data as 1')
        self.pushButton_exp2 = QPushButton('set data as 2')
        self.pushButton_cla = QPushButton('Reset all data')
        
        self.gridLayout.addWidget(self.pushButton_exp0,2,0)
        self.gridLayout.addWidget(self.pushButton_exp1,3,0)
        self.gridLayout.addWidget(self.pushButton_exp2,4,0)
        self.gridLayout.addWidget(self.pushButton_cla,4,1)
                
        self.pushButton_2.clicked.connect(self.clear_data)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.pushButton_exp0.clicked.connect(lambda : self.add_data(0))
        self.pushButton_exp1.clicked.connect(lambda : self.add_data(1))
        self.pushButton_exp2.clicked.connect(lambda : self.add_data(2))
        self.pushButton_cla.clicked.connect(self.cla)
        
        self.plainTextEdit.setPlainText(self.data)
        self.label.setText(label_text)   

    def cla(self,):
        self.exp_data_dict= {} # dictionary  
        return 
        
    def add_data(self,index=0):
        self.take_data() 
        self.exp_data_dict[index] = self.get_numpy_array() 
        return 
    
    def take_data(self) :
        # make the text ends with \n 
        mytext = (self.plainTextEdit.toPlainText()).strip()+'\n'
        self.data = mytext
        return 

    def clear_data(self):
        self.plainTextEdit.clear()
        return 
        
    def get_numpy_array(self,):
        """
        convert text into numpy array 
        """
        ff = StringIO(self.data)
        numpy_array = np.loadtxt(ff)
        return numpy_array

class DialogFresco_GUI(QDialog,form_omp):
    def __init__(self,data=None):
        super().__init__()
        self.setupUi(self)
        
        self.fresco_in = '' # input string 
        self.expdata = None # dictionay 
        self.fresco_results=None #output dictionary 
        #-- layout 
        self.fig = Figure() 
        self.canvas = FigureCanvas(self.fig)
        #self.plot_layout.addWidget(self.canvas)
        #--default plot 
        self.ax = self.fig.add_subplot(111)
        self.ax.text(0.5,0.5,'Plot Area' )
        self.canvas.draw() 
        self.plot_layout.addWidget(self.canvas)
        toolbar = NavigationToolbar(self.canvas, self)
        self.plot_layout.addWidget(toolbar)
        
        #-----data dialog
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.pushButton_load.clicked.connect(self.open_file)
        self.pushButton_save.clicked.connect(self.save_file)
        self.pushButton_run.clicked.connect(self.run_fresco)
        self.pushButton_plot.clicked.connect(self.plot)
        self.pushButton_expdata.clicked.connect(self.get_expdata)
        
        # replace existing widget in ui to new class     
    def get_expdata(self,):
        dlg = addData_Dialog()
        check = dlg.exec_()
        if check:
            self.expdata = dlg.exp_data_dict 
        return 
        
    def run_fresco(self,):
        self.Widget_text.clear()
        self.Widget_text.append('Running Fresco. Please wait.\n' )
        fresco_input_txt = self.Text_fresco_in.toPlainText().strip()+'\n'
        self.fresco_in = fresco_input_txt
        #---remove previous results------------
        
        out = run_fresco_v2.run_fresco_from_input_txt(
                    fresco_input_txt,
                    fresco_path='fresco.exe',
                    fresco_input_path='_test.in',
                    fresco_output_path='_test.out',
                    verbose=True)

        # check result
        if run_fresco_v2.chck_fresco_out(fname='_test.out')==0:
            ff=open('_test.out','r')
            outtxt = ff.read()
            ff.close() 
            self.Widget_text.append(outtxt)
            self.Widget_text.append('Fresco run finished.\n' )
        else:
            self.Widget_text.clear()  
            self.Widget_text.append('Error in Fresco run.\n' )
            ff=open('_test.out','r')
            ll=ff.readlines()
            ff.close()
            # show last few lines of output
            self.Widget_text.append(ll[-3]+ll[-2]+ll[-1])
            self.fresco_results=None
            return
        #---------------------------
        # Get results
        #---------------------------
        elastic_out = run_fresco_v2.get_elastic_result_from_fresco_out(
                fname='_test.out')
        angle = elastic_out[:,0]
        xs = elastic_out[:,1]
        ratio = elastic_out[:,2]
        out[0] = elastic_out # replace elastic part
        # neutral scattering case neutral_case=0
        # and np.sum(ratio) = 0.
        #--------
        # Store results for later plotting
        #--------
        self.fresco_results = out
        return 
    
    def plot(self,):
        self.ax.cla()
        if self.fresco_results is not None:
            #---get current status 
            elastic_opt = self.comboBox_type.currentIndex()
            channel =  self.spinBox_channel.value()
            y_scale = self.comboBox_scale.currentIndex()
            data_num = self.spinBox_dataN.value()
            
            #--set scale 
            if y_scale==0:
                self.ax.set_yscale('linear')
            else:
                self.ax.set_yscale('log')
            #---try plot     
            try:
                data_a = self.fresco_results[channel]
                if channel==0: #elastic case 
                    if elastic_opt==0:
                        self.ax.plot(data_a[:,0],data_a[:,2])
                    else:
                        self.ax.plot(data_a[:,0],data_a[:,1])
                else:
                    self.ax.plot(data_a[:,0],data_a[:,1])                    
            except:
                self.Widget_text.append('\n Plotting Error!\n')
                pass 
            try:    
                if self.expdata is not None:
                    exp_data = self.expdata[data_num] 
                    self.ax.errorbar(exp_data[:,0],exp_data[:,1],fmt='*',yerr=exp_data[:,2],label='exp.')
            except: 
                pass 
            self.canvas.draw() 

        else:
            return 
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
                ff= open(fileName,'r')
                self.fresco_in = ff.read()
                ff.close()
                self.Text_fresco_in.setPlainText(self.fresco_in)
                return 
            else:
                print('No file is chosen')
                return
    def save_file(self,):
            options = QFileDialog.Options()
            fileName, _filter = QFileDialog.getSaveFileName(self,
                        "Save file",
                        "",
                        "All Files (*)",
                        options=options)
            if fileName :
                print('Save Results to {}'.format(fileName))
                #---write to fileName
                ff=open(fileName,'w')
                txt = self.Text_fresco_in.toPlainText().strip()+'\n'
                ff.write(txt)
                ff.close()
                if self.expdata is not None:
                    ff=open(fileName+'_dat','w') 
                    for item in self.expdata.keys():
                        data_array = self.expdata[item]
                        print('save data') 
                    ff.close() 
                
                return 
            else:
                print('No file is chosen')
                return                                
#=============================================================================

#===================================================================================================       
if __name__ == "__main__":
    #---for test ---
    def run_app():
        """
        launcher of Qt in Spyder 
        to avoid error 
        """
        if not QApplication.instance():
            app = QApplication(sys.argv)
        else: 
            app = QApplication.instance() 
            
        myWindow = DialogFresco_GUI()
        myWindow.show()
        app.exec_() 
        return myWindow  
    m = run_app() 
    #sys.exit(app.exec_()) # when run in console
    #app.exec_()           # when run in cmd