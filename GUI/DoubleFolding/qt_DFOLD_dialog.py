# -*- coding: utf-8 -*-
"""
Created on Wed May 13 13:00:38 2020

@author: Y.-H. song
"""
import sys
import os 
#---current path where this file resides 
try:
    here = os.path.dirname(os.path.realpath(__file__))
except:
    here = '.'
sys.path.insert(0,here) # to import DFPot3 in this directory 


from PyQt5.QtWidgets import (QApplication,QMainWindow,QDialog,QFileDialog
                             ,QWidget,QComboBox,QVBoxLayout,QLineEdit,QFormLayout)
 
from PyQt5 import uic

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
        FigureCanvasQTAgg as FigureCanvas,
        NavigationToolbar2QT as NavigationToolbar)
        
from subprocess import (call,Popen)
import numpy as np

import DFPot3 

form_omp = uic.loadUiType(here+'/'+"dialog_DFOLD.ui")[0]

element_names = ["n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl",
		 "Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
		 "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb",
		 "Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er",
		 "Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At",
		 "Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No",
		 "Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og",
		 "119","120","121","122","123","124","125","126","127","128","129","130"]

def Sao_Paulo_density(A,Z):          
    """
    return WS potential parameter of density
    based on Sao-Paulo parametrization
    """
    N = A - Z
    # rho0=0.091 # fm^{-3} 
    R_p = 1.81*Z**(1./3.)-1.12
    a_p = 0.47-0.00083*Z 
    R_n = 1.49*N**(1./3.)-0.79
    a_n = 0.47+0.00046*N 
    return (R_p,a_p,R_n,a_n)


class ChooseDensity(QWidget):
    def __init__(self,A,Z,charge,
                 item_list = ['HFB-14', 'Gaussian', 'Woods-Saxon','Sao-Paulo-density' ]):
        super().__init__()
        self.A = A
        self.Z = Z 
        self.pn_choice = charge 
        self.density_function = None 
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        
        self.comboBox = QComboBox()
        self.layout.addWidget(self.comboBox)
        self.comboBox.addItems(item_list) 
        self.comboBox.currentIndexChanged.connect(self.change_form) 
      
    def change_form(self,):
        try: 
            self.layout.removeWidget(self.subWidget)
            self.subWidget.close() 
        except: 
            pass 
        # For other forms of density 
        self.subWidget = QWidget() 
        self.sublayout = QFormLayout() 
        self.subWidget.setLayout(self.sublayout)  
        
        if self.comboBox.currentText() == 'Gaussian':
            self.layout.addWidget(self.subWidget)
            self.lineEdit_alpha = QLineEdit()  
            self.sublayout.addRow('alpha',self.lineEdit_alpha) 
        elif self.comboBox.currentText() == 'Woods-Saxon': 
            self.layout.addWidget(self.subWidget)
            self.lineEdit_R0 = QLineEdit() 
            self.lineEdit_a0 = QLineEdit() 
            self.sublayout.addRow('R0',self.lineEdit_R0)
            self.sublayout.addRow('a0',self.lineEdit_a0)
        elif self.comboBox.currentText() == 'Sao-Paulo-density':     
            self.layout.addWidget(self.subWidget)
            self.lineEdit_R0 = QLineEdit() 
            self.lineEdit_a0 = QLineEdit() 
            self.sublayout.addRow('R0',self.lineEdit_R0)
            self.sublayout.addRow('a0',self.lineEdit_a0)
            (R_p,a_p,R_n,a_n) = Sao_Paulo_density(self.A,self.Z)
            if self.pn_choice == 0 :  #neutron
                self.lineEdit_R0.setText('{:.4}'.format( R_n) )
                self.lineEdit_a0.setText('{:.4}'.format( a_n) )
            elif self.pn_choice == 1 :  #neutron
                self.lineEdit_R0.setText('{:.4}'.format( R_p) )
                self.lineEdit_a0.setText('{:.4}'.format( a_p) )    
            
    def get_density(self,):
        """
        return (interpolating) function of density 
        
        in case of HFB-14 --> interpolate density 
        in case of Gaussian and others 
           --> need to take additional input parameters 
               to construct function also need to normalize 
        """
        if self.pn_choice == 0 : norm = self.A - self.Z 
        if self.pn_choice == 1 : norm = self.Z 
        
        if self.comboBox.currentIndex() == 0 : # HFB-14 
            try: 
               (r,rhon,rhop) = DFPot3.read_density(self.Z,self.A)
            except:    
                return (None, -1)
        
            if self.pn_choice==0: #neutron 
                self.density_function = DFPot3.interpolating_function(r, rhon)
            elif self.pn_choice==1 : #proton 
                self.density_function = DFPot3.interpolating_function(r, rhop)
            return (self.density_function, 0)     
        elif self.comboBox.currentIndex() == 1 : #Gaussian 
            alpha = float(self.lineEdit_alpha.text() ) 
            f_rho = lambda r : DFPot3.shape_Gaussian(r,alpha=alpha)
            C = DFPot3.normalize_density(f_rho, norm ,r_range=(0.,20.),num_quad=70)
            self.density_function = lambda r : DFPot3.shape_Gaussian(r,alpha=alpha,C=C) 
            return (self.density_function, 0)     
        elif self.comboBox.currentIndex() in [2,3] : #WS or Sao-Paulo
            R0 = float(self.lineEdit_R0.text() ) 
            a0 = float(self.lineEdit_a0.text() ) 
            f_rho = lambda r : DFPot3.shape_WS(r,R0=R0,a0=a0,V0=1.0)
            C = DFPot3.normalize_density(f_rho, norm ,r_range=(0.,20.),num_quad=70)
            self.density_function = lambda r : DFPot3.shape_WS(r,R0=R0,a0=a0,V0=C)
            return (self.density_function, 0)     
        else : # not-available yet 
            return (self.density_function,-1)     
          
class DialogDFOLD_GUI(QDialog,form_omp):
    def __init__(self,ap,zp,at,zt,e_a,data_path='./DoubleFolding/density-hfb14'):
        super().__init__()
        self.setupUi(self)
        #if fresco_input_object:
        #    self.fresco_input = fresco_input_object
        #else: 
        #    raise ValueError('Error, one have to enter fresco_input_object')
            
        self.input= {'ap':ap,'at':at,'zp':zp,'zt':zt,'E/A':e_a }
        self.densities = {'proton-projectile': None, 
                          'neutron-projectile':None,
                          'proton-target': None, 
                          'neutron-target': None } 
        self.DFpot = {'R':None, 'Coulomb': None, 
                      'Isoscalar': None, 
                      'Isovector': None }   
        self.status = True                       
        #--button action
        self.buttonBox.accepted.connect(self.take_data)
        self.buttonBox.rejected.connect(self.reject)
        #--title text
        text = ' {}{} + {}{} at Elab={} MeV'.format(
            ap,element_names[zp],at,element_names[zt],e_a*ap)
        self.label_Main.setText(text)
        
        # replace existing widget in ui to new class 
        pot_list =['M3Y_Paris_ZR', 'M3Y_Reid_ZR',
                   'M3Y_Paris_FR', 'M3Y_Reid_FR',
                   'DDM3Y_Reid',
                'DDM3Y_Paris' ,
                'BDM3Y_Reid' ,
                'BDM3Y_Paris',
                'CDM3Y1_Paris',
                'CDM3Y2_Paris',
                'CDM3Y3_Paris',
                'CDM3Y4_Paris',
                'CDM3Y5_Paris',
                'CDM3Y6_Paris',
                'to_be_added']
        self.comboBox_potential.addItems(pot_list) 
        
        self.comboBox_proton_projectile_density = ChooseDensity(ap,zp,1)
        self.gridLayout.addWidget( self.comboBox_proton_projectile_density,1,1)
        self.comboBox_neutron_projectile_density = ChooseDensity(ap,zp,0)
        self.gridLayout.addWidget( self.comboBox_neutron_projectile_density,1,2)
        self.comboBox_proton_target_density = ChooseDensity(at,zt,1)
        self.gridLayout.addWidget( self.comboBox_proton_target_density,1,3)
        self.comboBox_neutron_target_density = ChooseDensity(at,zt,0)
        self.gridLayout.addWidget( self.comboBox_neutron_target_density,1,4)
                
        self.comboBox_potential.currentIndexChanged.connect(self.pot_change) 
        self.pushButton.clicked.connect(self.update)  #<--compute !
        
        self.canvas = FigureCanvas(Figure() )
        self.plot_layout.addWidget(self.canvas)
        self.toolbar = NavigationToolbar(self.canvas,
                self.plot_Widget, coordinates =True)
        self.plot_layout.addWidget(self.toolbar)
                
    def take_data(self,):
        # actions to store data and update parent window
        # self.DFpot stores the folding potential
        # but whether to use it or not is determined in the parent window. 
        self.accept()
        
    def pot_change(self,):
        # potential changed and change density as such 
        if self.comboBox_potential.currentText() in ['to_be_added']:
            self.comboBox_proton_projectile_density.setEnabled(False)
            self.comboBox_neutron_projectile_density.setEnabled(False)
            self.comboBox_proton_target_density.setEnabled(False)
            self.comboBox_neutron_target_density.setEnabled(False)
        else: 
            self.comboBox_proton_projectile_density.setEnabled(True)
            self.comboBox_neutron_projectile_density.setEnabled(True)
            self.comboBox_proton_target_density.setEnabled(True)
            self.comboBox_neutron_target_density.setEnabled(True)
  
    def update(self,):
        # update info , densities and compute folding potential and plot 
        self.status = True 
        self.rm_plot()
        self.label_status.setText('')
        
        pot_choice = self.comboBox_potential.currentText()
        #-------prepare density--------------------------
        if pot_choice in ['to_be_added']:
            # no density function is necessary 
            pass 
        else:     
            # prepare density       
            den_sts=[0,0,0,0]
            den_names=['prot.-proj.','neut.-proj.','prot.-targ.','neut.-targ.']
            (self.densities['proton-projectile'],den_sts[0]) = self.comboBox_proton_projectile_density.get_density()
            (self.densities['neutron-projectile'],den_sts[1])= self.comboBox_neutron_projectile_density.get_density()
            (self.densities['proton-target'],den_sts[2]) = self.comboBox_proton_target_density.get_density() 
            (self.densities['neutron-target'],den_sts[3]) = self.comboBox_neutron_target_density.get_density() 
            text=""
            for i in range(4):
                if den_sts[i]<0: 
                    text += 'Error: {} density is not available!\n'.format(den_names[i]) 
                    self.status = False
            self.label_status.setText(text)            
                
        #------- compute DF potential -------------------    
        if self.status: # density is prepared  
            if 'M3Y' in pot_choice:
                (R,UC,US) = DFPot3.DF_M3Y_all(self.densities['neutron-projectile'],
                                  self.densities['proton-projectile'],
                                  self.densities['neutron-target'],
                                  self.densities['proton-target'],
                                  type = pot_choice ,
                                  E_A=self.input['E/A'],
                                  R_points=np.arange(0.1,18.0,0.1),
                                  k_range = (0.0,5.0), num_quadrature= 96,
                                  r_range =(0.0,15.0), iter_max=30)
                # --store data
                (self.DFpot['R'], self.DFpot['Coulomb'], self.DFpot['Isoscalar']) = (R, UC, US)
                self.label_status.setText('All prepared.') 
                # plot graph --------------
                fig1 = Figure()
                ax1f1 = fig1.add_subplot(1,2,1)
                ax1f1.plot(R,self.densities['neutron-projectile'](R)+self.densities['proton-projectile'](R),
                           label='projectile density')
                ax1f1.plot(R,self.densities['neutron-target'](R)+self.densities['proton-target'](R),
                           label='target density')           
                ax1f1.set_xlabel('fm')
                ax1f1.set_ylabel(r'$fm^{-3}$')
                ax1f1.legend()           
                ax2f1 = fig1.add_subplot(1,2,2) 
                ax2f1.plot(R,US,label='iso-scalar U') 
                ax2f1.set_xlabel('fm')
                ax2f1.set_ylabel('MeV')
                ax2f1.legend()            
                self.add_plot(fig1) 
            
            else:     
                self.label_status.setText('Error: potential option is not available yet!')   
                self.status = False # error 

            
    def add_plot(self,fig):
        # add canvas to layout of widget
        self.canvas = FigureCanvas(fig)
        self.plot_layout.addWidget(self.canvas) # add a canvas widget into layout
        self.canvas.draw()
        # add toolbar to widget
        self.toolbar = NavigationToolbar(self.canvas,
                self.plot_Widget, coordinates =True)
        self.plot_layout.addWidget(self.toolbar)

    def rm_plot(self,):
        self.plot_layout.removeWidget(self.canvas)
        self.canvas.close()
        self.plot_layout.removeWidget(self.toolbar)
        self.toolbar.close()
#=============================================================================

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
            
        ap=40;zp=20;at=80;zt=40;e_a=10.0; #default
        myWindow = DialogDFOLD_GUI(ap,zp,at,zt,e_a)
        myWindow.show()
        # app.exec_() 
        return myWindow  
    m = run_app() 
    #sys.exit(app.exec_()) # when run in console
    #app.exec_()           # when run in cmd