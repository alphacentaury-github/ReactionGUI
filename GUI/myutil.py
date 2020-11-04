#! /usr/bin/python3
#
"""
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
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import scipy.optimize 
import re 
import os 
from io import StringIO
import pickle 

hc=197.3269 # MeV fm
amu=931.494 #MeV
mp=938.27
mn=939.5653
mHe6=4*mn+2*mp-4.878*6
mHe4=2*mn+2*mp-7.073*4
mH3=2*mn+mp-2.8272*3

#----------------------------------------------------------------
class all_global_variables() :
    """ store global values which are used 
        in other part of codes
        as dictionary. 
        
        For the moment, it is to store path,
        but can be used for other purpose. 
    """
    # global over class
    savefile = '.gui_path'
    class_property = {'current_path': os.getcwd() , 
                      'fresco_filename': 'fresco.exe',
                      'sfresco_filename': 'sfresco.exe',
                      'frin_filename': '_test.in',
                      'frout_filename': '_test.out',
                      'minunit_filename': '_test.min',
                      'search_filename': '_test.search',
                      'omegt_filename' : os.getcwd()+'/omget_RIPL3/omget.exe'
                          }    
    def __init__(self,):
        self.object_property = 0 
    def set_global_value(self, **kwarg ):
        # add/change values 
        for item,value in kwarg.items():
            all_global_variables.class_property[item] = value
    def del_global_variables(self):
        # empty dictionary 
        all_global_variables.class_property ={} 
    def save_to_file(self):
        # save dictionary to file 
        ff= open(all_global_variables.savefile,'wb')
        pickle.dump(all_global_variables.class_property,ff)
        ff.close() 
    def load_from_file(self):
        #load dictionary from file 
        ff= open(all_global_variables.savefile,'rb')         
        para_dict = pickle.load(ff)
        ff.close() 
        all_global_variables.class_property = para_dict 
        return para_dict 
                
#---------------------------------------------------------------------
def gauleg(a,b,n):
    """
    Gauss Legendre quadrature points/weights in (a,b) of degree 2^n-1

    Parameters
    ----------
    a : float
        lower bound 
    b : float 
        upper bound 
    n : integer
        degree 
    Returns
    -------
    x : ndarray
       sample points 
    w : ndarray 
       weights 
    """
    x, y = np.polynomial.legendre.leggauss(n) 
    x = ((b-a)*x+b+a)*0.5 
    w = (b-a)*0.5* y 
    return (x,w) 

#
# Replace part of file strings
#
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def replace_line_infile(ff, line_num, strings,ffo=None):
   """ replace entire lines with new one in file 
       line number is counted from zero !
   """
   f=open(ff,'r')
   lines=f.readlines()
   if (strings[-1]=='\n'):
     lines[line_num]=strings
   else:
     lines[line_num]=strings+'\n' # new line
   f.close()
   if not ffo:  # ffo is empty
     ffo=ff
   f=open(ffo,'w')
   for i in lines:
     f.write(i)
   f.close()

def replace_line_column_infile(ff,line_num,col_num,strings,ffo=None):
  """ replace column in line with new one in file
      line number and column number counts from zero"""
  f=open(ff,'r')
  lines=f.readlines()
  words=lines[line_num].split()
  words[col_num]=strings
  # new replaced line
  ss=''
  for i in words:
    ss=ss+i+'  '  
  lines[line_num]=ss+'\n'
  f.close()
  if not ffo: # ffo is empty
    ffo=ff
  f=open(ffo,'w')
  for i in lines:
    f.write(i)
  f.close()

def tabulate(inx,iny,xnew,fill_value=None):
  """
    convert array (x,y) into a regular table of
    [xnew,ynew] in equal steps, 
    fill_value= (y_min, y_max) : set y=y_min if x< inx_min
                                     y=y_max if x> inx_max
                'extrapolate'  : use extrapolation if x is out of range of inx         
    return [xnew,ynew] array                    
  """
  if fill_value:
    f=interp1d(inx,iny,bounds_error=False,fill_value=fill_value)
  else :
    f=interp1d(inx,iny)
  out=np.column_stack([xnew,f(xnew)]) #numpy array 
  return out

def convert_table(fname,x_range,fill_value=None):
  """
   load data from file and convert it into a table in equal step size
   f_name : file name string
   x_range: array [xini,xfinal,xstep ]
   fill_value: option for extrapolation (y< value, y> value) or 'extrapolate'
  """
  xnew=np.arange(x_range[0],x_range[1],x_range[2])
  dat=np.loadtxt(fname)
  out=tabulate(dat[:,0],dat[:,1],xnew,fill_value=fill_value)
  return out 

def txt_to_array(text):
    # text table into array 
    ff= StringIO(text)
    array = np.loadtxt(ff)
    return array 

def array_to_txt(array):
    ff = StringIO()
    np.savetxt(ff,array)
    return ff.getvalue() 

def add_error_to_datatext(text):
    """
    text is assumed to contain
    
    columns of 
    angle, x-s, error 
    
    if there is no error in the text,
    add small error to the text 

    """
    data = txt_to_array(text)
    if data.shape[0] == 0: # empty text 
        return text  
    if data.shape[1] == 2: # no error is given 
        new_text=''
        for line in text.split('\n'): 
            line = line.strip()
            if len(line)==0:
                continue 
            if line[0]!='#' : # data line 
                new_text = new_text+line+' 0.001\n' 
            else:
                new_text = new_text+line+'\n'
        text = new_text         
    return text 


def fresco_input_formfactor(dat,x_range=None,fname='fort.4',
                            fill_value=None,
                            comment='',num_type='R',
                            write='w'):
  """
    Generate input potential form factor format 
    write to 'fname' file
    (or use write ='a'  to append )
    
    num_type='R','I','C' corresponds to shape=7,8,9 

    input data is an array shape (:,2)
    when 'R','I', it is a real-valued for real/imaginary part of potential
    if 'C', it is a complex number  
  """
  xnew=np.arange(x_range[0],x_range[1],x_range[2])
  npoints=len(xnew) 
  rstep=x_range[2]
  rfirst=xnew[0]
  out=tabulate(dat[:,0],dat[:,1],xnew,fill_value=fill_value)  
  ff=open(fname,write)
  ff.write(comment+'\n')
  ff.write('%i  %f  %f\n'%(npoints,rstep,rfirst))
  for i in out[:,1]:
      if num_type in ['R','I']: # i must be real-valued 
          ff.write('%f\n'%(np.real(i)))
      elif num_type=='C':
          ff.write('%f %f\n'%(np.real(i),np.imag(i)))
  ff.close()   
  return 

def read_namelist_input(text,upper_case_item=False):
    """
    asuume text is in written in namelist format
    
    item1 = value1 item2=value2  
    
    convert it into a dictionary 
    
    upper_case_item = True : change the all items into upper case 
    
    note: the values are still string. 
          
    """
    
    out_dict = {}
    text = re.sub('[,:{}\t\n"]',' ',text)
    text = re.sub('=',' = ',text) # always space around '=' 
    ww = text.split() 
    for j in range(len(ww)): 
        if ww[j]=='=':
            if upper_case_item:
                item = ww[j-1].upper()
                value = ww[j+1]
            else :
                item = ww[j-1]
                value = ww[j+1]
            out_dict[ item ] = value 
    return out_dict    
    
def write_namelist_form(dict_in,upper_case_item=False,separator=' '):
    """
    write dictionary into a namelist form text
    """
    out_text = ''
    for item in dict_in.keys():
        if upper_case_item :
            item_ = item.upper() 
        else:
            item_ = item 
        out_text += '{}={}{}'.format(item_, dict_in[item],separator)
    return out_text     
    
def to_lab_energy(mp,mt,ecm):
   """ convert ecm to elab
      for projectile mass mp
          target mass mt
      all are in MeV
   """
   mu=mp*mt/(mp+mt*1.0)
   return mp/mu*ecm

def to_cm_energy(mp,mt,elab):
   mu=mp*mt/(mp+mt*1.0)
   return mu/mp*elab


def BW_form(x,er,gamma):
   """ Breit-Wigner form for resonance
   """
   return 0.25*gamma**2/((x-er)**2+0.25*gamma**2)


def find_resonance(ecm,phase,guess=1.0,guess2=2.0):
  """ find the position of resonance and width
    from energy and phase array (degree)
    Assume the case with Breit-Wigner form 
    background phase shift=0.
    units are degree.
    center of mass energy in MeV

    resonance energy is defined by passing 90 degree phase shift
    width is guessed from 45 degree phase shift
    But, then cross section is fitted with BW-form  

  """
  # first make the phase shift positive, continuous 
  for i in range(len(phase)) :
    if  phase[i]< 0 : 
       phase[i] = phase[i]+180. 
  # now find the position of 90 degree
  f=interp1d(ecm,phase-90.,'cubic')
  # brentq method for solution
  #er=brentq(f,guess,guess2)  
  #  root finding
  sol=scipy.optimize.root(f,guess) 
  er=sol.x  # resonance energy
  # now search width 
  #  first guess width from result
  f=interp1d(ecm,phase-45.,'cubic')
  sol=scipy.optimize.root(f,er)
  gam=abs(sol.x-er)/2.
  # now curve fit with BW_form
  f=interp1d(ecm,phase,'cubic')
  xn=np.arange(er-gam,er+gam,2*gam/100.)
  ydata=np.sin(f(xn)*np.pi/180.)**2
  popt,pcov=scipy.optimize.curve_fit(BW_form,xn,ydata,p0=(er,gam))
  
  xn=np.arange(ecm[0],ecm[-1],(ecm[-1]-ecm[0])/200.)
  ydata=np.sin(f(xn)*np.pi/180.)**2
  gdata=BW_form(xn,popt[0],popt[1])
  plt.plot(xn,ydata,xn,gdata)
  plt.title('E=%f +i %f'%(popt[0],popt[1]))
  plt.savefig('phase.png')
  return popt[0],popt[1] 

def clean_comm(fname):
  # remove all comments @ from the file 
  #  and add blank line for END or & 
  f=open(fname,'r')
  lines=f.readlines()
  f.close() 
  out_fname=fname+'x' # default 
  f=open(out_fname,'w')
  for i in lines:
     if '@' in i:
        i='#'+i
     if ('END'in i) or ( '&' in i ):
        i='#'+i+'\n\n'
     f.write(i)
  f.close()
  return 


def read_fresco_res(fname):
  """ read text files
      assume different data are separated by a blank line
      (and '#' '!' are used for commenting )
      
      return dictionary of list of strings       
  """
  ff=open(fname,'r')
  lines=ff.readlines()
  ff.close()
  out={}
  j=0
  ll=[]
  for i in lines:
    w=i.split()
    if len(w)!=0 and (w[0][0] in ['#','!']):
      continue
    if len(w)!=0 :
      ll.append([ float(k) for k in w])
    if len(w)==0 and len(ll)==0 : # continuos blank
      continue
    if (len(w)==0 and len(ll)!=0) or (i ==lines[-1]): #if met blank, it means end of data
      out[j]=ll[:]
      j=j+1
      ll=[]
  return out 


def filon_integral1(x_array,f_val,omega,osc='sin'):
    """
    osc='sin' or 'cos' for integral , g(z)->osc(z) 
    
    I = \int_x0^x_{2n} f(x) g(omega x) dx 
          
    x_array should be equidistanced 2n+1 number of points x_0..x_{2n} 
    f_val is a value of f(x) at those points.
    
    omega is a frequency.
    
    REF: Abramowitz and Stegun p.890  
    """
    # check input 
    if (len(x_array)%2 ==0 ): # even number of points 
        print('Error: number of sample points should be odd! ')
        return None 
    if osc=='sin':
        g  = np.sin(omega*x_array) # g_1(omega x)
        gp = np.cos(omega*x_array) 
    elif osc=='cos':
        g = np.cos(omega*x_array) # g_1(omega x)
        gp = -np.sin(omega*x_array) 
    else: 
        print('Error: oscillating function should be sin or cos!')
        return None 
    #--compute alpha , beta, gamma 
    h = x_array[1]-x_array[0] 
    theta = omega*h  
    #---need to avoid theta=0 case... 
    alpha = 1.0/theta+np.sin(2*theta)/(2*theta**2)-2.*np.sin(theta)**2/theta**3
    beta  = 2.0*( (1.0+np.cos(theta)**2)/theta**2 -np.sin(2*theta)/theta**3 )
    gamma =4.0*(np.sin(theta)/theta**3 -np.cos(theta)/theta**2 )
    # compute even sum and odd sum of f(i)*g(omega*i). 
    term1 = f_val[0]*gp[0]-f_val[-1]*gp[-1] 
    even_sum = np.sum(f_val[0::2]*g[0::2])-0.5*(f_val[0]*g[0]+f_val[-1]*g[-1]) 
    odd_sum = np.sum(f_val[1::2]*g[1::2])
    total = h*(alpha*term1 +beta*even_sum+gamma*odd_sum) # negelect S' or C' terms 
    return total      

def spherical_bessel_decompose(L,z):
    """
    return coefficient a and b for 
    sin and cos decomposition of spherical bessel function 
    
    j_L(z) = a_L *sin(z) + b_L *cos(z) 
    
    REF: Abramowitz and Stegun, p.438 
    
    Note: a_L, b_L is singular at z=0.
         Thus, not reliable for small z
           
    """
    f_pos=[z**(-1),z**(-2)] # f_pos[n] = f_n 
    f_neg=[z**(-1),0.0] # f_neg[n] = f_{-n}
    for i in range(1,L+2):
        f_plus = (2*i+1)*f_pos[i]/z-f_pos[i-1] # f_{i+1}
        f_mins  = (2*(-i)+1)*f_neg[i]/z-f_neg[i-1] # f_{-i-1}
        f_pos.append(f_plus) 
        f_neg.append(f_mins)  
    return (f_pos[L], (-1)**(L+1)*f_neg[L+1])    
    
def SBT_filon(x_array,f_val,L,k):
    """
    compute spherical bessel transformation 
    by using Filon's method 
    
    integral= int_{x_0}^{x_2n} dx f(x) j_L(k x)

    Parameters
    ----------
    x_array : numpy array
         x[i]= x_0+ i* h 
         x_0 must be non-zero. 
    f_val : numpy array 
        sample values of f(x) at x_array 
    L : integder
        order of spherical bessel function
    k : float
        frequency of j_L(k x)

    Returns
    -------
    integral= int_{x_0}^{x_2n} dx f(x) j_L(k x) 

    Note
    -------
    By the nature of decomposition of spherical bessel into 
    sin and cos , code becomes singular at k*x=0. 
    Thus,  x_0 > 0. 
    
    """    
    # for test 
    # x_array =np.linspace(0.001,15.0,101); f_val= x_array;k = 1.0; L=1 
    # shift the points a little bit to avoid x=0. 
    if x_array[0] < 1.0e-20:
        x_array = x_array +1.0e-20
    z = k*x_array 
    if z[0] < 1.0e-20 : 
        z[0] = z[0]+1.0e-20 
    # decompose spherical bessel function into 
    (a,b) = spherical_bessel_decompose(L,z)
    # define new integrand f*a and f*b 
    f_s = f_val*a
    f_c = f_val*b 
    sin_integral = filon_integral1(x_array,f_s ,k ,osc='sin')
    cos_integral = filon_integral1(x_array,f_c ,k ,osc='cos')
    integral = sin_integral+cos_integral 
    return integral 

def frac2float(fnumber):   
    """
    convert 'string expression' of fraction into float
    """
    fnumber = fnumber.strip() 
    if '/' in fnumber:
       ww  = fnumber.split('/')
       ff = float(ww[0])/float(ww[1])
    else :
       ff = float(fnumber)
       
    return ff


#=====================MAIN===================================== 

""" This routine only works for a specific cutoff value """
if __name__ == '__main__':
   import sys 
   from subprocess import call
  # remove all comments @ from the file
   if (sys.argv[1]=='res'):
     clean_comm(str(sys.argv[2]))
     fname=sys.argv[2]+'x'
     out=np.loadtxt(fname)
     en=out[:,0];phase=out[:,3]
     ecm=to_cm_energy(2,4,en)
     er,gam=find_resonance(ecm,phase,guess=0.8,guess2=2.0)
   elif (sys.argv[1]=='clean'):
     clean_comm(str(sys.argv[2]))
   else :
   #---default case
   # usage: test.py [infile] [inf]
   #        convert fort.16 into .res file
     infile=str(sys.argv[1])
     print(infile)
     inf = infile.split('.')[0]
     if (len(sys.argv)==3):
        inf = str(sys.argv[2])
     call("fresco < "+infile+" > "+inf+".out",shell=True)
     clean_comm('fort.16')
     call('mv fort.16x '+inf+'.res',shell=True)
     call('rm fort.*',shell=True) # ! delete except fort.4 ??? 
     # call fresco
     # cleanup file
     # rename file 
