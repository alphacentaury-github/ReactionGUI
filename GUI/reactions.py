#! /usr/bin/python
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

#
# Utility programs for reaction calculations
#
# () leggau : Gauss-Legendre Quadrature in arbitrary range
#
# ()  
import numpy as np
#import scipy
#import matplotlib.pyplot as plt 
import sys
import math
#from scipy import interpolate
#from scipy import integrate
#from scipy.interpolate import interp1d

#----import personal library
#  sys.path.insert(0, '\home\yhsong\lib','\home\yhsong\bin')
#  import test


f_mass='mass16.txt'
f_nubase='nubase2016.txt'

hbarc = 197.326968 
amu = 931.4940954 # MeV
mp  = 938.272    
mn  = 939.5653
m_N = (mp+mn)/2.0
alpha = 1.0/137.03599

element_names = ["n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl",
		 "Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
		 "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb",
		 "Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er",
		 "Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At",
		 "Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No",
		 "Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og",
		 "119","120","121","122","123","124","125","126","127","128","129","130"]
#-----------------------------------------------
def leggau(rmin,rmax,degree):
  """ get quadradure points and weights
  """  
  import numpy.polynomial.legendre as lg
  x , w = lg.leggauss(degree) # in range (-1,1) 
  x = (rmax+rmin)/2.+x*(rmax-rmin)/2.
  w = (rmax-rmin)/2.*w
  return x, w

def rmcomm(text,num_arg=1):
  """
    read num_arg numbers of data from file
    ignoring comments in texts
    only takes values prior '#' or '!' 
  """
  text=text.replace('!','#')
  ww=text.split('#')
  out=ww[0].split()[0:num_arg]
  if num_arg==1:
    out=out[0]
  return out

#------------------------------------------------
def nuclei_name(a,z):
    if z < 119:
        name = '{}{}'.format(a,element_names[z])
    else :
        name = '{}-{}'.format(a,element_names[z])
    return name 


def read_mass(a,z,fname=f_mass,subtract_electron=False):
   """ Read the amotic mass table
       and return the bare nuclear mass of A,Z  
       by subtracting electron mass
       in units of amu.

   Note: (atomic mass)=(nuclear mass)+(elcectron mass) in amu
         (mass excess)=( (atomic mass)-(mass number) ) in MeV
                       thus mass excess includes electron mass

         (binding energy)=(atomic mass)-(neutron mass)
                          -(proton mass)-(electron mass) in MeV  
                         =(atomic mass)-(neutron atomic mass)
                          -(proton atomic mass)    in MeV   
                       thus binding energy do not inlcude electron
   """
   import sys
   electron_mass = 0.510998 # MeV 
   ff=open(fname,'r')
   lines=ff.readlines()
   ff.close
   for ll in lines[39:]: # start from 40th line
     zz=int(ll[9:14]) 
     aa=int(ll[14:19])
     mm=float(ll[96:99]+'.'+ll[100:106])
     if (a==aa and z==zz): 
       mass=mm
       if subtract_electron:
           bmass=mm-electron_mass*zz/amu  # remove electron mass 
       else:
           bmass=mm
       return bmass
   raise ValueError('#### ERROR: no mass found for a= %i z= %i'%(a,z ))

def read_nuclei(a,z,fname=f_nubase,subtract_electron=False):
   """
     return basic info of a nucleus
     mass, spin , parity 
   """
   ff=open(fname,'r')
   ll=ff.readlines()
   ff.close()
   #bmass=read_mass(a,z)
   out={}
   #out['mass']=bmass
   out['A']=a
   out['Z']=z
   for ii in ll:
      aa=int(ii[0:3])
      zz=int(ii[4:7])
      if ii[7] !='0' :  # not a ground state
         continue      
      if (a==aa and z==zz) :
         out['name']=ii[11:16].strip() 
         mass_excess=float(ii[18:29].replace('#','.'))/1000. #mass excess in MeV 
         out['mass_excess'] = mass_excess 
         if subtract_electron:
             mass=(mass_excess/amu+aa)-0.510998/amu*zz   # amu unit
         else:
             mass=(mass_excess/amu+aa)
         out['mass']=mass #in amu unit 
         ww=ii[79:87]
         if ('(' in ww) or ('#' in ww)    : # not confirmed or measured
           out['J']=None
           out['P']=None
           out['J_P']=None 
           return out
         if '/' in ww : #half integer spin
           jj=int(ww.split('/')[0])
           out['J_P']=ww.strip() 
           out['J']=jj/2.0
           if '+' in ww:
             out['P']=1
           if '-' in ww:
             out['P']=-1             
           return out
         else : #integer spin
           out['J_P']=ww.strip() 
           if '+' in ww:
             out['J']=int(ww.split('+')[0])
             out['P']=1
             return out
           if '-' in ww:
             out['J']=int(ww.split('-')[0])
             out['P']=-1
             return out     
   raise ValueError('#### ERROR: not found a= %i z= %i'%(a,z ))

def interp_nuclei_name(nuc_name):
    """
    convert nuclei name like 16O or A-Z with numbers  
    into A,Z
    """
    #---- a special case
    nuc_name = nuc_name.strip()  # if 16-O form is used     
    if '-' in nuc_name: #this case are numbers like 12-6
        ww = nuc_name.split('-')
        nuc_a = int(ww[0])
        nuc_z = int(ww[1])
        return (nuc_a,nuc_z)
    else:
        head = nuc_name.lstrip('0123456789') # remove numbers 
        nuc_a = int( nuc_name[:-len(head)] ) 
        
        if head == 'n': # special case to distinguish 'N' and 'n'
            nuc_z = 0
            return (nuc_a,nuc_z)
        else: 
            head = head.upper() 
        # normal case 
        for (z,el) in enumerate(element_names[1:]):
            el = el.upper() 
            if el == head : # problem when H in He... 
                nuc_z = z + 1 #because neutron is removed. 
                return (nuc_a,nuc_z)
    # faild to find nuclei         
    raise ValueError('Error in name form: use 16O or A-Z form for Z>120') 
    
def interp_spin_name(spin_name):
    """
    interpret spin_name like 1/2+, 3- , 0.5+
    """    
    spin_name = spin_name.strip() 
    parity = spin_name[-1]
    if '/' in spin_name[:-1]:
       ww  = spin_name[:-1].split('/')
       j = float(ww[0])/float(ww[1])
    else :
       j = float(spin_name[:-1])
    if parity=='+':
        p = 1
    elif parity == '-':
        p =-1
    return (j,p)     
       

def get_Qvalue(ap,zp,ex_p,at,zt,ex_t,ax,zx,ex_x,ar,zr,ex_r):
    """
    It is better to use mass excess directly,
    instead of converting mass into amu unit
    """
    out = read_nuclei(ap,zp); m1_excess = out['mass_excess'] #MeV 
    out = read_nuclei(at,zt); m2_excess = out['mass_excess']
    out = read_nuclei(ax,zx); m3_excess = out['mass_excess']
    out = read_nuclei(ar,zr); m4_excess = out['mass_excess']
    m1 = m1_excess+ex_p 
    m2 = m2_excess+ex_t 
    m3 = m3_excess+ex_x 
    m4 = m4_excess+ex_r 
    Q = (m1+m2-m3-m4 ) # MeV 
    return Q 
            
#------------------------------------------------
def kin():
  """ Interactive interface for kinematics
  """
  print('#=== Kinematics helper ===')
  print('# For reaction T(P,X)R      ' )
  print('# T(arget) P(rojectile) X(projectile like)')
  print('#     (AT,ZT)+(AP,ZP) -> (AX,ZX)+(AR,ZR) ')
  inp=input('Enter AT,ZT, AP, ZP (optional AX,ZX,) :\n')
  #ll = rmcomm(inp)
  #ww = ll.split() 
  ww = rmcomm(inp,num_arg=6) 
  AT,ZT,AP,ZP = [ int(i)  for i in ww[0:4]]
  if (len(ww)>4):
     AX,ZX = [ int(i)  for i in ww[4:6]]
     AR=AT+AP-AX
     ZR=ZT+ZP-ZX
  else :
     AX,ZX,AR,ZR = AP,ZP,AT,ZT
  inp=input(' for exact Q-value enter Ex for AT,AP,AX,AR in amu or blank:\n')
 # ll = rmcomm(inp)
 # ww = ll.split()[0:4]
  ww = rmcomm(inp,num_arg=4)
  if ww!=[]:
    ww = [float(i) for i in ww]
    EXT,EXP,EXX,EXR=(ww[0],ww[1],ww[2],ww[3]) 
  else :
    EXT,EXP,EXX,EXR=(0,0,0,0)
 
  inp=input('Enter energy(MeV) (optional : lab,cm,ea=E/A):\n')
  #ll = rmcomm(inp)
  #ww =ll.split()[0:2]  
  ww = rmcomm(inp,num_arg=2)
  if (len(ww)>1):
      EN= float(ww[0])
      if (ww[1]=='lab'):
          EN_type=0
      elif (ww[1]=='cm'):
          EN_type=2
      elif (ww[1]=='ea'):
          EN_type=1 
  else : #default
      EN = float(ww[0])
      EN_type = 0 # Elab 
  
  txt = kin2(AP,ZP,AT,ZT,EN,AX,ZX,AR,ZR,
         EXP,EXT,EXX,EXR,EN_type,mass_opt='NoAME16')
  print(txt) 
  return 


def kin2_simplified(AP,ZP,AT,ZT,EN, AX=None,ZX=None,AR=None,ZR=None,
          EXP=0,EXT=0,EXX=0,EXR=0,EN_type=0,mass_opt='AME16'    ):
    """ Interactive interface for kinematics
  
      EN_type=0 Lab 
              1 E/A
              2 CM  
    """
    amu=931.4940954 #MeV
    if (AX and ZX): # AX and ZX are given 
         AR=AT+AP-AX
         ZR=ZT+ZP-ZX
    else :          
        AX,ZX,AR,ZR = AP,ZP,AT,ZT
     
    try:  
          M_T=read_mass(AT,ZT)+EXT/amu # in amu units
          M_P=read_mass(AP,ZP)+EXP/amu
          M_X=read_mass(AX,ZX)+EXX/amu
          M_R=read_mass(AR,ZR)+EXR/amu  
          #Qval =(M_P+M_T-M_X-M_R)*amu # in MeV units # This includes Excitation energies
          # or directly use mass excess 
          Qval = get_Qvalue(AP,ZP,EXP,AT,ZT,EXT,AX,ZX,EXX,AR,ZR,EXR)
    except: 
          return 'mass is not known for nuclei' 
    if mass_opt=='NoAME16': # use mass as integer values. But, Qvalue still requires exact one 
          M_T=AT+EXT/amu # in amu units
          M_P=AP+EXP/amu
          M_X=AX+EXX/amu
          M_R=AR+EXR/amu  
          
    MUI=M_P*M_T/(M_P+M_T) # in amu unit
    MUF=M_X*M_R/(M_X+M_R)
   
    if (EN_type==0): # default lab energy chosen
      ELAB=EN 
      EpA=ELAB/AP    # E/A
      ECM=ELAB*MUI/M_P
    elif (EN_type==1): # E/A given 
      EpA=EN
      ELAB=EN*AP  
      ECM=ELAB*MUI/M_P  
    elif (EN_type==2): # CM given     
      ECM=EN
      ELAB=ECM*M_P/MUI
      EpA=ELAB/AP      
    # equivalent X+R reaction     
    ECMf=ECM+Qval
    ELABf=ECMf*M_X/MUF
    EpAf=ELABf/AX  # E/A for x
    return (Qval,ELAB,EpA,ECM,ELABf,EpAf,ECMf)
    
def kin2(AP,ZP,AT,ZT,EN,AX=None,ZX=None,AR=None,ZR=None,
         EXP=0,EXT=0,EXX=0,EXR=0,EN_type=0,mass_opt='AME16' ):
  """ Interactive interface for kinematics
  
      EN_type=0 Lab 
              1 E/A
              2 CM  
  """
  amu=931.4940954 #MeV
  if (AX and ZX): # AX and ZX are given 
     AR=AT+AP-AX
     ZR=ZT+ZP-ZX
  else :          
     AX,ZX,AR,ZR = AP,ZP,AT,ZT
     
  try:  
          M_T=read_mass(AT,ZT)+EXT/amu # in amu units
          M_P=read_mass(AP,ZP)+EXP/amu
          M_X=read_mass(AX,ZX)+EXX/amu
          M_R=read_mass(AR,ZR)+EXR/amu  
          #Qval =(M_P+M_T-M_X-M_R)*amu # in MeV units # This includes Excitation energies
          # or directly use mass excess 
          Qval = get_Qvalue(AP,ZP,EXP,AT,ZT,EXT,AX,ZX,EXX,AR,ZR,EXR)
  except: 
          return 'mass is not known for nuclei' 
  if mass_opt=='NoAME16': # use mass as integer values. But, Qvalue still requires exact one 
          M_T=AT+EXT/amu # in amu units
          M_P=AP+EXP/amu
          M_X=AX+EXX/amu
          M_R=AR+EXR/amu  
          
  MUI=M_P*M_T/(M_P+M_T) # in amu unit
  MUF=M_X*M_R/(M_X+M_R)
   
  if (EN_type==0): # default lab energy chosen
      ELAB=EN 
      EpA=ELAB/AP    # E/A
      ECM=ELAB*MUI/M_P
  elif (EN_type==1): # E/A given 
      EpA=EN
      ELAB=EN*AP  
      ECM=ELAB*MUI/M_P  
  elif (EN_type==2): # CM given     
      ECM=EN
      ELAB=ECM*M_P/MUI
      EpA=ELAB/AP
      
  #kinematic conversion factor rho
  rho_kinematic = np.sqrt(M_P*M_X/(M_T*M_R)*ECM/(ECM+Qval))    
  kcm_i=np.sqrt(2.*MUI*amu*ECM)
  L_i = kcm_i/hbarc*(1.2*AT**(1./3.)) #Grazing orbital 
  E_threshold = - Qval*(M_R+ M_X)/(M_R+M_X-M_P) 
  ECMf=ECM+Qval
  kcm_f=np.sqrt(2.*MUF*amu*ECMf)
  L_f = kcm_f/hbarc*(1.2*AR**(1./3.))
  ELABf=ECMf*M_X/MUF
  EpAf=ELABf/AX  # E/A for x
  # Coulomb scale k_c and Sommerfeld parameter 
  k_C=ZP*ZT*alpha*(MUI*amu)
  eta_C = k_C/kcm_i 
  #--multiplication factor for S-factor  
  # S(E)=E*exp(2*pi*eta)*sigma(E) 
  Sfactor_convert = ECM*np.exp(2*np.pi*eta_C) # MeV 
  # Allowed angles and energies 
  s = (ELAB*(M_R-M_P)+M_R*Qval)/(M_X+M_R)
  rtilde = np.sqrt(M_P*M_X*ELAB)/(M_X+M_R) 
  theta_lab_X_min = 0.0 
  E_X_max = (rtilde+np.sqrt(rtilde**2+s) )**2
  if (s > 0.0):
     theta_lab_X_max = np.pi 
     E_X_min = (max(0,-rtilde+np.sqrt(rtilde**2+s) ))**2
  elif (s == 0.0): # or use very small number     
     theta_lab_X_max =  np.pi/2.0 
     E_X_min = 0.0 
  else : # s<0 
     theta_lab_X_max = np.arccos(np.sqrt(-s)/rtilde) 
     E_X_min = (max(0,rtilde-np.sqrt(rtilde**2+s) ) )**2 
  
  #---for target-like particle 
  E_R_min = Qval+ELAB-E_X_max 
  E_R_max = Qval+ELAB-E_X_min  
  s_R = (Qval*M_X  - ELAB*(M_P-M_X))/(M_X+M_R)
  rtilde_R = np.sqrt(M_P*M_R*ELAB)/(M_X+M_R)
  theta_lab_R_min= 0.0 
  E_R_max = (rtilde_R+np.sqrt(rtilde_R**2+s_R) )**2
  if (s_R > 0.0) :
     theta_lab_R_max = np.pi
     E_R_min = (max(0,-rtilde_R+np.sqrt(rtilde_R**2+s_R) ))**2
  elif (s_R== 0.0) :
     theta_lab_R_max =  np.pi/2.0 
     E_R_min = 0.0
  else :  
     theta_lab_R_max = np.arccos(np.sqrt(-s_R)/rtilde_R) 
     E_R_min = (max(0,rtilde_R-np.sqrt(rtilde_R**2+s_R) ) )**2 
  
  # equivalent kinematics
  text =('#===kinematics =====================\n' 
       + '#     with amu=931.4940954 MeV \n '
       + '**** For T(P,X)R reaction \n'
       + '(1) Masses from (AME16+Excitation energies) \n'
       + '   mP/amu    mT/amu      mX/amu    mR/amu   \n'  
       + ' %.6f %.6f %.6f %.6f \n'%(M_P,M_T,M_X,M_R)
       + '(2) Q value for T(P,X)R and reduced masses \n'
       + ' Qval(MeV) red_in/amu red_out/amu \n'
       + ' %8.3f %8.3f %8.3f \n' %(Qval,MUI,MUF)
       + ' E_threshold = %8.3f \n '%(E_threshold) 
       + '(3) Additional informations  \n'
       + '    Eplab    E/A   Tcm_in   Tcm_out     kcm(MeV)  kcm_out   L_gr \n'
       + ' %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n '%(ELAB,EpA,ECM,ECMf,kcm_i,kcm_f,L_i)
       + '    rho_kinematic  \n'
       + ' %.3e \n'%(rho_kinematic)  
       + '**** for S-factor \n'
       + '  k_C(MeV)     eta_C       C(MeV) = S(E)[MeV-mb=keV.b]/sigma(E)[mb] \n'
       + ' %8.3f  %8.3f  %.4e \n'%(k_C, eta_C, Sfactor_convert) 
       + '**** R(X,P)T equivalent \n'
       + '    Exlab    E/A   Tcm_out  kcm(MeV)    1/k(fm)   Lgr \n'
       + ' %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n'%(ELABf,EpAf,ECMf,kcm_f,hbarc/kcm_f,L_f) 
       )
  
  if (ELAB < E_threshold):
      text = text + '\n\n !! Energy is below threshold for the reaction !!\n'
      
  text = text +'**** In Lab frame: \n'
  text = text +' %8.3f <= theta_X <= %8.3f '%(theta_lab_X_min*180/np.pi
                             ,theta_lab_X_max*180/np.pi)
  text = text +' %8.3f <= E_X <= %8.3f \n'%(E_X_min,E_X_max) 
  text = text +' %8.3f <= theta_R <= %8.3f '%(theta_lab_R_min*180/np.pi
                             ,theta_lab_R_max*180/np.pi)
  text = text +' %8.3f <= E_R <= %8.3f \n'%(E_R_min, E_R_max) 
  
  return text 

def frame_convert(data,Tcm,Q,m1,m2,m3,m4,mode=1):
   """ convert angle and differential cross section
       between lab frame and cm frame 

       Note that this convert cross sections between   
       normal kinematics cross section in CM 
         <-> normal/inverse kinematics cross section in CM.
       Thus, it is not usual CM<->LAB convertor!!
        
       !!! mode= -2 can be wrong because of ambiguity. 
       
       input data are numpy array
       theta, sigma, dsigma = data[:,0], data[:,1], data[:,2]
  
       In C.M.,            
       theta_cm is always defined as an angle between 
         m1 and m3 (projectile like fragment)
       sigma(theta_cm) is a cross section for 

       2(1,3)4 reaction 
         masses in input are assumed for this configuration

       mode=1 : normal kinematics 2(1,3)4 
                sigma(theta_cm) -> sigma(theta_lab)
       mode=2 : inverse kinematics, 1(2,4)3 
                heavy in -> heavy detected  
                sigma(theta_{cm})->sigma(theta_heavy) 
       mode=3 : inverse kinematics, 1(2,3)4 
                heavy in -> light detected
                sigma(theta_{cm})->sigma(theta_light) 

       mode=-1 : convert lab to cm 
                 sigma(theta_lab)->sigma(theta_cm)
       mode=-2 : inverse kinematics
                 sigma(theta_heavy) -> sigma(theta_cm)
       mode=-3 : inverse kinematics
                 sigma(theta_light) -> sigma(theta_cm)                      
   """
   theta=data[:,0]
   sigma=data[:,1]
   if abs(mode)==1 :   # m1(m2,m3)m4
       rho=np.sqrt(m1*m3/m2/m4*Tcm/(Q+Tcm))
       if mode==1:
         theta_cm=theta*np.pi/180. # input is sigma(theta_cm)
         theta_lab=[math.atan2(math.sin(i),rho+math.cos(i)) for i in theta_cm]
         theta_lab=np.array(theta_lab)
         theta_p=theta_lab*180./np.pi # returning angle
       if mode==-1:
         theta_lab=theta*np.pi/180. 
         theta_cm=theta_lab+np.arcsin(rho*np.sin(theta_lab))
         theta_p=theta_cm*180./np.pi #returning angle 
       fac=np.sqrt((1.0+2.0*rho*np.cos(theta_cm)+rho**2)**3)/np.abs(1.0+rho*np.cos(theta_cm))
       if mode==-1:
         fac=1./fac
   if abs(mode)==2 :   # m2(m1,m4)m3
       rho=np.sqrt(m2*m4/m1/m3*Tcm/(Q+Tcm))
       if mode==2:
         theta_cm=theta*np.pi/180. # input is sigma(theta_cm)
         theta_lab=[math.atan2(math.sin(i),rho+math.cos(i)) for i in theta_cm]
         theta_lab=np.array(theta_lab)
         theta_p=theta_lab*180./np.pi
       if mode==-2:
         theta_lab=theta*np.pi/180. 
         theta_cm=theta_lab+np.arcsin(rho*np.sin(theta_lab))
         theta_p=theta_cm*180./np.pi 
       fac=np.sqrt((1.0+2.0*rho*np.cos(theta_cm)+rho**2)**3)/np.abs(1.0+rho*np.cos(theta_cm))
       if mode==-2:
         fac=1./fac
   if abs(mode)==3 :   # m2(m1,m3)m4
       rho=np.sqrt(m2*m3/m1/m4*Tcm/(Q+Tcm))
       if mode==3:
         theta_cm=theta*np.pi/180. # input is sigma(theta_cm)
         theta_lab=[math.atan2(math.sin(i),rho-math.cos(i)) for i in theta_cm]
         theta_lab=np.array(theta_lab)
         theta_p=theta_lab*180./np.pi # degree
       if mode==-3:
         theta_lab=theta*np.pi/180. 
         theta_cm= np.pi-(theta_lab+np.arcsin(rho*np.sin(theta_lab)))
         theta_p= theta_cm*180./np.pi 
       fac=np.sqrt((1.0-2.0*rho*np.cos(theta_cm)+rho**2)**3)/np.abs(1.0-rho*np.cos(theta_cm))
       if mode==-3:
         fac=1./fac
   try:
     dsigma=data[:,2]
     output=np.column_stack((theta_p, fac*sigma, fac*dsigma))
   except:
     output=np.column_stack((theta_p, fac*sigma))

   return output

def xs_Lab_Cm(data,rho,opt='cm_to_lab'):
    """
    For given data array ,
      data[:,0] are angles 
      data[:,1] are cross sections 
      data[:,2] are errors
    (the shape of data array must be (:,3)! ) 
        
    data can be in lab frame or cm frame 
    
    rho is the kinematic conversion factor 
     rho = sqrt(ma*mb/mA/mB*Tcm/(Tcm+Q) )
    """
    if not data.shape[1]==3:
        raise ValueError('data shape is wrong. shape must be (:,3)!')
    text = ""
    if rho > 1 and opt=='lab_to_cm': 
        text += "--------  WARNING -------------------------------------- \n"
        text += "  LAB -> CM conversion is not unique for this reaction. \n " 
        text += "  Two different CM angle can correspond to one LAB angle.\n"
        text += "  Two possibilities are both listed .\n"
        text += "--------------------------------------------------------  \n"
    if opt=='cm_to_lab':
        th_cm = data[:,0]*np.pi/180.
        th_lab = np.arctan2(np.sin(th_cm),(rho+np.cos(th_cm)))*180/np.pi  
        fac = (1.0+2.0*rho*np.cos(th_cm)+rho**2)**(1.5)
        fac = fac/np.abs(1.0+rho*np.cos(th_cm)) 
        
        text += "  cm angle    lab angle    factor    xs     err\n"
        for i in range(len(data)):
            text +=' {:6.2f} {:6.2f} {:.3e} {:.3e} {:.3e}\n'.format(data[i,0],
                                            th_lab[i],fac[i],fac[i]*data[i,1],
                                            fac[i]*data[i,2])
    if opt=='lab_to_cm': 
        th_lab = data[:,0]*np.pi/180.
        th_cm1 = th_lab+np.arcsin(rho*np.sin(th_lab))
        fac1 = (1.0+2.0*rho*np.cos(th_cm1)+rho**2)**(1.5)
        fac1 = fac1/np.abs(1.0+rho*np.cos(th_cm1)) 
        fac1 = 1.0/fac1 
        if rho > 1: 
            th_cm2 = th_lab+np.pi-np.arcsin(rho*np.sin(th_lab))
            fac2 = (1.0+2.0*rho*np.cos(th_cm2)+rho**2)**(1.5)
            fac2 = fac2/np.abs(1.0+rho*np.cos(th_cm2))            
            fac2 = 1.0/fac2 
            
        text +=" lab angle  cm angle  factor  xs   err \n"     
        for i in range(len(data)):
            text +=' {:6.2f} {:6.2f} {:.3e} {:.3e} {:.3e}\n'.format(data[i,0],
                                th_cm1[i]*180/np.pi,fac1[i],fac1[i]*data[i,1],
                                        fac1[i]*data[i,2]) 
            if rho > 1:
                text +=' {:6.2f} {:6.2f} {:.3e} {:.3e} {:.3e}\n'.format(
                            data[i,0],th_cm2[i]*180/np.pi,
                            fac2[i],fac2[i]*data[i,1],fac2[i]*data[i,2] ) 
    return text 


def frame_convert2(data, m_a, m_A,m_b, m_B, E_a, en_opt=0 , mode=0 ):
    """
    For a reaction  A(a,b)B,
    convert angle and cross section between lab and cm frame       
    input: 
       masses are in amu units, 
       energy in MeV unit  can be 'lab','cm','E/A'

      en_opt= 0 Lab 
              1 E/A
              2 CM  
              
       mode=0 converts
        data=(theta_lab, sigma_lab) into (theta_cm, sigma_cm) 
        
       mode=1 converts
        data=(theta_cm, sigma_cm) into (theta_lab, sigma_lab) 
        
     Note CM energy is total relative energy in CM not only projectile energy.
     
     Note if rho > 1. as inverse kinematics. 
          lab angle -> cm angle have ambiguity. Thus, can be wrong. 
    """
    amu = 931.4940954 #MeV
    if en_opt==0:
        E_a = E_a
        E_cm = m_A/(m_a+m_A)*E_a 
    elif en_opt==1 :
        E_a = E_a*int(round(m_a)) # nearest integer 
        E_cm = m_A/(m_a+m_A)*E_a     
    elif en_opt==2 :
        E_cm = E_a
        E_a = E_cm*(m_a+m_A)/m_A           
    Qval =(m_a+m_A-m_b-m_B)*amu # MeV 
    rho = np.sqrt( m_a*m_b/(m_A*m_B)*E_cm/(Qval+E_cm)   )
    #---treat input    
    data = np.array(data)  
    data_type = data.shape[0] # number of input category 1= angle, 2= angle and sigma
    if data.shape == (data_type,) : # only one number 
        angles = data[0]
        if data_type > 1:
            sigma = data[1]
        if data_type > 2:
            error = data[2]
    else : # multiple angles 
        angles = data[:,0]
        if data_type > 1:
            sigma = data[:,1]
        if data_type > 2:
            error = data[:,2]
    #-- conversion         
    if mode==0 : # from lab to cm 
        if (rho>1.0):
            print('WARNING: rho>1 case lab->cm is ambiguous and can be wrong!\n')
        theta_lab = angles*np.pi/180. # randian  
        # This expression actually can be not unique.. 
        theta_cm = np.arccos(np.cos(theta_lab)*(rho*np.cos(theta_lab)
                                +np.sqrt(1.0-rho**2*np.sin(theta_lab)**2))-rho)
        #--conversion factor sigma_lab = fac * sigma_cm 
        fac = (1.0+2.0*rho*np.cos(theta_cm)+rho**2)**1.5/np.abs(1.0+rho*np.cos(theta_cm))
        fac = 1.0/fac # from lab to cm 
        sigma_lab = sigma 
        sigma_cm = fac*sigma_lab 
    if mode==1 : # from cm to lab 
        theta_cm = angles*np.pi/180.  # radian 
        theta_lab = np.arctan2( np.sin(theta_cm), (rho+np.cos(theta_cm)) ) #radian 
        #--conversion factor sigma_lab = fac * sigma_cm 
        fac = (1.0+2.0*rho*np.cos(theta_cm)+rho**2)**1.5/np.abs(1.0+rho*np.cos(theta_cm))
        sigma_cm = sigma 
        sigma_lab = fac*sigma_cm 
    
    if data_type==1: 
        output=np.column_stack((theta_lab*180/np.pi,theta_cm*180/np.pi))
    elif data_type==2:
        output=np.column_stack((theta_lab*180/np.pi,theta_cm*180/np.pi,fac*sigma))
    elif data_type==3:
        output=np.column_stack((theta_lab*180/np.pi,theta_cm*180/np.pi,fac*sigma,fac*error))    
    return output     

#-----------------------------------------------------------------------------
def WLH_global_optical_para(zp,at,zt,elab ):
    """
    Whitehead-Lim-Holt global optical potential 
    arxiv: 2009.08436

    Parameters
    ----------
    zp: int
        charge of nucleon either 0 or 1.
    at : integer
        mass number of target 
    zt : int
        charge number of target
    elab : float
        energy of nucleon
    Returns
    -------
    global optical potential parameters 

    """
    # parameter sets # eq.(4)
    # para[0] os for neutron ,para[1] is for proton 
    para_uV = [ [53.459,-0.2356,-0.00133,1.317e-5,-2.88e-8,-20.58,0.317,
         -0.00158,3.49e-6,-10.96,-0.0155],
        [54.154,-0.252,-0.0011,1.19e-5,-2.6e-8,20.87,-0.306,
         0.00172,-4.46e-6,-21.92,-0.01035]   ]  
    para_rV = [ [1.298,-5.41e-4,1.98e-6,-0.397],
        [1.310,-5.09e-4,1.98e-6,-0.391] ]
    para_aV = [
        [0.699,0.0023,-3.77e-5,2.38e-7,-5.14e-10,-0.171,0.625],
        [0.773,-1.28e-4,-6.1e-6,5.39e-8,-1.67e-10,-0.746,0.522]
        ]
    para_uW = [
        [2.23,0.262,-5.85e-4,-7.84,-0.046],
        [3.80,0.237,-4.67e-4,11.9,-0.077]
        ]
    para_rW = [
        [0.450,86.82,0.830,84.17,0.755,2.34e-6],
        [0.543,48.78,0.708,54.23,0.499,9.93e-7]
        ]
    para_aW =[
        [0.546,-0.287,-17.1,0.376,-0.00187],
        [0.407,-0.394,-9.39,0.0599,-0.853]
        ]
    para_uS =[
        [1.866,-0.0546,1.41e-4,-3.70,0.4218,-0.00842],
        [0.774,-0.0216,-1.281,0.0316,0.0,0.0]
        ]    
    para_rS =[
        [1.238,0.00275,-2.14e-4,-1.266],
        [1.134,-0.0120,-3.25e-8,-0.836]
        ]
    para_aS =[
        [0.578,0.0169,-3.82e-4,5.47e-4,3.91e-6],
        [0.146,0.00362,-3.47e-4,0.019,-1.32e-4]
        ]
    para_uSO = [
        [8.852,0.0127,-2.02e-4,5.57e-7],
        [8.852,0.0127,-2.02e-4,5.57e-7]
        ]
    para_rSO = [[1.260,-0.827],
                [1.260,-0.827]]
    para_aSO =[[0.663,0.0032,-2.83e-5,6.58e-8],
               [0.663,0.0032,-2.83e-5,6.58e-8]]
    # compute each parameters for E,A,Z eq.(11-14)
    AT = at
    ZT = zt
    delta =  (AT-2.0*ZT)/AT # (N-Z)/A
    U_V = ( para_uV[zp][0]+ para_uV[zp][1]*elab+ para_uV[zp][2]*elab**2
         + para_uV[zp][3]*elab**3+para_uV[zp][4]*elab**4
         +(para_uV[zp][5]+para_uV[zp][6]*elab
           +para_uV[zp][7]*elab**2+para_uV[zp][8]*elab**3)*delta
         +para_uV[zp][9]*np.exp(para_uV[zp][10]*elab)*delta**2 )
    r_V = ( para_rV[zp][0]+ para_rV[zp][1]*elab+para_rV[zp][2]*elab**2
          +para_rV[zp][3]*AT**(-1./3.)    )
    
    
    optical_pot_para={} 
    return optical_pot_para 


#=======================================================================
""" This routine only works for a specific cutoff value """
if __name__ == '__main__':
   #import matplotlib.pyplot as plt
   import numpy as np
   import sys 
   from subprocess import call
   #  
   method=None
   while (method!='st'):
     print('#== Reaction helper ====')
     print('') 
     print('#   kin   =   kinematics of two-body reaction ' ) 
     print('#   frame =  conversion of cross section between frame')
     print('#   TO BE ADDED                               ' )
     print('#   st    =   exit shell                      ' )
     print('')
     method=input('*** Enter the choice of method :\n')
     method=rmcomm(method,num_arg=1)
     if (method=='kin'):
       kin()
     if (method=='frame'):
       print('#--------------------------------------')
       print('# assuming m2(m1,m3)m4 normal kinematics')
       print('#  i.e. light particle 1 in -> light particle 3 out')
       print('#--------------------------------------')
       inp=input('Enter m1,m2,m3,m4 :\n')
       ww = rmcomm(inp,num_arg=4)
       m1,m2,m3,m4 = [ float(i)  for i in ww[0:4]]
       inp=input('# enter Tcm_in and reaction  Q-value (MeV):\n')
       ww = rmcomm(inp,num_arg=2)
       Tcm ,Q = [float(i) for i in ww[0:2]]
       print('#--------------------------')
       print('# mode>0 : sigma_cm -> sigma_lab ')
       print('# mode<0 : sigma_lab -> sigma_cm ')
       print('#  ')
       print('# |mode|=1 : normal kinematics m2(m1,m3)m4 ')
       print('# |mode|=2 : inverse kinematics  ')
       print('#            heavy in, heavy out , m1(m2,m4)m3 ')
       print('# |mode|=3 : inverse kinematics ')
       print('#            heavy in, light out , m1(m2,m3)m4 ')
       print('#-------------------------')
       inp=input('# choose the mode :\n')
       ww=rmcomm(inp,num_arg=1)
       mode=int(ww) 
       inp=input('# enter input data file name:\n')
       ww=rmcomm(inp,num_arg=1)
       dat=np.loadtxt(ww)
       out=frame_convert(dat,Tcm,Q,m1,m2,m3,m4,mode=mode)       
       inp=input('# enter output file name:\n')
       ll=rmcomm(inp)
       if (len(out[0,:])==2):
          np.savetxt(ll,out,fmt=['%15.6E','%15.6E'])
       if (len(out[0,:])==3):
          np.savetxt(ll,out,fmt=['%15.6E','%15.6E','%15.6E'])
   print('Bye~~') 
