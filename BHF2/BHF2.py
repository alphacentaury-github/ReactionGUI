# -*- coding: utf-8 -*-
"""
Last modified June,15th 2022
@author: Y.-H. Song

Program to solve Brueckner HF equation 
      in symmetric nuclear matter  
      neutron matter calculation is not complete. 
"""
#===imports===========================================
import numpy as np
import matplotlib.pyplot as plt 
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
import json 
import time    
#import scipy 
#import scipy.integrate
#import scipy.linalg 
#import scipy.interpolate
#from scipy import interpolate
#import scipy.optimize 
#import cProfile
import bonn_py

#===Global Parameters================================
hbarc =197.3269788 # MeV.fm
alpha_em =1./137.035999 # fine structure constrant 
amu = 931.4940954 # MeV/c^2 
mpi = 138.03898   # MeV/c^2 isoscalar pion mass 
mN  = 938.9187       # MeV/c^2 nucleon mass
#===Global Options===================================

neutron_matter = 0 # 0 : symmetric nuclear matter 
                   # 1 : pure neutron matter , incomplete  
                   
pot_choice = 1 # 0 : 1S0 separable 
               # 1 : bonn potential (exact type bonn depends on the Bonn input file)
               # 2 : potential library by C. Elster 
bonn_choice = 4 
		#  0: Bonn A ps , LS 
		#  1: Bonn B ps , LS
		#  2: Bonn C ps , LS
		#  3: Bonn A pv , LS
		#  4: Bonn B pv , LS 
		#  5: Bonn C pv , LS    
        # Note that one cannot change bonn_choice once loaded  
        # unless starting a new kernel !             

scp_choice = 1 # 0 : standard choice + effective mass approximation
               # 1 : continuous choice      

effective_mass= [0.8,80.0] # initial guess global effective mass variable
effective_mass_def = 1 # 0 or 1       
                
Jmax = 1 # channels to include 
Nk = 12  # number of quadrature for the integration to get U(p)  
Nq=40    # number of quadrature for the integration to get g(q,k)  
iter_max= 6 # number of iteration to compute U(p) 

#====Utility function=======================================
def make_channel_list(Jmax,Jmin=0,save=False):
    """
    construct list of NN channels up to Jmax 
    for isospin (T,T_z)
    
    return a dictionary 
    dict ={'J': array
           'S': array
           'L': array 
           'T': array
           'icoup': array}
    
    when two i and i+1 channels are icoup icoup becomes non-zero 
    icoup[i] = i+1    icoup[i+1] = -i
    |icoup| represent icoup other channel index  
    """
    #--1S0 and 1P0
    if Jmin==0:
        J_list=[0,0]
        S_list=[0,1]
        L_list=[0,1]
        T_list=[1,1]
        coupl_list=[0,0]
        indx=1
        Jmin=1 # start from J=1 
    else:
        J_list=[]
        S_list=[]
        L_list=[]
        T_list=[]
        coupl_list=[]
        indx=0
    for J in range(Jmin,Jmax+1):
        for S in range(1+1):
            for L in range(J,J-S-1,-1):
                if ((L+S)%2)==0: #even number 
                    T=1
                else:
                    T=0
                if L==J:
                    indx = indx +1 
                    J_list.append(J)
                    S_list.append(S)
                    L_list.append(L)
                    T_list.append(T)
                    coupl_list.append(0)
                if L==J-1: #icoup channel case
                    # add L=J-1 channel
                    indx = indx+1 
                    J_list.append(J)
                    S_list.append(S)
                    L_list.append(J-1)
                    T_list.append(T)
                    coupl_list.append(indx+1) 
                    # add L=J-1 channel 
                    J_list.append(J)
                    S_list.append(S)
                    L_list.append(J+1)
                    T_list.append(T)
                    coupl_list.append(-indx)
                    indx = indx+1 
    channels = {'J':J_list, 'S':S_list,'L':L_list,'T':T_list
                ,'icoup': coupl_list  }
    
    if save:
        with open('channels.txt','w') as ff:
            ff.write(json.dumps(channels))
    return channels 
                    
# gauss lengendre quadrature (-1,1)
# leggau(n) function can be used
def leggau(n,xmin=-1,xmax=+1):
    """ Gaussian quadrature in range (xmin,xmax)
        from quadrature (-1,1)
    """ 
    x,w = np.polynomial.legendre.leggauss(n)
    xp=(xmax-xmin)*0.5*(x+1.0)+xmin
    wp=(xmax-xmin)*0.5*w 
    return xp,wp

def read_output_Amos(fname='nntgmat_Loc.out'):
    """
    Read scp potential output
    separated by 'Iteration'
    """
    ff = open(fname,'r')
    ll = ff.readlines()
    ff.close()    
    # find positions
    positions=[]
    for i in range(len(ll)):
        if ('Initial input' in ll[i]) or ('Iteration ' in ll[i]):
            positions.append(i) 
    length = positions[1]-positions[0]        
    # convert into dictionary     
    out={}
    for i in range(len(positions)):
        out[i] = [] 
        for j in range(positions[i]+1,positions[i]+length):
            ww=ll[j].split() 
            out[i].append([float(k) for k in ww])
        out[i]=np.array(out[i])
    return out     

#=====NN interaction======================================

def NN_interaction(k,kp,L=0,Lp=0,S=0,J=0,T=0,Tz=0,opt=0):
    """
    potential in momentum space, 
    < Lp,kp|V^{S,J}|L, k> 
    (spin conserving interaction at the moment)
    k and kp are momentum in fm^-1 units. 
    <V> is returned in fm power
    !
    ! Bonn potential case 
    ! v(1) = V_{J,J}^{S=0}(kp,k)
    ! v(2) = V_{J,J}^{S=1}(kp,k)
    ! v(3) = V_{Lp=J+1,L=J+1}^{S=1}(kp,k)
    ! v(4) = V_{Lp=J-1,L=J-1}^{S=1}(kp,k)
    ! v(5) = V_{Lp=J+1,L=J-1}^{S=1}(kp,k)
    ! v(6) = V_{Lp=J-1,L=J+1}^{S=1}(kp,k)
    !
    ! 
    !Need to check again for the ordering!!
    """
    if opt==0: # simple separable potential eq.(3.15) in Haftel-Tabakin
        if (L==0 and Lp==0 and S==0 and J==0): #1S0 only 
            alph=3.6; a=1.2;
            g0=alph/(k**2+a**2)
            g0p=alph/(kp**2+a**2)
            return -g0*g0p 
        else:
            return 0.0 
    elif opt==1: # bonn potential 
        xmev = kp*hbarc 
        ymev = k*hbarc 
        vv = bonn_py.bonn_wrap(xmev,ymev,J,bonn_choice) #vv is MeV^-2 unit 
        # 0v, 1v, v++, v--, v+-, v-+ (lsj formalism) order 
        #--unit conversion to fm , (m*pi/2)
        vv = vv*(mN*np.pi/2*hbarc) 
        if (L==J and Lp==J and S==0): # singlet 
            return vv[0]
        elif (L==J and Lp==J and S==1): # triplet uncoupled, careful for python array index  
            return vv[1]
        elif (L==J-1 and Lp==J-1 and S==1): 
            return vv[3]
        elif (L==J-1 and Lp==J+1 and S==1):
            return vv[4]
        elif (L==J+1 and Lp==J-1 and S==1):
            return vv[5]
        elif (L==J+1 and Lp==J+1 and S==1):
            return vv[2]
        else:
            raise ValueError('unphysical partial waves in interaction')
    else :
        raise ValueError('ERROR: interaction choice is not available ')

def NN_interaction2(k,kp,J=0,T=0,Tz=0,opt=0):
    """
    Similar to NN_interaction 
    but, instead of return one potential matrix element,
    return all 6 potential matrix elements for given total angular momentum J 
    just like Bonn potential.     
    ! Bonn potential case 
    ! v(1) = V_{J,J}^{S=0}(kp,k)
    ! v(2) = V_{J,J}^{S=1}(kp,k)
    ! v(3) = V_{Lp=J+1,L=J+1}^{S=1}(kp,k)
    ! v(4) = V_{Lp=J-1,L=J-1}^{S=1}(kp,k)
    ! v(5) = V_{Lp=J+1,L=J-1}^{S=1}(kp,k)
    ! v(6) = V_{Lp=J-1,L=J+1}^{S=1}(kp,k)
    # 0v, 1v, v++, v--, v+-, v-+ (lsj formalism) order 
    """
    vv = np.zeros(6)
    if opt==0: # simple separable potential eq.(3.15) in Haftel-Tabakin
        if J==0:
            alph=3.6; a=1.2;
            g0=alph/(k**2+a**2)
            g0p=alph/(kp**2+a**2)
            vv[0] = -g0*g0p 
    elif opt==1: # bonn potential 
        xmev = kp*hbarc 
        ymev = k*hbarc 
        vv = bonn_py.bonn_wrap(xmev,ymev,J,bonn_choice) #vv is MeV^-2 unit 
        #--unit conversion to fm , (m*pi/2)
        vv = vv*(mN*np.pi/2*hbarc) 
    else :
        raise ValueError('ERROR: interaction choice is not available ')
    return vv 
#======potential matrix preparation=========================================
def make_Vmat(xq_ext,ich,chan_dic,Vopt=0):
    """
    create potential matrix for g-matrix equation. 
    
    xq_ext : extended momentums in fm^-1  
            quadratures+on-shell momentum 
            [q_i, k0] N_q+1 dimension
    ich : channel index         
    chan_dic : dictionary of channel informations 
    Vopt : choice of potential  
        
    NN_interaction is used with potential choice opt  
    """
    J = chan_dic['J'][ich]
    L = chan_dic['L'][ich]
    S = chan_dic['S'][ich]
    T = chan_dic['T'][ich]
    Nq = len(xq_ext)-1
    if (chan_dic['icoup'][ich]==0): # uncoupled channel
        #---construct potential matrix 
        Vmat = np.zeros( (Nq+1 ,Nq+1 ) )+1j*0.0
        for i in range(Nq+1):
            for j in range(i,Nq+1):
                 Vmat[j,i] = NN_interaction(xq_ext[i],xq_ext[j]
                                         ,L=L,Lp=L,S=S,J=J,T=T,opt=Vopt)
                 Vmat[i,j] = Vmat[j,i] 
    elif (chan_dic['icoup'][ich] > 0): # icoup channel 
        #---extend 2Nq+2 dim matrix  
        xq_ext = np.concatenate( (xq_ext,xq_ext) )
        Vmat = np.zeros( (2*Nq+2 ,2*Nq+2 ) )+1j*0.0
        for i in range(Nq+1):
            for j in range(i,Nq+1):
                Vmat[j,i] = NN_interaction(xq_ext[i],xq_ext[j]
                                   ,L=L,Lp=L,S=S,J=J,opt=Vopt)
                Vmat[i,j] = Vmat[j,i] 
                #=== V_{Lp,L}(kp,k)
                #  careful that NN_interaction is defined with arguments 
                #               (k,kp,L,Lp)
                Vmat[j+Nq+1,i] = NN_interaction(xq_ext[i],xq_ext[j]
                                   ,L=L,Lp=L+2,S=S,J=J,T=T,opt=Vopt)
                Vmat[i,j+Nq+1] = Vmat[j+Nq+1,i]
                
                Vmat[j,i+Nq+1] = NN_interaction(xq_ext[i],xq_ext[j]
                                   ,L=L+2,Lp=L,S=S,J=J,T=T,opt=Vopt)
                Vmat[i+Nq+1,j] = Vmat[j,i+Nq+1]
             
                Vmat[j+Nq+1,i+Nq+1] = NN_interaction(xq_ext[i],xq_ext[j]
                                   ,L=L+2,Lp=L+2,S=S,J=J,T=T,opt=Vopt)
                Vmat[i+Nq+1,j+Nq+1] = Vmat[j+Nq+1,i+Nq+1]     
    elif (chan_dic['icoup'][ich] < 0 ):
        raise ValueError('Vmat is requested for icoup <0.')
    return Vmat     

def make_Vmat_dic(xq_ext,chan_dic,Vopt=0):
    """
    construct potential dictionary 
    
    Vmat_dic[(p0,k0,ich)] = Vmat  

    k0 = xq_ext[-1] 
    
    Size of Vmat is len(xq_ext) for uncoupled channel
                   2*len(xq_ext) for coupled channel
    """
    Vmat_dic = {} 
    for ich in range(len(chan_dic['J'])):
        if chan_dic['icoup'][ich]>=0:
            Vmat_dic[ich] = make_Vmat(xq_ext,ich,chan_dic,Vopt=Vopt)
    return Vmat_dic 

def make_Vmat_dic2(xq_ext,chan_dic,Vopt=0):
    """
    Instead of iterating over channel index, 
    iterate over J and use NN_interaction2 
    
    ! Bonn potential case 
    ! v(1) = V_{J,J}^{S=0}(kp,k)
    ! v(2) = V_{J,J}^{S=1}(kp,k)
    ! v(3) = V_{Lp=J+1,L=J+1}^{S=1}(kp,k)
    ! v(4) = V_{Lp=J-1,L=J-1}^{S=1}(kp,k)
    ! v(5) = V_{Lp=J+1,L=J-1}^{S=1}(kp,k)
    ! v(6) = V_{Lp=J-1,L=J+1}^{S=1}(kp,k)
    !
    ! Be careful that 3P0 channel corresponds to v(4) not v(2) !
    """
    Vmat_dic = {} 
    J_list = np.unique(chan_dic['J'])
    Nq = len(xq_ext)-1
    #---initialize Vmats 
    for ich in range(len(chan_dic['J'])):
        if chan_dic['icoup'][ich]==0:
            Vmat_dic[ich] = np.zeros( (Nq+1 ,Nq+1 ) )+1j*0.0
        elif chan_dic['icoup'][ich]> 0:    
            Vmat_dic[ich] = np.zeros( (2*Nq+2 ,2*Nq+2 ) )+1j*0.0
    #---iterate over J, k ,kp 
    for J in J_list: 
        for i in range(Nq+1):
            for j in range(i,Nq+1):
                # <kp|V|k> 
                kk = xq_ext[j]
                kk_p = xq_ext[i]
                vv = NN_interaction2(kk,kk_p,J,opt=Vopt) 
                #---assign vv(6) to Vmat 
                #---Note! the particular structure of chan_dic is used. Not general. 
                ich = chan_dic['J'].index(J) # find first channel index for J 
                #---uncoupled channels 
                if J==0: 
                    #special case 3P0 is triplet uncoupled,however, with L=J-1 not L=J. 
                    Vmat_dic[ich][i,j]= vv[0] 
                    Vmat_dic[ich+1][i,j] = vv[2] 
                    
                    Vmat_dic[ich][j,i]=  vv[0]
                    Vmat_dic[ich+1][j,i] = vv[2]
                else:
                    Vmat_dic[ich][i,j]= vv[0] 
                    Vmat_dic[ich+1][i,j] = vv[1]
                    
                    Vmat_dic[ich][j,i]=  vv[0]
                    Vmat_dic[ich+1][j,i] = vv[1]
                #---coupled channels 
                if (J>0): 
                    Vmat_dic[ich+2][i,j] = vv[3] #< J-1,1|V(kp,k)| J-1,1>
                    Vmat_dic[ich+2][i+Nq+1,j+Nq+1] = vv[2] #< J+1,1|V(kp,k)| J+1,1>
                    Vmat_dic[ich+2][i+Nq+1,j] = vv[4] #< J+1,1|V(kp,k)| J-1,1>
                    Vmat_dic[ich+2][i,j+Nq+1] = vv[5] #< J-1,1|V(kp,k)| J+1,1>
                    
                    Vmat_dic[ich+2][j,i] = vv[3]
                    Vmat_dic[ich+2][j+Nq+1,i+Nq+1] = vv[2]
                    Vmat_dic[ich+2][j+Nq+1,i] = vv[5]
                    Vmat_dic[ich+2][j,i+Nq+1] = vv[4] 
    return Vmat_dic             


#=======BHF calculation=====================================================
def pauli(Kcm,kp,kF):
    # angle averaged Q 
    #eq.(3.9) in REF
    #
    #input K,kp,kF are in fm**-1 unit
    # kF is float
    # kp , Kcm can be array but with the same shape 
    #         
    kp = np.atleast_1d(kp) #vectorized for kp 
    if (kF < 1.e-6):
        return np.ones(len(kp))
    qbar = (Kcm**2+kp**2-kF**2)/(2.0*Kcm*kp)    
    # conditional new assignment 
    qbar[kp**2 <= kF**2-Kcm**2] = 0.0
    qbar[(kp-Kcm)**2 >= kF**2] = 1.0
    if np.any(qbar >1):
        raise ValueError('qbar >1 !')
    return qbar     

def interpolate_U(x,p0_array,U_array,xq):
    """
    interpolate U(p) at U(x) 
    from U(p_list) 
    return U(x) and U'(x) 

    important!!
    
    How the extrapolation is done is crucial 
    for the calculation. 
    """
    if x <= p0_array[0] : #extrapolation  x < p0_min 
        x0 = p0_array[0]     #linear approximation 
        x1 = p0_array[1]
        y0 = U_array[0]
        y1 = U_array[1]
        dUdx = (y0-y1)/(x0-x1)
        Ux = dUdx*(x-x0)+y0          
    elif x >= p0_array[-1] : # extrapolation x > p0_max 
        #----method 1
        # assume U(xq[-1]) = 0 and linear extrapolation. 
        dUdx = U_array[-1]/(p0_array[-1] - xq[-1]) 
        Ux = U_array[-1]+dUdx*(x-p0_array[-1])  
        #----method 2
        #Ux = 0.0   # take zero   
        #dUdx = 0.0  
        #----method 3 cubic spline           
        #interp_f = CubicSpline(p0_array,U_array )
        #Ux = interp_f(x)
        #dUdx = interp_f(x,1)        
    else : # spline interpolation  
        interp_f = CubicSpline(p0_array,U_array )
        Ux = interp_f(x)
        dUdx = interp_f(x,1)
    return Ux, dUdx     
            
def get_effective_mass(k_F,p0_array,U_array,method=0,step_size=0.1):
    """
    determine ratio (= m*/m) and U0 from scp U(p) 
    
    p^2/2m+U(p) = p^2/(2 m*) -U0 

    How to fix these two ?
    possible methods are 
    (0) match at the U(k_F) and U(k_F-step_size)   
    (1) fit ratio and U0 in range of (p0_array[0],k_F)          
    """
    U_func = CubicSpline(p0_array,np.real(U_array) )
    if method==0:
        p1 = k_F 
        p2 = (k_F-step_size)
        U1 = U_func(p1)
        U2 = U_func(p2) 
        U0 = (p2**2*U1-p1**2*U2)/(p1**2-p2**2)   
        ratio = 1.0/(1+2*mN/(p1*hbarc)**2*(U1+U0))
    elif method==1:
        def U_fit(p,ratio,U0):
            return (p*hbarc)**2/(2*mN)*(1./ratio-1.)-U0 # fm^-1 unit 
        p = np.arange(p0_array[0],k_F,0.1) # only up to k_F 
        popt, pcov = curve_fit(U_fit,p,U_func(p)) #fit 
        ratio = popt[0]
        U0 = popt[1]
    else:
        raise ValueError('Unknown option for effective mass')
    return [ratio, U0 ]

def energy_2N_U(Kcm,q, kF,p0_array,U_array,xq):
    """
    two nucleon energy with auxiliary potential U
    
    E(krel,Kcm) = krel^2+ (m/hbacc^2)*( U(|k+Kcm|)+U(|k-Kcm|) )
    where |k+/-Kcm| is obtained Brueckner angle averaged prescription.
    (REF: Amos paper )
    
    * Ignore Kcm part which always cancel. 
    * nucleon mass is multiplied 
      E_{2N} is fm^-1 originally. 
      However, following the H-T convention, 
      one have to use mN*E_{2N} instead of E_{2N}
      Thus, the final output here is in unit of fm^-2.
    * Only real part of  U(p) for p < k_F 
       
    * Kav,q,kF are in fm^-1 units
    * U is in MeV units 
    * Note : q can be array ,but Kav must be a number 
    """      
    interp_U = lambda x : interpolate_U(x,p0_array,U_array,xq)[0] 
    q = np.atleast_1d(q)
    temp = q**2 #kinetic energy term 
    if (kF < 1.e-6): # free particle case 
            return temp
        #--nuclear matter angle average approximation 
    Qbar = pauli(Kcm,q,kF)
    q_sq_plus = Kcm**2+q**2+2./np.sqrt(3.)*q*Kcm*Qbar**(1.5) #|k+K|^2
    q_sq_mnus = Kcm**2+q**2-2./np.sqrt(3.)*q*Kcm*Qbar**(1.5) #|k-K|^2
    if np.any(q_sq_mnus < 0.0):
        print('ERROR ! q_sq_mnus <0 !')
        for i in range(len(q)):
            print('(Kav,q,q_mnus,Qbar^1.5/sqrt(3))=({},{},{},{})'.format(
            Kcm,q[i],q_sq_mnus[i],Qbar[i]**(1.5)/np.sqrt(3.0) ) )       
        raise ValueError(' q_sq_mnus <0 !')
        
    q_plus = np.sqrt(q_sq_plus)
    q_mnus = np.sqrt(q_sq_mnus)
    U_plus = np.array([ interp_U(i) for i in q_plus])   
    U_mnus = np.array([ interp_U(i) for i in q_mnus])
    
    #------ remove imaginary part when p <k_F
    U_plus[q_plus <= kF] = np.real(U_plus[q_plus <= kF] )
    U_mnus[q_mnus <= kF] = np.real(U_mnus[q_mnus <= kF] )
    # continuous definition of U(p) 
    E_2N = temp + (U_plus + U_mnus)*(mN/hbarc**2) #  HT convention fm^-2 unit 
    return E_2N  # HT convention, fm^-2 unit  

def construct_k0_quadrature(p_0,k_F,nleg=40,opt=0):
    """
    prepare quadratures for k0 integral 
    in calculation of auxiliary potential from g-matrix 

    p_0 : incident nucleon momentum, fm^-1 unit  
          not array  
    k_F : Fermi momentum , fm^-1 unit 
          not array 
    nleg : number of quadratures for the integration       
          use even number for simplicity  
          
    return       
      quadrature xk(:), wXk(:) such that 
      
    \int_0^{(k_F+p_0)/2} dk k^2 X(k) f(k) = sum_{i} wXk(i) *f(xk(i))              
    
      also avergaed cm momentum at quadratures K_av(xk,p_0,k_F)
    """
    k_max = 0.5*(k_F+p_0)
    k_mid = 0.5*np.abs(k_F-p_0) # for both p_0> k_F and p_0< k_F !
    
    N1 = int(nleg*(k_mid)/(k_max))      
    if (p_0< 1.e-5):
        N1=nleg
    if (np.abs(k_F-p_0) <1.e-5):
        N1=0
    if (p_0 > k_F):
        N1=0       
    N2 = nleg - N1
    
    if (N1 > 0):
        if opt==0: # get N1 quadrature from (0,A)
            xk1, wk1 = leggau(N1,xmin=0.0,xmax=k_mid)
        else : # get N1 quadrature from 2N of (-A,A)  
            xk1, wk1 = leggau(2*N1,xmin=-k_mid,xmax=k_mid)
            xk1 = xk1[N1:2*N1]
            wk1 = wk1[N1:2*N1]
            
        wXk1 = 1.0*xk1**2*wk1 # X(k)=1.0 
        Kav1 = np.sqrt( p_0**2+ xk1**2)
    if (N2 > 0):
        if opt==0:
            xk2, wk2 = leggau(N2,xmin=k_mid,xmax=k_max)
        else: 
            xk2, wk2 = leggau(2*N2,xmin=-1,xmax=1) #doubled
            xk2 = xk2[N2:2*N2]*(k_max-k_mid)+k_mid #returned to N2  
            wk2 = wk2[N2:2*N2]*(k_max-k_mid)
            
        wXk2 = (0.25*(k_F**2-p_0**2)-xk2*(xk2-p_0))/(2.*p_0)*xk2*wk2        
        Kav2 = np.sqrt( p_0**2+xk2**2-0.25*(2*xk2+p_0-k_F)*(2*xk2+p_0+k_F) )       
    if (N1==nleg):
        xk = xk1 
        wXk = wXk1 
        Kav = Kav1 
    elif (N2==nleg):
        xk = xk2
        wXk = wXk2 
        Kav = Kav2 
    elif (N1 > 0 and N2 >0):
        xk = np.concatenate((xk1,xk2)) # concatenate xk1 and xk2 
        wXk = np.concatenate((wXk1,wXk2)) #concatenate Xwk1 and Xwk2
        Kav = np.concatenate((Kav1,Kav2))       
    return xk, wXk, Kav       

def K_av(k0,p_0,k_F):
    # average cm momentum 
    k_mid = np.abs(k_F-p_0)/2.0
    k_max = (k_F+p_0)/2.0 
    if (p_0 > k_F):
        if (k0 >=  k_mid and k0 <= k_max ):
            Kav = np.sqrt( p_0**2+k0**2-0.25*(2*k0+p_0-k_F)*(2*k0+p_0+k_F) )
        else:
            raise ValueError('k0 is not in range (k_mid,k_max)')
    else: 
        if (k0 < k_mid ):
            Kav = np.sqrt(p_0**2+k0**2)
        elif (k0 >=  k_mid and k0 <= k_max ):
            Kav = np.sqrt( p_0**2+k0**2-0.25*(2*k0+p_0-k_F)*(2*k0+p_0+k_F) )
        else:
            raise ValueError('k0 is not in range (0,k_max)')   
    return Kav 


def construct_q_quadrature(nleg=40):
    """
    prepare q integration for the g-matrix 
    
    use change of variable for infinite range integral 

    \int_{0}^\infty dq f(q)
     \int_0^{pi/2} dz C/(cos^2 z) f(q= C tan(z) )
    """
    xk, wk = leggau(nleg,xmin=0.0,xmax=np.pi/2.0)
    Ck = 3.0 # somewhat arbitrary 
    qq = Ck*np.tan(xk) 
    ww = Ck*wk/np.cos(xk)**2
    return qq, ww 

def construct_q_quadrature2(Kcm,k_F,nleg=40,
                            n_1=0,n_2=0,
                            q_1=2.0,q_2=20.0):
    """
    same as construct_q_quadrature
    
    except the range of integration have q_min 
    and whole integration range is splitted by three region 
    
    \int_{qmin}^\infty dq 
     = \int_{qmin}^{q_1}dq + \int_{q_1}^{q_2} dq + \int_{q_2}^\infty dq 
    
    each region have (n_1,n_2,ntot-n1-n2) number of quadratures
    
    The last region can be obtaine by change of variable 
    \int_{q_min}^\infty dq
    = \int_0^{pi/2} dz C/(cos^2 z)  with q= C tan(z)+q_min 
    = \int_{-1}^1 dx (pi/2/2) C/(cos^2 z)  with z= pi/2*(x+1)/2 
    
    if n1=0,(or n1+n2=0), set q_2=q_min 
    if n1>0 and n2>0, use input q_2 value  
    
    """
    if (k_F**2-Kcm**2 > 0.0):
        q_min = np.sqrt(k_F**2-Kcm**2)
    else : 
        q_min = 0.0 
    q_1= k_F+Kcm  # input override !  
    if (n_1 > 0 and n_2 > 0): # separate region case 
        xk1, wk1 = leggau(n_1,xmin=q_min,xmax=q_1) #first region
        xk2, wk2 = leggau(n_2,xmin=q_1,xmax=q_2)   #second region
        #---third region 
        xk, wk = leggau(nleg-n_1-n_2,xmin=0.0,xmax=np.pi/2.0)
        Ck = 3.0 # somewhat arbitrary 
        xk3 = Ck*np.tan(xk) + q_2  
        wk3 = Ck*wk/np.cos(xk)**2
        #--combine 
        qq = np.concatenate((xk1,xk2,xk3))
        ww = np.concatenate((wk1,wk2,wk3))
    else:  # no separation of region 
        xk, wk = leggau(nleg,xmin=0.0,xmax=np.pi/2.0)
        Ck = 3.0 # somewhat arbitrary 
        qq = Ck*np.tan(xk) + q_min 
        ww = Ck*wk/np.cos(xk)**2
    return qq, ww 


        
def make_uq_array(xq,wq,k0,Kav0,k_F,p0_array,U_array):
    """
    from xq[1:Nq], wq[1:Nq] and k0
    construct Nq+1 dimensional xq_ext and uq_ext   

    Do we need to be careful for the extrapolation of U? 
    also,  set 
    im_U(p) = 0 and im_dU(p)=0 if p < kF 
    in the energy calculation ?       
    
    Here, xq,wq are assumed to be independent of p0,k0.                      
    """ 
    #--prepare extended quadrature, weights
    xq_ext = np.concatenate( (xq,np.atleast_1d(k0)) ) # Nq+1
    #--constrct u_q
    temp1 = -2./np.pi*wq*xq**2*pauli(Kav0,xq,k_F)
    if scp_choice==0: #standard choice, effective mass approximation
        ratio = effective_mass[0] 
        U0 = effective_mass[1] 
        temp2 = (xq**2+Kav0**2  
             -( (k0**2+Kav0**2)/ratio-2*mN*U0/hbarc**2) )
    elif scp_choice==1: #continuous choice 
        temp2 = ( energy_2N_U(Kav0,xq,k_F,p0_array,U_array,xq) 
             -energy_2N_U(Kav0,k0,k_F,p0_array,U_array,xq) )
    #--- u_q(1:Nq) = temp1/temp2 
    #-----now for the last term 
    temp3 = 2./np.pi*np.sum( wq/temp2 )*k0**2*pauli(Kav0,k0,k_F)
    #-----imaginary part treatment 
    if scp_choice==0: # standard choice, no imaginary part 
        temp4 = temp3  
    elif scp_choice==1: # continuous choice
        #------
        # E'(k0) using p0_array and U_array 
        #
        # dE/dk0 is obtained by using cubic spline 
        #-----
        choice = 1 #note!
        if (choice==0):
            #----similar to Amos code ?               
            Qbar_k0 = pauli(Kav0,k0,k_F)
            k0_sq_plus = Kav0**2+k0**2+2./np.sqrt(3.)*k0*Kav0*Qbar_k0**(1.5) #|k+K|^2
            k0_sq_mnus = Kav0**2+k0**2-2./np.sqrt(3.)*k0*Kav0*Qbar_k0**(1.5) #|k-K|^2        
            k0_plus = np.sqrt(k0_sq_plus)
            k0_mnus = np.sqrt(k0_sq_mnus)
            U_plus, dU_plus = interpolate_U(k0_plus, p0_array, U_array, xq)   
            U_mnus, dU_mnus = interpolate_U(k0_mnus, p0_array, U_array, xq)          
            dEk0 = (2*k0/(mN/hbarc)+dU_plus+dU_mnus)*(mN/hbarc) # careful for units        
        elif (choice==1):
            #---------alternative numerical derivative
            fcs= CubicSpline(xq ,temp2) # interpolation of E(k)-E(k0)
            dEk0= fcs(k0,1) # interp/extrap of derivative dE/dk 
        #--------------------------------------------------------------          
        temp4 = temp3 - 2.*1j*k0**2*pauli(Kav0,k0,k_F)/np.abs(dEk0) 
        #temp4 = 0.0 #test for no principal value trick and no complex pole
    #----construct uq_ext      
    uq_ext = np.concatenate( (temp1/temp2,np.atleast_1d(temp4))) # Nq+1 dimension array 
    return xq_ext, uq_ext 

def solve_gmatrix(k0,Kav0,k_F,
                  xq_ext,uq_ext,ich,chan_dic,Vmat):
    """
    solve g-matrix 
    """
    #----construct potential matrix 
    #and extend xq_ext for 2Nq+2 if icoup channel      
    if (chan_dic['icoup'][ich]==0): # unicoup channel
        # contruct uV matrix 
        uVmat = np.zeros( (Nq+1 ,Nq+1 ) )+1j*0.0 
        #--test shows ui*(mN/hbarc)=-uq_ext and vv= vmat 
        for i in range(Nq+1):
            uVmat[i,:] = Vmat[i,:]*uq_ext[:]
        #---solve g-matrix     
        gmat = np.matmul( np.linalg.inv(np.eye(len(xq_ext)) - uVmat) , Vmat  ) 
    elif (chan_dic['icoup'][ich] > 0): # icoup channel         
        #---extend 2Nq+2 dim matrix  
        xq_ext = np.concatenate( (xq_ext,xq_ext) )
        uq_ext = np.concatenate( (uq_ext,uq_ext) )              
        # contruct uV matrix 
        uVmat = np.zeros( (2*Nq+2 ,2*Nq+2 ) )+1j*0.0 
        for i in range(2*Nq+2):
            uVmat[i,:] = Vmat[i,:]*uq_ext[:]
        #---solve g-matrix     
        gmat = np.matmul( np.linalg.inv(np.eye(len(xq_ext)) - uVmat) , Vmat  ) 
    
    return gmat

def get_U_from_g(chan_dic,gmat_dic,wXk,chan_choice=None):
    # From the g_matrix results computed for p=p0
    # compute U(p=p0) 
    # (also need integration weights wXk)
    #
    #---gmatrix for all channels, all k0 are obtained          
    # construct U(p_0)
    #
    # chan_choice can be used to only consider particular channels 
    factor = 8./np.pi/(mN/hbarc)
    U_p0 = 0.0
    # U_p0_chns=[]
    if chan_choice is not None: 
        list_of_chan_index=chan_choice 
    else:
        list_of_chan_index= range(len(chan_dic['J']))
        
    for ich in list_of_chan_index:
        J = chan_dic['J'][ich]
        L = chan_dic['L'][ich]
        S = chan_dic['S'][ich]
        T = chan_dic['T'][ich]
        sums = 0.0 
        if ( chan_dic['icoup'][ich] >= 0):
            on_shell_idx = Nq   # careful for python convention
        elif ( chan_dic['icoup'][ich] < 0):
            on_shell_idx = 2*Nq +1 # careful for python convention
        
        for jj in range(len(wXk)):
            sums = sums +wXk[jj]*gmat_dic[ich,jj][on_shell_idx,on_shell_idx]
            
        if neutron_matter==0:   # nuclear matter 
            U_p0 = U_p0 + factor*(2*J+1)*(2*T+1)*sums  
        elif neutron_matter==1: # neutron matter 
            if np.mod(L+S,2)==0: # L+S+1 is odd, only T=1,Tz=1 should be summed  
                U_p0 = U_p0 + factor*(2*J+1)*sums  
            else:
                pass                 
        # U_p0_chns.append(factor*(2*J+1)*(2*T+1)*sums) # separate channel contribution   
    return U_p0*hbarc # to MeV unit  


def get_gmat_dic_for_p0(it,ip,k_F,chan_dic,Nq,Nk,p0_list,U_list,Vmat_dic):    
    """
    compute gmat_dic[ich,ik] for given ip 
    
    it : index of iteration 
    ip : index for p0 momentum 
    ik : index of k0 quadrature 

    """
    gmat_dic={} #g-matrix at p_0 value
    p_0 = p0_list[ip]
    xk, wXk, Kav = construct_k0_quadrature(p_0,k_F,nleg=Nk,opt=1)
    for ik in range(len(xk)):
        #-----now construct the g-matrix equation    
        k0 = xk[ik]
        Kav0 = Kav[ik]                
        xq, wq = construct_q_quadrature2(Kav0, k_F,nleg=Nq,
                            n_1=18,n_2=18,q_1=2.0,q_2=20.0) #depends on K_av 
        #----construct Nq+1 arrays 
        xq_ext, uq_ext = make_uq_array(xq,wq,k0,Kav0,k_F,p0_list,U_list)
        
        #----construct/load potential table
        if it==1: # only at the first iteration 
            #Vmat_dic[ip,ik] = make_Vmat_dic(xq_ext,chan_dic,Vopt=pot_choice)
            Vmat_dic[ip,ik] = make_Vmat_dic2(xq_ext,chan_dic,Vopt=pot_choice)
            #print('To do: save_Vmat_dic(p0,k0)')                      
        else:
            pass
            #print('To do: load_Vmat_dic(p0,k0)')                    
        for ich in range(len(chan_dic['J'])):                     
            if (chan_dic['icoup'][ich]==0):
                Vmat = Vmat_dic[ip,ik][ich]
                gmat = solve_gmatrix(k0,Kav0,k_F
                          ,xq_ext,uq_ext,ich,chan_dic,Vmat)
                gmat_dic[ich,ik]=gmat
            if (chan_dic['icoup'][ich] >0):
                Vmat = Vmat_dic[ip,ik][ich]
                gmat = solve_gmatrix(k0,Kav0,k_F
                          ,xq_ext,uq_ext,ich,chan_dic,Vmat)
                gmat_dic[ich,ik]=gmat # gmat is stored in both icoup channels
                gmat_dic[ich+1,ik]=gmat
    
    return gmat_dic ,wXk  ,Vmat_dic, xq_ext

def get_gmat_dic_nuclear_matter(k_F,chan_dic,Nq,Nk,
                                p0_list,U_list,Vmat_dic={},it=1):
    """
    Same as get_gmat_dic_for_p0 
    except using tilde{K}_{av} for nuclear matter energy calculation. 
    
    And xk quadrature in range (0,k_F) independent of p0
    
    assume U(p) is already obtained self consistently 
    
    Thus,no iteration is necessary. 
    
    it = 1 -> Vmat_dic is re-calculated 
         other wise Vmat_dic is taken as input 
    """    
    gmat_dic={} 
    xk, wk = leggau(Nk,xmin=0,xmax=k_F)
    ratio = xk/k_F
    wXk = wk*xk**2*(1.0-1.5*ratio+0.5*ratio**3)   # for energy integration  
    tildeKav= np.sqrt(3.0/5.0*k_F**2*(1.0-ratio)
                      *(1.0+ratio**2/(3.0*(2+ratio))))
    for ik in range(len(xk)):
        k0 = xk[ik]
        Kav0 = tildeKav[ik]
        xq, wq = construct_q_quadrature2(Kav0, k_F,nleg=Nq,
                            n_1=18,n_2=18,q_1=2.0,q_2=20.0)
        #xq, wq = construct_q_quadrature(Nq)
        xq_ext, uq_ext = make_uq_array(xq,wq,k0,Kav0,k_F,p0_list,U_list)
        #----construct/load potential table
        if it==1: # only at the first iteration 
            Vmat_dic[ik] = make_Vmat_dic2(xq_ext,chan_dic,Vopt=pot_choice)
            #print('To do: save_Vmat_dic(p0,k0)')                      
        else:
            pass
            #print('To do: load_Vmat_dic(p0,k0)')                                
        for ich in range(len(chan_dic['J'])):                     
            if (chan_dic['icoup'][ich]==0):
                Vmat = Vmat_dic[ik][ich]
                gmat = solve_gmatrix(k0,Kav0,k_F
                          ,xq_ext,uq_ext,ich,chan_dic,Vmat)
                gmat_dic[ich,ik]=gmat
            if (chan_dic['icoup'][ich] >0):
                Vmat = Vmat_dic[ik][ich]
                gmat = solve_gmatrix(k0,Kav0,k_F
                          ,xq_ext,uq_ext,ich,chan_dic,Vmat)
                gmat_dic[ich,ik]=gmat # gmat is stored in both icoup channels
                gmat_dic[ich+1,ik]=gmat        
    return gmat_dic, wXk,Vmat_dic
    
def Energy_nuclear_matter(k_F,chan_dic,Nq,Nk,p0_list,U_list):    
    """
    Compute energy of Nuclear Matter in BHF 
    
    Note that in fact, p0_list and U_list depends on k_F 
    """
    gmat_dic, wXk, Vmat_dic = get_gmat_dic_nuclear_matter(
                                  k_F,chan_dic,Nq,Nk,p0_list,U_list,it=1)
    
    E_over_A = 0.0  
    
    factor = 4.0/np.pi/(mN/hbarc)
    for ich in range(len(chan_dic['J'])):
        J = chan_dic['J'][ich]
        T = chan_dic['T'][ich]
        if ( chan_dic['icoup'][ich] >= 0):
            on_shell_idx = Nq   # careful for python convention
        elif ( chan_dic['icoup'][ich] < 0):
            on_shell_idx = 2*Nq +1 # careful for python convention
            
        sums = 0.0    
        for jj in range(len(wXk)):
            sums = sums +wXk[jj]*np.real(gmat_dic[ich,jj][on_shell_idx,on_shell_idx])  
        E_over_A = E_over_A + factor*(2*J+1)*(2*T+1)*sums 
    #----add free gas energy  
    E_over_A = 3.0/5.0*k_F**2/(2*mN/hbarc) + E_over_A     
    return E_over_A *(hbarc) # in MeV 

def Energy_nuclear_matter_another(k_F,Np,p0_list,U_list):
    """
    Instead of using G-matrix, use self consistent potential U
    for the calculation of energy. 
    
    E/A = 3/5*k_F^2/m +3/(2 k_F^3)*\int_0^{k_F} dp p^2 U(p) 
    
    """
    interp_U = CubicSpline(p0_list, np.real(U_list)) # interpolate Re(U)(p) 
    xp, wp = leggau(Np,xmin=0,xmax=k_F)
    integral = np.sum( wp*xp**2*interp_U(xp)/hbarc) # integration over p
    E_over_A = 3.0/5.0*k_F**2/(2*mN/hbarc) + 3.0/(2*k_F**3)*integral
    return E_over_A*(hbarc)  
        
#============================================================
if __name__ == '__main__':    
    start_time = time.time()    
    # prepare channels
    chan_dic = make_channel_list(Jmax,save=True)        
    #======for test ==============
    # To simulate Amos code results, one can choose particular channels 
    # which contributes to U(p) 
    chan_choice = None #default 
    #chan_choice = [0,2,6,1,3,7,4,5,8,9] # Amos channel order 
    #chan_choice = [0]
    #prepare quadratures 
    
    #-----k_F values for E/A plot  
    # In case of kF<=0, solve free space t-matrix 
    
    k_F_list = np.arange(0.2,2.4,0.2)
    #k_F_list = np.array([0])
    EA_list = []
    save_gmat_dic={'channels': chan_dic,
                   'k_F_list': k_F_list }
    for k_F in k_F_list:
        print('Computing U(p) function at k_F=',k_F)
        #----p values to compute U(p) 
        #--default choice of p0_list 
        if k_F > 0. :
            if scp_choice==0:
                p0_list = np.linspace(1.e-4,k_F,10)
            elif scp_choice==1:
                p0_list = np.linspace(1.e-4,3*k_F,20)
        else: 
            #--free space calculation 
            p0max = 200. 
            p0_list = np.linspace(1.e-4,p0max,40)
        #--manual choice of p0_list     
        # p0_list = np.linspace(1.e-4,3*k_F,20)
        # p0_list = np.array([0.1e-3,0.2,0.4,0.6,0.8,1.0,1.2,1.4]) 
        # p0_list= np.array([0.1e-3,0.6,1.2,1.8,2.4,3.0,3.6,4.2,4.8,5.4,6.0,6.6,7.2])
        # p0_list = np.linspace(1.e-4,k_F,10)
        #---check p0_list consistency 
        if scp_choice==0:
            if ( p0_list-k_F >0).any():
                raise ValueError('p>k_F for strandard choice')        
        elif scp_choice==1:
            if (np.max(p0_list) < 2*k_F):
                print('Warning: p_max < 2*k_F.'+
                      ' extrapolation U(p) may not be reliable.')
        #--initial guess U(p) -----------------------------------------------
        #interp_U = lambda x: -100./(np.exp(10.*x**2)+1)  #WS form 
        interp_U = lambda x: 0.0  #initial guess         # zero 
        U_list = np.array([ interp_U(i) for i in p0_list])     
        #--iterations for U(p)------------------------------------------
        U_results={}
        U_results[0] = U_list #initialize 
        
        Vmat_dic={} 
        if k_F <=0. :
            iter_it = 1 # free space 
        else :
            iter_it = iter_max +1 
            
        for it in range(1,iter_it):
            print('iteration for U(p) at it=',it)
            U_new = []         
            for ip in range(len(p0_list)) :        
                #print('running iter={}, ip={}'.format(it,ip))
                gmat_dic, wXk, Vmat_dic, xq_ext = get_gmat_dic_for_p0(it,ip,k_F,chan_dic,Nq,Nk,
                                                    p0_list,U_list,Vmat_dic) 
                Up0 = get_U_from_g(chan_dic,gmat_dic,wXk,chan_choice=chan_choice)           
                U_new.append(Up0)
            U_new = np.array(U_new)
            #----convergence check /update U(p)  
            if scp_choice==0:
                effective_mass_new = get_effective_mass(k_F, p0_list, U_new,
                                                        method=effective_mass_def)
                diff=np.array(effective_mass)-np.array(effective_mass_new) 
                if ( np.abs(diff[0])<1.0e-3 and np.abs(diff[1]< 1.0e-3)):
                    print('effective mass converged. ratio={},U0={}'.format(
                        effective_mass[0],effective_mass[1])) 
                else:
                    effective_mass = effective_mass_new 
                    
            elif scp_choice==1:
                diff = U_list - U_new 
                diff_per_point = np.sqrt(np.sum(diff**2))/len(diff)
                if diff_per_point < 1.0e-3:
                    print('auxiliary potential converged')  
            #---update auxiliary potential        
            U_list=np.array(U_new) # unit conversion to MeV 
            #----store----------------------
            U_results[it]= U_list    
        #------plot---------------------------------------------------------------    
        final_time = time.time() 
        print('Execution time(sec)=',final_time-start_time,'for k_F=',k_F)    
        #=====save gmatrix=================================
          
        
        #------compare with amos_code_results 
        # amos_dic = read_output_Amos('nntgmat_Loc_ch_1_to_10.out')
        # # Re U
        # for i in range(len(U_results)):    
        #      plt.plot(p0_list,np.real(U_results[i]),label='My {}'.format(i)  )   
        # for i in range(len(amos_dic)):    
        #      plt.plot(amos_dic[i][:,0],amos_dic[i][:,1],'*',label='Amos {}'.format(i) )         
        # plt.legend()   
        #=====compute nuclear matter energy 
        E_over_A = Energy_nuclear_matter_another(k_F,20,p0_list,U_list)
        print('E/A={} at k_F={}'.format( E_over_A, k_F ) )
        EA_list.append(E_over_A) 
    #---plot E/A graph------------------    
    if len(k_F_list) > 3 :
        plt.figure()     
        plt.plot(k_F_list, EA_list)
        plt.title('E/A from BHF calculation')
        plt.xlabel(r'$k_F(fm^{-1})$')
        plt.ylabel('E/A (MeV)')
        plt.savefig('Energy_Neutron.png')
        EA_list = np.array(EA_list)
        data = np.zeros( (len(k_F_list),2)   ) 
        data[:,0]=k_F_list 
        data[:,1]=EA_list 
        np.savetxt('EA_BonnB_neutron.txt',data,header='k_F  E/A')
    
        rho = k_F_list**3/3/np.pi**2 #neutron matter 
        plt.plot(rho, EA_list)
        plt.xlim(0,0.26)
        plt.ylim(0,30)
    
    
    