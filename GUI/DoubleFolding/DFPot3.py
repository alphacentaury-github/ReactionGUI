# -*- coding: utf-8 -*-
"""
Double Folding Potential code
for spherical density

Created on Tue Jun  2 17:20:48 2020

@author: Y.-H. Song
"""

import numpy as np
import scipy
import scipy.misc
from scipy.special import spherical_jn
import matplotlib.pyplot as plt

import scipy.integrate
import scipy.interpolate
import os

#---current path where this file resides
try:
    here = os.path.dirname(os.path.realpath(__file__))
except:
    here = '.'

# M3Y potential Zero range or finite range exchange potential
M3Y_core_para= {'M3Y_Reid': [ [7999.0, -2134.0, 0.0],[ 4.0 , 2.5, 0.707],
                            [4631.4,-1787.1,-7.847], [4.0, 2.5, 0.707],
                            [0.002,-276.0]],
                'M3Y_Paris': [[11061.625,-2537.5, 0.0 ], [4.0,2.5,0.707],
                             [-1524.25, -518.75, -7.847], [4.0,2.5,0.707],
                             [0.003, -592.0]]}

# CDD,alphaDD,betaDD,gammaDD of density dependent M3Y pot
M3Y_DD_para = { 'M3Y_Reid' : (1.0,0.0,0.0,0.0) ,
                'M3Y_Paris': (1.0,0.0,0.0,0.0) ,
                'DDM3Y_Reid': (0.2845, 3.6391, 2.9605, 0.0) ,
                'DDM3Y_Paris': (0.2963,3.7231,3.7384,0.0) ,
                'BDM3Y_Reid' : (1.2253, 0.0, 0.0, 1.5124),
                'BDM3Y_Paris': (1.2521, 0.0, 0.0, 1.7452),
                'CDM3Y1_Paris': (0.3429, 3.0232, 3.5512, 0.5),
                'CDM3Y2_Paris': (0.3346, 3.0357, 3.0685, 1.0),
                'CDM3Y3_Paris': (0.2985, 3.4528, 2.6388, 1.5  ),
                'CDM3Y4_Paris': (0.3052, 3.2998, 2.3180, 2.0  ),
                'CDM3Y5_Paris': (0.2728, 3.7367, 1.8294, 3.0  ),
                'CDM3Y6_Paris': (0.2658, 3.8033, 1.4099, 4.0) }

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

# Shape forms

def shape_Gaussian(r,alpha=1.0,pwr=1.0,C=1.0):
    return  C*r**pwr*np.exp(-(r/alpha)**2)

def shape_Yukawa(r,beta=1.0,V=1.0):
    return V*np.exp(-beta*r)/(beta*r)

def shape_3p_Fermi(r,RR=1.0,AR=1.0,W=1.0,C=1.0):
    return C*(1.0+W*r**2)/(1.0+np.exp((r-RR)/AR))

def shape_WS(r,R0=1.0,a0=1.0,V0=1.0,PWR=1.0,IDER=0):
    """
    return Woods-Saxon shape or derivative

    IDER==0 case:
        f(r)=V0/(1+exp((r-R0)/a0) )^PWR

    IDER!=0 case:
        df(r)/dr

    """
    EXPR = np.exp((r-R0)/a0)
    if IDER ==0:
        return V0/(1.0+EXPR)**PWR
    else:
        return -PWR*V0*EXPR/(a0*(1.0+ EXPR)**(PWR+1))

def shape_delta(r,V0=1.0,small=1.0e-15):
    """
    delta function V0*delta(r)

    defined only for very small r
    """
    if r < small:
        return V0
    else:
        return 0.0

def moment_integral(r_array,rho_array,k=0):
    """
    moment integral of rho(x) using samples

    (4*pi) \int rho*r^(2+k) dr

    Parameters
    ----------
    r_array :  array of sample points
             length should be 2**(integer)+1
    rho_array :  array of sample values of rho(x)
        length should be the same as r_array
    k : integer, optional
        order of moment. The default is 0.

    Returns
    -------
    volume :  float
        result of k-moment intgeral
        (4*pi) \int rho*r^(2+k) dr
    rms :  float
        ratio of (k+2)-moment and k-moment
        ( (4*pi) \int rho*r^(4+k) dr) /(volume)

    Example
    -------

    >>> romb(xxx)

    """
    step = r_array[1]-r_array[0]
    volume  = scipy.integrate.romb(rho_array*r_array**(2+k))*step*4*np.pi
    rms     = scipy.integrate.romb(rho_array*r_array**(4+k))*(
                step*4*np.pi/volume)
    return (volume ,rms )
#-------------------------------------------------------------------
def read_density(Z,A,directory_path=here+'/density-hfb14'):
    """
    Parameters
    ----------
    Z : integer
        charge number in range Z=8-110
    A : integer
        mass number
    directory_path : string, optional
        path of the directory which store files.
        The default is 'density-hfb14'.

    Returns
    -------
    (r, rhon, rhop)

    r: array of radial points
    rhon : neutron density
    rhop : proton density
    """
    filename =directory_path+'/z{:03d}.dat'.format(Z)
    ff = open(filename,'r')
    ll = ff.readlines()
    ff.close()
    found= False
    for (i,line) in enumerate(ll):
        if  'Z= {:3d} A= {:3d}'.format(Z,A) in line:
            found = True
            start_i = i
    if not found:
        return (-1,-1,-1) # error
    # found
    ww = ll[start_i].replace('=',' ')
    ww = ww.split()
    #beta2 = float(ww[6]) # beta2 is not necessary for the moment
    nrho  = int(ww[8])
    r = []
    rhon = []
    rhop = []
    for i in range(start_i+2,start_i+nrho+2):
         ww= ll[i].split()
         r.append( float( ww[0] ))
         rhon.append( float( ww[1] ))
         rhop.append( float( ww[2] ))
    r = np.array(r)
    rhon= np.array(rhon)
    rhop= np.array(rhop)
    return (r,rhon,rhop)

def interpolating_function(points,values,fill_value='extrapolate'):
    """
    return interpolating/extrapolating function
    from given points and values f(r)

    """
    from scipy import interpolate
    f = interpolate.interp1d(points,values
                            ,kind='cubic'
                            ,bounds_error=False
                            ,fill_value = fill_value )
    return f

#----------------------------------------------------------------------

def SBT(k,r_i,w_i,f_i,L=0):
    """
    Compute spherical bessel transformation
    using quadrature r_i and weights w_i
    and function values at r_i

    g(k) = \int_rmin^rmax dr j_0( k r )*f(r)
         => sum_i w_i* j_0(k*r_i)*f(r_i)

    Note that r^2 factor should be included in the function f(r)
    """
    jL = spherical_jn(L,k*r_i)
    return np.sum(jL*w_i*f_i)


#-------------------------------------------------------------------
def M3YZR_pot(r,type='M3Y_Reid'):
    """
    return M3Y potential value in MeV units
    ( v00(r),v01(r) ) for given r

    Zero-range version for exchange part.
    Thus, the delta potential need to be treated separately

    """

    if type in M3Y_core_para.keys():
        VD_00 = M3Y_core_para[type][0]
        MD_00 = M3Y_core_para[type][1]
        VE_00_delta = M3Y_core_para[type][4][1]
    else:
        raise ValueError('ERROR: not available option in M3Y.\n'
           +' It should be among: {}'.format(M3Y_core_para.keys()))

    if r < 1e-4:
        # in fact, this is not correct
        # the potential diverges ar zero range
        return (VE_00_delta,None)
    else:
        V00 = 0.0
        #V01 = 0.0
        for i in range(3):
            V00 = V00+VD_00[i]*np.exp(-MD_00[i]*r)/(MD_00[i]*r)
            #V01 = V01+VD_01[i]*np.exp(-MD_01[i]*r)/(MD_01[i]*r)
        return (V00,None)

def M3YZR_pot_k(k,E_A=0.0,type='M3Y_Reid'):
    """
    return F.T. of M3Y potential
    including exchange term in zero-range approximation

    k in fm^{-1} unit

    Energy dependence appears in exchange delta contribution
    as (1-0.002*E/A)

    output is in MeV.fm^3 unit
    """
    if type in M3Y_core_para.keys():
        VD_00 = M3Y_core_para[type][0]
        MD_00 = M3Y_core_para[type][1]
        VE_00_delta = M3Y_core_para[type][4][1]
        en_fac_00 = M3Y_core_para[type][4][0]
        VE_00_delta = VE_00_delta*(1-en_fac_00*E_A)
    else:
        raise ValueError('ERROR: not available option in M3Y.\n'
           +' It should be among: {}'.format(M3Y_core_para.keys()))

    V00D = 0.0
    #V01 = 0.0
    for i in range(3):
        term = 4*np.pi*VD_00[i]/(MD_00[i]*(MD_00[i]**2+k**2))
        V00D = V00D + term
        #term = 4*np.pi*VD_01[i]/(MD_00[i]*(MD_01[i]**2+k**2))
        #V01 = V01 + term
    # add exchange correction term
    V00E = VE_00_delta
    V00  = V00D + V00E
    return (V00D,V00E,V00,None,None,None)

#============================================================================
def DF_M3YZR(f_rhon_p,f_rhop_p,f_rhon_t,f_rhop_t,
           E_A=0.0, type='M3Y_Reid',
           R_points= np.arange(0.1,30.0,0.1),
           k_range = (0.0,3.0), num_quadrature=70,
           r_range = (0.0,15.0) ):
    """
    Compute Double Folding potential using
    input density functions rho(r)
    and bare M3Y interaction potential.
    ('M3Y-Reid')

    Direct folding and zero-range exchange only

    Parameters
    ----------
    f_rhon_p, f_rhop_p : function
        neurton and proton density function of projectile

    f_rhon_t, f_rhop_t : function
        neurton and proton density function of target
    E_A : float
        incident energy in MeV per nucleon
    R_points : array, optional
        points of interest for optical potential.
        The default is np.arange(0.1,30.0,0.1).

    Returns
    -------
    R and Coulomb U_C and isoscalar folding V_D(R;E)

    """
    # get integral quadrature weight
    (xr, wr) = gauleg(r_range[0],r_range[1],num_quadrature)
    (xk, wk) = gauleg(k_range[0],k_range[1],num_quadrature)
    rhon_1 = f_rhon_p(xr) # density values at quasrature
    rhop_1 = f_rhop_p(xr)
    rhon_2 = f_rhon_t(xr)
    rhop_2 = f_rhop_t(xr)

    # compute spherical Bessel transforms
    # M3Y potential V_00(k)*k^2 and V_01(k)*k^2
    alpha_em = 197.3269/137.03599 # in MeV.fm unit
    Vk_out = M3YZR_pot_k(xk,E_A=E_A,type=type) # Vk_out=[VD,VE,VD+VE,...]
    tilde_v00 = Vk_out[2] # VD+VE
    #(tilde_v00 , tilde_v01 ) =  M3YZR_pot_k(xk,E_A=E_A,type=type)
    tilde_vksq_00 = tilde_v00*xk**2  # v_00(k)*k^2

    # SBT of target projectile density
    tilde_rhon_1 = np.array([4.0*np.pi*SBT(k,xr,wr,rhon_1*xr**2) for k in xk])
    tilde_rhop_1 = np.array([4.0*np.pi*SBT(k,xr,wr,rhop_1*xr**2) for k in xk])
    tilde_rhon_2 = np.array([4.0*np.pi*SBT(k,xr,wr,rhon_2*xr**2) for k in xk])
    tilde_rhop_2 = np.array([4.0*np.pi*SBT(k,xr,wr,rhop_2*xr**2) for k in xk])

    # combined one # k^2 factors are already included.
    tilde_UC = tilde_rhop_1*tilde_rhop_2 # Coulomb
    tilde_U_isoscalar = ( (tilde_rhon_1+tilde_rhop_1)
                         *tilde_vksq_00*(tilde_rhon_2+tilde_rhop_2) )
    # Get isoscalar U and isovector U
    UC = []
    U_isos = []
    #U_isov = []
    for Ri in R_points:
        UC.append( alpha_em*2.0/np.pi*SBT(Ri,xk,wk,tilde_UC) )
        U_isos.append( 1./(2*np.pi**2)*SBT(Ri,xk,wk,tilde_U_isoscalar) )
        #U_isov.append(1./(2*np.pi**2)*SBT(Ri,xk,wk,tilde_U_isovector) )
    return (R_points, np.array(UC), np.array(U_isos) )
#--------------------------------------------------------------
#M3Y finite range version
#-------------------------------------------------------------
def M3YFR_pot(r,E_A=0.0,type='M3Y_Reid'):
    """
    Finite range version M3Y potential  in the exchange part

    at the moment iso-scalar part only v_00

    r can be array.
    """
    # M3YFR parameters
    if type in M3Y_DD_para.keys():
        (C,alpha,beta,gamma) = M3Y_DD_para[type]
    else :
         raise ValueError('ERROR: not available option in M3Y.\n'
                          +' It should be among: {}'.format(M3Y_DD_para.keys()))

    ww = type.split('_')

    if (ww[1]=='Reid'):
        #iso-scalar
        VD_00 = [7999.0, -2134.0, 0.0] # MeV
        MD_00 = [ 4.0 , 2.5, 0.707] # fm^-1
        VE_00 = [4631.4,-1787.1,-7.847]
        ME_00 = [4.0, 2.5, 0.707]
        g_E = 1.0-0.002*E_A # 0.002 or 0.005 ?
        #iso-vector
        #VD_01 = [-4886.0, 1176.0, 0.0] # MeV
        #MD_01 = [4.0,2.5,1.0] # fm^{-1}
        #    en_fac_01 = 1.0-0.005*E_A # this is not clear
        # VE_01_delta = 217.0*en_fac_01
    elif (ww[1]=='Paris'):
        #iso-scalar
        VD_00 = [11061.625,-2537.5, 0.0 ]
        MD_00 = [4.0,2.5,0.707]
        VE_00 = [-1524.25, -518.75, -7.847]
        ME_00 = [4.0,2.5,0.707]
        g_E = 1.0-0.003*E_A
        #iso-vector is not certain
    else :
        print('ERROR: option is not available!')
        return (None,None)
    VD_00_r = 0.0
    VE_00_r = 0.0
    # V_01 = 0.0
    for i in range(3) :
        VD_00_r = VD_00_r+VD_00[i]*np.exp(-MD_00[i]*r)/(MD_00[i]*r)
        VE_00_r = VE_00_r+VE_00[i]*np.exp(-ME_00[i]*r)/(ME_00[i]*r)
    return (g_E*VD_00_r,g_E*VE_00_r, C,alpha,beta,gamma )

def M3YFR_pot_k(k,E_A=0.0,type='M3Y_Reid'):
    """
    Momentum space expression of M3Y-finite range potential

    However, note that the exchange term does not use this form
    because of additional phase factor

    k can be array in fm^-1 units
    """
    # M3YFR parameters
    if type in M3Y_DD_para.keys():
        (C,alpha,beta,gamma) = M3Y_DD_para[type]
    else :
         raise ValueError('ERROR: not available option in M3Y.\n'
                          +' It should be among: {}'.format(M3Y_DD_para.keys()))
    ww = type.split('_')
    if (ww[1]=='Reid'):
        #iso-scalar
        VD_00 = [7999.0, -2134.0, 0.0] # MeV
        MD_00 = [ 4.0 , 2.5, 0.707] # fm^-1
        VE_00 = [4631.4,-1787.1,-7.847]
        ME_00 = [4.0, 2.5, 0.707]
        g_E = 1.0-0.002*E_A # 0.002 or 0.005 ?
        #iso-vector
        #VD_01 = [-4886.0, 1176.0, 0.0] # MeV
        #MD_01 = [4.0,2.5,1.0] # fm^{-1}
        #    en_fac_01 = 1.0-0.005*E_A # this is not clear
        # VE_01_delta = 217.0*en_fac_01
    elif (ww[1]=='Paris'):
        #iso-scalar
        VD_00 = [11061.625,-2537.5, 0.0 ]
        MD_00 = [4.0,2.5,0.707]
        VE_00 = [-1524.25, -518.75, -7.847]
        ME_00 = [4.0,2.5,0.707]
        g_E = 1.0-0.003*E_A
        #iso-vector is not certain
    else :
        print('ERROR: option is not available!')
        return (None,None)
    VD_00_k = 0.0
    VE_00_k = 0.0
    for i in range(3):
        term = 4*np.pi*VD_00[i]/(MD_00[i]*(MD_00[i]**2+k**2))
        VD_00_k = VD_00_k + term
        term = 4*np.pi*VE_00[i]/(ME_00[i]*(ME_00[i]**2+k**2))
        VE_00_k = VE_00_k + term
    return (VD_00_k*g_E,VE_00_k*g_E, C,alpha,beta,gamma)


#============================================================================
def DF_M3YFR(f_rhon_p,f_rhop_p,f_rhon_t,f_rhop_t,
           E_A=0.0, A_p=1, A_t=1,
           type='M3Y_Reid',
           R_points= np.arange(0.1,30.0,0.1),
           k_range = (0.0,3.0), num_quadrature=70,
           r_range = (0.0,15.0),  iter_max = 20  ):
    """
    Compute Double Folding potential using
    input density functions rho(r)
    and bare M3Y interaction potential. (v_{D,E}(s) of M3Y)
    (density dependence is treated seperately
     v_{D,E}(E,rho,s) = F(rho) * v_{D,E}(s) )
    Direct folding and Finite-range exchange

    Reference : D.T. Khoa, G.R. Satchler, Nuclear Physics A 668(2000),3-41

    Parameters
    ----------
    f_rhon_p, f_rhop_p : function
        neurton and proton density function of projectile

    f_rhon_t, f_rhop_t : function
        neurton and proton density function of target
    E_A : float
        incident energy in MeV per nucleon
    R_points : array, optional
        points of interest for optical potential.
        The default is np.arange(0.1,30.0,0.1).

    """
    # get integral quadrature weight
    (xr, wr) = gauleg(r_range[0],r_range[1],num_quadrature)
    (xk, wk) = gauleg(k_range[0],k_range[1],num_quadrature)
    s_range = (0.0,30.0) # what should be the range ...
    (xs, ws) = gauleg(s_range[0],s_range[1],num_quadrature)

    # M3Y potential in momentum space
    alpha_em = 197.3269/137.03599 # Coulomb in MeV.fm unit
    amu=931.4940954      # fm^-1 unit
    hbarc = 197.326968  # MeV.fm

    (v_D_k,v_E_k,CDD,alphaDD,betaDD,gammaDD) = M3YFR_pot_k(xk,E_A,type=type)
    (v_D_s,v_E_s,CDD,alphaDD,betaDD,gammaDD) = M3YFR_pot(xs,E_A,type=type)
    #----First obtain Direct potential term
    # eq. (A.1)
    fourpi_rsq = 4.0*np.pi*xr**2
    rhon_1 = f_rhon_p(xr) # density values at quadrature
    rhop_1 = f_rhop_p(xr)
    rhon_2 = f_rhon_t(xr)
    rhop_2 = f_rhop_t(xr)

    # # derivative of densities ....
    # drhon_1 = scipy.misc.derivative(f_rhon_p,xr,dx=0.01,n=1)
    # drhop_1 = scipy.misc.derivative(f_rhop_p,xr,dx=0.01,n=1)
    # drhon_2 = scipy.misc.derivative(f_rhon_t,xr,dx=0.01,n=1)
    # drhop_2 = scipy.misc.derivative(f_rhop_t,xr,dx=0.01,n=1)
    # ddrhon_1 = scipy.misc.derivative(f_rhon_p,xr,dx=0.01,n=2)
    # ddrhop_1 = scipy.misc.derivative(f_rhop_p,xr,dx=0.01,n=2)
    # ddrhon_2 = scipy.misc.derivative(f_rhon_t,xr,dx=0.01,n=2)
    # ddrhop_2 = scipy.misc.derivative(f_rhop_t,xr,dx=0.01,n=2)
    # # laplacian
    # lap_rhon_1 = ddrhon_1+2.0*drhon_1/xr
    # lap_rhop_1 = ddrhop_1+2.0*drhop_1/xr
    # lap_rhon_2 = ddrhon_2+2.0*drhon_2/xr
    # lap_rhop_2 = ddrhop_2+2.0*drhop_2/xr

    # iso-scalar densities
    rho_1 = rhon_1 + rhop_1 # matter density rho = rhop+rhon
    rho_2 = rhon_2 + rhop_2
    brho_1 = rho_1*np.exp(-betaDD*rho_1) # bar_rho(r)
    brho_2 = rho_2*np.exp(-betaDD*rho_2)
    trho_1 = rho_1**2 # for tilde{rho}(r)
    trho_2 = rho_2**2

    # mass and charge could be obtained from density integral.
    # mass number
    A_1 = np.sum(wr*rho_1*np.pi*4.0*xr**2)
    A_2 = np.sum(wr*rho_2*np.pi*4.0*xr**2)
    # charge number
    Z_1 = np.sum(wr*rhop_1*np.pi*4.0*xr**2)
    Z_2 = np.sum(wr*rhop_2*np.pi*4.0*xr**2)

    # #-----Fermi momentum
    # # At the moment no surface corrections are included...
    #C_S =1.0/36.0  # 1/36 or 1/4 ?
    #kF_1  = np.sqrt(  (3./2.*np.pi**2*rho_1)**(2./3.)
    #                  +5.0*C_S/3.0*(drhop_1+drhon_1)**2/rho_1**2
    #                  +5.0/36.0*(lap_rhop_1+lap_rhon_1)/rho_1 )
    kF_1  = np.sqrt(  (3./2.*np.pi**2*rho_1)**(2./3.) )

    #kF_2  = np.sqrt(  (3./2.*np.pi**2*rho_2)**(2./3.)
    #                  +5.0*C_S/3.0*(drhop_2+drhon_2)**2/rho_2**2
    #                  +5.0/36.0*(lap_rhop_2+lap_rhon_2)/rho_2 )
    kF_2  = np.sqrt(  (3./2.*np.pi**2*rho_2)**(2./3.) )

    rho_1_k = np.zeros(len(xk))
    rho_2_k = np.zeros(len(xk))
    rhop_1_k = np.zeros(len(xk))
    rhop_2_k = np.zeros(len(xk))
    brho_1_k = np.zeros(len(xk))
    brho_2_k = np.zeros(len(xk))
    trho_1_k = np.zeros(len(xk))
    trho_2_k = np.zeros(len(xk))
    for ki in range(len(xk)):
        rho_1_k[ki] = SBT(xk[ki],xr,wr,rho_1*fourpi_rsq,L=0) # rho(k) of (A.2)
        rho_2_k[ki] = SBT(xk[ki],xr,wr,rho_2*fourpi_rsq,L=0)
        rhop_1_k[ki] = SBT(xk[ki],xr,wr,rhop_1*fourpi_rsq,L=0) #for Coulomb
        rhop_2_k[ki] = SBT(xk[ki],xr,wr,rhop_2*fourpi_rsq,L=0)
        brho_1_k[ki] = SBT(xk[ki],xr,wr,brho_1*fourpi_rsq,L=0) #bar{rho}(k) of (A.2)
        brho_2_k[ki] = SBT(xk[ki],xr,wr,brho_2*fourpi_rsq,L=0)
        trho_1_k[ki] = SBT(xk[ki],xr,wr,trho_1*fourpi_rsq,L=0) #tilde{rho}(k) of (A.2)
        trho_2_k[ki] = SBT(xk[ki],xr,wr,trho_2*fourpi_rsq,L=0)

    # U_ND is eq.(A.1)
    U_CD = np.zeros(len(R_points))
    U_ND = np.zeros(len(R_points))
    for Ri in range(len(R_points)):
        term1 = SBT( R_points[Ri],xk,wk, rho_1_k*rho_2_k*v_D_k*xk**2,L=0)
        term2 = alphaDD*SBT( R_points[Ri],xk,wk, brho_1_k*brho_2_k*v_D_k*xk**2,L=0)
        term3 = -gammaDD*SBT( R_points[Ri],xk,wk,
                      (trho_1_k*rho_2_k +rho_1_k*trho_2_k)*v_D_k*xk**2,L=0)
        U_ND[Ri]  = CDD/(2*np.pi**2)*( term1+term2+term3) # g(E) is inlcuded in v_D
        U_CD[Ri]  = alpha_em*2.0/np.pi*SBT(R_points[Ri],xk,wk
                                         ,rhop_1_k*rhop_2_k )

    # #----compute exchange part of DF potential
    # #    For Coulomb interaction, one only use direct term. U_CD
    # #    No Coulomb exchange potential is considered for the moment.
    U_NE = np.zeros(len(R_points))
    # # needs f(k,s), bf(k,s), tf(k,s) , kF(r)
    f_1_ks = np.zeros( (len(xk),len(xs)) ) # for f_1(r,s)
    f_2_ks = np.zeros( (len(xk),len(xs)) ) # for f_2(r,s)
    bf_1_ks = np.zeros( (len(xk),len(xs)) ) # for bar{f}_1(r,s)
    bf_2_ks = np.zeros( (len(xk),len(xs)) ) # for bar{f}_2(r,s)
    tf_1_ks = np.zeros( (len(xk),len(xs)) ) # for tilde{f}_1(r,s)
    tf_2_ks = np.zeros( (len(xk),len(xs)) ) # for tilde{f}_2(r,s)
    #fp_1_ks = np.zeros( (len(xk),len(xs)) )
    #fp_2_ks = np.zeros( (len(xk),len(xs)) )

    for si in range(len(xs)):
        xx1 = kF_1*xs[si]
        hatj_1 = 3.0*(np.sin(xx1)-xx1*np.cos(xx1))/xx1**3  # j_1(kF r)
        xx2 = kF_2*xs[si]
        hatj_2 = 3.0*(np.sin(xx2)-xx2*np.cos(xx2))/xx2**3  # j_1(kF r)
        #xx1 = kFp_1*xs[si]
        #hatjp_1 = 3.0*(np.sin(xx1)-xx1*np.cos(xx1))/xx1**3
        #xx2 = kFp_2*xs[si]
        #hatjp_2 = 3.0*(np.sin(xx2)-xx2*np.cos(xx2))/xx2**3
        for ki in range(len(xk)):
            f_1_ks[ki,si] = SBT(xk[ki],xr,wr,rho_1*fourpi_rsq*hatj_1,L=0 ) #f(k,s) in (A.9)
            f_2_ks[ki,si] = SBT(xk[ki],xr,wr,rho_2*fourpi_rsq*hatj_2,L=0 )
            #fp_1_ks[ki,si] = SBT(xk[ki],xr,wr,rhop_1*fourpi_rsq*hatjp_1,L=0 )
            #fp_2_ks[ki,si] = SBT(xk[ki],xr,wr,rhop_2*fourpi_rsq*hatjp_2,L=0 )
            bf_1_ks[ki,si] = SBT(xk[ki],xr,wr,brho_1*fourpi_rsq*hatj_1,L=0 ) # bar{f}(k,s) in (A.9)
            bf_2_ks[ki,si] = SBT(xk[ki],xr,wr,brho_2*fourpi_rsq*hatj_2,L=0 )
            tf_1_ks[ki,si] = SBT(xk[ki],xr,wr,trho_1*fourpi_rsq*hatj_1,L=0 ) #tilde{f}(k,s) in (A.9)
            tf_2_ks[ki,si] = SBT(xk[ki],xr,wr,trho_2*fourpi_rsq*hatj_2,L=0 )
    #construct G0(R,s) of eq.(A.8) of REF
    G0_Rs = np.zeros( (len(R_points),len(xs)))
    for si in range(len(xs)):
        term1 = f_1_ks[:,si]*f_2_ks[:,si]*xk**2
        term2 = alphaDD*bf_1_ks[:,si]*bf_2_ks[:,si]*xk**2
        term3 = -gammaDD*(f_1_ks[:,si]*tf_2_ks[:,si] + tf_1_ks[:,si]*f_2_ks[:,si])*xk**2
        for Ri in range(len(R_points)):
            G0_Rs[Ri,si] = SBT(R_points[Ri],xk,wk,term1+term2+term3 ,L=0)/(2.0*np.pi**2)

    M = A_1*A_2/(A_1+A_2)
    if (E_A >= 0.0001 ):
        E_cm = M*E_A
    else:    #default
        E_cm = Z_1*Z_2/(A_1**(1./3.)+A_2**(1./3.))*alpha_em

    # iteration for search K(R)
    # what is the convergence condition?
    for Ri in range(len(R_points)-1,-1,-1 ):
        U_NE[Ri]= 0.0
        UNDFP = U_ND[Ri]+U_NE[Ri]+U_CD[Ri]
        converged = False
        for iter in range(iter_max):
            # is it okay to use absolute value???
            KR = np.sqrt(2*amu*M*np.abs(E_cm-UNDFP))/hbarc #fm^-1

            j0_KR = spherical_jn(0,KR*xs/M)
            temp = 4.0*np.pi*CDD*G0_Rs[Ri ,:]*j0_KR*v_E_s*xs**2
            U_NE[Ri] = np.sum( ws*temp)  # eq. (A.7) nuclear folding
            UNDFP_old = UNDFP
            UNDFP = U_ND[Ri]+U_NE[Ri]+U_CD[Ri] # total direct+exchange+coulomb
            if np.abs(UNDFP-UNDFP_old)/(UNDFP+UNDFP_old+1.e-5) < 1.e-5:
                converged = True
                break
        if not converged :
            raise ValueError('not converged after max iteration at Ri={}'.format(Ri) )

    return (R_points, U_ND,U_NE, U_CD)

def normalize_density(f_rho,norm,r_range=(0.,20.),num_quad=70):
    """
    compute normalization factor 'C'
    so that the

    norm = 4pi\int dr r^2 f_rho(r)*C

    f_rho : density function of r
    norm  : normalization factor
    """
    (xr,wr) = gauleg(r_range[0],r_range[1],num_quad)
    rho = f_rho(xr)
    vol = np.sum(wr*rho*xr**2)*4.0*np.pi
    norm_factor = norm/vol
    return norm_factor

def DF_M3Y_all(f_rhon_p,f_rhop_p,f_rhon_t,f_rhop_t,
               type = 'DDM3Y_Paris_ZR' ,
               E_A= 0.0 ,
               R_points=np.arange(0.1,18.0,0.1),
               k_range = (0.0,5.0), num_quadrature= 96,
               r_range =(0.0,15.0), iter_max=30):
    """
    Combinding Zero range and finite range version of M3Y.

    if the type name ends with ZR -> use ZR version
    if the type name ends with FR or stats with B,C,D -> use FR version
    """
    list_of_names = ['M3Y_Paris_ZR', 'M3Y_Reid_ZR',
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
                'CDM3Y6_Paris']
    if not (type in list_of_names):
        raise ValueError('Error: unavailable M3Y potential . ')

    if type[-2:]=='ZR':
        name = type[:-3]
        (R, UC, US) = DF_M3YZR(f_rhon_p,f_rhop_p,f_rhon_t,f_rhop_t,
                            E_A= E_A, type= name,
                            R_points= R_points,
                            k_range = k_range , num_quadrature=num_quadrature,
                            r_range = r_range )
    elif type[-2:]=='FR':
        name = type[:-3]
        (R, UD, UE, UC) = DF_M3YFR(f_rhon_p,f_rhop_p,f_rhon_t,f_rhop_t,
                            type = name ,
                            E_A= E_A ,
                            R_points= R_points,
                            k_range = k_range, num_quadrature= num_quadrature,
                            r_range = r_range, iter_max= iter_max)
        US = UD+UE
    else :
        name = type
        (R, UD, UE, UC) = DF_M3YFR(f_rhon_p,f_rhop_p,f_rhon_t,f_rhop_t,
                            type = name ,
                            E_A= E_A ,
                            R_points= R_points,
                            k_range = k_range, num_quadrature= num_quadrature,
                            r_range = r_range, iter_max= iter_max)
        US = UD+UE
    return (R,UC,US)

#------------------test M3Y-------------------------------------------------
if __name__ == '__main__':

    #-----density from HFB case-----------------------------------------
    zp=10;ap=20;zt=20;at=40;
    try:
        r_p, rhon_p, rhop_p = read_density(zp,ap)
        r_t, rhon_t, rhop_t = read_density(zt,at)
            # normalized densities
        f_rhon_p = interpolating_function(r_p, rhon_p)
        f_rhop_p = interpolating_function(r_p, rhop_p)
        f_rhon_t = interpolating_function(r_t, rhon_t)
        f_rhop_t = interpolating_function(r_t, rhop_t)
    except:
        print('Eroror: projectile/target density is not available ')

    #-----density of alpha +40Ca  in simple form ------------------------
    # f_rhon_p = lambda r : shape_Gaussian(r,alpha=1.19,pwr=0.0,C=1.0)*0.426279/2.0
    # f_rhop_p = lambda r : shape_Gaussian(r,alpha=1.19,pwr=0.0,C=1.0)*0.426279/2.0
    # f_rhon_t = lambda r : shape_WS(r,R0=3.572,a0=0.55,V0=0.084895)
    # f_rhop_t = lambda r : shape_WS(r,R0=3.588,a0=0.55,V0=0.083906)

    # #---case of DFMSPH 40Ca+124Sn of WS type--------------------------
    # zp=20;ap=40;zt=50;at=124  ;
    # f_rhon_p = lambda r : shape_WS(r,R0=3.776,a0=0.5860)
    # f_rhop_p = lambda r : shape_WS(r,R0=3.776,a0=0.5860)
    # f_rhon_t = lambda r : shape_WS(r,R0=5.490,a0=0.4920)
    # f_rhop_t = lambda r : shape_WS(r,R0=5.490,a0=0.4920)

    # factor_p  = normalize_density(f_rhon_p, ap)
    # factor_t  = normalize_density(f_rhon_t, at)

    # f_rhon_p = lambda r : shape_WS(r,R0=3.776,a0=0.5860)*(ap-zp)/ap*factor_p
    # f_rhop_p = lambda r : shape_WS(r,R0=3.776,a0=0.5860)*zp/ap*factor_p
    # f_rhon_t = lambda r : shape_WS(r,R0=5.490,a0=0.4920)*(at-zt)/at*factor_t
    # f_rhop_t = lambda r : shape_WS(r,R0=5.490,a0=0.4920)*zt/at*factor_t
    #--------------------------------------------------------------------
    # -- test of Zero range exchange case-----------------------------------
    # (R, UC,US) = DF_M3YZR(f_rhon_p,f_rhop_p,f_rhon_t,f_rhop_t,
    #                        E_A=0.0, type='M3Y_Paris',
    #                        R_points= np.arange(0.1,30.0,0.1),
    #                        k_range = (0.0,3.0), num_quadrature=70,
    #                        r_range = (0.0,15.0) )

    # # compare with DFMSPH
    # dat = np.loadtxt('DFMSPH/OUT1_DFMSPH.dat')
    # # Coulomb potential
    # plt.plot(dat[:,0],dat[:,4],R,UC,'.');plt.xlim([8,13])
    # # Direct nuclear+ZR exchange
    # plt.plot(dat[:,0],dat[:,7],R,US,'.');plt.xlim([8,13])

    # plt.plot(dat[:,0],dat[:,4],R,UC1,'.');plt.xlim([8,13])
    # plt.plot(dat[:,0],dat[:,5],R,UD1,'.');plt.xlim([8,13]) # direct term only


    #--test of Finite Range exchange case-------------------------------------
    (R, U_ND, UE, UC1) = DF_M3YFR(f_rhon_p,f_rhop_p,f_rhon_t,f_rhop_t,
                            type = 'M3Y_Paris' ,
                            E_A= 0.0 ,
                            R_points=np.arange(0.1,18.0,0.1),
                            k_range = (0.0,5.0), num_quadrature= 96,
                            r_range =(0.0,15.0), iter_max=30)
    (R2, U_ND2, UE2, UC2) = DF_M3YFR(f_rhon_p,f_rhop_p,f_rhon_t,f_rhop_t,
                            type = 'DDM3Y_Paris' ,
                            E_A= 0.0 ,
                            R_points=np.arange(0.1,18.0,0.1),
                            k_range = (0.0,5.0), num_quadrature= 96,
                            r_range =(0.0,15.0), iter_max=30)
    (R1, UC,US1) = DF_M3YZR(f_rhon_p,f_rhop_p,f_rhon_t,f_rhop_t,
                            E_A=0.0, type='M3Y_Paris',
                            R_points= np.arange(0.1,30.0,0.1),
                            k_range = (0.0,3.0), num_quadrature=70,
                            r_range = (0.0,15.0) )
    plt.plot(R,U_ND+UE,label='M3Y_Paris')
    plt.plot(R1,US1,label='M3Y_Paris_ZR')
    plt.plot(R2,U_ND2+UE2,label='DDM3Y1_Paris')
    plt.legend()

    #-----compare with DFMSPH results-----------------------
    # dat = np.loadtxt('DFMSPH/OUT1_DFMSPH.dat')
    # plt.plot(dat[:,0],dat[:,5],label='DFMSPH_UND')
    # plt.plot(dat[:,0],dat[:,6],label='DFMSPH_UNE')
    # plt.xlim([7,15]);plt.ylim([-200,200]);
    # plt.legend()
