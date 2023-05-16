# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 12:24:51 2022

@author: seand
"""

import numpy as np


def a_sph(m):
    top = m**2 - 1
    bottom = m**2 +2
    return top/bottom

def a_cde1(m):
    top = 2 * m**2 * np.log(m**2)
    bottom = m**2 - 1
    return top/bottom - 2

def a_cde2(m):
    '''
    calculates the polarizability per unit volume for a continuous distributions
    of ellipsoids where near-spherical particle shapes are most probable

    Parameters
    ----------
    m : complex
        Complex index of refraction. m = n + ik 

    Returns
    -------
    alpha: complex

    '''
    a = (m**6)/((m**2 - 1)**4)
    b = np.log(m**2)
    c = (m**2 - 1)**(-3)
    d = 2.5 * (m**2 - 1)**(-2)
    e = 11/(6 * (m**2 - 1))
    f = 0.25
    return -a*b + c + d + e + f








def a_cds_unlim(m):
    a = (1/3) * np.log(m**2)
    b = (4/3) * np.log((1/2)*(m**2 +1))
    return a+b


def a_cds_lim(m, a, b):
    '''
    Parameters
    ----------
    m : complex
        complex index of refraction of dust.
    a : float
        minimal spheroidal axis, units in micron
    b : float
        maximal spheroidal axis, units in micron

    Returns
    -------
    complex
        average polarizability per unit volume for a limited distribution of 
        spheroids.

    '''

    xi_min = a/b
    xi_max = b/a
    
    # need to calculate Lmax and Lmin
    #Lmin = problate
    #Lmax = oblate
    e2min = 1.0 - xi_min**2 
    e2max = 1.0 - xi_max**2
    Lmin1 = (1.0 - e2min)/(e2min)
    Lmin2 = (-1.0 + (1.0/(2.0 * np.sqrt(e2min))) * np.ln((1 + np.sqrt(e2min))/
                                                         (1 - np.sqrt(e2min))))
    Lmin = Lmin1 * Lmin2
    
    Lmax1 = 1.0/e2max
    Lmax2 = 1.0 - (np.sqrt(1 - e2max))/np.sqrt(e2max) * np.arcsin(np.sqrt(e2max))
    Lmax = Lmax1 * Lmax2    
    
    a1 = 1.0/(3.0 * (Lmax-Lmin))
    a2 = np.ln((1.0 + Lmax * (m**2 - 1.0)))
    a3 = 4.0 * np.ln((1.0 + m**2 + Lmin * (1.0 - m**2))/
                     (1.0 + m**2 + Lmax * (1.0 - m**2)))
    a = a1 * (a2 + a3)
    return a











def a_dhs(m):
    top = 6*m**2 + 3
    bottom = 2*m**2 - 2
    inner = 2*m**2 + 1
    innertop = m**2 + 2
    innerbot = 9*m**2
    return (top/bottom)* np.log(inner * (innertop/innerbot))









def absorption(a, lamda, volume, shape):
    C_abs = []
    if shape == 'spherical':
        for i in len(lamda):
            term1 = 4.0 * np.pi * (2 * np.pi / lamda[i])
            term2 = a[i].imag()
            C_abs.append(term1 * term2 * volume)
    else:
        for i in len(lamda):
            term1 = (4.0/3.0) * np.pi * (2.0 * np.pi / lamda[i])
            term2 = a[i].imag()
            C_abs.append(term1 * term2 * volume)
    return C_abs

def scattering(a, lamda, volume, shape, sigma):
    C_sca = []
    if shape == 'spherical':
        for i in len(volume):
            term1 = (8.0 / 3.0) * np.pi * ((2.0 * np.pi / lamda[i])**4.0)
            term2 = abs(a[i].imag())**2
            C_sca.append(term1 * term2)
    else:
        for i in len(volume):
            term1 = (4.0) * np.pi * (2.0 * np.pi / lamda[i])
            term2 = 3.0 * sigma[i]
            term3 = term1 / term2
            term4 = a[i].imag()
            C_sca.append(term3 * term4 * volume[i])
    return C_sca

def sigma(m, lamda):
    sig = []
    for i in lamda:
        term1 = (9.0 / (2.0 * (2.0 * np.pi / lamda[i])**3))
        term2 = (m[i]**2).imag()
        term3 = 1.0 / abs(m[i]**2 - 1)**2
        sig.append(term1 * term2 * term3)
    return sig










def spheroid_distributions(a, b):
    if a/b <= 1:        #problate
        #on this line
        e = np.sqrt(1- (a**2)/(b**2))
    elif a/b >1:        #oblate
        e = np.sqrt(1- (b**2)/(a**2))
    return e 

    

    
    
    
    
    
    
    
    
    
    
    
    
    



