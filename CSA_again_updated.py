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

def a_cde(m):
    top = 2 * m**2 * np.log(m**2)
    bottom = m**2 - 1
    return top/bottom - 2

def a_cds(m):
    a = (1/3) * np.log(m**2)
    b = (4/3) * np.log((1/2)*(m**2 +1))
    return a+b

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









