# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 13:30:58 2023

@author: seand
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spit



def volume_integrand_mrn(r, q):
    v = r**(-q)
    return v

def volume_integrand_kmh(r, q, a0=0.2):
    v = r**(-q) * np.exp(-r/a0)
    return v

rmin = 0.005
rmax = 0.25
q = 3.5

r_mrn = spit.quad(volume_integrand_mrn, rmin, rmax, args=q)
r_mrn_avg = ((1/(rmax - rmin)) * r_mrn[0])**(1/-q)
v_mrn_avg = (4./3.) * np.pi * r_mrn_avg**3


r_kmh = spit.quad(volume_integrand_kmh, rmin, rmax, args=q)
r_kmh_avg = ((1/(rmax - rmin)) * r_kmh[0])**(1/-q)
v_kmh_avg = (4./3.) * np.pi * r_mrn_avg**3

r_mrn_avg_5sf = "%.5f" %r_mrn_avg
r_kmh_avg_5sf = "%.5f" %r_kmh_avg
x = np.linspace(rmin, rmax)

y1 = volume_integrand_mrn(x, 3.5)

y2 = volume_integrand_kmh(x, 3.5)

fig,ax = plt.subplots()
title = "Size Distributions: MRN and KMH"
ax.set(yscale='log')
ax.plot(x, y1, label='MRN')
ax.plot(x, y2, label='KMH')
ax.set_xlabel(r'$\lambda (\mu m)$', fontsize=14)
ax.set_ylabel('n(r)', fontsize=14)
ax.legend()
text = f'Average R: MRN = {r_mrn_avg_5sf} \nAverage R: KMH = {r_kmh_avg_5sf}'
ax.text(0.15, 1e6, text)
ax.set_title(title)
plt.show()




def integrand_Li(s, a, b, c):
    '''
    

    Parameters
    ----------
    s : TYPE
        DESCRIPTION.
    a : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    c : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    return (a*b*c)/(2 * np.sqrt((s+a**2)**3 * (s+b**2) * (s+c**2)))
    
    
    
    
    
    
    
    
    
    
    
def shape_factors(a,b,c):
    '''
    Calculates the Geometric factor L_i for a spheroid of given axes a, b, c

    Parameters
    ----------
    a : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    c : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    L1 = spit.quad(integrand_Li, 0, np.inf, args=(a,b,c))
    L2 = spit.quad(integrand_Li, 0, np.inf, args=(b,c,a))
    L3 = spit.quad(integrand_Li, 0, np.inf, args=(c,b,a))
    return np.array([L1[0], L2[0], L3[0], L1[0]+L2[0]+L3[0]])








print('hello world')