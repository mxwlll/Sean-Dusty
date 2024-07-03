#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 10:23:54 2023
@author: Sean Dillon
This code will allow us to specify a probability shape distribution that can 
then be used to calculate the average absorption and scattering cross sections 
over the geometric factors L1 and L2
"""


import time
AA = time.time()
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spit
import os
# from CSA_math import *
# print('hello world')


dust_dir = ['/home/physics/Research/DUSTY/DUSTY/Lib_nk/', 
            'C:/Users/Max/Documents/DUSTY/dusty-master/data/Lib_nk']
# this is the possible locations of where dust can be


nk_path = dust_dir[1]               #where the dust is 
dust = 'sil-dlee.nk'                  #DUST NAME HERE #grf
rho = 3.33 #grams cm**-3            #density
pathy = os.path.join(nk_path, dust) #pipeline is open
wavelen, n_dust, k_dust = np.loadtxt(pathy, skiprows=12, unpack=True)
                                    #lamda, n, and k values are extracted
m = np.array([complex(n_dust[i], k_dust[i]) for i in range(len(wavelen))])
                                    #joins n, k values into complex number

##############################################################################
### Bringing in a second dust
dust2 = 'amC-hann.nk' #sil-dl
pathb = os.path.join(nk_path, dust2)
wavelen2, n_dust2, k_dust2 = np.loadtxt(pathb, skiprows=12, unpack=True)
                                    #lamda, n, and k values are extracted
wavelen2 = wavelen2**(-1) * 10000   #Convert wavelen2 to waveLENGTH from waveNUMBER
m2 = np.array([complex(n_dust2[i], k_dust2[i]) for i in range(len(wavelen2))])
                                    #joins n, k values into complex number


##############################################################################
### Bringing in a third dust
dust3 = 'gloliMg50.nk'
pathc = os.path.join(nk_path, dust3)
wavelen3, n_dust3, k_dust3 = np.loadtxt(pathc, skiprows=12, unpack=True)
                                    #lamda, n, and k values are extracted
wavelen3 = wavelen3**(-1) * 10000   #Convert wavelen2 to waveLENGTH from waveNUMBER
m3 = np.array([complex(n_dust3[i], k_dust3[i]) for i in range(len(wavelen3))])
                                    #joins n, k values into complex number


wt_a = 0.34
wt_b = 0.33
wt_c = 0.33
avg_wts=[wt_a,wt_b,wt_c]
# m_avg = np.array([np.average((m[j], m2[j], m3[j]),weights=avg_wts)for j in range(min(len(m),len(m2),len(m3)))])




#wavelen = wavelen**(-1) * 10000   #Convert wavelen2 to waveLENGTH from waveNUMBER



def probability(dis_name, l1, l2, lmin=0.05, m1=0, m2=0, d=0):
    '''
    This is the probability distribution as a function of L1 and L2, the 
    geometric parameters. This parameter gets inserted into the integral that 
    calculates the average polarizability per unit volume

    Parameters
    ----------
    dis_name : String
        This specifies the distribution we will be using. 
        'CDE' = Continuous Distribution of Ellipsoids
        'ERCDE' = Externally Restricted CDE, returns CDE if lmin=0
        'tCDE' = truncated CDE, REQUIRES MORE WORK
    l1 : Float
        Largest Geometric Constant, lmin<l1<1.0
    l2 : Float
        Second Largest Geometric Constant, lmin<l2<=l1
    lmin : Float, optional
        Minimum allowed geometric constant.  The default is 0.

    Returns
    -------
    Float
        Function dependent on l1 and l2

    '''
    l3 = 1 - l1 - l2
    if dis_name == 'CDE':
        return 2
    elif dis_name == 'CDE2':
        return 120 * l1 * l2 * l3
    elif dis_name == 'ERCDE':
        return 2/((1 - (3*lmin))**2)
    elif dis_name == 'tCDE':
        return 1/((1-d-m2)*(1-m1-m2-d) - 0.5*((1-d-m2)**2) - m1**2)

def sigma(m, lamda, v):
    sig = []
    for i in range(len(lamda)):
        k = (2.0 * np.pi)/lamda[i]
        term1 = (6.0*np.pi) / (v * (k**3))
        term2 = np.imag((m[i]**2))
        term3 = 1.0 / abs(m[i]**2 - 1)**2
        sig.append(term1 * term2 * term3)
    return sig


def bounds_l1():
    return [0,1]

def bounds_l2(l1):
    return [0,1-l1]


cabs_cde = []
for j in range(len(m)):
    def f(l1, l2, n=m[j], dis_name='CDE'):
        b = 1/(n**2 - 1)
        term1 = 1/3 * 1/(b + l1)
        term2 = 1/3 * 1/(b + l2)
        term3 = 1/3 * 1/(b + 1 - l1 - l2)
        # r = np.real((term1 + term2 + term3)*probability(dis_name, l1, l2))
        j = np.imag((term1 + term2 + term3)*probability(dis_name, l1, l2))
        return j
        # return np.real((term1 + term2 + term3)*probability(dis_name, l1, l2)) + np.imag((term1 + term2 + term3)*probability(dis_name, l1, l2))
    cabs_cde.append(spit.nquad(f, [bounds_l2, bounds_l1])[0])
kappa_cde = np.array((cabs_cde))
kappa_cde *= (2 * np.pi / (1e-4*wavelen)) / rho

cabs_cde2 = []
for j in range(len(m)):
    def f(l1, l2, n=m[j], dis_name='CDE2'):
        b = 1/(n**2 - 1)
        term1 = 1/3 * 1/(b + l1)
        term2 = 1/3 * 1/(b + l2)
        term3 = 1/3 * 1/(b + 1 - l1 - l2)
        # r = np.real((term1 + term2 + term3)*probability(dis_name, l1, l2))
        j = np.imag((term1 + term2 + term3)*probability(dis_name, l1, l2))
        return j
        # return np.real((term1 + term2 + term3)*probability(dis_name, l1, l2)) + np.imag((term1 + term2 + term3)*probability(dis_name, l1, l2))
    cabs_cde2.append(spit.nquad(f, [bounds_l2, bounds_l1])[0])

kappa_cde2 = np.array((cabs_cde2))
kappa_cde2 *= (2 * np.pi / (1e-4*wavelen)) / rho




cabs_ercde = []
for j in range(len(m)):
    def f(l1, l2, n=m[j], dis_name='ERCDE'):
        b = 1/(n**2 - 1)
        term1 = 1/3 * 1/(b + l1)
        term2 = 1/3 * 1/(b + l2)
        term3 = 1/3 * 1/(b + 1 - l1 - l2)
        # r = np.real((term1 + term2 + term3)*probability(dis_name, l1, l2))
        j = np.imag((term1 + term2 + term3)*probability(dis_name, l1, l2))
        return j
        # return np.real((term1 + term2 + term3)*probability(dis_name, l1, l2)) + np.imag((term1 + term2 + term3)*probability(dis_name, l1, l2))
    cabs_ercde.append(spit.nquad(f, [bounds_l2, bounds_l1])[0])

kappa_ercde = np.array((cabs_ercde))
kappa_ercde *= (2 * np.pi / (1e-4*wavelen)) / rho



# Determine the bounds
x_min = min(wavelen)
x_max = max(wavelen)
y_min = min(min(kappa_cde), min(kappa_cde2), min(kappa_ercde))
y_max = max(max(kappa_cde), max(kappa_cde2), max(kappa_ercde))


# Filter data for wavelengths below 50 microns
mask = wavelen < 70
wavelen_sub50 = wavelen[mask]
kappa_cde_sub50 = kappa_cde[mask]
kappa_cde2_sub50 = kappa_cde2[mask]
kappa_ercde_sub50 = kappa_ercde[mask]

# Determine the bounds for the second plot
x_min_sub50 = min(wavelen_sub50)
x_max_sub50 = max(wavelen_sub50)
y_min_sub50 = min(min(kappa_cde_sub50), min(kappa_cde2_sub50), min(kappa_ercde_sub50))
y_max_sub50 = max(max(kappa_cde_sub50), max(kappa_cde2_sub50), max(kappa_ercde_sub50))

# Create the first subplot for all wavelengths
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# First subplot
ax1.set(xscale='linear', yscale='log')
ax1.set_xlim(x_min, x_max)
ax1.set_ylim(y_min, y_max)
ax1.set_title('Effects of Different Distributions of Grain Shapes \n on Average Index of Refraction', fontsize=16)
ax1.set_xlabel(r'$\lambda (\mu m)$', fontsize=14)
ax1.set_ylabel(r'$<\kappa>$ cm$^{2}$ g$^{-1}$', fontsize=14)
ax1.plot(wavelen, kappa_cde, label='CDE')
ax1.plot(wavelen, kappa_cde2, label='CDE2')
ax1.plot(wavelen, kappa_ercde, label='ERCDE')
ax1.legend()

# 50 micron subplot
ax2.set(xscale='linear', yscale='log')
ax2.set_xlim(x_min_sub50, x_max_sub50)
ax2.set_ylim(y_min_sub50, y_max_sub50)
ax2.set_title('Effects of Different Distributions of Grain Shapes \n on Average Index of Refraction (λ < 50 μm)', fontsize=16)
ax2.set_xlabel(r'$\lambda (\mu m)$', fontsize=14)
ax2.set_ylabel(r'$<\kappa>$ cm$^{2}$ g$^{-1}$', fontsize=14)
ax2.plot(wavelen_sub50, kappa_cde_sub50, label='CDE')
ax2.plot(wavelen_sub50, kappa_cde2_sub50, label='CDE2')
ax2.plot(wavelen_sub50, kappa_ercde_sub50, label='ERCDE')
ax2.legend()

plt.tight_layout()

plt.show()

def volume_integrand_mrn(r, q):
    v = r**(-q)
    return v


rmin = 0.005
rmax = 0.25
q = 3.5

r_integral = spit.quad(volume_integrand_mrn, rmin, rmax, args=q)
r_average = ((1/(rmax - rmin)) * r_integral[0])**(1/-q)
v_avg = (4./3.) * np.pi * r_average**3


# Cabs_array = np.array((cabs))
# sig = np.array((sigma(m_avg, wavelen, v_avg)))
# Csca_array = Cabs_array/sig

# output = np.transpose((wavelen, Cabs_array, Csca_array))



##############################################################################
# This file must start with a three-line header of arbitrary text followed by a 
# three column tabulation of wavelength in um, absorption, and scattering cross 
# Sections

line0 = f'standard mrn mixture (a_min = {rmin}, amax = {rmax})\n'
line1 = f'dust, weight: {dust[:-3]}, {wt_a}; {dust2[:-3]}, {wt_b}; {dust3[:-3]}, {wt_c}\n'
line2 = ' lambda \t <C_abs>/v \t <C_sca>/v\n'

# f = open('look at this.dat', 'w')
# f.write(line0)
# f.write(line1)
# f.write(line2)
# # f.write(str(output))


# for i in range(len(output)):
#     f.write(f"{output[i,0]} \t {output[i,1]} \t {output[i,2]}\n")


# f.close()


BB=time.time()
print('This took {} seconds to run'.format(BB-AA))
