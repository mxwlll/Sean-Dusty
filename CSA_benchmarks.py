# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 11:42:15 2023

@author: seand


purpose of this file is to benchmark my code using Min et al 03, Fabian et al 
2001, Corman 2010, depew speck 06 as references. 

Dr Speck said that if this works, i can get a paper out of this
"""
import numpy as np
import matplotlib.pyplot as plt
from CSA_math import *
import os
plt.close('all')

dust_dir = ['/home/physics/Research/DUSTY/DUSTY/Lib_nk/', 
            "C:/UTSA/Research/DUSTY/DUSTY/Lib_nk/"]


nk_path = dust_dir[0]               #where the dust is 
dust = 'oliv_nk_y.nk'                  #DUST NAME HERE
rho = 3.33 #grams cm**-3            #density
pathy = os.path.join(nk_path, dust) #pipeline is open

wavelen, n_dust, k_dust = np.loadtxt(pathy, skiprows=12, unpack=True)
                                    #lamda, n, and k values are extracted

m = np.array([complex(n_dust[i], k_dust[i]) for i in range(len(wavelen))])
                                    #joins n, k values into complex number

wavelen = wavelen**(-1) * 10000   #Convert wavelen2 to waveLENGTH from waveNUMBER


##############################################################################
### Bringing in a second dust
dust2 = 'oliv_nk_x.nk'
pathb = os.path.join(nk_path, dust2)
wavelen2, n_dust2, k_dust2 = np.loadtxt(pathb, skiprows=12, unpack=True)
                                    #lamda, n, and k values are extracted
wavelen2 = wavelen2**(-1) * 10000   #Convert wavelen2 to waveLENGTH from waveNUMBER
m2 = np.array([complex(n_dust2[i], k_dust2[i]) for i in range(len(wavelen2))])
                                    #joins n, k values into complex number


##############################################################################
### Bringing in a third dust
dust3 = 'oliv_nk_z.nk'
pathc = os.path.join(nk_path, dust3)
wavelen3, n_dust3, k_dust3 = np.loadtxt(pathc, skiprows=12, unpack=True)
                                    #lamda, n, and k values are extracted
wavelen3 = wavelen3**(-1) * 10000   #Convert wavelen2 to waveLENGTH from waveNUMBER
m3 = np.array([complex(n_dust3[i], k_dust3[i]) for i in range(len(wavelen3))])
                                    #joins n, k values into complex number



##############################################################################
### Spheres

# asph = np.array(a_sph(m))           #spherical, dust1
# kap_sph = [((2 * np.pi / (1e-4*wavelen[i]))/rho)*np.imag(asph[i])for i in range(len(wavelen))]
# asph2 = np.array(a_sph(m2))         #spherical, dust2
# kap_sph2 = [((2 * np.pi / (1e-4*wavelen2[i]))/rho2)*np.imag(asph2[i])for i in range(len(wavelen2))]





#This calls the function from CSA_again_updated to solve Eqn 39 from Min 03

##############################################################################
### Ellipsoids
acde1 = np.array(a_cde1(m))           #CDE1, dust1
kap_cde1 = [((2 * np.pi / (1e-4*wavelen[i]))/rho)*np.imag(acde1[i])for i in range(len(wavelen))]

# acde2 = np.array(a_cde2(m))           #CDE2, dust1
# kap_cde2 = [((2 * np.pi / (1e-4*wavelen[i]))/rho)*np.imag(acde2[i])for i in range(len(wavelen))]

acde1_2 = np.array(a_cde1(m2))        #CDE1, dust2
kap_cde1_2 = [((2 * np.pi / (1e-4*wavelen2[i]))/rho)*np.imag(acde1_2[i])for i in range(len(wavelen2))]

# acde2_2 = np.array(a_cde2(m2))        #CDE2, dust2
# kap_cde2_2 = [((2 * 20 * np.pi / (1e-4*wavelen2[i]))/rho)*np.imag(acde2_2[i])for i in range(len(wavelen2))]


acde1_3 = np.array(a_cde1(m3))        #CDE1, dust3
kap_cde1_3 = [((2 * np.pi / (1e-4*wavelen3[i]))/rho)*np.imag(acde1_3[i])for i in range(len(wavelen3))]




###############################################################################
#Bring in benchmarking data
testdata = np.loadtxt("cde1_fab01_fig7_olivine.csv", delimiter=',', skiprows=1)

testlamda = np.array([testdata[i,0] for i in range(len(testdata))])
testkappa = np.array([testdata[i,1] for i in range(len(testdata))])

testdatashort = np.loadtxt("cde1_fab01_fig7_olivine_short.csv", delimiter=',', skiprows=1)
testlamda2 = np.array([testdatashort[i,0] for i in range(len(testdatashort))])
testkappa2 = np.array([testdatashort[i,1] for i in range(len(testdatashort))])


# ##############################################################################
### Comparing different shapes
fig1, ax1 = plt.subplots()
title = 'Benchmark - OlivineX - Fabian 2001'
ax1.set(xscale='linear', yscale='log', xlim=(8.5,44), ylim=(1, 10000))
ax1.set_title(title, fontsize=16)
ax1.set_xlabel(r'$\lambda (\mu m)$', fontsize=14)
ax1.set_ylabel(r'$<\kappa>$ cm$^{2}$ g$^{-1}$', fontsize=14)
ax1.plot(wavelen, kap_cde1, label='CDE1 - {}'.format(dust[:-3]))
# ax1.plot(wavelen2, kap_cde1_2, label='CDE1 - {}'.format(dust2[:-3]))
ax1.plot(wavelen3, kap_cde1_3, label='CDE1 - {}'.format(dust3[:-3]))

ax1.plot(testlamda, testkappa,'red', label='Benchmark')
ax1.plot(testlamda2, testkappa2, 'red')


ax1.legend()
fig1.savefig(title+'.png')





##############################################################################
#### Comparing different dusts
# fig2, ax2 = plt.subplots()
# ax2.set(title='Average Mass Absorption Coefficient '+dust, xlabel= r'$\lambda (\mu m)$',
#         ylabel= r'$<\kappa>$ g cm$^{-3}$', xscale='log', yscale='log', xlim=(8.5, 13.0))
# ax2.plot(wavelen, kap_cde, 'b', label=dust)
# ax2.plot(1e-4*wavelen2, kap_cde2, 'k', label=dust2)
# ax2.legend()
# fig2.savefig('Mass Absorption Coefficient Benchmark '+dust[:-3] + ' vs ' +dust2[:-3]+'.png')


# fig3, ax3 = plt.subplots()
# ax3.loglog(wavelen, n_dust, label=dust)
# ax3.loglog(wavelen2, n_dust2, label=dust2)
# ax3.legend()

##############################################################################
#### Cross Section Area



# fig4, ax4 = plt.subplots()
# ax4.set()
# ax4.plot






plt.show()





print('hello world')
