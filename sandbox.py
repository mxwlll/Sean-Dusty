# -*- coding: utf-8 -*-
"""
Created on Tue May 21 15:06:56 2024

@author: seand
"""


import numpy as np
import matplotlib.pyplot as plt
import os



# lam_dust1, cabs_dust1, csca_dust1 = np.loadtxt('grph1-dlspheres.dat', unpack=True, skiprows=1)
# lam_dust2, cabs_dust2, csca_dust2 = np.loadtxt('sil-dleespheres.dat', unpack=True, skiprows=1)



# terp = np.interp(lam_dust1, lam_dust2, cabs_dust2)

'''

fig, ax = plt.subplots()

ax.plot(lam_dust1, cabs_dust1, label='cabs1')
ax.plot(lam_dust2, cabs_dust2, label='cabs2')
ax.legend()
plt.show()

'''



dustlist = [('grph1-dl.nk', 'spheres'), ('sil-dlee.nk', 'CDE'), ('FeO.nk', 'CDE2')]
namelist = [dustlist[j][0][:-3]+dustlist[j][1]+'.dat' for j in range(len(dustlist))]

lam_max, cabs_max, csca_max = np.loadtxt(max(namelist, key=os.path.getsize), unpack=True)

total_array = np.ndarray((3,len(lam_max),len(dustlist)))
total_array[:,:,0] = lam_max
total_array[0,:,1] = cabs_max
total_array[0,:,2] = csca_max
# for h in range(len(namelist)):
    