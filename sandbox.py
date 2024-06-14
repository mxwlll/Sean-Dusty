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





dustlist = [('grph1-dl.nk', 'spheres'), ('sil-dlee.nk', 'CDE'), ('FeO.nk', 'CDE2')]
namelist = [dustlist[j][0][:-3]+dustlist[j][1]+'.dat' for j in range(len(dustlist))]
weightlist = [1.0, 200.0, 10.0]
# lam_max, cabs_max, csca_max = np.loadtxt(max(namelist, key=os.path.getsize), unpack=True)
lam_final = np.geomspace(0.2, 500, num=500)
# terp = np.interp(lam_final, lam_max, cabs_max)

total_array = np.ndarray((3,len(lam_final),len(dustlist)))
total_array[:,:,0] = lam_final
# total_array[0,:,1] = cabs_max
# total_array[0,:,2] = csca_max
# for h in range(len(namelist)):


# for k in range(len(namelist)):
#     lam, cabs, csca = np.loadtxt(namelist[k], unpack=True)
#     total_array[k,:,1] = np.interp(lam_final, lam, cabs)
#     total_array[k,:,2] = np.interp(lam_final, lam, csca)


avg_array = np.ndarray((len(lam_final),3))
avg_array[:,0] = lam_final
for j in range(len(lam_final)):
    avg_array[j,1] = np.average(total_array[:,j,1], weights=weightlist)
    avg_array[j,2] = np.average(total_array[:,j,2], weights=weightlist)

titlestring=''
for g in range(len(namelist)):
    titlestring += namelist[g][:3] + str(weightlist[g]).replace('.','')


'''
fig, ax = plt.subplots()

ax.scatter(lam, cabs, marker='o', label='original')
# ax.scatter(x, terp, marker='o',label='interpolated')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()
plt.show()

'''
