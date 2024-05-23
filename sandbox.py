# -*- coding: utf-8 -*-
"""
Created on Tue May 21 15:06:56 2024

@author: seand
"""


import numpy as np





lam_dust1, cabs_dust1, csca_dust1 = np.loadtxt('grph1-dlspheres.dat', unpack=True, skiprows=1)
lam_dust2, cabs_dust2, csca_dust2 = np.loadtxt('sil-dleespheres.dat', unpack=True, skiprows=1)
