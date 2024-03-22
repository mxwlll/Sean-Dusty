# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 13:27:28 2024

@author: seand
"""

import matplotlib.pyplot as plt
import numpy as np
print('hello world')



class Cabs(object):
    """Gives the cross section absorption
    
    attributes"""
    def __init__(self, m, l1=(0,1), l2=(0,1-l1), dis_name='CDE',lmin=0.05, m1=0, m2=0, d=0):
        self.m = m
        self.dis_name = dis_name 
        self.l1 = l1
        self.l2 = l2
        self.dis_name = dis_name
        self.lmin = lmin
        self.m1 = m1
        self.m2 = m2
        self.d = d

    def probability_cde(self):
        return 2
    
    def probability_cde2(self, l1, l2):
        l3 = 1 - l1 - l2
        return 120 * l1 * l2 * l3
    
    def probability_ercde(self, l1, l2, lmin):
        ercde = 2/((1 - (3*lmin))**2)
        return ercde
    
    def probability_tcde(self, m1, m2, d):
        tcde = 1/((1-d-m2)*(1-m1-m2-d) - 0.5*((1-d-m2)**2) - m1**2)
        return tcde












print('hello world')













