# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 13:32:01 2022

@author: nvandermaele
"""
import numpy as np


Prod_rate = 1.00
half_life_Be10 = 1387000.0  
L = 0.69314718056/half_life_Be10
def Surface_exp_age(conc):
    T_exp=-(1/L)*np.log((1-conc*L)/Prod_rate)
    return T_exp


print (Surface_exp_age(96000))