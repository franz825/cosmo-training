# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 11:37:53 2022

@author: nvandermaele
"""
# Keep erosion fixed (try e.g. 1, 10 and 100 m/Myr)
import matplotlib.pyplot as plt
import numpy as np 


def conc_surf(eros_rate,t_exp):     # erosion rate is in cm/yr, conversion factor to obtain m/Myr is *10^4 (1 cm/yr = 10000 m/Myr)
    Atn = 152                         # Attenuation length for neutrons, g/cm2

    # Production rate and production of 10Be (only neutrons)
    Prod_rate = 5.00                  # Prod rate for the exercice


    # Independent constants
    e = 2.718281828459                # neperian number
    rho = 1.7                        # bulk density for sands

    # Lives and decays for 10Be
    half_life_Be10 = 1387000.0  # in years
    L = 0.69314718056/half_life_Be10  # decay constant (lambda in most litterature)


    C_surface=((Prod_rate)/(L+(eros_rate*rho)/Atn))*(1-e**(-t_exp*(L+((eros_rate*rho)/Atn))))
    #print ("10Be concentration=",C_surface)
    return C_surface
    
    
age=[1000,10000,20000,30000, 50000, 100000,200000,300000,400000,500000,1000000]    
conc_plot=np.zeros(11)
    
for i in range(len(age)): 
    Be_concentration = conc_surf(0.001,age[i])
    conc_plot[i]=Be_concentration
    print (conc_plot)
for i in range(len(conc_plot)):
    conc_plot[i]=conc_plot[i]/1000
    
plt.plot(age,conc_plot, 'o')
plt.xlabel("Exposure age (Ma)")
plt.ylabel('Concentration (10$^3$ at/g)')    
    