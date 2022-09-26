# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 11:37:53 2022

@author: nvandermaele
"""
#%% Setup
# Keep erosion fixed (try e.g. 1, 10 and 100 m/Myr)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import display
import plotly.express as px

def conc_surf(erosion_rate, t_exp):     # erosion rate is in cm/yr, conversion factor to obtain m/Myr is *10^4 (1 cm/yr = 10000 m/Myr)

    Atn = 152                           # Attenuation length for neutrons, g/cm2

    # Production rate and production of 10Be (only neutrons)
    Prod_rate = 5.00                    # Prod rate for the exercice

    # Independent constants
    e = 2.718281828459                  # neperian number
    rho = 1.7                           # bulk density for sands

    # Lives and decays for 10Be
    half_life_Be10 = 1387000.0          # in years
    L = 0.69314718056/half_life_Be10    # decay constant (lambda in most litterature)

    # Surface concentration (at/g)
    surface_concentration = ((Prod_rate) / (L + (erosion_rate*rho)/Atn)) * (1 - e**(-t_exp * (L + ((erosion_rate * rho) / Atn))))
    return surface_concentration

# Constant erosion rate of the surface (cm/yr)
erosion_rates = [0.000001, 0.0001, 0.001, 0.03]

# Vector containing the exposure ages (yr) to process
ages = [1000, 10000, 20000, 30000, 50000, 100000, 200000, 300000, 400000, 500000, 1000000]    

#%% Loop

# Containers for 10Be concentrations results
crn_table = pd.DataFrame({"exposure_age": ages})
crn_plot = pd.DataFrame()

# Loop within erosion rate values
for i in range(len(erosion_rates)):

    # Get current erosion rate
    erosion_rate = erosion_rates[i]

    # Array to contain CRN concentration for each 
    crn_concentrations = np.zeros(len(ages))

    # Loop within exposure age values
    for i in range(len(ages)): 
        
        # Get current exposure age
        exposure_age = ages[i]

        # Compute 10Be concentration at the surface 
        Be_concentration = conc_surf(erosion_rate, exposure_age)

        # Divide 10Be concentration by 1000 (and store in array
        crn_concentrations[i] = Be_concentration

    # Convert array of concentrations in a named serie
    crn_concentrations = pd.Series(crn_concentrations, name = "crn_conc_erate_" + str(erosion_rate), dtype = 'Float64')
    # Include array of concentrations in a complete dataframe
    crn_df = pd.DataFrame({"exposure_age": ages, "crn_concentration": crn_concentrations, "erosion_rate": pd.Series(np.repeat(erosion_rate, len(ages)))})

    # Collect concentrations formatted for plots
    crn_plot = pd.concat([crn_plot, crn_df], axis = 0)
    # Collect concentrations into pandas dataframe for csv export
    crn_table = pd.concat([crn_table, crn_concentrations], axis = 1)

#%% Write ouputs to files

# Write results to csv file
crn_table.to_csv("exposure_ages.csv")

# Print the table containing 10Be concentration for each exposure age
display(crn_table)

#%%

# Coerce erosion rate as string and convert at/g to 10^3 at/g
crn_plot["erosion_rate"] = crn_plot["erosion_rate"].astype(str)
crn_plot["crn_concentration"] = crn_plot["crn_concentration"] / 1000

# Plot 10Be concentrations for various erosion rates 
plot = px.scatter(crn_plot, x = "exposure_age", y = "crn_concentration", color = "erosion_rate", labels = dict(exposure_age = "Exposure age (Ma)", crn_concentration = "Concentration (10<sup>3</sup> at/g)", erosion_rate = "Erosion rate (cm /yr)"), log_x = True, log_y = True)
plot.show()


#%% LEGACY
# for i in range(len(crn_concentrations)):
#     crn_concentrations[i] = crn_concentrations[i]/1000
    
# plt.plot(ages, crn_concentrations, 'o')
# plt.xlabel("Exposure age (Ma)")
# plt.ylabel('Concentration (10$^3$ at/g)')    
    
# %%
