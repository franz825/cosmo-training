def surface_concentration(eros_rate,t_exp):     # erosion rate is in cm/yr, conversion factor to obtain m/Myr is *10^4 (1 cm/yr = 10000 m/Myr)
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