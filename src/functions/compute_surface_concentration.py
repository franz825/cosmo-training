def compute_surface_concentration(erosion_rate, t_exp, parameters):     # erosion rate is in cm/yr, conversion factor to obtain m/Myr is *10^4 (1 cm/yr = 10000 m/Myr)

    # Import packages
    import math

    # Surface concentration (at/g)
    surface_concentration = ((parameters["prod_rate"]) / (parameters['L'] + (erosion_rate * parameters["rho"]) / parameters["atn"])) * (1 - math.exp(-t_exp * (parameters["L"] + ((erosion_rate * parameters["rho"]) / parameters["atn"]))))
    
    return surface_concentration