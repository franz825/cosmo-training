def compute_exposure_age(surface_concentration, erosion_rate, parameters):
    
    import math

    # Compute exposure age based on surface concentration and erosion rate
    exposure_age = -(1 / (parameters["L"] + ((parameters["rho"] * erosion_rate) / parameters["atn"]))) * math.log(1 - ((surface_concentration * (parameters["L"] + ((parameters["rho"] * erosion_rate) / parameters["atn"]))) / parameters["prod_rate"]))

    return exposure_age