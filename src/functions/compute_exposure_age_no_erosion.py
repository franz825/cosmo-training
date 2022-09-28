def surface_exposure_age_no_erosion(surface_concentration, parameters):

    import math

    # Compute exposure age based on surface concentration, considering no erosion
    exposure_age = - (1 / parameters["L"]) * math.log((1 - surface_concentration * parameters["L"]) / parameters["prod_rate"])

    return exposure_age