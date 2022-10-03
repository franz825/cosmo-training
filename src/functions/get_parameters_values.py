def get_parameters_values():
    atn = 160                           # Attenuation length for neutrons, g/cm2

    # Production rate and production of 10Be (only neutrons)
    prod_rate = 4.25                    # Prod rate for the exercice

    # Independent constants
    rho = 2.7                           # bulk density for sands

    # Lives and decays for 10Be
    half_life_Be10 = 1387000.0          # in years
    L = 0.69314718056/half_life_Be10    # decay constant (lambda in most litterature)

    parameters = {
        "atn": atn,
        "prod_rate": prod_rate,
        "rho": rho,
        "half_life_Be10": half_life_Be10,
        "L": L
    }

    return parameters
