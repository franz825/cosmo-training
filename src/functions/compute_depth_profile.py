def compute_nash_index(concentrations_observed, relevant_CRN_to_fit_reality_Be):  # computes the Nash Sutcliffe efficiency index
    
    import numpy as np

    # a list that will collect the difference between observed data and simulated ones (squared), numerator of NSE
    v = []
    # a list that will collect the difference between observed data and the mean observed data (squared), denominator of NSE
    w = []

    for i in range(len(concentrations_observed)):
        
        # variable taking difference value btn observed and simulated
        b = (concentrations_observed[i]-relevant_CRN_to_fit_reality_Be[i])**2
        # variable taking difference value btn observed and mean(observed)
        d = (concentrations_observed[i]-np.mean(concentrations_observed))**2
        # each difference (for each pair of points, one observed and one simulated), is transferred to list 'v'
        v.append(b)
        # each difference (between observed and mean of observations), is transferred to list 'w'
        w.append(d)

    nse_be = 1.0-(sum(v)/sum(w))  # Nash-Sutcliffe efficiency index formula

    return nse_be
    # print('NSE=', nse_be)

def be_accumulator(concentrations_observed, depths_observed, t_exp, erosion_cm, N_inh, parameters):

    import numpy as np

    relevant_CRN_to_fit_reality_Be = []
        
    # Column depth extent in cm, note that prod_depth_spal that distribute production over depth, must always be greater or equal in depth extent than the size of the column
    extent_profile = 5000

    # ATTENUATION LENGTHS OF PRODUCTION PATHWAYS
    Atn = 160                         # Attenuation length for neutrons, g/cm2
    Atnm = 1500                       # Attenuation length for negative muons, g/cm2
    Atfm = 4320                       # Attenuation length for fast muons, g/cm2

    # CONSTANTS FOR 10Be ACCUMULATION (MARGRETH et al. 2016)

    # Relative contribution to 10Be production
    Pn = 0.98039  # relative contribution of neutrons to Be production (unitless)
    Pnm = 0.01471  # relative contribution of negative muons for Be (unitless)
    Pfm = 0.00490  # relative contribution of fast muons for Be (unitless)

    # Production rate and production of 10Be
    Prod_rate = 1.00                   # Prod rate for the exercice
    # Production of neutrons (at/g/yr) via spallation
    Prod_n = Prod_rate*Pn
    Prod_nm = Prod_rate*Pnm            # Production of negative muons (at/g/yr)
    Prod_fm = Prod_rate*Pfm            # Production of fast muons (at/g/yr)

    # GENERATE PROD RATIO

    # depth range, in cm, goes down from 0 at the surface to 3500 cm at depth.
    depth = range(extent_profile)

    cntrlist_spal = []  # empty list where to store the concentration in CRN as a function of depth, from spallation effect
    # empty list where to store the concentration in CRN as a function of depth, from negative muon effect
    cntrlist_neg_muon = []
    # empty list where to store the concentration in CRN as a function of depth, from fast muon effect
    cntrlist_fast_muon = []

    # PRODUCTION OF CRN OVER DEPTH

    for z in depth:                 # for every depth from 0 cm to 3500 cm depth

        P_spal_depth = Prod_n * np.exp(-z * parameters["rho"] / Atn)         # production from spallation
        # production from negative muons
        P_neg_muon_depth = Prod_nm*np.exp(-z * parameters["rho"] / Atnm)
        P_fast_muon_depth = Prod_fm*np.exp(-z * parameters["rho"] / Atfm)  # producion from fast muons

        # gives a list of crns produced via spallation over depth
        cntrlist_spal.append(P_spal_depth)
        # gives a list of crns produced via negative muons over depth
        cntrlist_neg_muon.append(P_neg_muon_depth)
        # gives a list of crns produced via fast muons over depth
        cntrlist_fast_muon.append(P_fast_muon_depth)

    # change list to array vector (improving model performance)
    Prod_depth_spal = np.array(cntrlist_spal)
    Prod_depth_neg_muon = np.array(cntrlist_neg_muon)
    Prod_depth_fast_muon = np.array(cntrlist_fast_muon)

    # Be column
    # array vector building a 1500 cm depth space where to store the values of accumulated CRNs from spallation
    arr_cntrlistfinal_spal = np.zeros(extent_profile)
    # array vector building a 1500 cm depth space where to store the values of accumulated CRNs from negative muons
    arr_cntrlistfinal_neg_muon = np.zeros(extent_profile)
    # array vector building a 1500 cm depth space where to store the values of accumulated CRNs from fast muons
    arr_cntrlistfinal_fast_muon = np.zeros(extent_profile)
    # an array vector building the space to store the final results of the simulation, summing in the end effects from spallation, negative and fast muons
    arr_cntrlistfinal = np.zeros(extent_profile)

    # Calculation of Final Inheritance and adding to any layer of the depth column
    # N_inh_fin is the final inheritance that is the inheritance at t=0 that has decayed over the exposure age period (t_fin)
    N_inh_fin = N_inh*np.exp(-(parameters["L"]*t_exp))
    for i in range(len(arr_cntrlistfinal)):
        # we already store the inheritance in the vector collecting all CRN results at the end of the simulation
        arr_cntrlistfinal[i] = arr_cntrlistfinal[i]+N_inh_fin

    # Calculation of erosion rate
    eros_rate = erosion_cm/t_exp
    cm_eroded = erosion_cm  # amount of top centimeters removed due to erosion
    cm_eroded = int(cm_eroded)

    # removes the uppermost values of CRN's accumulated from spallation
    arr_cntrlistfinal_spal = np.delete(
        arr_cntrlistfinal_spal, range(cm_eroded))
    # removes the uppermost values of CRN's accumulated from negative muons
    arr_cntrlistfinal_neg_muon = np.delete(
        arr_cntrlistfinal_neg_muon, range(cm_eroded))
    arr_cntrlistfinal_fast_muon = np.delete(arr_cntrlistfinal_fast_muon, range(
        cm_eroded))    # removes the uppermost values of CRN's accumulated from fast muons
    # same for the vector aggregating the total CRN's concentration accumulated
    arr_cntrlistfinal = np.delete(arr_cntrlistfinal, range(cm_eroded))

    for i in range(len(arr_cntrlistfinal_spal)):
        arr_cntrlistfinal_spal[i] = arr_cntrlistfinal_spal[i]+((Prod_depth_spal[i])/(
            parameters["L"]+(eros_rate*parameters["rho"])/Atn))*(1-np.exp(-t_exp*(parameters["L"]+(eros_rate*parameters["rho"])/Atn)))

    for i in range(len(arr_cntrlistfinal_neg_muon)):
        arr_cntrlistfinal_neg_muon[i] = arr_cntrlistfinal_neg_muon[i]+(
            (Prod_depth_neg_muon[i])/(parameters["L"]+(eros_rate*parameters["rho"])/Atnm))*(1-np.exp(-t_exp*(parameters["L"]+(eros_rate*parameters["rho"])/Atnm)))

    for i in range(len(arr_cntrlistfinal_fast_muon)):
        arr_cntrlistfinal_fast_muon[i] = arr_cntrlistfinal_fast_muon[i]+(
            (Prod_depth_fast_muon[i])/(parameters["L"]+(eros_rate*parameters["rho"])/Atfm))*(1-np.exp(-t_exp*(parameters["L"]+(eros_rate*parameters["rho"])/Atfm)))

    # Sum up the different contributors to 10Be production

    for i in range(len(arr_cntrlistfinal)):
        arr_cntrlistfinal[i] = arr_cntrlistfinal[i]+arr_cntrlistfinal_spal[i] + \
            arr_cntrlistfinal_neg_muon[i]+arr_cntrlistfinal_fast_muon[i]

    for i in range(len(depths_observed)):
        # fulfills the list with simulated [CRN] to compare with observed [CRN], at respective same depth
        relevant_CRN_to_fit_reality_Be.append(
            arr_cntrlistfinal[(depths_observed[i])])

    nash_index = compute_nash_index(concentrations_observed, relevant_CRN_to_fit_reality_Be)
    # print("Exposure age =", t_exp)
    # print ("Erosion =", erosion_cm)

    result = {
        "exposure_age": t_exp,
        "erosion_rate": erosion_cm,
        "nash_index": nash_index,
        "crn_fitted": relevant_CRN_to_fit_reality_Be
    }

    return result