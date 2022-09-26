
# ONSET OF THE MODEL
#import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import matplotlib as plt


#import sys

#nb_simus = int(sys.argv[-2])
#name_output = str(sys.argv[-1])


def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)


# OBSERVED DATA IMPLEMENTATION

# write the observed 10Be concentrations here, from top to bottom, and their uncertainty at 1 standard deviation (sd)
observed_data = np.array([150000, 120000, 110000, 103000])
sd_observed_data = np.array([])

# Write the depth corresponding to each 10Be value mentioned above in cm
depth_of_data = [100,170,230,300]


# NASH SUTCLIFFE EFFICIENCY FUNCTION

# list defined to collect model outputs comparable (same depth) with observed [CRN]
relevant_CRN_to_fit_reality_Be = []



def Nash(observed_data, relevant_CRN_to_fit_reality_Be):  # computes the Nash Sutcliffe efficiency index
    # a list that will collect the difference between observed data and simulated ones (squared), numerator of NSE
    v = []
    # a list thet will collect the difference between observed data and the mean observed data (squared), denominator of NSE
    w = []
    for i in range(len(observed_data)):
        
        # variable taking difference value btn observed and simulated
        b = (observed_data[i]-relevant_CRN_to_fit_reality_Be[i])**2
        # variable taking difference value btn observed and mean(observed)
        d = (observed_data[i]-mean(observed_data))**2
        # each difference (for each pair of points, one observed and one simulated), is transferred to list 'v'
        v.append(b)
        # each difference (between observed and mean of observations), is transferred to list 'w'
        w.append(d)
    nse_be = 1.0-(sum(v)/sum(w))  # Nash-Sutcliffe efficiency index formula
    print('NSE=', nse_be)




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


# CONSTANTS FOR EQUATIONS OF THE MODEL (Laloye ET AL. 2019)

# Independent constants
e = 2.718281828459                # neperian number
rho = 2.7                        # bulk density for Granite

# Lives and decays for 10Be
half_life_Be10 = 1387000.0  # in years
L = 0.69314718056/half_life_Be10  # decay constant (lambda in most litterature)


# GENERATE PROD RATIO

# depth range, in cm, goes down from 0 at the surface to 3500 cm at depth.
depth = range(5000)

cntrlist_spal = []  # empty list where to store the concentration in CRN as a function of depth, from spallation effect
# empty list where to store the concentration in CRN as a function of depth, from negative muon effect
cntrlist_neg_muon = []
# empty list where to store the concentration in CRN as a function of depth, from fast muon effect
cntrlist_fast_muon = []


# PRODUCTION OF CRN OVER DEPTH

for z in depth:                 # for every depth from 0 cm to 3500 cm depth

    P_spal_depth = Prod_n*e**(-z*rho/Atn)         # production from spallation
    # production from negative muons
    P_neg_muon_depth = Prod_nm*e**(-z*rho/Atnm)
    P_fast_muon_depth = Prod_fm*e**(-z*rho/Atfm)  # producion from fast muons

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


# CRN ACCUMULATION OVER DEPTH AND OVER EXPOSURE HISTORY

def Be_accumulator(t_exp, Erosion_cm, N_inh):

    relevant_CRN_to_fit_reality_Be = []
    # Column depth extent in cm, note that prod_depth_spal that distribute production over depth, must always be greater or equal in depth extent than the size of the column
    extent = 5000

    # Be column
    # array vector building a 1500 cm depth space where to store the values of accumulated CRNs from spallation
    arr_cntrlistfinal_spal = np.zeros(extent)
    # array vector building a 1500 cm depth space where to store the values of accumulated CRNs from negative muons
    arr_cntrlistfinal_neg_muon = np.zeros(extent)
    # array vector building a 1500 cm depth space where to store the values of accumulated CRNs from fast muons
    arr_cntrlistfinal_fast_muon = np.zeros(extent)
    # an array vector building the space to store the final results of the simulation, summing in the end effects from spallation, negative and fast muons
    arr_cntrlistfinal = np.zeros(extent)

    # Calculation of Final Inheritance and adding to any layer of the depth column
    # N_inh_fin is the final inheritance that is the inheritance at t=0 that has decayed over the exposure age period (t_fin)
    N_inh_fin = N_inh*e**(-(L*t_exp))
    for i in range(len(arr_cntrlistfinal)):
        # we already store the inheritance in the vector collecting all CRN results at the end of the simulation
        arr_cntrlistfinal[i] = arr_cntrlistfinal[i]+N_inh_fin

    # Calculation of erosion rate
    eros_rate = Erosion_cm/t_exp
    cm_eroded = Erosion_cm  # amount of top centimeters removed due to erosion
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
            L+(eros_rate*rho)/Atn))*(1-e**(-t_exp*(L+(eros_rate*rho)/Atn)))

    for i in range(len(arr_cntrlistfinal_neg_muon)):
        arr_cntrlistfinal_neg_muon[i] = arr_cntrlistfinal_neg_muon[i]+(
            (Prod_depth_neg_muon[i])/(L+(eros_rate*rho)/Atnm))*(1-e**(-t_exp*(L+(eros_rate*rho)/Atnm)))

    for i in range(len(arr_cntrlistfinal_fast_muon)):
        arr_cntrlistfinal_fast_muon[i] = arr_cntrlistfinal_fast_muon[i]+(
            (Prod_depth_fast_muon[i])/(L+(eros_rate*rho)/Atfm))*(1-e**(-t_exp*(L+(eros_rate*rho)/Atfm)))

    # Sum up the different contributors to 10Be production

    for i in range(len(arr_cntrlistfinal)):
        arr_cntrlistfinal[i] = arr_cntrlistfinal[i]+arr_cntrlistfinal_spal[i] + \
            arr_cntrlistfinal_neg_muon[i]+arr_cntrlistfinal_fast_muon[i]

    for i in range(len(depth_of_data)):
        # fulfills the list with simulated [CRN] to compare with observed [CRN], at respective same depth
        relevant_CRN_to_fit_reality_Be.append(
            arr_cntrlistfinal[(depth_of_data[i])])

    Nash_result = Nash(observed_data, relevant_CRN_to_fit_reality_Be)
    print("Exposure age =", t_exp)
    print ("Erosion =", Erosion_cm)
    return relevant_CRN_to_fit_reality_Be    #









    # Plots
    arr_depth_of_data = np.array(depth_of_data)

    plt.plot(arr_cntrlistfinal_spal[0:1000],
             depth[0:1000], 'g-', linewidth=1.5)
    plt.plot(arr_cntrlistfinal_neg_muon[0:1000],
             depth[0:1000], 'b-', linewidth=1.5)
    plt.plot(arr_cntrlistfinal_fast_muon[0:1000],
             depth[0:1000], 'm-', linewidth=1.5)
    plt.gca().invert_yaxis()
    plt.gca().xaxis.set_ticks_position('top')
    plt.gca().xaxis.set_label_position('top')
    plt.xlabel('Concentration (at/g)')
    plt.ylabel('depth (cm)')
    plt.show()

    plt.plot(arr_cntrlistfinal_spal[400:1000],
             depth[400:1000], 'g-', linewidth=1.5)
    plt.plot(arr_cntrlistfinal_neg_muon[400:1000],
             depth[400:1000], 'b-', linewidth=1.5)
    plt.plot(arr_cntrlistfinal_fast_muon[400:1000],
             depth[400:1000], 'm-', linewidth=1.5)
    plt.gca().invert_yaxis()
    plt.gca().xaxis.set_ticks_position('top')
    plt.gca().xaxis.set_label_position('top')
    plt.xlabel('Concentration (at/g)')
    plt.ylabel('depth (cm)')
    plt.show()

    #plt.plot(arr_cntrlistfinal_spal[0:1000], depth[0:1000], 'g-', linewidth = 1.5)
    #plt.plot(arr_cntrlistfinal_neg_muon[0:1000], depth[0:1000], 'b-', linewidth = 1.5)
    #plt.plot(arr_cntrlistfinal_fast_muon[0:1000], depth[0:1000], 'm-', linewidth = 1.5)
    plt.plot(arr_cntrlistfinal[0:1000], depth[0:1000], 'k-', linewidth=2)
    plt.gca().invert_yaxis()
    plt.gca().xaxis.set_ticks_position('top')
    plt.gca().xaxis.set_label_position('top')
    plt.xlabel('Concentration (at/g)')
    plt.ylabel('depth (cm)')
    plt.show()

    # errorbar_list=[6934,6763,5791,5693,5727,5712,6387,6153,6267,7114,5665,9279,7336,7524]

    for i in range(len(arr_cntrlistfinal)):
        arr_cntrlistfinal[i] = arr_cntrlistfinal[i]/10000.0

    plt.plot(arr_cntrlistfinal[60:700], depth[60:700], 'k-', linewidth=2)
    plt.plot([18.0504, 16.8231, 12.6051, 14.3094, 14.0484, 13.1744, 12.9244, 12.5922, 12.7087, 14.0898, 12.1969,
             12.8295, 20.8074, 15.1297], [70, 87, 113, 157, 177, 197, 210, 240, 300, 370, 432, 504, 550, 586], 'ro')
    plt.errorbar([18.0504, 16.8231, 12.6051, 14.3094, 14.0484, 13.1744, 12.9244, 12.5922, 12.7087, 14.0898, 12.1969, 12.8295, 20.8074, 15.1297], [70, 87, 113, 157, 177, 197, 210, 240, 300, 370, 432, 504, 550, 586], yerr=[
                 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5], xerr=[0.6934, 0.6763, 0.5791, 0.5693, 0.5727, 0.5712, 0.6387, 0.6153, 0.6267, 0.7114, 0.5665, 0.9279, 0.7336, 0.7524], fmt='none', ecolor='red')
    plt.gca().invert_yaxis()
    plt.gca().xaxis.set_ticks_position('top')
    plt.gca().xaxis.set_label_position('top')
    plt.xlabel('Concentration (10$^4$ at/g)')
    plt.ylabel('depth (cm)')
    plt.savefig('scenario_simpler_800000.jpg', format='jpeg', dpi=1200)
    plt.show()

# Ã
    for i in range(len(observed_data)):
        observed_data[i] = observed_data[i]/1000.0
        arr_depth_of_data[i] = arr_depth_of_data[i]
    for i in range(len(arr_cntrlistfinal)):
        arr_cntrlistfinal[i] = arr_cntrlistfinal[i]*10.0
    xerr = [0.6934, 0.6763, 0.5791, 0.5693, 0.5727, 0.5712, 0.6387,
            0.6153, 0.6267, 0.7114, 0.5665, 0.9279, 0.7336, 0.7524, 0.1]
    for i in range(len(xerr)):
        xerr[i] = xerr[i]*10
    # yerr=[0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05]
    yerr = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
    plt.figure()
    ax = plt.gca()
    for i in range(len(observed_data)):
        observed_data[i] = observed_data[i]-xerr[i]
        arr_depth_of_data[i] = arr_depth_of_data[i]-yerr[i]

    for i in range(len(arr_depth_of_data)):
        ax.add_patch(
            plt.patches.Rectangle(
                # x, y (c'est le point inférieur gauche, donc en Y(depth) il faut retracté 5cm (1sigma) à la valeur depth de base, et en 10Beconc il faut rétracter aussi 1 sigma à l'a valeur propre)
                    (observed_data[i], arr_depth_of_data[i]),
                # width (il faut donc ajouter 2 sigmas pour en avoir un de chaque côté de la valeur centrale)
                xerr[i]*2.0,
                yerr[i]*2.0,        # height (donc il faut ajouter 5+5cm)
                color='darkblue', alpha=0.5
            ))

    # ax.set_aspect(1)
    ax.set_aspect('auto')
    plt.plot(arr_cntrlistfinal[60:700], depth[60:700], 'b', linewidth=2)
    plt.gca().invert_yaxis()
    plt.gca().xaxis.set_ticks_position('top')
    plt.gca().xaxis.set_label_position('top')
    plt.xlabel('Concentration (10$^3$ at/g)')
    plt.ylabel('depth (cm)')

    plt.xlim(80.0, 280.0)
    plt.ylim(700, 60)
    plt.savefig('scenario_box_sans_fit.jpg', format='jpeg', dpi=1200)
    plt.show()

 
#Be_accumulator(275000, 100, 20000)
erosion=[1,2,5,10,20,50]
age=[1000,10000,20000,50000,100000,200000]

for i in range(len(erosion)):
    for j in range(len(age)):
        Be_accumulator(age[j],erosion[i],20000)
        
        

