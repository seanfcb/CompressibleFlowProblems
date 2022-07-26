import numpy as np
import matplotlib.pyplot as plt
import datetime
from scipy.optimize import bisect
from CompressibleFlowFunctions.algos import *
from CompressibleFlowFunctions.Isentropic import *
from CompressibleFlowFunctions.Fanno import *
from CompressibleFlowFunctions.NSW import *
from CompressibleFlowFunctions.misc import *


from fannocvloss_algo_co2 import * ##### This script runs the input file fannocvloss_algo.py
from tabulate import tabulate
import sys
##==================================================================##
## Include new functions in CompressibleFlowFunctions.py when ready
##==================================================================##
print_statements = []
def column(matrix, i):
    return [row[i] for row in matrix]
def tabular_print(*args):
    print_statements.append(args)

def delta_mass(mdot,M,Po,Rs,To,gamma,A): #wrapping delta_mass_stag to iterate over mdot in g/s
    return delta_mass_stag(M,mdot,Po,Rs,To,gamma,A)



##==================================================================##
## 2-part algorithm for calculation
##==================================================================##

## PART 1: CALCULATE FRICTION LOSSES
# Step 1) M1 from eqn 2.2
# Step 2) Friction factor f from Colebrook-White equation
# Step 3) Check PHI(M1) > 4fL/D_h
# Step 4) Calculate PHI(M2) using eqn 2.5
# Step 5) Po2/Po1 using 2.3
# Step 6) Return Po2 & M2

## PART 2: CALCULATE PRESSURE DROP FROM CV ON FLOW DEVICE
##==================================================================##
## Define parameters
##==================================================================##

Dpipe      = 11/64
Dpipe      = Dpipe*0.0254 #(0.75-2*0.065)*0.0254
Apipe      = np.pi*Dpipe**2/4
##==================================================================##
## Define exit conditions
##==================================================================##
Po2 = 14.7
M2  = 1

##==================================================================##
## Calculate conditions at exit of check valve
##==================================================================##
def back_fanno(Po2, Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, SG,Cv_solen,L_hose,fluid):
    print_statements = []
    Po2_init = Po2
    mdot = bisect(delta_mass, 0.0001, 1000000, args=(M2,Po2*101325/14.7,Rs,To,gamma,Apipe))
    Q =  mdot_to_scfh(mdot,Rs,SG)
    P1 = p_from_pratio(Po2,gamma,1)
    ##==================================================================##
    ## Calculate conditions upstream of solenoid valve
    ##==================================================================##
    P2, Po2, M2 = valve_losses_backwards(P1,Cv_solen,SG,Q,mdot,Rs,To,gamma,Apipe)
    f, Re       = fanning_and_reynolds(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon,fluid)
    Lstar2      =  Lstar_fanno(f,Dpipe,M2,gamma)
    tabular_print("Before solenoid valve",round(P2,2),round(Po2,2),round(M2,4),Lstar2)


    ##==================================================================##
    ## Calculate inlet conditions
    ##==================================================================##

    M1, Po1, P1, Po2, P2, Lstar1, Lstar2 = fanno_losses_backwards(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon,L_hose,fluid)
    tabular_print("Bottle regulator exit",round(P1,2),round(Po1,2),round(M1,4),Lstar1)


    tabular_print("Location","P_static","P_total","Mach number","Lstar") #set headers
    tabular_print("Mass Flow Rate",mdot,"g/s") #set headers

    #print(tabulate(reversed(print_statements)))
    print_statements = []
    #print(tabulate(print_statements))
    return Po1, mdot
def fanno_iterator(Po2, Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, fluid, SG,Cv_solen,L_hose):
    Po1, mdot = back_fanno(Po2, Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, SG,Cv_solen,L_hose,fluid)
    return Po1-Po1_initial
#Pbottle = back_fanno(Po1_initial,Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, SG,Cv_ballv,Cv_nvalv,Cv_check,L_to_bval1,L_to_bval2,L_to_needle,L_to_bval3,L_to_check,L_thru_flange)

#
Cv_min           = 0.1
Cv_solen        = 1.2#np.linspace(Cv_min,Cv_nvalv,round((Cv_nvalv-Cv_min)/0.1)+1)#[1,2,3,4,5,6,7,8,9]
result           = []
result_nohead    = []

result.append(["Cv_solen","mdot (g/s)","Pbottle (psi)","Pexit (psi)"])
Pexit = bisect(fanno_iterator,0.1,Po1_initial,args=(Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, fluid, SG,Cv_solen,L_hose))
Pbottle, mdot = back_fanno(Pexit, Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, SG,Cv_solen,L_hose,fluid)
# t = newton(spacer_sizing,10,args=(Pexit*101325/14.7,To,Rs,mdot,gamma,specs))/0.0254*1000
result_nohead.append([Cv_solen,mdot,Pbottle,Pexit])
result.append([Cv_solen,mdot,Pbottle,Pexit])


# Pexit = bisect(fanno_iterator,0.1,Po1_initial,args=(Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, SG,Cv_ballv,Cv_needle,Cv_check,L_to_bval1,L_to_bval2,L_to_needle,L_to_bval3,L_to_check,L_thru_flange))
# Pbottle, mdot, nvalv_rat = back_fanno(Pexit, Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, SG,Cv_ballv,Cv_needle,Cv_check,L_to_bval1,L_to_bval2,L_to_needle,L_to_bval3,L_to_check,L_thru_flange)
# result.append([Cv_needle,mdot,nvalv_rat,Pbottle,Pexit])


print(tabulate(print_statements))
print(tabulate(result))

m = 101325*0.53/(Rs*To)

seconds = m/(mdot/1000) #Time in seconds to fill sandworm
print("Time to fill sandworm with CO2 (min): ")
time = str(datetime.timedelta(seconds=seconds))
print(time)

# Cv_solen   = column(result_nohead,0)
# pressrat = column(result_nohead,2)
# dotm     = column(result_nohead,1)
#
# figure, axis = plt.subplots(1,2)
#
# axis[0].plot(Cv_n,pressrat)
# axis[0].set_xlabel("Cv_needle valve")
# axis[0].set_ylabel("P2/P1")
#
# axis[1].plot(Cv_n,dotm)
# axis[1].set_xlabel("Cv_needle valve")
# axis[1].set_ylabel("mdot (g/s)")
#
# plt.show()
