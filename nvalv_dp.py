import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import *
from CompressibleFlowFunctions import *
from fannocvloss_algo import * ##### This script runs the input file fannocvloss_algo.py
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

Dpipe      = Dopipe-2*PipeT
Dpipe      = Dpipe*0.0254 #(0.75-2*0.065)*0.0254
Apipe      = np.pi*Dpipe**2/4
##==================================================================##
## Define exit conditions
##==================================================================##
Po2 = 294.2
M2  = 1

##==================================================================##
## Calculate conditions at exit of check valve
##==================================================================##
def back_fanno(Po2, Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, SG,Cv_ballv,Cv_nvalv,Cv_check,L_to_bval1,L_to_bval2,L_to_needle,L_to_bval3,L_to_check,L_thru_flange):
    print_statements = []

    mdot = newton(delta_mass,100, args=(M2,Po2*101325/14.7,Rs,To,gamma,Apipe))
    Q =  mdot_to_scfh(mdot,Rs,SG)

    M1, Po1, P1, Po2, P2, Lstar1, Lstar2 = fanno_losses_backwards(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon,L_thru_flange)

    tabular_print("Before engine",round(P2,2),round(Po2,2),round(M2,4),Lstar2)
    tabular_print("After needle valve",round(P1,2),round(Po1,2),round(M1,4),Lstar1)

    ##==================================================================##
    ## Calculate conditions upstream of needle valve
    ##==================================================================##
    P2, Po2, M2 = valve_losses_backwards(P1,Cv_nvalv,SG,Q,mdot,Rs,To,gamma,Apipe)
    f, Re       = fanning_and_reynolds(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon)
    Lstar2      =  Lstar_fanno(f,Dpipe,M2,gamma)
    tabular_print("Before needle valve",round(P2,2),round(Po2,2),round(M2,4),Lstar2)
    nvalv_rat   = P2/P1
    M1, Po1, P1, Po2, P2, Lstar1, Lstar2 = fanno_losses_backwards(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon,L_to_bval1)
    tabular_print("Bottle valve inlet conditions",round(P1,2),round(Po1,2),round(M1,4),Lstar1)



    tabular_print("Location","P_static","P_total","Mach number","Lstar") #set headers
    tabular_print("Mass Flow Rate",mdot,"g/s") #set headers

    #print(tabulate(reversed(print_statements)))
    print_statements = []
    #print(tabulate(print_statements))
    return Po1, mdot, nvalv_rat
def fanno_iterator(Po2, Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, SG,Cv_ballv,Cv_nvalv,Cv_check,L_to_bval1,L_to_bval2,L_to_needle,L_to_bval3,L_to_check,L_thru_flange):
    Po1, mdot, nvalv_rat = back_fanno(Po2, Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, SG,Cv_ballv,Cv_nvalv,Cv_check,L_to_bval1,L_to_bval2,L_to_needle,L_to_bval3,L_to_check,L_thru_flange)
    return Po1-Po1_initial
#Pbottle = back_fanno(Po1_initial,Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, SG,Cv_ballv,Cv_nvalv,Cv_check,L_to_bval1,L_to_bval2,L_to_needle,L_to_bval3,L_to_check,L_thru_flange)

#
Cv_needle        = np.linspace(0.1,Cv_nvalv,round((Cv_nvalv-0.1)/0.1)+1) #[1,2,3,4,5,6,7,8,9]
result           = []
result_nohead    = []

result.append(["Cv_needle","mdot (g/s)", "Pratio, needle valve","Pbottle (psi)","Pexit (psi)"])
for x in Cv_needle:
    Pexit = newton(fanno_iterator,0.1,args=(Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, SG,Cv_ballv,x,Cv_check,L_to_bval1,L_to_bval2,L_to_needle,L_to_bval3,L_to_check,L_thru_flange))
    Pbottle, mdot,nvalv_rat = back_fanno(Pexit, Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, SG,Cv_ballv,x,Cv_check,L_to_bval1,L_to_bval2,L_to_needle,L_to_bval3,L_to_check,L_thru_flange)
    result_nohead.append([round(x,2),round(mdot,2),round(nvalv_rat,2),round(Pbottle,2),round(Pexit,2)])
    result.append([round(x,2),round(mdot,2),round(nvalv_rat,2),round(Pbottle,2),round(Pexit,2)])
result.append(["Cv_needle","mdot (g/s)", "Pratio, needle valve","Pbottle (psi)","Pexit (psi)"])


# Pexit = bisect(fanno_iterator,0.1,Po1_initial,args=(Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, SG,Cv_ballv,Cv_needle,Cv_check,L_to_bval1,L_to_bval2,L_to_needle,L_to_bval3,L_to_check,L_thru_flange))
# Pbottle, mdot, nvalv_rat = back_fanno(Pexit, Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, SG,Cv_ballv,Cv_needle,Cv_check,L_to_bval1,L_to_bval2,L_to_needle,L_to_bval3,L_to_check,L_thru_flange)
# result.append([Cv_needle,mdot,nvalv_rat,Pbottle,Pexit])


print(tabulate(print_statements))
print(tabulate(result))

Cv_n     = column(result_nohead,0)
pressrat = column(result_nohead,2)
dotm     = column(result_nohead,1)

figure, axis = plt.subplots(1,2)

axis[0].plot(Cv_n,pressrat)
axis[0].set_xlabel("Cv_needle valve")
axis[0].set_ylabel("P2/P1")

axis[1].plot(Cv_n,dotm)
axis[1].set_xlabel("Cv_needle valve")
axis[1].set_ylabel("mdot (g/s)")

plt.show()
