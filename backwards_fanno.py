import numpy as np
from scipy.optimize import bisect
from CompressibleFlowFunctions import *
from fannocvloss_algo import * ##### This script runs the input file fannocvloss_algo.py
from tabulate import tabulate
#from test import * ##temp function
import sys
##==================================================================##
## Include new functions in CompressibleFlowFunctions.py when ready
##==================================================================##
print_statements = []

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

    mdot = bisect(delta_mass, 0.0001, 1000000, args=(M2,Po2*101325/14.7,Rs,To,gamma,Apipe))
    Q =  mdot_to_scfh(mdot,Rs,SG)

    M1, Po1, P1, Po2, P2, Lstar1, Lstar2 = fanno_losses_backwards(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon,L_thru_flange)

    tabular_print("Before engine",round(P2,2),round(Po2,2),round(M2,4),Lstar2)
    tabular_print("After check valve",round(P1,2),round(Po1,2),round(M1,4),Lstar1)

    ##==================================================================##
    ## Calculate conditions upstream of check valve
    ##==================================================================##
    P2, Po2, M2 = valve_losses_backwards(P1,Cv_check,SG,Q,mdot,Rs,To,gamma,Apipe)
    f, Re       = fanning_and_reynolds(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon)
    Lstar2      =  Lstar_fanno(f,Dpipe,M2,gamma)
    tabular_print("Before check valve",round(P2,2),round(Po2,2),round(M2,4),Lstar2)

    ##==================================================================##
    ## Calculate conditions downstream of run valve
    ##==================================================================##

    M1, Po1, P1, Po2, P2, Lstar1, Lstar2 = fanno_losses_backwards(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon,L_to_check)
    tabular_print("After run valve",round(P1,2),round(Po1,2),round(M1,4),Lstar1)
    ##==================================================================##

    ##==================================================================##
    ## Calculate conditions upstream of run valve
    ##==================================================================##
    P2, Po2, M2 = valve_losses_backwards(P1,Cv_ballv,SG,Q,mdot,Rs,To,gamma,Apipe)
    f, Re       = fanning_and_reynolds(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon)
    Lstar2      =  Lstar_fanno(f,Dpipe,M2,gamma)
    tabular_print("Before run valve",round(P2,2),round(Po2,2),round(M2,4),Lstar2)

    ##==================================================================##
    ## Calculate conditions downstream of needle valve
    ##==================================================================##

    M1, Po1, P1, Po2, P2, Lstar1, Lstar2 = fanno_losses_backwards(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon,L_to_bval3)
    tabular_print("After needle valve",round(P1,2),round(Po1,2),round(M1,4),Lstar1)
    ##==================================================================##

    ##==================================================================##
    ## Calculate conditions upstream of needle valve
    ##==================================================================##
    P2, Po2, M2 = valve_losses_backwards(P1,Cv_nvalv,SG,Q,mdot,Rs,To,gamma,Apipe)
    f, Re       = fanning_and_reynolds(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon)
    Lstar2      =  Lstar_fanno(f,Dpipe,M2,gamma)
    tabular_print("Before needle valve",round(P2,2),round(Po2,2),round(M2,4),Lstar2)

    ##==================================================================##
    ## Calculate conditions downstream of manual interlock valve
    ##==================================================================##

    M1, Po1, P1, Po2, P2, Lstar1, Lstar2 = fanno_losses_backwards(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon,L_to_needle)
    tabular_print("After manual interlock valve",round(P1,2),round(Po1,2),round(M1,4),Lstar1)

    ##==================================================================##
    ## Calculate conditions upstream of manual interlock valve
    ##==================================================================##
    P2, Po2, M2 = valve_losses_backwards(P1,Cv_ballv,SG,Q,mdot,Rs,To,gamma,Apipe)
    f, Re       = fanning_and_reynolds(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon)
    Lstar2      =  Lstar_fanno(f,Dpipe,M2,gamma)
    tabular_print("Before manual interlock  valve",round(P2,2),round(Po2,2),round(M2,4),Lstar2)

    ##==================================================================##
    ## Calculate conditions downstream of bottle valve
    ##==================================================================##

    M1, Po1, P1, Po2, P2, Lstar1, Lstar2 = fanno_losses_backwards(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon,L_to_bval2)
    tabular_print("After bottle valve",round(P1,2),round(Po1,2),round(M1,4),Lstar1)

    ##==================================================================##
    ## Calculate conditions upstream of bottle valve
    ##==================================================================##
    P2, Po2, M2 = valve_losses_backwards(P1,Cv_ballv,SG,Q,mdot,Rs,To,gamma,Apipe)
    f, Re       = fanning_and_reynolds(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon)
    Lstar2      =  Lstar_fanno(f,Dpipe,M2,gamma)
    tabular_print("Before bottle valve",round(P2,2),round(Po2,2),round(M2,4),Lstar2)

    ##==================================================================##
    ## Calculate inlet conditions
    ##==================================================================##

    M1, Po1, P1, Po2, P2, Lstar1, Lstar2 = fanno_losses_backwards(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon,L_to_bval1)
    tabular_print("Bottle valve inlet conditions",round(P1,2),round(Po1,2),round(M1,4),Lstar1)



    tabular_print("Location","P_static","P_total","Mach number","Lstar") #set headers
    tabular_print("Mass Flow Rate",mdot,"g/s") #set headers

    #print(tabulate(reversed(print_statements)))
    print_statements = []
    #print(tabulate(print_statements))
    return Po1-Po1_initial
#Pbottle = back_fanno(Po1_initial,Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, SG,Cv_ballv,Cv_nvalv,Cv_check,L_to_bval1,L_to_bval2,L_to_needle,L_to_bval3,L_to_check,L_thru_flange)

#
Pexit = bisect(back_fanno,14.7,Po1_initial,args=(Po1_initial, M2, Rs, To, gamma, Apipe, Dpipe, mu, epsilon, SG,Cv_ballv,Cv_nvalv,Cv_check,L_to_bval1,L_to_bval2,L_to_needle,L_to_bval3,L_to_check,L_thru_flange))
print(tabulate(print_statements))
