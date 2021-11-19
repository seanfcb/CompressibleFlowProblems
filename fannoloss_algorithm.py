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
#
#
# ##==================================================================##
# ## 2-part algorithm for calculation
# ##==================================================================##
#
# ## PART 1: CALCULATE FRICTION LOSSES
# # Step 1) M1 from eqn 2.2
# # Step 2) Friction factor f from Colebrook-White equation
# # Step 3) Check PHI(M1) > 4fL/D_h
# # Step 4) Calculate PHI(M2) using eqn 2.5
# # Step 5) Po2/Po1 using 2.3
# # Step 6) Return Po2 & M2
#
# ## PART 2: CALCULATE PRESSURE DROP FROM CV ON FLOW DEVICE
# ##==================================================================##
# ## Define parameters
# ##==================================================================##
#
mdot = 963.404
Q          = mdot_to_scfh(mdot,Rs,SG) #Volumetric Flow rate in scfh of air at 14.7 psia and 60F.
Dpipe      = Dopipe-2*PipeT
Dpipe      = Dpipe*0.0254 #(0.75-2*0.065)*0.0254
Apipe      = np.pi*Dpipe**2/4
P1, Po1, M1, Lstar, P2, Po2, M2, Re = fanno_losses(mdot,Rs,SG,Dpipe,Apipe,Po1_initial,Po1_initial*101325/14.7,To,gamma,mu,epsilon,L_to_bval1)
tabular_print("Location","P_static","P_total","Mach number","Lstar")
tabular_print("At the pipe inlet",round(P1,2),round(Po1,2),round(M1,4))
tabular_print("Before bottle valve",round(P2,2),round(Po2,2),round(M2,4),Lstar)

##==================================================================##
#Part 2: Calculate pressure loss through a flow device
#Return Po2 & M2
##==================================================================##
P1,M1,Po1 = valve_losses(P2,Cv_ballv,SG,Q,mdot,Rs,To,gamma,Apipe)
darcy      = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
fanning    = darcy/4
Lstar   = Lstar_fanno(fanning,Dpipe,M1,gamma)
tabular_print("After bottle valve",round(P1,2),round(Po1,2),round(M1,4),Lstar)

##==================================================================##
#Part 1
##==================================================================##
P1, Po1, M1, Lstar, P2, Po2, M2, Re = fanno_losses(mdot,Rs,SG,Dpipe,Apipe,Po1,Po1*101325/14.7,To,gamma,mu,epsilon,L_to_bval2)
tabular_print("Before manual interlock valve",round(P2,2),round(Po2,2),round(M2,4),Lstar)

##==================================================================##
#Part 2
##==================================================================##
P1,M1,Po1 = valve_losses(P2,Cv_ballv,SG,Q,mdot,Rs,To,gamma,Apipe)
darcy      = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
fanning    = darcy/4
Lstar   = Lstar_fanno(fanning,Dpipe,M1,gamma)
tabular_print("After manual interlock valve",round(P1,2),round(Po1,2),round(M1,4),Lstar)


##==================================================================##
#Part 1
##==================================================================##
P1, Po1, M1, Lstar, P2, Po2, M2, Re = fanno_losses(mdot,Rs,SG,Dpipe,Apipe,Po1,Po1*101325/14.7,To,gamma,mu,epsilon,L_to_needle)
tabular_print("Before needle valve",round(P2,2),round(Po2,2),round(M2,4),Lstar)

##==================================================================##
#Part 2
##==================================================================##
P1,M1,Po1 = valve_losses(P2,Cv_nvalv,SG,Q,mdot,Rs,To,gamma,Apipe)
darcy      = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
fanning    = darcy/4
Lstar   = Lstar_fanno(fanning,Dpipe,M1,gamma)
tabular_print("After needle valve",round(P1,2),round(Po1,2),round(M1,4),Lstar)

##==================================================================##
#Part 1
##==================================================================##
P1, Po1, M1, Lstar, P2, Po2, M2, Re = fanno_losses(mdot,Rs,SG,Dpipe,Apipe,Po1,Po1*101325/14.7,To,gamma,mu,epsilon,L_to_bval3)
tabular_print("Before run valve",round(P2,2),round(Po2,2),round(M2,4),Lstar)

##==================================================================##
#Part 2
##==================================================================##
P1,M1,Po1 = valve_losses(P2,Cv_ballv,SG,Q,mdot,Rs,To,gamma,Apipe)
darcy      = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
fanning    = darcy/4
Lstar   = Lstar_fanno(fanning,Dpipe,M1,gamma)
tabular_print("After run valve",round(P1,2),round(Po1,2),round(M1,4),Lstar)

##==================================================================##
#Part 1
##==================================================================##
P1, Po1, M1, Lstar, P2, Po2, M2, Re = fanno_losses(mdot,Rs,SG,Dpipe,Apipe,Po1,Po1*101325/14.7,To,gamma,mu,epsilon,L_to_check)
tabular_print("Before check valve",round(P2,2),round(Po2,2),round(M2,4),Lstar)

##==================================================================##
#Part 2
##==================================================================##
P1,M1,Po1 = valve_losses(P2,Cv_check,SG,Q,mdot,Rs,To,gamma,Apipe)
darcy      = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
fanning    = darcy/4
Lstar   = Lstar_fanno(fanning,Dpipe,M1,gamma)
tabular_print("After check valve",round(P1,2),round(Po1,2),round(M1,4),Lstar)

P1, Po1, M1, Lstar, P2, Po2, M2, Re = fanno_losses(mdot,Rs,SG,Dpipe,Apipe,Po1,Po1*101325/14.7,To,gamma,mu,epsilon,L_thru_flange)
tabular_print("Before engine",round(P2,2),round(Po2,2),round(M2,4),Lstar)


print(tabulate(print_statements))
