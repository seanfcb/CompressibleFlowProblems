import numpy as np
from scipy.optimize import bisect
from CompressibleFlowFunctions import *
from fannocvloss_algo import * ##### This script runs the input file fannocvloss_algo.py
from tabulate import tabulate
import sys
##==================================================================##
## Include new functions in CompressibleFlowFunctions.py when ready
##==================================================================##
print_statements = []
def tabular_print(*args):
    print_statements.append(args)


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

mdot       = 774.2 ##Uncomment to override value in fannocvloss.py
Dopipe     = 3/4 ##Uncomment to override value in fannocvloss.py
PipeT      = 0.065 ##Uncomment to override value in fannocvloss.py
Q          = mdot_to_scfh(mdot,Rs,SG) #Volumetric Flow rate in scfh of air at 14.7 psia and 60F.
Dpipe      = Dopipe-2*PipeT
Dpipe      = Dpipe*0.0254 #(0.75-2*0.065)*0.0254
Apipe      = np.pi*Dpipe**2/4
##==================================================================##
##============================PART 1================================##
##==================================================================##
#Calculate the Mach number at the inlet of the pipe.
##==================================================================##
M1  = bisect(delta_mass_stag,0.0001,0.99,args=(mdot,Po1_metric,Rs,To,gamma,Apipe))
##==================================================================##
#Calculate the Fanning friction factor
##==================================================================##
P1_initial = p_from_pratio(Po1_initial,gamma,M1)
Ti         = T_from_Tratio(To,gamma,M1)
rhoi       = P1_initial*(101325/14.7)/(Ti*Rs)
Re         = rhoi*M1*np.sqrt(gamma*Rs*Ti)*Dpipe/mu
darcy      = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
fanning    = darcy/4


tabular_print("Location","P_static","P_total","Mach number","Lstar")
tabular_print("At the dip-stick inlet",round(P1_initial,2),round(Po1_initial,2),round(M1,4))

##==================================================================##
#Check that PHI(M1) > 4fL/D and calculate PHI(M2) if possible
##==================================================================##
fanno_constant = 4*fanning*L_diptube/Dpipe
PHI1           = fanno_equation(M1,gamma)
if PHI1 < fanno_constant:
    print(tabulate(print_statements))
    sys.exit("This pipe will choke before the next flow device")
else:
    PHI2 = PHI1 - fanno_constant

##==================================================================##
#Calculate Po2/Po1 using 2.3
#Return Po2 & M2
##==================================================================##
Lstar   = Lstar_fanno(fanning,Dpipe,M1,gamma)
L_int   = Lstar - L_diptube
M2      = mach_fanno(L_int,fanning,Dpipe,gamma)
Poratf  = fanno_po_ratio(M1,gamma)
Postar  = Po1_initial/Poratf
Po2     = Postar*fanno_po_ratio(M2,gamma)
P2      = p_from_pratio(Po2,gamma,M2)
tabular_print("Before green valve",round(P2,2),round(Po2,2),round(M2,4),Lstar)


print(tabulate(print_statements))
