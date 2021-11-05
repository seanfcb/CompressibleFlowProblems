import numpy as np
from scipy.optimize import bisect
from CompressibleFlowFunctions import *
from fannocvloss import *
from tabulate import tabulate
import sys
##### This script runs the input file fannocvloss.py

##==================================================================##
## 2-part algorithm for calculation
## Include new functions in CompressibleFlowFunctions.py when ready
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

#mdot       = 450 ##Uncomment to override value in fannocvloss.py
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
print("Sequence:","P_static","P_total","Mach number")
print("At the dip-stick inlet",round(P1_initial,2),round(Po1_initial,2),round(M1,4))

##==================================================================##
#Check that PHI(M1) > 4fL/D and calculate PHI(M2) if possible
##==================================================================##
fanno_constant = 4*fanning*L/Dpipe
PHI1           = fanno_equation(M1,gamma)
if PHI1 < fanno_constant:
    sys.exit("This pipe will choke before the next flow device")
else:
    PHI2 = PHI1 - fanno_constant

##==================================================================##
#Calculate Po2/Po1 using 2.3
##==================================================================##
