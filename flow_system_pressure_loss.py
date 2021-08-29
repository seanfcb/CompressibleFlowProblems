import numpy as np
from scipy import optimize
from CompressibleFlowFunctions import *


def flowrates(P2,P1,Cv,Gf,Q):
    return 42.2*Cv*np.sqrt((P1-P2)*(P1+P2)/Gf) - Q ##From the Deltrol catalog, equation relating volume flow rate Q to

### Define parameters
mdot       = 150 #grams per second
Cv_green   = 7 #catalog Cv for green valve
Cv_black   = 1.5 #catalog Cv for black valve
Cv_solen   = 1.5 #catalog Cv for solenoid valve
Cv_nvalv   = 0.917 #catalog Cv for needle valve
Cv_check   = 1.8 #catalog Cv for check valve
Gf         = 1.1044 #Specific gravity of oxygen
P1_initial = 800 #Bottle pressure (psi)
Q          = 14840 #Volumetric Flow rate in scfh of air at 14.7 psia and 60F.

#Pressure drop through green valve
root = optimize.bisect(flowrates, 0, P1_initial,args=(P1_initial,Cv_green,Gf,Q))
P1 = root
#Pressure drop through black valve
root = optimize.bisect(flowrates, 0, P1,args=(P1,Cv_black,Gf,Q))
P1 = root
#Pressure drop through solenoid valve
root = optimize.bisect(flowrates, 0, P1,args=(P1,Cv_solen,Gf,Q))
P1 = root
#Pressure drop through needle valve
root = optimize.bisect(flowrates, 0, P1,args=(P1,Cv_nvalv,Gf,Q))
P1 = root
#Pressure drop through check valve
root = optimize.bisect(flowrates, 0, P1,args=(P1,Cv_check,Gf,Q))
P1 = root

print("The stagnation pressure after the check valve is {} psi".format(P1))

#Calculating the effective area
Aeff = area_from_mass(P1,300,8314.5/32,1.4,0.15) #assuming a choked line
print("The effective area of the oxygen line, assuming the line is choked, is {} sq. mm at {} g/s".format(Aeff*1e6,mdot))
