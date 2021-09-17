import numpy as np
from scipy.optimize import bisect
from CompressibleFlowFunctions import *
from flowsysloss import *
##### This script runs the input file flowsysloss.py

def flowrates(P2,P1,Cv,SG,Q):
    return 42.2*Cv*np.sqrt((P1-P2)*(P1+P2)/SG) - Q ##From the Deltrol catalog, equation relating volume flow rate Q to



### Define parameters
Dpipe = Dopipe-2*PipeT
Dpipe   = Dpipe*0.0254
Apipe      = np.pi*Dpipe*Dpipe/4
Q          = mdot_to_scfh(mdot,Rs,SG) #Volumetric Flow rate in scfh of air at 14.7 psia and 60F.

#Calculate the Mach number at the inlet of the pipe.
Po1_metric = Po1_initial*101325/14.7
M_initial  = mach_from_G(Po1_metric,Rs,To,gamma,mdot/1000,Dpipe,'subsonic')
P1_initial = p_from_pratio(Po1_initial,gamma,M_initial)


#Pressure drop through green valve
P1 = bisect(flowrates, 0, P1_initial,args=(P1_initial,Cv_green,SG,Q))
print(P1)
#Pressure drop through black valve
P1 = bisect(flowrates, 0, P1,args=(P1,Cv_black,SG,Q))
print(P1)
#Pressure drop through solenoid valve
P1 = bisect(flowrates, 0, P1,args=(P1,Cv_solen,SG,Q))
print(P1)
#Pressure drop through needle valve
P1 = bisect(flowrates, 0, P1,args=(P1,Cv_nvalv,SG,Q))
print(P1)
#Pressure drop through check valve
P1 = bisect(flowrates, 0, P1,args=(P1,Cv_check,SG,Q))
print(P1)

def delta_mass(M,mdot,P,Rs,To,gamma,A):
    return mdot/1000 - P*(1+(gamma-1)/2*M*M)**(gamma/(gamma-1))*A*np.sqrt(gamma/(Rs*To))*M*(1+(gamma-1)/2*M*M)**(-(gamma+1)/(2*(gamma-1)))

Mf  = bisect(delta_mass,0.0001,0.99,args=(mdot,P1*101325/14.7,Rs,To,gamma,Apipe))
Pof = P1*(1+((gamma-1)/2)*Mf**2)**((gamma)/(gamma-1))
print("The static pressure after the check valve is {} psi".format(P1))
print("The Mach number after the check valve is {} psi".format(Mf))
print("The stagnation pressure after the check valve is {} psi".format(Pof))

# #Calculating the effective area
# Aeff = area_from_mass(P1,300,8314.5/32,1.4,0.15) #assuming a choked line
# print("The effective area of the oxygen line, assuming the line is choked, is {} sq. mm at {} g/s".format(Aeff*1e6,mdot))
