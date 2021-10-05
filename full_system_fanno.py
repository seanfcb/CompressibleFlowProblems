import numpy as np
from scipy.optimize import bisect
from CompressibleFlowFunctions import *
from fannocvloss import *
##### This script runs the input file fannocvloss.py


### Define parameters
mdot       = 156 ## This value overrides mdot in the input script fannocvloss.py
Q          = mdot_to_scfh(mdot,Rs,SG) #Volumetric Flow rate in scfh of air at 14.7 psia and 60F.
Dpipe      = Dopipe-2*PipeT
Dpipe      = Dpipe*0.0254
Apipe      = np.pi*Dpipe**2/4
fanning    = 0.002241746878251434 ## Calculated using the find_fanning.py script
Lpipe      = L_diptube+L_to_black+L_to_solen+L_to_needle+L_to_check+L_thru_flange

#Calculate the Mach number at the inlet of the pipe.
Po1_metric = Po1_initial*101325/14.7
M_initial  = bisect(delta_mass_stag,0.0001,0.99,args=(mdot,Po1_metric,Rs,To,gamma,Apipe))
P1_initial = p_from_pratio(Po1_initial,gamma,M_initial)
##new block
# T = T_from_Tratio(To,gamma,M_initial)
# Re = P1_initial*101325/14.7*Dpipe*M_initial*np.sqrt(gamma/(Rs*T))/mu
# fanning    = bisect(colebrook_white,0,1,args=(Re,Dpipe,epsilon))/4
##end block
print("Sequence:","P_static","P_total","Mach number")
print("At the dip-stick inlet",round(P1_initial,2),round(Po1_initial,2),round(M_initial,4))


#Calculating Mach and pressure before the green valve.
Lstar   = Lstar_fanno(fanning,Dpipe,M_initial,gamma)
print("Lstar = ",Lstar)
print("Lpipe = ",Lpipe)
P1, Po_bval, M, L_int = fanno_losses(fanning,Po1_initial,Lstar,L_diptube,Dpipe,M_initial,gamma)
print("Before green valve",round(P1,2),round(Po_bval,2),round(M,4))

#Static and Stagnation Pressure drop through green valve
P1,M_aval,Po_aval = valve_losses(P1,Cv_green,SG,Q,mdot,Rs,To,gamma,Apipe)
print("After green valve",round(P1,2),round(Po_aval,2),round(M_aval,4))


#Friction losses between green and black valve
Lstar   = Lstar_fanno(fanning,Dpipe,M_aval,gamma)
P1, Po_bval, M, L_int = fanno_losses(fanning,Po_aval,Lstar,L_to_black,Dpipe,M_aval,gamma)
print("Before black valve",round(P1,2),round(Po_bval,2),round(M,4))

#Pressure drop through black valve
P1,M_aval,Po_aval = valve_losses(P1,Cv_black,SG,Q,mdot,Rs,To,gamma,Apipe)
print("After black valve",round(P1,2),round(Po_aval,2),round(M_aval,4))


#Friction losses between black and solenoid valve
Lstar   = Lstar_fanno(fanning,Dpipe,M_aval,gamma)
P1, Po_bval, M, L_int = fanno_losses(fanning,Po_aval,Lstar,L_to_solen,Dpipe,M_aval,gamma)
print("Before solenoid valve",round(P1,2),round(Po_bval,2),round(M,4))

#Pressure drop through solenoid valve
P1,M_aval,Po_aval = valve_losses(P1,Cv_solen,SG,Q,mdot,Rs,To,gamma,Apipe)
print("After solenoid valve",round(P1,2),round(Po_aval,2),round(M_aval,4))

#Friction losses between solenoid and needle valve
Lstar   = Lstar_fanno(fanning,Dpipe,M_aval,gamma)
P1, Po_bval, M, L_int = fanno_losses(fanning,Po_aval,Lstar,L_to_needle,Dpipe,M_aval,gamma)
print("Before needle valve",round(P1,2),round(Po_bval,2),round(M,4))

#Pressure drop through needle valve
P1,M_aval,Po_aval = valve_losses(P1,Cv_nvalv,SG,Q,mdot,Rs,To,gamma,Apipe)
print("After needle valve",round(P1,2),round(Po_aval,2),round(M_aval,4))

#Friction losses between needle and check valve
Lstar   = Lstar_fanno(fanning,Dpipe,M_aval,gamma)
P1, Po_bval, M, L_int = fanno_losses(fanning,Po_aval,Lstar,L_to_check,Dpipe,M_aval,gamma)
print("Before check valve",round(P1,2),round(Po_bval,2),round(M,4))


#Pressure drop through check valve
P1,M_aval,Po_aval = valve_losses(P1,Cv_check,SG,Q,mdot,Rs,To,gamma,Apipe)
print("After check valve",round(P1,2),round(Po_aval,2),round(M_aval,4))

#Friction losses into chamber
Lstar   = Lstar_fanno(fanning,Dpipe,M_aval,gamma)
print("Lstar = ",Lstar)
print("Lremaining = ",L_thru_flange)
P1, Po_bval, M, L_int = fanno_losses(fanning,Po_aval,Lstar,L_thru_flange,Dpipe,M_aval,gamma)
print("Open end: ",round(P1,2),round(Po_bval,2),round(M,4))
