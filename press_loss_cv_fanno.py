import numpy as np
from scipy.optimize import bisect
from CompressibleFlowFunctions import *
from fannocvloss import *
from tabulate import tabulate
##### This script runs the input file fannocvloss.py

print_statements = []
def tabular_print(*args):
    print_statements.append(args)


### Define parameters
mdot = 450
Q          = mdot_to_scfh(mdot,Rs,SG) #Volumetric Flow rate in scfh of air at 14.7 psia and 60F.
Dpipe      = Dopipe-2*PipeT
Dpipe      = (0.75-2*0.065)*0.0254#Dpipe*0.0254
Apipe      = np.pi*Dpipe**2/4



#Calculate the Mach number at the inlet of the pipe.
try:
    Po1_metric = Po1_initial*101325/14.7
    M_initial  = bisect(delta_mass_stag,0.0001,0.99,args=(mdot,Po1_metric,Rs,To,gamma,Apipe))
    P1_initial = p_from_pratio(Po1_initial,gamma,M_initial)
    Ti     = T_from_Tratio(To,gamma,M_initial)
    rhoi   = P1_initial*(101325/14.7)/(Ti*Rs)
    Re     = rhoi*M_initial*np.sqrt(gamma*Rs*Ti)*Dpipe/mu
    darcy = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
    fanning = darcy/4
    badmath = 1/0
except Exception as err:
    print("No.")
tabular_print("Sequence:","P_static","P_total","Mach number","Lstar")
tabular_print("At the dip-stick inlet",round(P1_initial,2),round(Po1_initial,2),round(M_initial,4))
print("Sequence:","P_static","P_total","Mach number","Lstar")
print("At the dip-stick inlet",round(P1_initial,2),round(Po1_initial,2),round(M_initial,4))


#Calculating Mach and pressure before the green valve.
Lstar   = Lstar_fanno(fanning,Dpipe,M_initial,gamma)
L_int   = Lstar - 0.4572
M       = mach_fanno(L_int,fanning,Dpipe,gamma)
Poratf  = fanno_po_ratio(M_initial,gamma)
Postar  = Po1_initial/Poratf
Po_bval = Postar*fanno_po_ratio(M,gamma)
P1      = p_from_pratio(Po_bval,gamma,M)
print("Before green valve",round(P1,2),round(Po_bval,2),round(M,4),Lstar)
tabular_print("Before green valve",round(P1,2),round(Po_bval,2),round(M,4),Lstar)


#Static and Stagnation Pressure drop through green valve
P1 = bisect(flowrates, 0, P1,args=(P1,Cv_green,SG,Q))
M_aval  = bisect(delta_mass_static,0.0001,0.99,args=(mdot,P1*101325/14.7,Rs,To,gamma,Apipe))
Po_aval = P1*(1+((gamma-1)/2)*M_aval**2)**((gamma)/(gamma-1))
print("After green valve",round(P1,2),round(Po_aval,2),round(M_aval,4))
tabular_print("After green valve",round(P1,2),round(Po_aval,2),round(M_aval,4))


#Friction losses between green and black valve
Po1_metric = Po_aval*101325/14.7
M_initial  = bisect(delta_mass_stag,0.0001,0.99,args=(mdot,Po1_metric,Rs,To,gamma,Apipe))
P1_initial = p_from_pratio(Po_aval,gamma,M_initial)
Ti     = T_from_Tratio(To,gamma,M_initial)
rhoi   = P1_initial*(101325/14.7)/(Ti*Rs)
Re     = rhoi*M_initial*np.sqrt(gamma*Rs*Ti)*Dpipe/mu
darcy = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
fanning = darcy/4

Lstar   = Lstar_fanno(fanning,Dpipe,M_initial,gamma)
L_int   = Lstar - 1.8
M       = mach_fanno(L_int,fanning,Dpipe,gamma)
Poratf  = fanno_po_ratio(M_initial,gamma)
Postar  = Po_aval/Poratf
Po_bval = Postar*fanno_po_ratio(M,gamma)
P1      = p_from_pratio(Po_bval,gamma,M)
print("Before black valve",round(P1,2),round(Po_bval,2),round(M,4),Lstar)
tabular_print("Before black valve",round(P1,2),round(Po_bval,2),round(M,4),Lstar)


#Pressure drop through black valve
P1 = bisect(flowrates, 0.001, P1,args=(P1,Cv_black,SG,Q))
M_aval  = bisect(delta_mass_static,0.0001,0.99,args=(mdot,P1*101325/14.7,Rs,To,gamma,Apipe))
Po_aval = P1*(1+((gamma-1)/2)*M_aval**2)**((gamma)/(gamma-1))
print("After black valve",round(P1,2),round(Po_aval,2),round(M_aval,4))
tabular_print("After black valve",round(P1,2),round(Po_aval,2),round(M_aval,4))


#Friction losses between black and solenoid valve
Po1_metric = Po_aval*101325/14.7
M_initial  = bisect(delta_mass_stag,0.0001,0.99,args=(mdot,Po1_metric,Rs,To,gamma,Apipe))
P1_initial = p_from_pratio(Po_aval,gamma,M_initial)
Ti     = T_from_Tratio(To,gamma,M_initial)
rhoi   = P1_initial*(101325/14.7)/(Ti*Rs)
Re     = rhoi*M_initial*np.sqrt(gamma*Rs*Ti)*Dpipe/mu
darcy = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
fanning = darcy/4

Lstar   = Lstar_fanno(fanning,Dpipe,M_initial,gamma)
L_int   = Lstar - 0.15
M       = mach_fanno(L_int,fanning,Dpipe,gamma)
Poratf  = fanno_po_ratio(M_initial,gamma)
Postar  = Po_aval/Poratf
Po_bval = Postar*fanno_po_ratio(M,gamma)
P1      = p_from_pratio(Po_bval,gamma,M)
print("Before solenoid valve",round(P1,2),round(Po_bval,2),round(M,4),Lstar)
tabular_print("Before solenoid valve",round(P1,2),round(Po_bval,2),round(M,4),Lstar)

#Pressure drop through solenoid valve
P1 = bisect(flowrates, 0.001, P1,args=(P1,Cv_solen,SG,Q))
M_aval  = bisect(delta_mass_static,0.0001,0.99,args=(mdot,P1*101325/14.7,Rs,To,gamma,Apipe))
Po_aval = P1*(1+((gamma-1)/2)*M_aval**2)**((gamma)/(gamma-1))
print("After solenoid valve",round(P1,2),round(Po_aval,2),round(M_aval,4))
tabular_print("After solenoid valve",round(P1,2),round(Po_aval,2),round(M_aval,4))

#Friction losses between solenoid and needle valve
Po1_metric = Po_aval*101325/14.7
M_initial  = bisect(delta_mass_stag,0.0001,0.99,args=(mdot,Po1_metric,Rs,To,gamma,Apipe))
P1_initial = p_from_pratio(Po_aval,gamma,M_initial)
Ti     = T_from_Tratio(To,gamma,M_initial)
rhoi   = P1_initial*(101325/14.7)/(Ti*Rs)
Re     = rhoi*M_initial*np.sqrt(gamma*Rs*Ti)*Dpipe/mu
darcy = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
fanning = darcy/4

Lstar   = Lstar_fanno(fanning,Dpipe,M_initial,gamma)
L_int   = Lstar - 0.3
M       = mach_fanno(L_int,fanning,Dpipe,gamma)
Poratf  = fanno_po_ratio(M_initial,gamma)
Postar  = Po_aval/Poratf
Po_bval = Postar*fanno_po_ratio(M,gamma)
P1      = p_from_pratio(Po_bval,gamma,M)
print("Before needle valve",round(P1,2),round(Po_bval,2),round(M,4),Lstar)
tabular_print("Before needle valve",round(P1,2),round(Po_bval,2),round(M,4),Lstar)

#Pressure drop through needle valve
P1 = bisect(flowrates, 0.001, P1,args=(P1,Cv_nvalv,SG,Q))
M_aval  = bisect(delta_mass_static,0.0001,0.99,args=(mdot,P1*101325/14.7,Rs,To,gamma,Apipe))
Po_aval = P1*(1+((gamma-1)/2)*M_aval**2)**((gamma)/(gamma-1))
print("After needle valve",round(P1,2),round(Po_aval,2),round(M_aval,4))
tabular_print("After needle valve",round(P1,2),round(Po_aval,2),round(M_aval,4))


#Friction losses between needle and check valve
Po1_metric = Po_aval*101325/14.7
M_initial  = bisect(delta_mass_stag,0.0001,0.99,args=(mdot,Po1_metric,Rs,To,gamma,Apipe))

P1_initial = p_from_pratio(Po_aval,gamma,M_initial)
Ti     = T_from_Tratio(To,gamma,M_initial)
rhoi   = P1_initial*(101325/14.7)/(Ti*Rs)
Re     = rhoi*M_initial*np.sqrt(gamma*Rs*Ti)*Dpipe/mu
darcy = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
fanning = darcy/4

Lstar   = Lstar_fanno(fanning,Dpipe,M_initial,gamma)
L_int   = Lstar - 1.5
M       = mach_fanno(L_int,fanning,Dpipe,gamma)
Poratf  = fanno_po_ratio(M_initial,gamma)
Postar  = Po_aval/Poratf
Po_bval = Postar*fanno_po_ratio(M,gamma)
P1      = p_from_pratio(Po_bval,gamma,M)
print("Before check valve",round(P1,2),round(Po_bval,2),round(M,4),Lstar)
tabular_print("Before check valve",round(P1,2),round(Po_bval,2),round(M,4),Lstar)


#Pressure drop through check valve
P1 = bisect(flowrates, 0.001, P1,args=(P1,Cv_check,SG,Q))
Mf  = bisect(delta_mass_static,0.0001,0.99,args=(mdot,P1*101325/14.7,Rs,To,gamma,Apipe))
Pof = P1*(1+((gamma-1)/2)*Mf**2)**((gamma)/(gamma-1))

##Final Lstar after check valve
Lstar   = Lstar_fanno(fanning,Dpipe,Mf,gamma)
print("After check valve",round(P1,2),round(Pof,2),round(Mf,4),Lstar)
tabular_print("After check valve",round(P1,2),round(Pof,2),round(Mf,4),Lstar)


print(tabulate(print_statements))
# print("The static pressure after the check valve is {} psi".format(P1))
# print("The Mach number after the check valve is {} ".format(Mf))
# print("The stagnation pressure after the check valve is {} psi".format(Pof))
