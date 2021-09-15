import numpy as np
from scipy.optimize import bisect
from CompressibleFlowFunctions import *
from fannocvloss import *
##### This script runs the input file fannocvloss.py

def flowrates(P2,P1,Cv,SG,Q):
    return 42.2*Cv*np.sqrt((P1-P2)*(P1+P2)/SG) - Q ##From the Deltrol catalog, equation relating volume flow rate Q to


### Define parameters
Q          = mdot_to_scfh(mdot,Rs,SG) #Volumetric Flow rate in scfh of air at 14.7 psia and 60F.
Dpipe      = Dopipe-2*PipeT
Dpipe      = Dpipe*0.0254
Apipe      = np.pi*Dpipe**2/4



#Calculate the Mach number at the inlet of the pipe.
Po1_metric = Po1_initial*101325/14.7
Astar_ini  = area_from_mass(Po1_metric,To,Rs,gamma,mdot/1000)
M_initial  = mach_from_aratio(Apipe,Astar_ini,gamma,'subsonic')
P1_initial = p_from_pratio(Po1_initial,gamma,M_initial)
Ti     = T_from_Tratio(To,gamma,M_initial)
rhoi   = P1_initial*(101325/14.7)/(Ti*Rs)
Re     = rhoi*M_initial*np.sqrt(gamma*Rs*Ti)*Dpipe/mu
darcy = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
fanning = darcy/4

#Calculating Mach and pressure before the green valve.
Lstar   = Lstar_fanno(fanning,Dpipe,M_initial,gamma)
L_int   = Lstar - 0.3
M       = mach_fanno(L_int,fanning,Dpipe,gamma)
Poratf  = (1/M_initial)*((2+(gamma-1)*M_initial**2)/(gamma+1))**((gamma+1)/(2*gamma-2))
Postar  = Po1_initial/Poratf
Po_bval = Postar*(1/M)*((2+(gamma-1)*M**2)/(gamma+1))**((gamma+1)/(2*gamma-2))
P1      = p_from_pratio(Po_bval,gamma,M)

#Static and Stagnation Pressure drop through green valve
P1 = bisect(flowrates, 0, P1,args=(P1,Cv_green,SG,Q))
print(P1)
M_aval  = bisect(delta_mass,0.0001,0.99,args=(mdot,P1*101325/14.7,Rs,To,gamma,Apipe))
Po_aval = P1*(1+((gamma-1)/2)*M_aval**2)**((gamma)/(gamma-1))

#Friction losses between green and black valve
Po1_metric = Po_aval*101325/14.7
Astar_ini  = area_from_mass(Po1_metric,To,Rs,gamma,mdot/1000)
M_initial  = mach_from_aratio(Apipe,Astar_ini,gamma,'subsonic')
P1_initial = p_from_pratio(Po_aval,gamma,M_initial)
Ti     = T_from_Tratio(To,gamma,M_initial)
rhoi   = P1_initial*(101325/14.7)/(Ti*Rs)
Re     = rhoi*M_initial*np.sqrt(gamma*Rs*Ti)*Dpipe/mu
darcy = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
fanning = darcy/4

Lstar   = Lstar_fanno(fanning,Dpipe,M_initial,gamma)
L_int   = Lstar - 1.8
M       = mach_fanno(L_int,fanning,Dpipe,gamma)
Poratf  = (1/M_initial)*((2+(gamma-1)*M_initial**2)/(gamma+1))**((gamma+1)/(2*gamma-2))
Postar  = Po1_initial/Poratf
Po_bval = Postar*(1/M)*((2+(gamma-1)*M**2)/(gamma+1))**((gamma+1)/(2*gamma-2))
P1      = p_from_pratio(Po_bval,gamma,M)

#Pressure drop through black valve
P1 = bisect(flowrates, 0.001, P1,args=(P1,Cv_black,SG,Q))
print(P1)
M_aval  = bisect(delta_mass,0.0001,0.99,args=(mdot,P1*101325/14.7,Rs,To,gamma,Apipe))
Po_aval = P1*(1+((gamma-1)/2)*M_aval**2)**((gamma)/(gamma-1))


#Friction losses between black and solenoid valve
Po1_metric = Po_aval*101325/14.7
Astar_ini  = area_from_mass(Po1_metric,To,Rs,gamma,mdot/1000)
M_initial  = mach_from_aratio(Apipe,Astar_ini,gamma,'subsonic')
P1_initial = p_from_pratio(Po_aval,gamma,M_initial)
Ti     = T_from_Tratio(To,gamma,M_initial)
rhoi   = P1_initial*(101325/14.7)/(Ti*Rs)
Re     = rhoi*M_initial*np.sqrt(gamma*Rs*Ti)*Dpipe/mu
darcy = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
fanning = darcy/4

Lstar   = Lstar_fanno(fanning,Dpipe,M_initial,gamma)
L_int   = Lstar - 0.15
M       = mach_fanno(L_int,fanning,Dpipe,gamma)
Poratf  = (1/M_initial)*((2+(gamma-1)*M_initial**2)/(gamma+1))**((gamma+1)/(2*gamma-2))
Postar  = Po1_initial/Poratf
Po_bval = Postar*(1/M)*((2+(gamma-1)*M**2)/(gamma+1))**((gamma+1)/(2*gamma-2))
P1      = p_from_pratio(Po_bval,gamma,M)

#Pressure drop through solenoid valve
P1 = bisect(flowrates, 0.001, P1,args=(P1,Cv_solen,SG,Q))
print(P1)
M_aval  = bisect(delta_mass,0.0001,0.99,args=(mdot,P1*101325/14.7,Rs,To,gamma,Apipe))
Po_aval = P1*(1+((gamma-1)/2)*M_aval**2)**((gamma)/(gamma-1))

#Friction losses between solenoid and needle valve
Po1_metric = Po_aval*101325/14.7
Astar_ini  = area_from_mass(Po1_metric,To,Rs,gamma,mdot/1000)
M_initial  = mach_from_aratio(Apipe,Astar_ini,gamma,'subsonic')
P1_initial = p_from_pratio(Po_aval,gamma,M_initial)
Ti     = T_from_Tratio(To,gamma,M_initial)
rhoi   = P1_initial*(101325/14.7)/(Ti*Rs)
Re     = rhoi*M_initial*np.sqrt(gamma*Rs*Ti)*Dpipe/mu
darcy = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
fanning = darcy/4

Lstar   = Lstar_fanno(fanning,Dpipe,M_initial,gamma)
L_int   = Lstar - 0.3
M       = mach_fanno(L_int,fanning,Dpipe,gamma)
Poratf  = (1/M_initial)*((2+(gamma-1)*M_initial**2)/(gamma+1))**((gamma+1)/(2*gamma-2))
Postar  = Po1_initial/Poratf
Po_bval = Postar*(1/M)*((2+(gamma-1)*M**2)/(gamma+1))**((gamma+1)/(2*gamma-2))
P1      = p_from_pratio(Po_bval,gamma,M)

#Pressure drop through needle valve
P1 = bisect(flowrates, 0.001, P1,args=(P1,Cv_nvalv,SG,Q))
print(P1)
M_aval  = bisect(delta_mass,0.0001,0.99,args=(mdot,P1*101325/14.7,Rs,To,gamma,Apipe))
Po_aval = P1*(1+((gamma-1)/2)*M_aval**2)**((gamma)/(gamma-1))

#Friction losses between needle and check valve
Po1_metric = Po_aval*101325/14.7
Astar_ini  = area_from_mass(Po1_metric,To,Rs,gamma,mdot/1000)
M_initial  = mach_from_aratio(Apipe,Astar_ini,gamma,'subsonic')
P1_initial = p_from_pratio(Po_aval,gamma,M_initial)
Ti     = T_from_Tratio(To,gamma,M_initial)
rhoi   = P1_initial*(101325/14.7)/(Ti*Rs)
Re     = rhoi*M_initial*np.sqrt(gamma*Rs*Ti)*Dpipe/mu
darcy = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
fanning = darcy/4

Lstar   = Lstar_fanno(fanning,Dpipe,M_initial,gamma)
L_int   = Lstar - 1.5
M       = mach_fanno(L_int,fanning,Dpipe,gamma)
Poratf  = (1/M_initial)*((2+(gamma-1)*M_initial**2)/(gamma+1))**((gamma+1)/(2*gamma-2))
Postar  = Po1_initial/Poratf
Po_bval = Postar*(1/M)*((2+(gamma-1)*M**2)/(gamma+1))**((gamma+1)/(2*gamma-2))
P1      = p_from_pratio(Po_bval,gamma,M)


#Pressure drop through check valve
P1 = bisect(flowrates, 0.001, P1,args=(P1,Cv_check,SG,Q))
print(P1)

def delta_mass(M,mdot,P,Rs,To,gamma,A):
    return mdot/1000 - P*(1+(gamma-1)/2*M*M)**(gamma/(gamma-1))*A*np.sqrt(gamma/(Rs*To))*M*(1+(gamma-1)/2*M*M)**(-(gamma+1)/(2*(gamma-1)))

Mf  = bisect(delta_mass,0.0001,0.99,args=(mdot,P1*101325/14.7,Rs,To,gamma,Apipe))
Pof = P1*(1+((gamma-1)/2)*Mf**2)**((gamma)/(gamma-1))
print("The static pressure after the check valve is {} psi".format(P1))
print("The Mach number after the check valve is {} ".format(Mf))
print("The stagnation pressure after the check valve is {} psi".format(Pof))

#Calculating the effective area\
m_star = mass_from_area(Mf,Pof*101325/14.7,To,Rs,gamma,Apipe)*1000
Aeff = area_from_mass(Pof*101325/14.7,To,Rs,gamma,m_star/1000) #assuming a choked line
print("The effective area of the oxygen line is {} sq. mm at {} g/s".format(Aeff*1e6,mdot))
