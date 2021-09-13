import numpy as np
from CompressibleFlowFunctions import *
from scipy.optimize import bisect


##Calculate the injection area required to choke the injector, and the needle valve effective diameter with the following flow values:
Po_ox     = 800*101325/14.7 #psi
Po_fuel   = 800*101325/14.7
To         = 300 #Kelvin
Dopipe     = 1/2 #pipe outer diameter in inches
PipeT      = 0.049 #wall thickness in inches
gamma_ox   = 1.4
gamma_fuel = 1.24
mdot_t     = 400 #g/s
Mm_ox      = 32 #molar mass
Mm_fuel    = 28 #molar mass
Ru         = 8314.5
Y_ox       = 96/124
Y_fuel     = 1-Y_ox
epsilon    = 0.25e-6


##Preliminary calculations
mdot_ox   = mdot_t*Y_ox/1000
mdot_fuel = mdot_t*Y_fuel/1000

Rs_ox   = Ru/Mm_ox
Rs_fuel = Ru/Mm_fuel

Dpipe = Dopipe-2*PipeT
Dpipe   = Dpipe*0.0254
Apipe   = np.pi*Dpipe**2/4

## Step 1: Calculate the Darcy friction factor using the Colebrook-White equation and assuming Re = inf
darcy = (-2*np.log10(epsilon/(3.7*Dpipe)))**(-2)

## Step 2: Calculate Astar from flow properties
Astar_ox   = area_from_mass(Po_ox,To,Rs_ox,gamma_ox,mdot_ox)
Astar_fuel = area_from_mass(Po_fuel,To,Rs_fuel,gamma_fuel,mdot_fuel)

## Step 3: Calculate the inlet Mach number
Mi_ox   = mach_from_aratio(Apipe,Astar_ox,gamma_ox,'subsonic')
Mi_fuel = mach_from_aratio(Apipe,Astar_fuel,gamma_fuel,'subsonic')

## Step 4: Knowing the inlet Mach number, calculate inlet pressure
Pi_ox   = p_from_pratio(Po_ox,gamma_ox,Mi_ox)
Pi_fuel = p_from_pratio(Po_fuel,gamma_fuel,Mi_fuel)

## Step 5: Using the equation for friction flow, calcualte the length Lstar to choke the pipe
def fanno(Lstar,f,D,M,gamma): #Define the Fanno equation to iterate on
    return (1-M**2)/(gamma*M**2) + (gamma+1)/(2*gamma)*np.log(((gamma+1)*M**2)/(2*(1+(gamma-1)/2*M**2))) - 4*f*Lstar/D

Lstar_ox   = bisect(fanno,0,1e6,args=(darcy,Dpipe,Mi_ox,gamma_ox))
Lstar_fuel = bisect(fanno,0,1e6,args=(darcy,Dpipe,Mi_fuel,gamma_fuel))
