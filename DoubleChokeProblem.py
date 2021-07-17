import numpy as np
from CompressibleFlowFunctions import *

##The double choke problem##

##Calculate the injection area required to choke the injector, and the needle valve effective diameter with the following flow values:
Po1        = 800 #psi
Po2        = 200 #psi
To         = 300 #Kelvin
Dopipe     = 3/8 #pipe outer diameter in inches
PipeT      = 0.049 #wall thickness in inches
gamma_ox   = 1.4
gamma_fuel = 1.24
mdot_t     = 200 #g/s
Mm_ox      = 32 #molar mass
Mm_fuel    = 28 #molar mass
Ru         = 8314.5
Y_ox       = 96/124
Y_fuel     = 1-Y_ox
Dhole      = 1/32 #inches

## Part 1: Calculate the injector area knowing the stagnation conditions at state 2
Rs_ox      = Ru/Mm_ox
Rs_fuel    = Ru/Mm_fuel
Astar_ox   = area_from_mass(Po2, To, Rs_ox, gamma_ox, mdot_t*Y_ox/1000)
Astar_fuel = area_from_mass(Po2, To, Rs_fuel, gamma_fuel, mdot_t*Y_fuel/1000)

## Part 1.1: How many should each injector have knowing the hole size is 1/32"?
numholes_ox   = hole_numbers(Dhole,Astar_ox)
numholes_fuel = hole_numbers(Dhole,Astar_fuel)

##Part 2: Assuming there is a NSW at the exit of the needle valve, what needle valve area and effective diameter allow this?
Dpipe = Dopipe-2*PipeT

M_ox         = mach_from_pressure_ratio(Po1, Po2, gamma_ox)
Aneedle_ox, Dneedle_ox   = astar_all_else_known(Dpipe,M_ox,gamma_ox)

M_fuel       = mach_from_pressure_ratio(Po1, Po2, gamma_fuel)
Aneedle_fuel, Dneedle_fuel = astar_all_else_known(Dpipe,M_fuel,gamma_fuel)

Dneedle_ox = Dneedle_ox/0.0254
Dneedle_fuel = Dneedle_fuel/0.0254
