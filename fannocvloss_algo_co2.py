######  Input file for flow_system_pressure_loss.py
import numpy as np

##Quantities specific to Oxygen (O2)
fluid      = 'oxygen'
Mm         = 32 #Molar Mass of oxygen
gamma      = 1.4 #gamma of oxygen
To         = 300
Rs         = 8314.5/Mm # calculate the specific gas constant
SG         = 1.1044 #specific gravity of O2 vs Air
mdot       = 150 #grams per second
mu         = 12e-6
gammarat   = (gamma+1)/(gamma-1)
specs      = [2.862626*0.0254,30*np.pi/180] #Format for spacer_sizing.[dia,alpha_inj]

# ##Quantities specific to Ethylene (C2H4)
# fluid      = 'ethylene'
# Mm         = 28 #Molar Mass of ethylene
# gamma      = 1.24 #gamma of ethylene
# To         = 300
# Rs         = 8314.5/Mm # calculate the specific gas constant
# SG         = 0.9686 #specific gravity of C2H4 vs Air
# #mdot       = 155 #grams per second
# mu         = 90e-6
# gammarat   = (gamma+1)/(gamma-1)
# specs      = [1.892*0.0254,30*np.pi/180] #Format for spacer_sizing.[dia,alpha_inj]


# ##Quantities specific to Hydrogen (H2)
# fluid      = 'hydrogen'
# Mm         = 2 #Molar Mass of hydrogen
# gamma      = 1.4 #gamma of hydrogen
# To         = 300
# Rs         = 8314.5/Mm # calculate the specific gas constant
# SG         = 0.0696 #specific gravity of H2 vs Air
# mdot       = 155 #grams per second
# mu         = TBD
# gammarat   = (gamma+1)/(gamma-1)
# specs      = [1.892*0.0254,30*np.pi/180] #Format for spacer_sizing.[dia,alpha_inj]

# ##Quantities specific to Hydrogen (H2)
# fluid      = 'co2'
# Mm         = 44 #Molar Mass of hydrogen
# gamma      = 1.28 #gamma of hydrogen
# To         = 300
# Rs         = 8314.5/Mm # calculate the specific gas constant
# SG         = 1.5189 #specific gravity of H2 vs Air
# mdot       = 155 #grams per second
# mu         = TBD
# gammarat   = (gamma+1)/(gamma-1)
# specs      = [1.892*0.0254,30*np.pi/180] #Format for spacer_sizing.[dia,alpha_inj]


## Flow system specific quantities
epsilon     = 0.5e-6
Dpipe       = 11/64 #inches
Cv_solen    = 1.4 #catalog Cv for solenoid valve
Po1_initial = 60 #Bottle pressure (psi)
Po1_metric  = Po1_initial*101325/14.7

## Pipe lengths (currently approx)
L_hose    = 10/3.281
