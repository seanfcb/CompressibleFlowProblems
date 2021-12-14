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



## Flow system specific quantities
epsilon     = 0.25e-6
Dopipe      = 3/4 #inches
PipeT       = 0.065 #wall thickness in inches
Cv_ballv    = 13.6
Cv_nvalv    = 9.6 #catalog Cv for needle valve
Cv_check    = 4.7 #catalog Cv for check valve
Po1_initial = 800 #Bottle pressure (psi)
Po1_metric  = Po1_initial*101325/14.7

## Pipe lengths (currently approx)
L_to_bval1    = 1.5/3.281
L_to_bval2    = 6/3.281
L_to_needle   = 1/3.281
L_to_bval3    = 1/3.281
L_to_check    = 1.66667/3.281#(2/12)/3.281
L_thru_flange = 5/3.281
