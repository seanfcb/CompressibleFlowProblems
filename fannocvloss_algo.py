######  Input file for flow_system_pressure_loss.py

##Quantities specific to Oxygen (O2)
Mm         = 32 #Molar Mass of oxygen
gamma      = 1.4 #gamma of oxygen
To         = 300
Rs         = 8314.5/Mm # calculate the specific gas constant
SG         = 1.1044 #specific gravity of O2 vs Air
mdot       = 150 #grams per second
mu         = 12e-6
gammarat   = (gamma+1)/(gamma-1)


# ##Quantities specific to Ethylene (C2H4)
# Mm         = 28 #Molar Mass of ethylene
# gamma      = 1.24 #gamma of ethylene
# To         = 300
# Rs         = 8314.5/Mm # calculate the specific gas constant
# SG         = 0.9686 #specific gravity of C2H4 vs Air
# mdot       = 155 #grams per second
# mu         = 90e-6


# ##Quantities specific to Hydrogen (H2)
# Mm         = 2 #Molar Mass of hydrogen
# gamma      = 1.4 #gamma of hydrogen
# To         = 300
# Rs         = 8314.5/Mm # calculate the specific gas constant
# SG         = 0.0696 #specific gravity of H2 vs Air
# mdot       = 155 #grams per second
# mu         = TBD



## Flow system specific quantities
epsilon     = 0.25e-6
Dopipe      = 3/8 #inches
PipeT       = 0.049 #wall thickness in inches
Cv_green    = 7 #catalog Cv for green valve
Cv_black    = 1.5 #catalog Cv for black valve
Cv_solen    = 1.5 #catalog Cv for solenoid valve
Cv_nvalv    = 0.917 #catalog Cv for needle valve
Cv_check    = 1.8 #catalog Cv for check valve
Po1_initial = 800 #Bottle pressure (psi)
Po1_metric  = Po1_initial*101325/14.7

## Pipe lengths (currently approx)
L_diptube     = 1.5/3.281
L_to_black    = 6/3.281
L_to_solen    = 0.5/3.281
L_to_needle   = 1/3.281
L_to_check    = 5/3.281
L_thru_flange = 5/3.281
