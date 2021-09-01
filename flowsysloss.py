######  Input file for flow_system_pressure_loss.py

# ##Quantities specific to Oxygen (O2)
# Mm         = 32 #Molar Mass of oxygen
# gamma      = 1.4 #gamma of oxygen
# To         = 300
# Rs         = 8314.5/Mm # calculate the specific gas constant
# SG         = 1.1044 #specific gravity of O2 vs Air

##Quantities specific to Ethylene (C2H4)
Mm         = 28 #Molar Mass of oxygen
gamma      = 1.24 #gamma of oxygen
To         = 300
Rs         = 8314.5/Mm # calculate the specific gas constant
SG         = 0.9686 #specific gravity of O2 vs Air


Dpipe      = 0.277*0.0254 #inches to meters
mdot       = 45 #grams per second
Cv_green   = 7 #catalog Cv for green valve
Cv_black   = 1.5 #catalog Cv for black valve
Cv_solen   = 1.5 #catalog Cv for solenoid valve
Cv_nvalv   = 0.917 #catalog Cv for needle valve
Cv_check   = 1.8 #catalog Cv for check valve
Po1_initial = 800 #Bottle pressure (psi)
