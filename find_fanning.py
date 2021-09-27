import numpy as np
from scipy.optimize import bisect
from CompressibleFlowFunctions import *

def fanning_search(fanning,Lstar_measured,L_diptube,L_flex1,L_flex2,L_aflange,Dpipe,gamma,Po1_initial,Cv_green,Cv_solen,mdot,Rs,To,Q,SG):
    ##Calculate P, Po, M before green valve

    L_int   = Lstar_measured - L_diptube
    M       = mach_fanno(L_int,fanning,Dpipe,gamma)
    Poratf  = fanno_po_ratio(M_inlet,gamma)
    Postar  = Po1_initial/Poratf
    Po_bval = Postar*fanno_po_ratio(M,gamma)
    P1      = p_from_pratio(Po_bval,gamma,M)
    print("Begin iteration")
    print("Sequence:","P_static","P_total","Mach number")
    print("Before green valve",round(P1,2),round(Po_bval,2),round(M,4))

    #Static and Stagnation Pressure drop through green valve
    P1 = bisect(flowrates, 0, P1,args=(P1,Cv_green,SG,Q))
    M_aval  = bisect(delta_mass_static,0.0001,0.99,args=(mdot,P1*101325/14.7,Rs,To,gamma,Apipe))
    Po_aval = P1*(1+((gamma-1)/2)*M_aval**2)**((gamma)/(gamma-1))
    print("After green valve",round(P1,2),round(Po_aval,2),round(M_aval,4))

    ##Calculate P, Po, M before solenoid valve
    Lstar_measured = L_int
    L_int   = Lstar_measured - L_flex1
    M       = mach_fanno(L_int,fanning,Dpipe,gamma)
    Poratf  = fanno_po_ratio(M_aval,gamma)
    Postar  = Po_aval/Poratf
    Po_bval = Postar*fanno_po_ratio(M,gamma)
    P1      = p_from_pratio(Po_bval,gamma,M)
    print("Before solenoid valve",round(P1,2),round(Po_bval,2),round(M,4))

    #Static and Stagnation Pressure drop through solenoid valve
    P1 = bisect(flowrates, 0, P1,args=(P1,Cv_solen,SG,Q))
    M_aval  = bisect(delta_mass_static,0.0001,0.99,args=(mdot,P1*101325/14.7,Rs,To,gamma,Apipe))
    Po_aval = P1*(1+((gamma-1)/2)*M_aval**2)**((gamma)/(gamma-1))
    print("After solenoid valve",round(P1,2),round(Po_aval,2),round(M_aval,4))

    ##Calculate P, Po, M at open end of pipe in chamber
    Lstar_measured = L_int
    L_int   = Lstar_measured - (L_flex2+L_aflange)
    M_out   = mach_fanno(L_int,fanning,Dpipe,gamma)
    Poratf  = fanno_po_ratio(M_inlet,gamma)
    Postar  = Po_aval/Poratf
    Po_bval = Postar*fanno_po_ratio(M,gamma)
    P1      = p_from_pratio(Po_bval,gamma,M)
    Aeff    = area_from_mass(Po_bval*101325/14.7,To,Rs,gamma,mdot/1000)*1e6
    print("Open end",round(P1,2),round(Po_bval,2),round(M_out,4))
    print("effective area: ",Aeff)
    print("effective area - pipe area: ",round(Aeff-Apipe*1e6,6))
    print("fanning friction factor: ",fanning)
    print("End iteration")
    return Apipe*1e6-Aeff


## Known quantities
Dopipe      = 3/8
PipeT       = 0.049
Mm          = 32 #Molar Mass of oxygen
gamma       = 1.4 #gamma of oxygen
Po1_initial = 800
Po1_metric  = Po1_initial*101325/14.7
To          = 300
Rs          = 8314.5/Mm # calculate the specific gas constant
SG          = 1.1044 #specific gravity of O2 vs Air
mu          = 12e-6
mdot        = 288 ## g/s O2
Cv_green    = 7 #catalog Cv for green valve
Cv_solen    = 1.5 #catalog Cv for solenoid valve

## Length scales
L_diptube = 1.5/3.281
L_flex1    = 3/3.281
L_flex2    = 3/3.281
L_aflange  = 2/3.281
Lstar_measured = L_diptube + L_flex1 + L_flex2 + L_aflange ## Total length of tubing: Dip stick=18". Green valve to solenoid=3'. Solenoid to flange=3'. Flange to open tube=2'.

##Calculate pipe area and pipe diameter in SI units
Dpipe      = (Dopipe - 2*PipeT)*0.0254
Apipe      = np.pi*Dpipe**2/4

##Calculate the volume flow rate Q
Q         = mdot_to_scfh(mdot,Rs,SG) #Volumetric Flow rate in scfh of air at 14.7 psia and 60F.

##Mach number entering the diptube
M_inlet  = bisect(delta_mass_stag,0.0001,0.99,args=(mdot,Po1_initial*101325/14.7,Rs,To,gamma,Apipe))

##First guess for fanning number
fanning = 0.003
P1_initial = p_from_pratio(Po1_initial,gamma,M_inlet)

print("Sequence:","P_static","P_total","Mach number")
print("At the dip-stick inlet",round(P1_initial,2),round(Po1_initial,2),round(M_inlet,4))


fanning = bisect(fanning_search,0.0025,0.003,args=(Lstar_measured,L_diptube,L_flex1,L_flex2,L_aflange,Dpipe,gamma,Po1_initial,Cv_green,Cv_solen,mdot,Rs,To,Q,SG))#,full_output=True)
