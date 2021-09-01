import numpy as np
import sys
from scipy.optimize import bisect

###All functions take as an input: pressure in PSI, Temperature in Kelvin, Pipe diameters in inches
###All functions output answers in SI units

def area_from_mass(Po,To,Rs,gamma,mdot):
    ##Function calculates the choking area using the compressible area ratio
    Gstar = Po*np.sqrt(gamma/Rs/To)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))##We call Gstar the ratio mdot/Astar
    Astar = mdot/Gstar
    return Astar


def mass_from_area(M,Po,To,Rs,gamma,Dpipe):
    ##Function calculates the mass flow rate using the compressible area ratio
    Astar = np.pi*Dpipe*Dpipe/4
    #Gstar = Po*np.sqrt(gamma/Rs/To)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))##We call Gstar the ratio mdot/Astar
    Gstar = Po*np.sqrt(gamma/Rs/To)*M*(1+(gamma-2)/2*M*M)**(-(gamma+1)/(2*(gamma-1)))##We call Gstar the ratio mdot/Astar
    mdot  = Gstar*Astar
    return mdot

def prat_from_mach(gamma,M):
    ##For a known pre-shock mach number, what is the stagnation pressure ratio
    pratio = (((gamma+1)*M*M)/((gamma-1)*M*M+2))**(gamma/(gamma-1))*((gamma+1)/(2*gamma*M*M-(gamma-1)))**(1/(gamma-1))
    return pratio

def mach_from_pressure_ratio(Po1,Po2,gamma):
    ##For a desired stagnation pressure ratio, this function calculates the Mach number before a NSW
    Por  = Po2/Po1 ##Desired pressure ratio
    def Prat(Mi,gamma,Por):
        return (((gamma+1)*Mi*Mi)/((gamma-1)*Mi*Mi+2))**(gamma/(gamma-1))*((gamma+1)/(2*gamma*Mi*Mi-(gamma-1)))**(1/(gamma-1)) - Por
    M = bisect(Prat,1,100,args=(gamma,Por))
    return M

def mach_after_shock(M1,gamma):
    ##Calculates the Mach number after a NSW knowing the pre-shock Mach number
    M2 = np.sqrt(((gamma-1)*M1*M1+2)/(2*gamma*M1*M1-(gamma-1)))
    return M2

def pstatic_after_shock(M,gamma,P):
    ##Calculates the static pressure after a NSW knowing the pre-choke Mach number, gamma, and static pressure
    P2 = P*(2*gamma*M*M-(gamma-1))/(gamma+1)
    return P2

def pstag_after_shock(M,gamma,Po1):
    ##Calculates the stagnation pressure after a NSW knowing the pre-choke Mach number, gamma
    Po2 = Po1*(((gamma+1)*M*M)/((gamma-1)*M*M+2))**(gamma/(gamma-1))*((gamma+1)/(2*gamma*M*M-(gamma-1)))**(1/(gamma-1)) ##Solution from the NSW stagnation pressure ratio equation.
    return Po2

def astar_all_else_known(Dpipe,M,gamma):
    ##Function calculates the choking area using the compressible area ratio knowing all other properties
    Apipe  = np.pi*Dpipe**2/4
    Aratio = aratio_from_mach(M,gamma)
    Astar  = Apipe/Aratio
    Dstar  = np.sqrt(Astar*4/np.pi)
    return Astar, Dstar

def mach_from_G(Po,Rs,To,gamma,mdot,Dpipe,subsuper):
    Apipe = np.pi*Dpipe*Dpipe/4
    def delta_G(M,Po,Rs,To,gamma,mdot,Apipe):
        return mdot/Apipe - Po*np.sqrt(gamma/Rs/To)*M*(1+(gamma-2)/2*M*M)**(-(gamma+1)/(2*(gamma-1)))
    if subsuper == 'subsonic':
        M = bisect(delta_G,0.00001,0.99,args=(Po,Rs,To,gamma,mdot,Apipe))
    elif subsuper == 'supersonic':
        M = bisect(delta_G,1,99,args=(Po,Rs,To,gamma,mdot,Apipe))
    else:
        sys.exit('Please specify whether you want to resolve to the "subsonic" or "supersonic" branch when calling mach_from_G')

    return M

def mach_from_aratio(Aexit,Astar,gamma,subsuper):
    ##Function calculates the Mach number at a given location using the compressible area ratio knowing all other properties
    Apipe   = Aexit
    def arat_delta(M,gamma,Apipe,Astar):
        return Apipe/Astar - aratio_from_mach(M,gamma)
    if subsuper == 'subsonic':
        M = bisect(arat_delta,0.00001,0.99,args=(gamma,Apipe,Astar))
    elif subsuper == 'supersonic':
        M = bisect(arat_delta,1,99,args=(gamma,Apipe,Astar))
    else:
        sys.exit('Please specify whether you want to resolve to the "subsonic" or "supersonic" branch when calling mach_from_aratio')
    return M

def mach_from_massflow(Aexit,mdot,Po,To,Rs,gamma,subsuper):
    ##Function calculates the Mach number at a given location using the compressible area ratio knowing all other properties
    Apipe   = Aexit
    Dpipe   = np.sqrt(4*Apipe/np.pi)
    def g_delta(M,Po,To,Rs,gamma,Dpipe,Apipe,mdot):
        return mdot/Apipe - mass_from_area(M,Po,To,Rs,gamma,Dpipe)
    if subsuper == 'subsonic':
        M = bisect(g_delta,0.00001,0.99,args(Po,To,Rs,gamma,Dpipe,Apipe,mdot))
    elif subsuper == 'supersonic':
        M = bisect(g_delta,1,99,args(Po,To,Rs,gamma,Dpipe,Apipe,mdot))
    else:
        sys.exit('Please specify whether you want to resolve to the "subsonic" or "supersonic" branch when calling mach_from_massflow')

    return M

def aratio_from_mach(M,gamma):
    #Function calculates the compressible area ratio A/Astar knowing the Mach number and gamma
    Aratio = ((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))*(1+(gamma-1)/2*M*M)**((gamma+1)/(2*(gamma-1)))/M
    return Aratio

def p_from_pratio(Po,gamma,M):
    ##Function calculates the static pressure knowing the gas properties, Mach number, and stagnation pressure using the isentropic pressure ratio equation P/Po
    P_static = Po*(1+((gamma-1)/2)*M**2)**(-(gamma)/(gamma-1))
    return P_static

def T_from_Tratio(To,gamma,M):
    ##Function calculates the static pressure knowing the gas properties, Mach number, and stagnation pressure using the isentropic pressure ratio equation P/Po
    T_static = To/(1+((gamma-1)/2)*M**2)
    return T_static

def mdot_to_scfh(mdot,Rs,G):
    mdot = mdot*60*2.205/1000 #g/s to lbm/min
    rho  = 101325/(Rs*288.7)*2.205/(3.28084**3) #calculate the stp density, converting to lb/cu.ft
    Q    = mdot/rho
    scfh = Q/np.sqrt(1/G)*60
    return scfh

def hole_numbers(Dhole,Astar):
    numholes = 4*Astar/np.pi/Dhole/Dhole
    return numholes
