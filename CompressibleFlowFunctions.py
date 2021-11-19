import numpy as np
import sys
from scipy.optimize import bisect

###All functions take as an input: pressure in PSI, Temperature in Kelvin, Pipe diameters in inches
###All functions output answers in SI units

def fanno_losses_backwards(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon,L): #function to be added to CompressibleFlowFunctions.py
    '''
    Function calculates initial conditions in a friction pipe knowing the exit conditions
    Expected inputs:
    Po2      : Exit Stagnation pressure, PSI
    To       : Stagnation temperature, K
    gamma    : Ratio of specific heats
    M2       : Exit Mach number
    Rs       : Specific gas constant, J/kgK (double check units)
    Dpipe    : Pipe diameter, meters
    mu       : Dynamic viscosity
    epsilon  : Surface roughness
    L        : Pipe length, meters
    '''
    PHI2           = fanno_equation(M2,gamma)
    f, Re          = fanning_and_reynolds(Po2,To,gamma,M2,Rs,Dpipe,mu,epsilon)
    Lstar2         = Lstar_fanno(f,Dpipe,M2,gamma)
    fanno_constant = 4*f*L/Dpipe
    PHI1           = fanno_constant + PHI2
    Lstar1         = Lstar2 + L
    M1             = bisect(delta_fanno,0.001,0.9999,args=(Lstar1,f,Dpipe,gamma))
    Poratf  = fanno_po_ratio(M2,gamma)
    Postar  = Po2/Poratf
    Po1     = Postar*fanno_po_ratio(M1,gamma)
    P1      = p_from_pratio(Po1,gamma,M1)
    P2      = p_from_pratio(Po2,gamma,M2)

    return M1, Po1, P1, Po2, P2, Lstar1, Lstar2

def flowrates_backwards(P1,P2,Cv,SG,Q):
    '''
    Function simply wraps the flowrates() function to iterate on inlet pressure. Returns inlet pressure
    Expected inputs:
    P1 and P2: Pressures upstream and downstream, PSI
    Cv       : Flow coefficient
    SG       : Specific gravity w.r.t. air
    Q        : Volumetric flow rate, SCFH (see mdot_to_scfh)
    '''
    return flowrates(P2,P1,Cv,SG,Q)

def valve_losses_backwards(P1,Cv,SG,Q,mdot,Rs,To,gamma,Apipe):
    P_bval  = bisect(flowrates_backwards,P1, 10000, args=(P1,Cv,SG,Q))
    M_bval  = bisect(delta_mass_static,0.0000001,0.99999999,args=(mdot,P_bval*101325/14.7,Rs,To,gamma,Apipe))
    Po_bval = P_bval/(1+((gamma-1)/2)*M_bval**2)**(-(gamma)/(gamma-1))

    return P_bval, Po_bval, M_bval


def fanning_and_reynolds(Po1,To,gamma,M,Rs,Dpipe,mu,epsilon):
    P1         = p_from_pratio(Po1,gamma,M)
    T1         = T_from_Tratio(To,gamma,M)
    rhoi       = P1*(101325/14.7)/(T1*Rs)
    Re         = rhoi*M*np.sqrt(gamma*Rs*T1)*Dpipe/mu
    darcy      = bisect(colebrook_white,1e-6,1,args=(Re,Dpipe,epsilon))
    fanning    = darcy/4

    return fanning, Re

def fanno_losses(mdot,Rs,SG,Dpipe,Apipe,Po1,Po1_metric,To,gamma,mu,epsilon,L):
    ##==================================================================##
    ##============================PART 1================================##
    ##==================================================================##
    #Calculate the Mach number at the inlet of the pipe.
    ##==================================================================##
    M1  = bisect(delta_mass_stag,0.0001,0.99,args=(mdot,Po1_metric,Rs,To,gamma,Apipe))
    ##==================================================================##
    #Calculate the Fanning friction factor
    ##==================================================================##
    P1         = p_from_pratio(Po1,gamma,M1)
    fanning, Re = fanning_and_reynolds(Po1,To,gamma,M1,Rs,Dpipe,mu,epsilon)



    ##==================================================================##
    #Check that PHI(M1) > 4fL/D and calculate PHI(M2) if possible
    ##==================================================================##
    fanno_constant = 4*fanning*L/Dpipe
    PHI1           = fanno_equation(M1,gamma)
    if PHI1 < fanno_constant:
        sys.exit("This pipe will choke before the next flow device")
    else:
        PHI2 = PHI1 - fanno_constant

    ##==================================================================##
    #Calculate Po2/Po1 using 2.3
    #Return Po2 & M2
    ##==================================================================##
    Lstar1   = Lstar_fanno(fanning,Dpipe,M1,gamma)
    L_int   = Lstar1 - L
    M2      = mach_fanno(L_int,fanning,Dpipe,gamma)
    Poratf  = fanno_po_ratio(M1,gamma)
    Postar  = Po1/Poratf
    Po2     = Postar*fanno_po_ratio(M2,gamma)
    P2      = p_from_pratio(Po2,gamma,M2)

    return P1, Po1, M1, Lstar1, P2, Po2, M2, Re


def valve_losses(P1,Cv,SG,Q,mdot,Rs,To,gamma,Apipe):
    #P2 = bisect(flowrates, 0, P1,args=(P1,Cv,SG,Q))
    P2 = bisect(flowrates,0,P1,args=(P1,Cv,SG,Q))
    M_aval  = bisect(delta_mass_static,0.0001,0.99,args=(mdot,P2*101325/14.7,Rs,To,gamma,Apipe))
    Po_aval = P2/(1+((gamma-1)/2)*M_aval**2)**(-(gamma)/(gamma-1))
    return P2,M_aval,Po_aval

def flowrates(P2,P1,Cv,SG,Q):
    '''
    Calculates the static pressure drop through a flow device rated by Cv
    Expected inputs:
    P1 and P2: Pressures upstream and downstream, PSI
    Cv       : Flow coefficient
    SG       : Specific gravity w.r.t. air
    Q        : Volumetric flow rate, SCFH (see mdot_to_scfh)
    '''
    return 42.2*Cv*np.sqrt((P1-P2)*(P1+P2))/np.sqrt(SG) - Q ##From the Deltrol catalog, equation relating volume flow rate Q to

def area_from_mass(Po,To,Rs,gamma,mdot):
    '''
    Function calculates the choking area using the compressible area ratio
    Expected inputs:
    Po       : Stagnation pressure, Pa
    To       : Stagnation temperature, K
    Rs       : Specific gas constant, J/kgK (double check units)
    gamma    : Ratio of specific heats
    mdot     : Mass flow rate, kg/s
    '''
    Gstar = Po*np.sqrt(gamma/Rs/To)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))##We call Gstar the ratio mdot/Astar
    Astar = mdot/Gstar
    return Astar

def mass_from_area(M,Po,To,Rs,gamma,Area):
    '''
    Function calculates the mass flow rate using the compressible area ratio
    Expected inputs:
    M        : Mach number
    Po       : Stagnation pressure, Pa
    To       : Stagnation temperature, K
    Rs       : Specific gas constant, J/kgK (double check units)
    gamma    : Ratio of specific heats
    Area     : Pipe area, sq. m
    '''
    #Astar = np.pi*Dpipe*Dpipe/4
    #Gstar = Po*np.sqrt(gamma/Rs/To)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))##We call Gstar the ratio mdot/Astar
    Gstar = Po*np.sqrt(gamma/Rs/To)*M*(1+(gamma-2)/2*M*M)**(-(gamma+1)/(2*(gamma-1)))##We call Gstar the ratio mdot/Astar
    mdot  = Gstar*Area
    return mdot

def prat_from_mach(gamma,M):
    '''
    For a known pre-shock mach number, what is the stagnation pressure ratio
    Expected inputs:
    M        : Mach number
    gamma    : Ratio of specific heats
    '''
    pratio = (((gamma+1)*M*M)/((gamma-1)*M*M+2))**(gamma/(gamma-1))*((gamma+1)/(2*gamma*M*M-(gamma-1)))**(1/(gamma-1))
    return pratio

def mach_from_pressure_ratio(Po1,Po2,gamma):
    '''
    For a desired stagnation pressure ratio, this function calculates the Mach number before a NSW
        Expected inputs:
        Po1      : Stagnation pressure before a normal shock wave, units same as Po2
        Po2      : Stagnation pressure after a normal shock wave, units same as Po1
        gamma    : Ratio of specific heats
    '''
    Por  = Po2/Po1 ##Desired pressure ratio
    def Prat(Mi,gamma,Por):
        return (((gamma+1)*Mi*Mi)/((gamma-1)*Mi*Mi+2))**(gamma/(gamma-1))*((gamma+1)/(2*gamma*Mi*Mi-(gamma-1)))**(1/(gamma-1)) - Por
    M = bisect(Prat,1,100,args=(gamma,Por))
    return M

def mach_after_shock(M1,gamma):
    '''
    Calculates the Mach number after a NSW knowing the pre-shock Mach number
    Expected inputs:
    M1       : Mach number
    gamma    : Ratio of specific heats
    '''
    M2 = np.sqrt(((gamma-1)*M1*M1+2)/(2*gamma*M1*M1-(gamma-1)))
    return M2

def pstatic_after_shock(M,gamma,P):
    '''
    Calculates the static pressure after a NSW knowing the pre-choke Mach number, gamma, and static pressure
    Expected inputs:
    M        : Mach number
    gamma    : Ratio of specific heats
    P        : Static pressure before a normal shock wave
    '''
    P2 = P*(2*gamma*M*M-(gamma-1))/(gamma+1)
    return P2

def pstag_after_shock(M,gamma,Po1):
    '''
    Calculates the stagnation pressure after a NSW knowing the pre-choke Mach number, gamma
    Expected inputs:
    M        : Mach number
    gamma    : Ratio of specific heats
    Po1      : Stagnation pressure before a normal shock wave
    '''
    Po2 = Po1*(((gamma+1)*M*M)/((gamma-1)*M*M+2))**(gamma/(gamma-1))*((gamma+1)/(2*gamma*M*M-(gamma-1)))**(1/(gamma-1)) ##Solution from the NSW stagnation pressure ratio equation.
    return Po2

def astar_all_else_known(Dpipe,M,gamma):
    '''
    Function calculates the choking area using the compressible area ratio knowing all other properties
    Expected inputs:
    Dpipe    : Diameter of pipe, meters
    M        : Mach number
    gamma    : Ratio of specific heats

    '''
    Apipe  = np.pi*Dpipe**2/4
    Aratio = aratio_from_mach(M,gamma)
    Astar  = Apipe/Aratio
    Dstar  = np.sqrt(Astar*4/np.pi)
    return Astar, Dstar

def mach_from_G(Po,Rs,To,gamma,mdot,Dpipe,subsuper):
    '''
    Calculates the Mach number knowing all other flow properties. This function allows the user to specify whether to resolve to the subsonic or supersonic branch
    Expected inputs:
    Po       : Stagnation pressure, Pa
    Rs       : Specific gas constant, J/kgK (double check units)
    To       : Stagnation temperature, K
    gamma    : Ratio of specific heats
    mdot     : Mass flow rate, kg/s
    Dpipe    : Diameter of pipe, meters
    subsuper : Specify either 'subsonic' or 'supersonic'
    '''
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

def mach_from_aratio(Apipe,Astar,gamma,subsuper):
    '''
    Function calculates the Mach number at a given location using the compressible area ratio knowing all other properties
    Expected inputs:
    Apipe    : Pipe area, sq. m
    Astar    : Choking area, sq. m
    gamma    : Ratio of specific heats
    subsuper : Specify either 'subsonic' or 'supersonic'
    '''
    def arat_delta(M,gamma,Apipe,Astar):
        return Apipe/Astar - aratio_from_mach(M,gamma)
    if subsuper == 'subsonic':
        M = bisect(arat_delta,0.00001,0.99,args=(gamma,Apipe,Astar))
    elif subsuper == 'supersonic':
        M = bisect(arat_delta,1,99,args=(gamma,Apipe,Astar))
    else:
        sys.exit('Please specify whether you want to resolve to the "subsonic" or "supersonic" branch when calling mach_from_aratio')
    return M

def mach_from_massflow(Apipe,mdot,Po,To,Rs,gamma,subsuper):
    '''
    Function calculates the Mach number at a given location using the compressible area ratio knowing all other properties
    Expected inputs:
    Apipe    : Pipe area, sq. m
    mdot     : Mass flow rate, kg/s
    Po       : Stagnation pressure, Pa
    To       : Stagnation temperature, K
    Rs       : Specific gas constant, J/kgK (double check units)
    gamma    : Ratio of specific heats
    subsuper : Specify either 'subsonic' or 'supersonic'
    '''
    Dpipe   = np.sqrt(4*Apipe/np.pi)
    def g_delta(M,Po,To,Rs,gamma,Dpipe,Apipe,mdot):
        return mdot/Apipe - mass_from_area(M,Po,To,Rs,gamma,Dpipe)
    if subsuper == 'subsonic':
        M = bisect(g_delta,0.00001,0.99,args=(Po,To,Rs,gamma,Dpipe,Apipe,mdot))
    elif subsuper == 'supersonic':
        M = bisect(g_delta,1,99,args=(Po,To,Rs,gamma,Dpipe,Apipe,mdot))
    else:
        sys.exit('Please specify whether you want to resolve to the "subsonic" or "supersonic" branch when calling mach_from_massflow')

    return M

def aratio_from_mach(M,gamma):
    '''
    Function calculates the compressible area ratio A/Astar knowing the Mach number and gamma
    Expected inputs:
    M        : Mach number
    gamma    : Ratio of specific heats
    '''
    Aratio = ((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))*(1+(gamma-1)/2*M*M)**((gamma+1)/(2*(gamma-1)))/M
    return Aratio

def p_from_pratio(Po,gamma,M):
    '''
    Function calculates the static pressure knowing the gas properties, Mach number, and stagnation pressure using the isentropic pressure ratio equation P/Po
    Expected inputs:
    Po       : Stagnation pressure, any units can be used. Static pressure will be returned in the same units provided for stagnation pressure
    gamma    : Ratio of specific heats
    M        : Mach number
    '''
    P_static = Po*(1+((gamma-1)/2)*M**2)**(-(gamma)/(gamma-1))
    return P_static

def T_from_Tratio(To,gamma,M):
    '''
    Function calculates the static pressure knowing the gas properties, Mach number, and stagnation pressure using the isentropic pressure ratio equation P/Po
    Expected inputs:
    To       : Stagnation temperature, K
    gamma    : Ratio of specific heats
    M        : Mach number

    '''
    T_static = To/(1+((gamma-1)/2)*M**2)
    return T_static

def mdot_to_scfh(mdot,Rs,G):
    '''
    Function calculates the volumetric flow rate in standard cubic feet per hour (scfh) of Nitrogen.
    Expected inputs:
    mdot     : Mass flow rate, in g/s
    Rs       : Specific gas constant, J/kgK (double check units)
    G        : Specific gravity of the studied fluid
    '''
    mdot = mdot*60*2.205/1000 #g/s to lbm/min
    rho  = 101325/(Rs*288.7)*2.205/(3.28084**3) #calculate the stp density, converting to lb/cu.ft
    Q    = mdot/rho
    scfh = Q/np.sqrt(1/G)*60
    return scfh

def hole_numbers(Dhole,Astar):
    '''
    Function calculates the number of holes on an injector knowing the choking area Astar and the drill diameter.
    Expected inputs:
    Dhole    : Hole diameter. Units compatible with Astar
    Astar    : Choking area. Units compatible with Dhole
    '''
    numholes = 4*Astar/np.pi/Dhole/Dhole
    return numholes

def colebrook_white(f,Re,D,epsilon):
    '''
    Subtracts both sides of the Colebrook-White equation to calculate the Darcy friction factor.
    Divide the result by 4 for the Fanning friction factor.
    Expected inputs:
    f       : Darcy friction factor
    Re      : Reynolds Number
    D       : Pipe diameter
    epsilon : Surface roughness in micrometers
    '''
    return 1/np.sqrt(f) - (-2)*np.log10(epsilon/(3.7*D) + 2.51/(Re*np.sqrt(f)))

def fanno_equation(M,gamma):
    return ((1-M**2)/(gamma*M**2) + (gamma+1)/(2*gamma)*np.log(((gamma+1)*M**2)/(2*(1+(gamma-1)/2*M**2))))


def delta_fanno(M,L,f,D,gamma):
    '''
    Function returns the sum of the left hand side and right hand side of the Fanno equation.
    Expected inputs:
    M       : Inlet Mach number
    L       : Choking pipe length
    f       : Fanning friction factor
    D       : Pipe diameter
    gamma   : Ratio of specific heats
    '''
    return ((1-M**2)/(gamma*M**2) + (gamma+1)/(2*gamma)*np.log(((gamma+1)*M**2)/(2*(1+(gamma-1)/2*M**2))))-4*f*L/D



def Lstar_fanno(f,D,M,gamma): #Define the Fanno equation to iterate on
    '''
    Function directly calculates Lstar in the Fanno equation.
    Expected inputs:
    f       : Fanning friction factor
    D       : Pipe diameter
    M       : Inlet Mach number
    gamma   : Ratio of specific heats
    '''

    return ((1-M**2)/(gamma*M**2) + (gamma+1)/(2*gamma)*np.log(((gamma+1)*M**2)/(2*(1+(gamma-1)/2*M**2))))*D/(4*f)

def mach_fanno(L,f,D,gamma): #Define the Fanno equation to iterate on
    '''
    Wraps the delta_fanno function to calculate a Mach number
    Expected inputs:
    L       : Choking pipe length
    f       : Fanning friction factor
    D       : Pipe diameter
    gamma   : Ratio of specific heats
    '''
    M = bisect(delta_fanno,0.001,0.99,args=(L,f,D,gamma))
    return M

def delta_mass_static(M,mdot,P,Rs,To,gamma,A):
    '''
    Using the mdot over astar equation for choked flow combined with the P/Po equation, provides an equation to iterate on knowing all other parameters.
    Expected inputs:
    M        : Mach number of the flow
    mdot     : Mass flow rate, g/s
    P        : Static pressure, Pa
    Rs       : Specific gas constant, J/kgK (double check units)
    To       : Stagnation temperature, K
    gamma    : Ratio of specific heats
    A        : Cross-sectional area of the pipe in sq. m
    '''

    return mdot/1000 - P*(1+(gamma-1)/2*M*M)**(gamma/(gamma-1))*A*np.sqrt(gamma/(Rs*To))*M*(1+(gamma-1)/2*M*M)**(-(gamma+1)/(2*(gamma-1)))

def delta_mass_stag(M,mdot,Po,Rs,To,gamma,A):
    '''
    Using the mdot over astar equation for choked flow, provides an equation to iterate on knowing all other parameters.
    Expected inputs:
    M        : Mach number of the flow
    mdot     : Mass flow rate, g/s
    Po       : Stagnation pressure, Pa
    Rs       : Specific gas constant, J/kgK (double check units)
    To       : Stagnation temperature, K
    gamma    : Ratio of specific heats
    A        : Cross-sectional area of the pipe in sq. m
    '''
    return mdot/1000 - A*Po*np.sqrt(gamma/(Rs*To))*M*(1+(gamma-1)/2*M*M)**((-(gamma+1))/(2*(gamma-1)))

def fanno_po_ratio(M,gamma):
    return (1/M)*((2+(gamma-1)*M**2)/(gamma+1))**((gamma+1)/(2*(gamma-1)))
