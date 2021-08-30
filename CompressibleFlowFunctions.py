import numpy as np
import CoolProp.CoolProp as CP
###All functions take as an input: pressure in PSI, Temperature in Kelvin, Pipe diameters in inches
###All functions output answers in SI units
###In all functions, fluid is a string that refers to a name available in CoolProp

def area_from_mass(Po,To,Rs,gamma,mdot, fluid):
    ##Function calculates the choking area using the compressible area ratio
    Po = Po*101325/14.7
    Gstar = Po*np.sqrt(gamma/Rs/To)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))##We call Gstar the ratio mdot/Astar
    Astar = mdot/Gstar
    return Astar


def mass_from_area(Po,To,Rs,gamma,Dpipe, fluid):
    ##Function calculates the mass flow rate using the compressible area ratio
    Gstar = Po*np.sqrt(gamma/Rs/To)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))##We call Gstar the ratio mdot/Astar
    mdot = Gstar*Astar
    return mdot

def prat_from_mach(gamma,M):
    ##For a known pre-shock mach number, what is the stagnation pressure ratio
    pratio = (((gamma+1)*M*M)/((gamma-1)*M*M+2))**(gamma/(gamma-1))*((gamma+1)/(2*gamma*M*M-(gamma-1)))**(1/(gamma-1))
    return pratio

def mach_from_pressure_ratio(Po1,Po2,gamma, fluid):
    ##For a desired stagnation pressure ratio, this function calculates the Mach number before a NSW
    tol  = 1e-9 #Tolerance for convergence method
    Mi   = 20 #First guess for Mach number
    Por  = Po2/Po1 ##Desired pressure ratio
    Prat = (((gamma+1)*Mi*Mi)/((gamma-1)*Mi*Mi+2))**(gamma/(gamma-1))*((gamma+1)/(2*gamma*Mi*Mi-(gamma-1)))**(1/(gamma-1)) ##Solution from the NSW stagnation pressure ratio equation.
    zero = Por - Prat
    a = 1
    b = Mi
    M = (a+b)/2
    while abs(zero) > tol:
        Prat = (((gamma+1)*M*M)/((gamma-1)*M*M+2))**(gamma/(gamma-1))*((gamma+1)/(2*gamma*M*M-(gamma-1)))**(1/(gamma-1)) ##Solution from the NSW stagnation pressure ratio equation.
        zero = Por-Prat
        if zero<0:
            a = M
            M = (a+b)/2
        elif zero > 0:
            b = M
            M = (a+b)/2
    return M

def mach_after_shock(M1,gamma, fluid):
    ##Calculates the Mach number after a NSW knowing the pre-shock Mach number
    M2 = np.sqrt(((gamma-1)*M1*M1+2)/(2*gamma*M1*M1-(gamma-1)))
    return M2

def pstatic_after_shock(M,gamma,P, fluid):
    ##Calculates the static pressure after a NSW knowing the pre-choke Mach number, gamma, and static pressure
    P2 = P*(2*gamma*M*M-(gamma-1))/(gamma+1)
    return P2

def pstag_after_shock(M,gamma,Po1, fluid):
    ##Calculates the stagnation pressure after a NSW knowing the pre-choke Mach number, gamma
    Po2 = Po1*(((gamma+1)*M*M)/((gamma-1)*M*M+2))**(gamma/(gamma-1))*((gamma+1)/(2*gamma*M*M-(gamma-1)))**(1/(gamma-1)) ##Solution from the NSW stagnation pressure ratio equation.
    return Po2

def astar_all_else_known(Dpipe,M,gamma, fluid):
    ##Function calculates the choking area using the compressible area ratio knowing all other properties
    Dpipe  = Dpipe*0.0254
    Apipe  = np.pi*Dpipe**2/4
    Aratio = aratio_from_mach(M,gamma)
    Astar  = Apipe/Aratio
    Dstar  = np.sqrt(Astar*4/np.pi)
    return Astar, Dstar

def mach_from_aratio_subsonic(Aexit,Astar,gamma, fluid):
    ##Function calculates the Mach number at a given location using the compressible area ratio knowing all other properties
    tol     = 1e-9 #Tolerance for convergence method
    Mguess  = 0.99
    Apipe   = Aexit
    Aratio1 = Apipe/Astar
    Aratio2 = aratio_from_mach(Mguess,gamma)
    zero    = Aratio1 - Aratio2
    a = 0
    b = Mguess
    M = (a + b)/2
    while abs(zero) > tol:
        Aratio = aratio_from_mach(M,gamma)
        zero   = Aratio1 - Aratio
        if zero < 0:
            a = M
            M = (a+b)/2
        elif zero > 0:
            b = M
            M = (a+b)/2
    return M

def mach_from_massflow_subsonic(Aexit,mdot,Po,To,Rs,gamma):
    ##Function calculates the Mach number at a given location using the compressible area ratio knowing all other properties
    tol     = 1e-9 #Tolerance for convergence method
    Mguess  = 0.99
    Apipe   = Aexit
    Dpipe   = np.sqrt(4*Apipe/np.pi)
    G1 = mdot/Apipe
    G2 = mass_from_area(Po,To,Rs,gamma,Dpipe)
    a = 0
    b = Mguess
    M = (a + b)/2
    zero = 1
    while abs(zero) > tol:
        Gcalc = mass_from_area(Po,To,Rs,gamma,Dpipe)
        zero   = G1 - Gcalc
        if zero < 0:
            a = M
            M = (a+b)/2
        elif zero > 0:
            b = M
            M = (a+b)/2
    return M
    
def mach_from_aratio_supersonic(Aexit,Astar,gamma, fluid):
    ##Function calculates the Mach number at a given location using the compressible area ratio knowing all other properties
    tol     = 1e-9 #Tolerance for convergence method
    Mguess  = 100
    Aratio1 = Aexit/Astar
    Aratio2 = aratio_from_mach(Mguess,gamma)
    zero    = Aratio1 - Aratio2
    a = 1
    b = Mguess
    M = Mguess
    while abs(zero) > tol:
        Aratio = aratio_from_mach(M,gamma)
        zero   = Aratio1 - Aratio
        if zero > 0:
            a = M
            M = (a+b)/2
        elif zero < 0:
            b = M
            M = (a+b)/2
    return M

def aratio_from_mach(M,gamma, fluid):
    #Function calculates the compressible area ratio A/Astar knowing the Mach number and gamma
    Aratio = ((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))*(1+(gamma-1)/2*M*M)**((gamma+1)/(2*(gamma-1)))/M
    return Aratio

def p_from_pratio(To,Po,gamma,M, fluid):
    ##Function calculates the static pressure knowing the gas properties, Mach number, and stagnation pressure using the isentropic pressure ratio equation P/Po
    Po = Po*101325/14.7
    ho = CP.PropsSI('H','P',Po,'T',To,fluid) #ho = h(Po,To)
    so = CP.PropsSI('S','P',Po,'T',To,fluid) #s = so = s(Po,To)
    Tsol = T_from_Tratio(To,Po,gamma,M,fluid)
    T_static = Tsol[0]
    P_static = CP.PropsSI('P','T',T_static,'S',so,fluid)
    return P_static

def T_from_Tratio(To,Po,gamma,M, fluid):
    ##Function calculates the static pressure knowing the gas properties, Mach number, and stagnation pressure using the isentropic pressure ratio equation P/Po
    Po = Po*101325/14.7
    ho = CP.PropsSI('H','P',Po,'T',To,fluid) #ho = h(Po,To)
    so = CP.PropsSI('S','P',Po,'T',To,fluid) #s = so = s(Po,To)
    Ttriple = CP.PropsSI('T_Triple',fluid)
    #fixed point iteration
    #initial guess, h = h0
    tol = 1e-9
    i = 0

    gamma_o = CP.PropsSI('isentropic_expansion_coefficient','H',ho,'S',so,fluid)
    T      = To / (1 + 0.5*(gamma_o-1)*M**2)
    h      = CP.PropsSI('H','T',T,'S',so,fluid)
    hold   = 1e3*h
    while abs((hold - h)/h) > tol:
      #print(i)
      hold = h
      c  = CP.PropsSI('speed_of_sound','H',h,'S',so,fluid)
      V  = M*c
      h = ho - (M*c)**2/2
      print(hold,h)
      i = i + 1
        
    T_static = CP.PropsSI('T','H',h,'S',so,fluid)
    phase = CP.PhaseSI('H',h,'S',so,fluid)
    if phase[-3:] != 'gas':
      print('WARNING: Fluid is no longer in a gas phase')
    return (T_static, phase)


def hole_numbers(Dhole,Astar):
    Dhole = Dhole*0.0254
    numholes = 4*Astar/np.pi/Dhole/Dhole
    return numholes
