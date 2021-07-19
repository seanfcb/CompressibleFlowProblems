import numpy as np
###All functions take as an input: pressure in PSI, Temperature in Kelvin, Pipe diameters in inches
###All functions output answers in SI units

def area_from_mass(Po,To,Rs,gamma,mdot):
    ##Function calculates the choking area using the compressible area ratio
    Po = Po*101325/14.7
    Gstar = Po*np.sqrt(gamma/Rs/To)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))##We call Gstar the ratio mdot/Astar
    Astar = mdot/Gstar
    return Astar


def mass_from_area(Po,To,Rs,gamma,Dpipe):
    ##Function calculates the mass flow rate using the compressible area ratio
    Gstar = Po*np.sqrt(gamma/Rs/To)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))##We call Gstar the ratio mdot/Astar
    mdot = Gstar*Astar
    return mdot

def mach_from_pressure_ratio(Po1,Po2,gamma):
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
    Dpipe  = Dpipe*0.0254
    Apipe  = np.pi*Dpipe**2/4
    Aratio = ((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))*(1+(gamma-1)/2*M*M)**((gamma+1)/(2*(gamma-1)))/M
    Astar  = Apipe/Aratio
    Dstar  = np.sqrt(Astar*4/np.pi)
    return Astar, Dstar

def mach_from_aratio_subsonic(Aexit,Astar,gamma):
    ##Function calculates the Mach number at a given location using the compressible area ratio knowing all other properties
    tol     = 1e-9 #Tolerance for convergence method
    Mguess  = 0.99
    Apipe   = Aexit
    Aratio1 = Apipe/Astar
    Aratio2 = ((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))*(1+(gamma-1)/2*Mguess*Mguess)**((gamma+1)/(2*(gamma-1)))/Mguess
    zero    = Aratio1 - Aratio2
    a = 0
    b = Mguess
    M = (a + b)/2
    while abs(zero) > tol:
        Aratio = ((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))*(1+(gamma-1)/2*M*M)**((gamma+1)/(2*(gamma-1)))/M
        zero   = Aratio1 - Aratio
        if zero < 0:
            a = M
            M = (a+b)/2
        elif zero > 0:
            b = M
            M = (a+b)/2
    return M

def mach_from_aratio_supersonic(Aexit,Astar,gamma):
    ##Function calculates the Mach number at a given location using the compressible area ratio knowing all other properties
    tol     = 1e-9 #Tolerance for convergence method
    Mguess  = 100
    Aratio1 = Aexit/Astar
    Aratio2 = ((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))*(1+(gamma-1)/2*Mguess*Mguess)**((gamma+1)/(2*(gamma-1)))/Mguess
    zero    = Aratio1 - Aratio2
    a = 1
    b = Mguess
    M = Mguess
    while abs(zero) > tol:
        Aratio = ((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))*(1+(gamma-1)/2*M*M)**((gamma+1)/(2*(gamma-1)))/M
        zero   = Aratio1 - Aratio
        if zero > 0:
            a = M
            M = (a+b)/2
        elif zero < 0:
            b = M
            M = (a+b)/2
    return M

def p_from_pratio(Po,gamma,M):
    ##Function calculates the static pressure knowing the gas properties, Mach number, and stagnation pressure using the isentropic pressure ratio equation P/Po
    Po = Po*101325/14.7
    P_static = Po*(1+((gamma-1)/2)*M**2)**(-(gamma)/(gamma-1))
    return P_static

def hole_numbers(Dhole,Astar):
    Dhole = Dhole*0.0254
    numholes = 4*Astar/np.pi/Dhole/Dhole
    return numholes
