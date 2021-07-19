# ChokedInjectors
Here is a list of functions available in CompressibleFlowFunctions.py:

area_from_mass(Po,To,Rs,gamma,mdot): 

    ##Function calculates the choking area Astar using the compressible area ratio

mass_from_area(Po,To,Rs,gamma,Dpipe): 

    ##Function calculates the mass flow rate using the compressible area ratio m_dot/Astar

mach_from_pressure_ratio(Po1,Po2,gamma):
    ##For a desired stagnation pressure ratio, this function calculates the Mach number before a NSW using the normal shock equations

mach_after_shock(M1,gamma):

    ##Calculates the Mach number after a NSW knowing the pre-shock Mach number

pstatic_after_shock(M,gamma,P):

    ##Calculates the static pressure after a NSW knowing the pre-choke Mach number, gamma, and static pressure

pstag_after_shock(M,gamma,Po1):

    ##Calculates the stagnation pressure after a NSW knowing the pre-choke Mach number, gamma

astar_all_else_known(Dpipe,M,gamma):

    ##Function calculates the choking area using the compressible area ratio knowing all other properties
    
mach_from_aratio_subsonic(Aexit,Astar,gamma):

    ##Function calculates the Mach number at a given location using the compressible area ratio A/Astar knowing all other properties (areas, gamma)

mach_from_aratio_supersonic(Aexit,Astar,gamma):

    ##Function calculates the Mach number at a given location using the compressible area ratio knowing all other properties

aratio_from_mach(M,gamma):

    ##Function calculates the compressible area ratio A/Astar knowing the Mach number and gamma

p_from_pratio(Po,gamma,M):

    ##Function calculates the static pressure knowing the gas properties, Mach number, and stagnation pressure using the isentropic pressure ratio equation P/Po

hole_numbers(Dhole,Astar):

    ## Not part of the compressible flow solver per se, this function is used in injector design. 
    ## Once the choking area at the injector is known and using a specific hole diameter (drill bit size, etc.), this function spits out the number of holes required to choke.
