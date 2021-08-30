import numpy as np
from CompressibleFlowFunctions import *
from DoubleChokeProblem import Astar_ox, Astar_fuel, Apipe

To         = 300
gamma_ox   = 1.4
gamma_fuel = 1.4
Rs_ox      = 8314.5/32
Rs_fuel    = 8314.5/2
Y_ox       = 8/9
Y_fuel     = 1/9
m_dot      = 200/1000
mdot_ox    = Y_ox*m_dot
mdot_fuel  = Y_fuel*m_dot


mdot_astar_ox   = mdot_ox/Astar_ox
mdot_astar_fuel = mdot_fuel/Astar_fuel
f_gamma_ox      = np.sqrt(gamma_ox/Rs_ox/To)*((gamma_ox+1)/2)**(-(gamma_ox+1)/(2*(gamma_ox-1)))
f_gamma_fuel    = np.sqrt(gamma_fuel/Rs_fuel/To)*((gamma_fuel+1)/2)**(-(gamma_fuel+1)/(2*(gamma_fuel-1)))

Po2_ox   = mdot_astar_ox/f_gamma_ox*14.7/101325
Po2_fuel = mdot_astar_fuel/f_gamma_fuel*14.7/101325

print('Po2_ox = %f psi'%Po2_ox)
print('Po2_fuel = %f psi'%Po2_fuel)


##What Po1 allows for this
M2_ox    = mach_from_aratio_subsonic(Apipe,Astar_ox,gamma_ox)
M2_fuel  = mach_from_aratio_subsonic(Apipe/4,Astar_fuel,gamma_fuel)

M1_ox    = np.sqrt((2+M2_ox*M2_ox*(gamma_ox-1))/(M2_ox*M2_ox*2*gamma_ox-(gamma_ox-1)))
M1_fuel  = np.sqrt((2+M2_fuel*M2_fuel*(gamma_fuel-1))/(M2_fuel*M2_fuel*2*gamma_fuel-(gamma_fuel-1)))

Po1_ox   = prat_from_mach(gamma_ox,M1_ox)
Po1_fuel = prat_from_mach(gamma_fuel,M1_fuel)

print('Po1_ox = %f psi'%Po1_ox)
print('Po1_fuel = %f psi'%Po1_fuel)
