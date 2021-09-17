import numpy as np
from CompressibleFlowFunctions import *


To         = 300
gamma_ox   = 1.4
gamma_fuel = 1.4
Rs_ox      = 8314.5/32
Rs_fuel    = 8314.5/2
# Y_ox       = 8/9
# Y_fuel     = 1/9
# m_dot      = 200/1000
# mdot_ox    = Y_ox*m_dot
# mdot_fuel  = Y_fuel*m_dot
A_fuel       = 62*np.pi*(1/32*0.0254)**2/4
A_ox     = 42*(0.03125*0.0254)*(0.0117*0.0254)

Pof_ox     = 781.3*101325/14.7
Pof_fuel   = 783.4*101325/14.7

mdot_ox    = mass_from_area(1,Pof_ox,To,Rs_ox,gamma_ox,A_ox)
mdot_fuel  = mass_from_area(1,Pof_fuel,To,Rs_fuel,gamma_fuel,A_fuel)
