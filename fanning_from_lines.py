import numpy as np
from scipy.optimize import bisect
from CompressibleFlowFunctions import *
from fannocvloss import *

mdot   = 288 ## g/s O2
Astarf = 21.338e-6 #21 square mm to sq. m
Dpipe  = (Dopipe - 2*PipeT)*0.0254
Apipe  = np.pi*Dpipe**2/4
Lstar  = 7/3.2808 #approximate conversion for feet to m is /3.2808


M_inlet  = bisect(delta_mass_stag,0.0001,0.99,args=(mdot,Po1_initial*101325/14.7,Rs,To,gamma,Apipe))
fanning  = ((1-M_inlet**2)/(gamma*M_inlet**2) + (gamma+1)/(2*gamma)*np.log(((gamma+1)*M_inlet**2)/(2*(1+(gamma-1)/2*M_inlet**2))))*Dpipe/(4*Lstar)
