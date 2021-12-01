import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import bisect
from CompressibleFlowFunctions.algos import *
from CompressibleFlowFunctions.Isentropic import *
from CompressibleFlowFunctions.Fanno import *
from CompressibleFlowFunctions.NSW import *
from CompressibleFlowFunctions.misc import *


from fannocvloss_algo import * ##### This script runs the input file fannocvloss_algo.py
from tabulate import tabulate
import sys
mdot = 957
Q = mdot_to_scfh(mdot,Rs,SG)

def delmass(Po,M,mdot,Rs,To,gamma,A):
    return delta_mass_stag(M,mdot,Po,Rs,To,gamma,A)

def find_M_bval(M1,P2,Cv,mdot,Rs,To,gamma,Apipe,Q1):
    print(M1)
    Po1  = newton(delmass,500*101325/14.7,args=(M1,mdot,Rs,To,gamma,Apipe))
    P1   = p_from_pratio(Po1,gamma,M1)
    T1   = T_from_Tratio(To,gamma,M1)
    Q2   = flowrates_swagelok(T1,P2,P1,Cv,SG)
    print(Q2)
    delQ = Q1-Q2
    return delQ
# def find_P_bval(P1,P2,Cv,SG,Q,mdot,Rs,To,gamma,Apipe):
#     T1   = bisect(flowrates_swagelok,1,500,args=(P2,P1,Cv_nvalv,SG,Q))#newton(flowrates_swagelok,To-10,args=(P2,P1,Cv,SG,Q))
#     M1   = bisect(mach_from_Tratio,0.01,0.99,args=(To,T1,gamma))
#     Po1a = po_from_pratio(P1,gamma,M1)
#     Po1b = newton(delmass,Po1a*101325/14.7,args=(M1,mdot,Rs,To,gamma,Apipe))
#     return Po1a - Po1b
#
# P_bval = newton(find_P_bval,P2+10,args=(P2,Cv,SG,Q,mdot,Rs,To,gamma,Apipe))


def valve_losses_backwards(P2,Cv,SG,Q,mdot,Rs,To,gamma,Apipe): ##Based on the Swagelok equation
    '''
    Does this display anything

    '''

    return P_bval#, Po_bval, M_bval

def flowrates_swagelok(T1,P2,P1,Cv,SG):
    '''
    Calculates the static pressure drop through a flow device rated by Cv
    Expected inputs:
    P1 and P2: Pressures upstream and downstream, PSI
    Cv       : Flow coefficient
    SG       : Specific gravity w.r.t. air
    Q        : Volumetric flow rate, SCFH (see mdot_to_scfh)
    '''
    conv = 60*22.67*np.sqrt(5/9)#conversion constant
    delP = P1-P2
    Q = conv*Cv*P1*(1-(2/3)*(delP/P1))*np.sqrt(delP/(P1*SG*T1))
    return Q

Dpipe      = Dopipe-2*PipeT
Dpipe      = Dpipe*0.0254 #(0.75-2*0.065)*0.0254
Apipe      = np.pi*Dpipe**2/4

M1 = bisect(find_M_bval,0.01,0.999,args=(635,Cv_nvalv,mdot,Rs,To,gamma,Apipe,Q))

#P_bval = valve_losses_backwards(635,Cv_nvalv,SG,Q,mdot,Rs,To,gamma,Apipe)
