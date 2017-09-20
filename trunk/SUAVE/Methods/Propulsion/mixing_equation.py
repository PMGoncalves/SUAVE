## @ingroup Methods-Propulsion
# fm_solver.py
# 
# Created:  Sep 2017, P Goncalves

import numpy as np

from scipy.optimize import fsolve

# ----------------------------------------------------------------------
#  fm_solver
# ----------------------------------------------------------------------

## @ingroup Methods-Propulsion


def mixing_equation(psi3,gamma):
    """
    Function that takes in an area ratio and a Mach number associated to
    one of the areas and outputs the missing Mach number
    
    Inputs:
    M       [dimensionless]
    gamma   [dimensionless]
    Aratio  [dimensionless]
    
    Outputs:
    M1      [dimensionless]
    
    """
    func = lambda M_out: (M_out**2*(1+(gamma-1)/2*M_out**2)/(1+gamma*M_out**2)**2) -psi3
    
    #Initializing the array
    M_guess = 0.001*psi3/psi3
    
    # Creating index
    i_psi = psi3 < 1.0
    
    #--Supersonic solution
    M_guess[i_psi]= 0.001
 
    M_out= fsolve(func,M_guess)
    
    return M_out
