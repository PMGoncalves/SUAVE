## @ingroup Components-Energy-Converters
# Compressor.py
#
# Created:  Jul 2014, A. Variyar
# Modified: Jan 2016, T. MacDonald

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# SUAVE imports

import SUAVE

from SUAVE.Core import Units

# package imports
import numpy as np
from scipy.optimize import fsolve

from SUAVE.Components.Energy.Energy_Component import Energy_Component
from SUAVE.Methods.Propulsion.mixing_equation import mixing_equation


# ----------------------------------------------------------------------
#  Compressor Component
# ----------------------------------------------------------------------
## @ingroup Components-Energy-Converters
class Mixer(Energy_Component):
    """This is a mixer component, which takes in two streams with their own set
    of properties and computes the outcome properties of the mixing

    Source:
    NACA, "Effects oof parallel-jet mixing on downstream Mach number and stagnation
    pressure with application to engine testing in supersonic tunnels",
    Harry Bernstein
    """
    
    def __defaults__(self):
        """This sets the default values for the component to function.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        None
        """          
        #set the default values
        self.tag                             = 'Compressor'
        self.polytropic_efficiency           = 1.0
        self.pressure_ratio                  = 1.0
        self.inputs.stagnation_temperature   = 0.
        self.inputs.stagnation_pressure      = 0.
        self.outputs.stagnation_temperature  = 0.
        self.outputs.stagnation_pressure     = 0.
        self.outputs.stagnation_enthalpy     = 0.
    

    def compute(self,conditions):
        """ This computes the output values from the input values according to
        equations from the source.

        Assumptions:
        Constant polytropic efficiency and pressure ratio

        Source:
        https://web.stanford.edu/~cantwell/AA283_Course_Material/AA283_Course_Notes/

        Inputs:
        conditions.freestream.
          isentropic_expansion_factor         [-]
          specific_heat_at_constant_pressure  [J/(kg K)]
        self.inputs.
          stagnation_temperature              [K]
          stagnation_pressure                 [Pa]

        Outputs:
        self.outputs.
          stagnation_temperature              [K]  
          stagnation_pressure                 [Pa]
          stagnation_enthalpy                 [J/kg]
          work_done                           [J/kg]

        Properties Used:
        self.
          pressure_ratio                      [-]
          polytropic_efficiency               [-]
        """          
        #unpack the values
        
        #unpack from conditions
        gamma    = conditions.freestream.isentropic_expansion_factor
        
        #unpack from inputs
        Tt_in    = self.inputs.stagnation_temperature
        Pt_in    = self.inputs.stagnation_pressure
        M_in     = self.inputs.mach
        Tt_in_2  = self.inputs.stagnation_temperature_2
        Pt_in_2  = self.inputs.stagnation_pressure_2
#        M_in_2   = self.inputs.mach_2
        
        #unpack from self
        A3_A1    = self.inputs.area_ratio
        
        
        # Compute static properties coming from component 1
        P_in     = ((1+(gamma-1)/2*M_in**2)**(-gamma/(gamma-1)))*Pt_in
        T_in     = ((1+(gamma-1)/2*M_in**2)**-1)*Tt_in
        
        # Compute static properties coming from component 2
        
        # If M_in_2 is specified
#        P_in_2   = ((1+(gamma-1)/2*M_in_2**2)**(-gamma/(gamma-1)))*Pt_in_2
#        T_in_2   = (1+(gamma-1)/2*M_in_2**2)*Tt_in_2
#        
        P_in_2 = 1.0*Pt_in/Pt_in
        i_static = P_in <= Pt_in_2
        i_total  = P_in > Pt_in_2
        
        P_in_2[i_static] = P_in[i_static]
        P_in_2[i_total] = Pt_in_2[i_total]
        
        M_in_2 = np.sqrt((((Pt_in_2/P_in_2)**((gamma-1)/gamma))-1)*2/(gamma-1))
        T_in_2   = ((1+(gamma-1)/2*M_in_2**2)**-1)*Tt_in_2

        # Compute intermediate variables for mixing equation
        A2_A1    = A3_A1 - 1
        alpha    = A2_A1*P_in_2/P_in
        ksi      = M_in_2/M_in
        theta    = T_in_2/T_in   
        X        = (1 + gamma*ksi**2*M_in**2)/(1+gamma*M_in**2)
        Y        = (1+(gamma-1)/2*ksi**2*M_in**2)/(1+(gamma-1)/2*M_in**2)   
        psi      = M_in**2*(1+(gamma-1)/2*M_in**2)/(1+gamma*M_in**2)**2  
        psi3     = psi*(((1+alpha*ksi*Y*np.sqrt(theta))*(1+(alpha*ksi/np.sqrt(theta))))/(1+alpha*X)**2)
        
        q = 7
#        print 'Pressure: 1 ', P_in[q], 'and 2', P_in_2[q]
#        print 'Temperature: 1', T_in[q], 'and 2', T_in_2[q]
#        print 'Mach : 1', M_in[q], 'and 2', M_in_2[q]
#        print 'alpha', alpha[q]
#        print 'ksi', ksi[q]
#        print 'theta', theta[q]
#        print 'y', Y[q]
#        print 'x', X[q]
#        print 'psi', psi[q]
#        print 'psi' , psi3[q]
        # Compute M_out
            #-- Initialize arrays
        M_out  = 1*psi3/psi3
            #--Make i_psi the size of output arrays
        i_psi = psi3 < 10.0
            #-- Mixing equation
        M_out[i_psi]    = mixing_equation(psi3[i_psi],gamma)  
   
        # Compute output properties
        #-- Momentum conservation
        P_out   = (P_in*(1+gamma*M_in**2)+P_in_2*(1+gamma*M_in_2**2)*A2_A1)/((1+A2_A1)*(1+gamma*M_out**2))
        #-- Continuity equation
        T_out   = ((np.sqrt(T_in*T_in_2)*M_out*P_out*(1+A2_A1))/(P_in*M_in*np.sqrt(T_in_2) + (P_in_2*M_in_2*np.sqrt(T_in)*A2_A1)))**2
        Pt_out  = P_out*(1+(gamma-1)/2*M_out**2)**(gamma/(gamma-1))
        Tt_out  = T_out*(1+(gamma-1)/2*M_out**2)
        
        #pack computed quantities into the outputs
        self.outputs.stagnation_temperature  = Tt_out
        self.outputs.stagnation_pressure     = Pt_out
        self.outputs.mach_number             = M_out
    
    
    __call__ = compute
