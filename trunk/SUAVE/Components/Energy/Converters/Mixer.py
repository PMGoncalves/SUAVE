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

from SUAVE.Components.Energy.Energy_Component import Energy_Component

# ----------------------------------------------------------------------
#  Compressor Component
# ----------------------------------------------------------------------
## @ingroup Components-Energy-Converters
class Compressor(Energy_Component):
    """This is a compressor component typically used in a turbofan.
    Calling this class calls the compute function.
    
    Assumptions:
    Pressure ratio and efficiency do not change with varying conditions.

    Source:
    https://web.stanford.edu/~cantwell/AA283_Course_Material/AA283_Course_Notes/
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
        Cp       = conditions.freestream.specific_heat_at_constant_pressure
        R        = conditions.freestream.universal_gas_constant
        
        #unpack from inputs
        Tt_in    = self.inputs.stagnation_temperature_1
        Tt_in_2    = self.inputs.stagnation_temperature_2
        Pt_in    = self.inputs.stagnation_pressure_1
        Pt_in_2    = self.inputs.stagnation_pressure_2
        M_in     = self.inputs.mach_1
        
        # Compute static properties coming from one component
        P_in     = ((1+(gamma-1)/2*M_in**2)**(-gamma/(gamma-1)))*Pt_in
        T_in     = (1+(gamma-1)/2*M_in**2)*Tt_in
        rho_in   = P_in/(R*T_in)
        u_in     = np.sqrt(gamma*R*T_in)*M_in
        
        # Assume pressure must be the same for the mixing to occur
        P_in_2     = P_in
        
        # Assume small variable nozzle to allow for pressure match
        M_in_2     = np.sqrt((((Pt_in_2/P_in_2)**((gamma-1)/gamma))-1)*2/(gamma-1))
        T_in_2     = (1+(gamma-1)/2*M_in_2**2)*Tt_in_2
        rho_in_2   = P_in_2/(R*T_in_2)*M_in_2
        u_in_2    = np.sqrt(gamma*R*T_in_2)*M_in_2
        
        # Apply mass conservation 
        

        
        
        
        #unpack from self
        pid      = self.pressure_ratio
        etapold  = self.polytropic_efficiency
        
        #
        
        #Method to compute compressor properties
        
        #Compute the output stagnation quantities based on the pressure ratio of the component
        ht_in     = Cp*Tt_in
        Pt_out    = Pt_in*pid
        Tt_out    = Tt_in*pid**((gamma-1)/(gamma*etapold))
        ht_out    = Cp*Tt_out
        
        #compute the work done by the compressor(for matching with the turbine)
        work_done = ht_out- ht_in
        
        #pack computed quantities into the outputs
        self.outputs.stagnation_temperature  = Tt_out
        self.outputs.stagnation_pressure     = Pt_out
        self.outputs.stagnation_enthalpy     = ht_out
        self.outputs.work_done               = work_done
    
    
    __call__ = compute
