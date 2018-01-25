## @ingroup Components-Energy-Converters
# Combustor.py
#
# Created:  Oct 2014, A. Variyar
# Modified: Jan 2016, T. MacDonald

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE
from SUAVE.Core import Data
from SUAVE.Components.Energy.Energy_Component import Energy_Component


# ----------------------------------------------------------------------
#  Combustor Component
# ----------------------------------------------------------------------
## @ingroup Components-Energy-Converters
class Combustor(Energy_Component):
    """This is provides output values for a combustor
    Calling this class calls the compute function.
    
    Assumptions:
    None
    
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
        
        self.tag = 'Combustor'
        
        #-----setting the default values for the different components
        self.fuel_data                      = SUAVE.Attributes.Propellants.Jet_A()
        self.alphac                         = 0.0
        self.turbine_inlet_temperature      = 1.0
        self.inputs.stagnation_temperature  = 1.0
        self.inputs.stagnation_pressure     = 1.0
        self.outputs.stagnation_temperature = 1.0
        self.outputs.stagnation_pressure    = 1.0
        self.outputs.stagnation_enthalpy    = 1.0
        self.outputs.fuel_to_air_ratio      = 1.0
        self.fuel_data                      = Data()
    
    
    
    def compute(self,conditions):
        """ This computes the output values from the input values according to
        equations from the source.

        Assumptions:
        Constant efficiency and pressure ratio

        Source:
        https://web.stanford.edu/~cantwell/AA283_Course_Material/AA283_Course_Notes/

        Inputs:
        conditions.freestream.
          isentropic_expansion_factor         [-]
          specific_heat_at_constant_pressure  [J/(kg K)]
          temperature                         [K]
          stagnation_temperature              [K]
        self.inputs.
          stagnation_temperature              [K]
          stagnation_pressure                 [Pa]

        Outputs:
        self.outputs.
          stagnation_temperature              [K]  
          stagnation_pressure                 [Pa]
          stagnation_enthalpy                 [J/kg]
          fuel_to_air_ratio                   [-]

        Properties Used:
        self.
          turbine_inlet_temperature           [K]
          pressure_ratio                      [-]
          efficiency                          [-]
        """         
        # unpack the values
        
        # unpacking the values from conditions
        gamma  = conditions.freestream.isentropic_expansion_factor 
        Cp     = conditions.freestream.specific_heat_at_constant_pressure
        To     = conditions.freestream.temperature
        Tto    = conditions.freestream.stagnation_temperature
        
        # unpacking the values form inputs
        Tt_in  = self.inputs.stagnation_temperature
        Pt_in  = self.inputs.stagnation_pressure
        Tt4    = self.turbine_inlet_temperature
        pib    = self.pressure_ratio
        eta_b  = self.efficiency
        
        # unpacking values from self
        htf    = self.fuel_data.specific_energy        

        # method to compute combustor properties

        # method - computing the stagnation enthalpies from stagnation temperatures
        ht4     = Cp*Tt4
        ho      = Cp*To
        ht_in   = Cp*Tt_in

        # Using the Turbine exit temperature, the fuel properties and freestream temperature to compute the fuel to air ratio f
        f       = (ht4 - ht_in)/(eta_b*htf-ht4)

        # Computing the exit static and stagnation conditions
        ht_out  = Cp*Tt4
        Pt_out  = Pt_in*pib
        
        # pack computed quantities into outputs
        self.outputs.stagnation_temperature  = Tt4
        self.outputs.stagnation_pressure     = Pt_out
        self.outputs.stagnation_enthalpy     = ht_out
        self.outputs.fuel_to_air_ratio       = f 
    
    
    def compute_rocket(self,conditions):
        """ This computes the output values from the input values according to
        equations from the source.

        Assumptions:

        Source:

        Inputs:

        """         
        # -- unpack the values
        
        # unpacking the values from conditions
        To       = conditions.freestream.temperature
        
        # unpacking values from self
#        fuel     = self.fuel_data
#        oxidizer = self.oxidizer_data
        scale    = self.scaling_factor
        Pt_comb  = self.stagnation_pressure
        OF       = self.oxidizer_fuel_ratio
        
        OF       = 6.0
        # method to compute combustor properties

        # Initialize arrays
        Tt_ad_flame  = 1.0 * To / To  # Adiabatic flame temperature
        gamma        = 1.20438044 * To / To  # gas gamma
        Mw           = 13.60087 * To / To  # gas molar weight
        OF_opt       = 6 * To / To   # optimum OF
        Tt_comb      = 3629.8 * To / To  # Adiabatic flame temperature
        
        Tt_ad_flame = 3330. #function
        gamma       = 1.4   # function
        Mw          = 13.3  # function
        OF_opt      = 6.0   # function
        

        # Scaling factor: combustion temperature is always lower than adiabatic flame temperature
        scale = .8
        Tt_comb = scale * Tt_ad_flame
        
        # pack outputs
        Rm  = 8134./Mw
        Cp  = gamma/(gamma-1)*Rm
        
        # pack computed quantities into outputs
        self.outputs.stagnation_temperature             = Tt_comb
        self.outputs.stagnation_pressure                = Pt_comb
        self.outputs.oxidizer_fuel_ratio                = OF 
        self.outputs.isentropic_expansion_factor        = gamma
        self.outputs.specific_gas_constant              = Rm
        self.outputs.specific_heat_constant_pressure    = Cp
    
        
    __call__ = compute
