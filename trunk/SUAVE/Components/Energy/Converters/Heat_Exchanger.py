## @ingroup Components-Energy-Converters
# Compression_Nozzle.py
#
# Created:  Jul 2014, A. Variyar
# Modified: Jan 2016, T. MacDonald

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE

# python imports
from warnings import warn

# package imports
import numpy as np

from SUAVE.Components.Energy.Energy_Component import Energy_Component
from SUAVE.Core import Data

# ----------------------------------------------------------------------
#  Compression Nozzle Component
# ----------------------------------------------------------------------
## @ingroup Components-Energy-Converters
class Heat_Exchanger(Energy_Component):
    """This is a heat exchanger component inteded for use in a precooled engine
    Calling this class calls the compute function.
    
    Assumptions:
    
    Source:
        
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
        #setting the default values 
        self.tag = 'heat_exchanger'
        self.efficiency                      = 1.0
        self.pressure_ratio_A                = 1.0
        self.pressure_ratio_B                = 1.0
        self.inputs.stagnation_temperature_A = 0.
        self.inputs.stagnation_pressure_A    = 0.
        self.inputs.stagnation_temperature_B = 0.
        self.inputs.stagnation_pressure_B    = 35.
        self.refrigerant_data                = Data()
        self.heat_capacity_ratio             = 1.


    

    def compute(self,conditions):
        """ This computes the output values from the input values according to
        equations from the source.

        Assumptions:
        Efficiency based on delta T parameter:
            100% efficiency => delta T = 0.
            0%   efficiency => delta T = Tt_out - Tt_in
            
        Liquid H2 used for refrigerant (seen in value for Cp)

        Source:
        https://web.stanford.edu/~cantwell/AA283_Course_Material/AA283_Course_Notes/
        https://arc.aiaa.org/doi/abs/10.2514/6.IAC-06-D2.P.2.07
        

        Inputs:
        
        
        conditions.freestream.
          specific_heat_at_constant_pressure    [J/(kg K)]
        
        self.inputs.
          self.inputs.stagnation_temperature_A  [K]
          self.inputs.stagnation_temperature_B  [K]
          self.inputs.stagnation_pressure_A     [Pa]
          self.inputs.stagnation_pressure_B     [Pa]
          
        Outputs:
        self.outputs.
          stagnation_temperature                [K]  
          stagnation_pressure                   [Pa]
          stagnation_enthalpy                   [J/kg]
          mach_number                           [-]
          static_temperature                    [K]
          static_enthalpy                       [J/kg]
          velocity                              [m/s]

        Properties Used:
        self.
          pressure_ratio_A                      [-]
          pressure_ratio_B                      [-]
          efficiency                            [-]
          heat_capacity_ratio                   [-]
        """  

        #unpack the values
        
        #unpack from conditions
        Cp      = conditions.freestream.specific_heat_at_constant_pressure

        
        #unpack from inpust
        TtA_in  = self.inputs.stagnation_temperature_A
        TtB_in  = self.inputs.stagnation_temperature_B
        PtA_in  = self.inputs.stagnation_pressure_A
        PtB_in  = self.inputs.stagnation_pressure_B
        
        #unpack from self
        pid_A                   =  self.pressure_ratio_A
        pid_B                   =  self.pressure_ratio_B
        K                       =  self.heat_capacity_ratio
        eta                     =  self.efficiency
    
        # Method to compute the output variables
        
        #--Convert efficiency into delta Temperature
        dT = (1 - eta)*abs(TtA_in-TtB_in)
           
        #--Determinte hot and cold flow   
        #----Initialize array
        TtH_in  = 1.0*TtA_in/TtA_in
        TtC_in  = 1.0*TtA_in/TtA_in
        
        #----Determine conditions
        i_hot   = TtA_in > TtB_in
        i_cold  = TtA_in <= TtB_in
        
        TtH_in[i_hot]   = TtA_in[i_hot]
        TtH_in[i_cold]  = TtB_in[i_cold]
        
        TtC_in[i_hot]   = TtB_in[i_hot]
        TtC_in[i_cold]  = TtA_in[i_cold]
        
        #--Determine common output stagnation temperature
        Tt_out = (TtH_in + TtC_in*K + dT*(K-1)) / (1+K)
        
        #--Add dT to the final output stagnation temperatures of each flow
        TtH_out = Tt_out + dT
        TtC_out = Tt_out - dT
        
        #--Compute final stagnation pressure
        PtA_out = PtA_in*pid_A
        PtB_out = PtB_in*pid_B
        TtA_out = TtH_out
        TtB_out = TtC_out
        
        CpB = 34283.28
        CpA = Cp
        
        f       = K*CpA/CpB
        

        #pack computed quantities into outputs
        self.outputs.stagnation_temperature_A  = TtA_out
        self.outputs.stagnation_temperature_B  = TtB_out
        self.outputs.stagnation_pressure_A     = PtA_out
        self.outputs.stagnation_pressure_B     = PtB_out
        self.outputs.refrigerant_to_air_ratio  = f

    

    __call__ = compute