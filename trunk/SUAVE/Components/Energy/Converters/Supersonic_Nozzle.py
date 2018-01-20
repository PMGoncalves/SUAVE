## @ingroup Components-Energy-Converters
# Supersonic_Nozzle.py
#
# Created:  May 2015, T. MacDonald
# Modified: Jan 2016, T. MacDonald

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# SUAVE imports

import SUAVE

from SUAVE.Core import Units

# package imports
import numpy as np
import scipy as sp
from scipy.optimize import fsolve
from warnings import warn

from SUAVE.Components.Energy.Energy_Component import Energy_Component
from SUAVE.Methods.Propulsion.fm_id import fm_id

# ----------------------------------------------------------------------
#  Expansion Nozzle Component
# ----------------------------------------------------------------------

## @ingroup Components-Energy-Converters
class Supersonic_Nozzle(Energy_Component):
    """This is a nozzle component that allows for supersonic outflow.
    Calling this class calls the compute function.
    
    Assumptions:
    Pressure ratio and efficiency do not change with varying conditions.
    
    Source:
    https://web.stanford.edu/~cantwell/AA283_Course_Material/AA283_Course_Notes/
    """
    
    def __defaults__(self):
        """ This sets the default values for the component to function.
        
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
        
        #set the defaults
        self.tag = 'Nozzle'
        self.polytropic_efficiency           = 1.0
        self.pressure_ratio                  = 1.0
        self.inputs.stagnation_temperature   = 0.
        self.inputs.stagnation_pressure      = 0.
        self.outputs.stagnation_temperature  = 0.
        self.outputs.stagnation_pressure     = 0.
        self.outputs.stagnation_enthalpy     = 0.
    
    
    
    def compute(self,conditions):
        """This computes the output values from the input values according to
        equations from the source.
        
        Assumptions:
        Constant polytropic efficiency and pressure ratio
        
        Source:
        https://web.stanford.edu/~cantwell/AA283_Course_Material/AA283_Course_Notes/
        
        Inputs:
        conditions.freestream.
          isentropic_expansion_factor         [-]
          specific_heat_at_constant_pressure  [J/(kg K)]
          pressure                            [Pa]
          stagnation_pressure                 [Pa]
          stagnation_temperature              [K]
          universal_gas_constant              [J/(kg K)] (this is misnamed - actually refers to the gas specific constant)
          mach_number                         [-]
        self.inputs.
          stagnation_temperature              [K]
          stagnation_pressure                 [Pa]
                   
        Outputs:
        self.outputs.
          stagnation_temperature              [K]  
          stagnation_pressure                 [Pa]
          stagnation_enthalpy                 [J/kg]
          mach_number                         [-]
          static_temperature                  [K]
          static_enthalpy                     [J/kg]
          velocity                            [m/s]
          static_pressure                     [Pa]
          area_ratio                          [-]
                
        Properties Used:
        self.
          pressure_ratio                      [-]
          polytropic_efficiency               [-]
        """           
        
        #unpack the values
        
        #unpack from conditions
        gamma    = conditions.freestream.isentropic_expansion_factor
        Cp       = conditions.freestream.specific_heat_at_constant_pressure
        Po       = conditions.freestream.pressure
        Pto      = conditions.freestream.stagnation_pressure
        Tto      = conditions.freestream.stagnation_temperature
        R        = conditions.freestream.universal_gas_constant
        Mo       = conditions.freestream.mach_number
        
        #unpack from inputs
        Tt_in    = self.inputs.stagnation_temperature
        Pt_in    = self.inputs.stagnation_pressure
        
        #unpack from self
        pid      = self.pressure_ratio
        etapold  = self.polytropic_efficiency
        
        
        #Method for computing the nozzle properties
        
        #--Getting the output stagnation quantities
        Pt_out   = Pt_in*pid
        Tt_out   = Tt_in*pid**((gamma-1)/(gamma)*etapold)
        ht_out   = Cp*Tt_out
        
        
        #compute the output Mach number, static quantities and the output velocity
        Mach          = np.sqrt((((Pt_out/Po)**((gamma-1)/gamma))-1)*2/(gamma-1))
        
        #Remove check on mach numbers from expansion nozzle
        i_low         = Mach < 10.0
        
        #initializing the Pout array
        P_out         = 1.0 *Mach/Mach
        
        #Computing output pressure and Mach number for the case Mach <1.0
        P_out[i_low]  = Po[i_low]
        Mach[i_low]   = np.sqrt((((Pt_out[i_low]/Po[i_low])**((gamma-1)/gamma))-1)*2/(gamma-1))
        
        #Computing the output temperature,enthalpy, velocity and density
        T_out         = Tt_out/(1+(gamma-1)/2*Mach*Mach)
        h_out         = Cp*T_out
        u_out         = np.sqrt(2*(ht_out-h_out))
        rho_out       = P_out/(R*T_out)
        
        #Computing the freestream to nozzle area ratio (mainly from thrust computation)
        area_ratio    = (fm_id(Mo)/fm_id(Mach)*(1/(Pt_out/Pto))*(np.sqrt(Tt_out/Tto)))
        
        #pack computed quantities into outputs
        self.outputs.stagnation_temperature  = Tt_out
        self.outputs.stagnation_pressure     = Pt_out
        self.outputs.stagnation_enthalpy     = ht_out
        self.outputs.mach_number             = Mach
        self.outputs.static_temperature      = T_out
        self.outputs.density                 = rho_out
        self.outputs.static_enthalpy         = h_out
        self.outputs.velocity                = u_out
        self.outputs.static_pressure         = P_out
        self.outputs.area_ratio              = area_ratio
            
        
        
    def compute_rocket(self,conditions):
        """This computes the output values from the input values according to
        equations from the source.
        
        Assumptions:
        
        
        Source:
        
        
        Inputs:
        """           
        
        #unpack the values
        
        #unpack from conditions
        Po       = conditions.freestream.pressure
        R        = conditions.freestream.universal_gas_constant
        
        #unpack from inputs
        Tt_in    = self.inputs.stagnation_temperature
        Pt_in    = self.inputs.stagnation_pressure
        gamma    = self.inputs.isentropic_expansion_factor
        Rm       = self.inputs.specific_gas_constant
        Cp       = self.inputs.specific_heat_constant_pressure
        
        #unpack from self
        exp_ratio  = self.expansion_ratio
#        pid        = self.pressure_ratio
#        etapold    = self.polytropic_efficiency
        
        # -- Calculating flow properties
        
        # Vandenkerckhove function
        a = gamma
        b = ((1+gamma)/2)**((1+gamma)/(1-gamma))        
        Gf = np.sqrt(a*b)
        
        # P_out calculation
        a = Gf
        b = 2*gamma/(gamma-1)
        
        
        func = lambda P_out : exp_ratio - (  a / ( np.sqrt(b * (P_out/Pt_in)**(2/gamma) * (1-(P_out/Pt_in)**((gamma-1)/gamma)))))
        P_out = fsolve(func,Po,factor = 0.01)

        # in case pressures go too low
        if np.any(P_out<0.4*Po):
            warn('P_out goes too low',RuntimeWarning)
            P_out[P_out<0.4*Po] = 0.4*Po[P_out<0.4*Po]

        # Calculate other flow properties
        Pt_out = Pt_in #possibly add pressure ratio
        Tt_out = Tt_in #possibly add adiabatic efficiency
        u_out  = np.sqrt((2/(gamma-1))*(Rm)*Tt_out * (1 - (P_out/Pt_out)**((gamma-1)/gamma)))
        
        T_out  = Tt_out - u_out*u_out/(2*Cp)       
        Mach   = u_out / np.sqrt(gamma*Rm*T_out)

        
        #pack computed quantities into outputs
        self.outputs.stagnation_temperature  = Tt_out
        self.outputs.stagnation_pressure     = Pt_out
        self.outputs.velocity                = u_out
        self.outputs.static_pressure         = P_out
        self.outputs.mach_number             = Mach


    __call__ = compute