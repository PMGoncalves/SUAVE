## @ingroup Components-Energy-Converters
# Combustor.py
#
# Created:  Oct 2014, A. Variyar
# Modified: Jan 2016, T. MacDonald
#           Sep 2017, P. Goncalves
#           Jan 2018, W. Maier

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE

import numpy as np
from scipy.optimize import fsolve
from warnings import warn

from SUAVE.Core import Data
from SUAVE.Components.Energy.Energy_Component import Energy_Component
from SUAVE.Methods.Propulsion.rayleigh import rayleigh
from SUAVE.Methods.Propulsion.fm_solver import fm_solver

# ----------------------------------------------------------------------
#  Combustor Component
# ----------------------------------------------------------------------
## @ingroup Components-Energy-Converters
class Combustor(Energy_Component):
    """This provides output values for a combustor
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
        self.fuel_data                          = SUAVE.Attributes.Propellants.Jet_A()
        self.alphac                             = 0.0
        self.turbine_inlet_temperature          = 1.0
        self.inputs.stagnation_temperature      = 1.0
        self.inputs.stagnation_pressure         = 1.0
        self.inputs.static_pressure             = 1.0
        self.inputs.mach_number                 = 0.1
        self.outputs.stagnation_temperature     = 1.0
        self.outputs.stagnation_pressure        = 1.0
        self.outputs.static_pressure            = 1.0
        self.outputs.stagnation_enthalpy        = 1.0
        self.outputs.fuel_to_air_ratio          = 1.0
        self.fuel_data                          = Data()
        self.area_ratio                         = 1.0
        self.specific_heat_at_constant_pressure = 1510.
        self.isentropic_expansion_factor        = 1.238
        self.axial_fuel_velocity_ratio          = 0.5
        self.fuel_velocity_ratio                = 0.5
        self.burner_drag_coefficient            = 0.1
        self.temperature_reference              = 222.
        self.hf                                 = 0.


    
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
          area_ratio                          [-]
          fuel_data.specific_energy           [J/kg]
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
        ar     = self.area_ratio
        
        # compute pressure
        Pt_out = Pt_in*pib


        # method to compute combustor properties

        # method - computing the stagnation enthalpies from stagnation temperatures
        ht4     = Cp*Tt4
        ho      = Cp*To
        ht_in   = Cp*Tt_in
        
        # Using the Turbine exit temperature, the fuel properties and freestream temperature to compute the fuel to air ratio f
        f       = (ht4 - ht_in)/(eta_b*htf-ht4)

        # Computing the exit static and stagnation conditions
        ht_out  = Cp*Tt4
        
        # pack computed quantities into outputs
        self.outputs.stagnation_temperature  = Tt4
        self.outputs.stagnation_pressure     = Pt_out
        self.outputs.stagnation_enthalpy     = ht_out
        self.outputs.fuel_to_air_ratio       = f 


    def compute_rayleigh(self,conditions):
        """ This combutes the temperature and pressure change across the
        the combustor using Rayleigh Line flow; it checks for themal choking.

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
          area_ratio                          [-]
          fuel_data.specific_energy           [J/kg]
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
        Mach = self.inputs.mach_number
        Tt4    = self.turbine_inlet_temperature
        pib    = self.pressure_ratio
        eta_b  = self.efficiency
        
        # unpacking values from self
        htf    = self.fuel_data.specific_energy
        ar     = self.area_ratio
        
        # Rayleigh flow analysis, constant pressure burner
            
        # Initialize arrays
        M_out  = np.ones_like(Pt_in)
        M_aux  = np.ones_like(Pt_in)
        Ptr    = np.ones_like(Pt_in)
        Tt4_ray= np.ones_like(Pt_in)        
        
        # Isentropic decceleration through divergent nozzle     
        M_aux   = fm_solver(ar,Mach[:,0],gamma[:,0])  
        Mach[:,0] = M_aux
        
        # Determine max stagnation temperature to thermally choke flow                                     
        Tt4_ray = Tt_in*(1.+gamma*Mach*Mach)**2./((2.*(1.+gamma)*Mach*Mach)*(1.+(gamma-1.)/2.*Mach*Mach))    
        
#        Ttray_aux[:,0] = Tt4_ray[:,0]

        # Rayleigh limitations define Tt4, taking max temperature before choking
        Tt4 = Tt4 * np.ones_like(Tt4_ray)


        Tt4[Tt4_ray <= Tt4] = Tt4_ray[Tt4_ray <= Tt4]
        
        
        #Rayleigh calculations
        M_out[:,0], Ptr[:,0] = rayleigh(gamma[:,0],Mach[:,0],Tt4[:,0]/Tt_in[:,0]) 

        Pt_out     = Ptr*Pt_in
            
        # method to compute combustor properties

        # method - computing the stagnation enthalpies from stagnation temperatures
        ht4     = Cp*Tt4
        ho      = Cp*To
        ht_in   = Cp*Tt_in
        
        # Using the Turbine exit temperature, the fuel properties and freestream temperature to compute the fuel to air ratio f
        f       = (ht4 - ht_in)/(eta_b*htf-ht4)

        # Computing the exit static and stagnation conditions
        ht_out  = Cp*Tt4   #May be double counting here.....no need (maybe)

        
        # pack computed quantities into outputs
        self.outputs.stagnation_temperature  = Tt4
        self.outputs.stagnation_pressure     = Pt_out
        self.outputs.stagnation_enthalpy     = ht_out
        self.outputs.fuel_to_air_ratio       = f    
        self.outputs.mach_number             = M_out
  

    def compute_scramjet(self,conditions):
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
          axial_fuel_velocity_ratio           [-]
          fuel_velocity_ratio                 [-]
          burner_drag_coefficient             [-]
          temperature_reference               [K]
          Cpb                                 [-]
          hf
        """         
        # unpack the values
        
        # unpacking the values from conditions
        gamma  = conditions.freestream.isentropic_expansion_factor 
        Cp     = conditions.freestream.specific_heat_at_constant_pressure
        To     = conditions.freestream.temperature
        Tto    = conditions.freestream.stagnation_temperature
        R      = conditions.freestream.gas_specific_constant
        
        # unpacking the values form inputs
        Tt_in  = self.inputs.stagnation_temperature
        Pt_in  = self.inputs.stagnation_pressure
        Tt_4   = self.turbine_inlet_temperature
        M_in   = self.inputs.mach_number
        T_in   = self.inputs.static_temperature
        V_in   = self.inputs.velocity
        P_in   = self.inputs.static_pressure
        
        # unpacking values from self    
            #-- Fuel
        htf    = self.fuel_data.specific_energy
        f_st   = self.fuel_data.stoichiometric_air
        
            #-- Combustion process
        eta_b  = self.efficiency
        Vfx_V3 = self.axial_fuel_velocity_ratio
        Vf_V3  = self.fuel_velocity_ratio
        CfAwA3 = self.burner_drag_coefficient
        Tref   = self.temperature_reference
        hf     = self.hf
        
            #-- Combustor losses
        pid    = self.pressure_ratio
        
            #-- Gas properties
                #-- New flow properties --not computed, these values should change considerably
        Cpb    = self.specific_heat_at_constant_pressure
        gamma  = self.isentropic_expansion_factor 
        
        
        # Initialize
        f       = f_st*np.ones_like(Tt_in)
           
        V_out   = V_in*(((1+f*Vfx_V3)/(1+f))-(CfAwA3/(2*(1+f))))
        Tt_out  = (T_in/(1+f))*(1+(1/(Cpb*T_in))*(eta_b*f*htf+f*hf+f*Cpb*Tref+(1+f*(Vf_V3)**2)*V_in**2/2))
        T_out   = Tt_out - V_out**2/(2*Cpb)   
        M_out   = V_out/np.sqrt(gamma*R*T_out)
        
#        i_merda = M_out < 1.0
#        M = conditions.freestream.mach_number * Pt_in/Pt_in
#        print 'M erro :', M[i_merda]
#        
        if np.any(M_out<1.0):
            warn('Subsonic combustion, higher Mach required',RuntimeWarning)
            M_out[M_out<1.0] = 0.01
            f[M_out < 1.0] = 0.0

        # Computing the exit static and stagnation conditions
        ht_out  = Cp*Tt_out
        P_out   = P_in
        Pt_out  = P_out*(1+(gamma-1)/2*M_out**2)**(gamma/(gamma-1))
        
        # pack computed quantities into outputs
        self.outputs.stagnation_temperature  = Tt_out
        self.outputs.stagnation_pressure     = Pt_out
        self.outputs.stagnation_enthalpy     = ht_out
        self.outputs.fuel_to_air_ratio       = f
        self.outputs.static_temperature      = T_out
        self.outputs.static_pressure         = P_out
        self.outputs.velocity                = V_out
        self.outputs.mach_number             = M_out


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
        
        OF       = np.ones_like(To)*6.0
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
    
