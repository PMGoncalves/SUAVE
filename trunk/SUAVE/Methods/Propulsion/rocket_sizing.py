# rocket_sizing.py
# 
# Created:  May 2015, T. MacDonald 
# Modified: Jan 2016, E. Botero        

# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import SUAVE
import numpy as np
from SUAVE.Core import Data

# ----------------------------------------------------------------------
#   Sizing
# ----------------------------------------------------------------------

def rocket_sizing(rocket,mach_number = None, altitude = None, delta_isa = 0, conditions = None):  
    """ create and evaluate a gas turbine network
    """    
    
    #Unpack components
    #check if altitude is passed or conditions is passed
    if(conditions):
        #use conditions
        pass
        
    else:
        #check if mach number and temperature are passed
        if(mach_number==None or altitude==None):
            
            #raise an error
            raise NameError('The sizing conditions require an altitude and a Mach number')
        
        else:
            #call the atmospheric model to get the conditions at the specified altitude
            atmosphere = SUAVE.Analyses.Atmospheric.US_Standard_1976()
            atmo_data = atmosphere.compute_values(altitude,delta_isa)

            p   = atmo_data.pressure          
            T   = atmo_data.temperature       
            rho = atmo_data.density          
            a   = atmo_data.speed_of_sound    
            mu  = atmo_data.dynamic_viscosity   
        
            # setup conditions
            conditions = SUAVE.Analyses.Mission.Segments.Conditions.Aerodynamics()            
        
            # freestream conditions    
            conditions.freestream.altitude           = np.atleast_1d(altitude)
            conditions.freestream.mach_number        = np.atleast_1d(mach_number)
            conditions.freestream.pressure           = np.atleast_1d(p)
            conditions.freestream.temperature        = np.atleast_1d(T)
            conditions.freestream.density            = np.atleast_1d(rho)
            conditions.freestream.dynamic_viscosity  = np.atleast_1d(mu)
            conditions.freestream.gravity            = np.atleast_1d(9.81)
            conditions.freestream.gamma              = np.atleast_1d(1.4)
            conditions.freestream.Cp                 = 1.4*(p/(rho*T))/(1.4-1)
            conditions.freestream.R                  = p/(rho*T)
            conditions.freestream.speed_of_sound     = np.atleast_1d(a)
            conditions.freestream.velocity           = np.atleast_1d(a*mach_number)
            
            # propulsion conditions
            conditions.propulsion.throttle           =  np.atleast_1d(1.0)
    

    combustor                 = rocket.combustor
    nozzle                    = rocket.nozzle
    thrust                    = rocket.thrust   

    number_of_engines         = rocket.number_of_engines
    
    #Creating the network by manually linking the different components
    #Creating the network by manually linking the different components
#    combustor.oxidizer_data                               = self.oxidizer
 #   combustor.fuel_data                                   = self.fuel
    
    #flow through combustor
    combustor.compute_rocket(conditions)

    #link the nozzle to the combustor
    nozzle.inputs.stagnation_temperature                = combustor.outputs.stagnation_temperature
    nozzle.inputs.stagnation_pressure                   = combustor.outputs.stagnation_pressure
    nozzle.inputs.isentropic_expansion_factor           = combustor.outputs.isentropic_expansion_factor
    nozzle.inputs.specific_gas_constant                 = combustor.outputs.specific_gas_constant 
    nozzle.inputs.specific_heat_constant_pressure       = combustor.outputs.specific_heat_constant_pressure    

    #flow through the core nozzle
    nozzle.compute_rocket(conditions)
    
    # compute the thrust using the thrust component
    #link the thrust component to the combustor
    
    thrust.inputs.oxidizer_fuel_ratio                   = combustor.outputs.oxidizer_fuel_ratio
    thrust.inputs.number_of_engines                     = number_of_engines
    thrust.inputs.core_nozzle                           = nozzle.outputs
    thrust.inputs.stagnation_temperature                = combustor.outputs.stagnation_temperature
    thrust.inputs.stagnation_pressure                   = combustor.outputs.stagnation_pressure
    thrust.inputs.isentropic_expansion_factor           = combustor.outputs.isentropic_expansion_factor
    thrust.inputs.specific_gas_constant                 = combustor.outputs.specific_gas_constant
    
    #compute the thrust
    thrust.size_rocket(conditions)
    
    #update the design thrust value
    rocket.design_thrust = thrust.total_design

    #compute the sls_thrust
    #call the atmospheric model to get the conditions at the specified altitude
    atmosphere_sls = SUAVE.Analyses.Atmospheric.US_Standard_1976()
    atmo_data = atmosphere_sls.compute_values(0.0,0.0)
    
    p   = atmo_data.pressure          
    T   = atmo_data.temperature       
    rho = atmo_data.density          
    a   = atmo_data.speed_of_sound    
    mu  = atmo_data.dynamic_viscosity   

    # setup conditions
    conditions_sls = SUAVE.Analyses.Mission.Segments.Conditions.Aerodynamics()            

    # freestream conditions    
    conditions_sls.freestream.altitude           = np.atleast_1d(0.)
    conditions_sls.freestream.mach_number        = np.atleast_1d(0.01)
    conditions_sls.freestream.pressure           = np.atleast_1d(p)
    conditions_sls.freestream.temperature        = np.atleast_1d(T)
    conditions_sls.freestream.density            = np.atleast_1d(rho)
    conditions_sls.freestream.dynamic_viscosity  = np.atleast_1d(mu)
    conditions_sls.freestream.gravity            = np.atleast_1d(9.81)
    conditions_sls.freestream.gamma              = np.atleast_1d(1.4)
    conditions_sls.freestream.Cp                 = 1.4*(p/(rho*T))/(1.4-1)
    conditions_sls.freestream.R                  = p/(rho*T)
    conditions_sls.freestream.speed_of_sound     = np.atleast_1d(a)
    conditions_sls.freestream.velocity           = np.atleast_1d(a*0.01)
    
    # propulsion conditions
    conditions_sls.propulsion.throttle           =  np.atleast_1d(1.0)    
    
    state_sls = Data()
    state_sls.numerics = Data()
    state_sls.conditions = conditions_sls   
    results_sls = rocket.evaluate_thrust(state_sls)
    rocket.sealevel_static_thrust = results_sls.thrust_force_vector[0,0] / number_of_engines