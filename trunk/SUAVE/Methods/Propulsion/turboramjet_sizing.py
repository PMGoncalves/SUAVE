# turboramjet_sizing.py
# 
# Created:  Mar 2015, A. Variyar 
# Modified: Feb 2016, M. Vegh
#           Jan 2016, E. Botero

# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import SUAVE
import numpy as np
from SUAVE.Core import Data

# ----------------------------------------------------------------------
#   Sizing
# ----------------------------------------------------------------------

def turboramjet_sizing(turboramjet,mach_number = None, altitude = None, delta_isa = 0, conditions = None):  
    """ create and evaluate a turboramjet network
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

        
    # Shared components
    ram                       = turboramjet.ram
    inlet_nozzle              = turboramjet.inlet_nozzle
    combustor_2               = turboramjet.combustor_2
    core_nozzle               = turboramjet.core_nozzle
    mixer                     = turboramjet.mixer

        
    #-- turbojet-only components
    low_pressure_compressor   = turboramjet.low_pressure_compressor
    high_pressure_compressor  = turboramjet.high_pressure_compressor
    combustor                 = turboramjet.combustor
    high_pressure_turbine     = turboramjet.high_pressure_turbine
    low_pressure_turbine      = turboramjet.low_pressure_turbine

        
    thrust                    = turboramjet.thrust
    number_of_engines         = turboramjet.number_of_engines   
    turbojet_mach             = turboramjet.turbojet_mach
    ramjet_mach               = turboramjet.ramjet_mach
        
    #Creating the network by manually linking the different components
        
    #set the working fluid to determine the fluid properties
    ram.inputs.working_fluid                               = turboramjet.working_fluid
        
    #Flow through the ram , this computes the necessary flow quantities and stores it into conditions
    ram(conditions) 
    
    #link inlet nozzle to ram 
    inlet_nozzle.inputs.stagnation_temperature             = ram.outputs.stagnation_temperature 
    inlet_nozzle.inputs.stagnation_pressure                = ram.outputs.stagnation_pressure
    
    #Flow through the inlet nozzle
    inlet_nozzle(conditions)        
        
    #Bypass system
    
    Mo      = conditions.freestream.mach_number
    #-- Defines turbojet operation
    i_tj    = Mo <= turbojet_mach
    
    #-- Defines dual operation mode
    i_mx    = np.logical_and(Mo > turbojet_mach, Mo < ramjet_mach)
    #-- Defines ramjet operation
    i_rj    = Mo >= ramjet_mach
    
    #-- Determine ramjet bypass ratio
    # Only used for Mo between max Mach for turbojet and min Mach for ramjet
    # ram_bypass = 1 for Mo = min Mach for ramjet
    # ram_bypass = 0 for Mo = max Mach for turbojet
    ram_bypass = 1/(ramjet_mach - turbojet_mach) * (Mo - turbojet_mach)

 
    #----------------------
    # TURBOJET-ONLY OPERATION
    #----------------------

    #--link low pressure compressor to the inlet nozzle
    low_pressure_compressor.inputs.stagnation_temperature  = inlet_nozzle.outputs.stagnation_temperature
    low_pressure_compressor.inputs.stagnation_pressure     = inlet_nozzle.outputs.stagnation_pressure
        
    #Flow through the low pressure compressor
    low_pressure_compressor(conditions)
    
            
    #link the high pressure compressor to the low pressure compressor
    high_pressure_compressor.inputs.stagnation_temperature = low_pressure_compressor.outputs.stagnation_temperature
    high_pressure_compressor.inputs.stagnation_pressure    = low_pressure_compressor.outputs.stagnation_pressure
            
    #Flow through the high pressure compressor
    high_pressure_compressor(conditions)
            
    #link the combustor to the high pressure compressor
    combustor.inputs.stagnation_temperature                = high_pressure_compressor.outputs.stagnation_temperature
    combustor.inputs.stagnation_pressure                   = high_pressure_compressor.outputs.stagnation_pressure
        
    #flow through the high pressor comprresor
    combustor(conditions)
                
    #link the high pressure turbine to the combustor
    high_pressure_turbine.inputs.stagnation_temperature    = combustor.outputs.stagnation_temperature
    high_pressure_turbine.inputs.stagnation_pressure       = combustor.outputs.stagnation_pressure
    high_pressure_turbine.inputs.fuel_to_air_ratio         = combustor.outputs.fuel_to_air_ratio
    
    #link the high pressuer turbine to the high pressure compressor
    high_pressure_turbine.inputs.compressor                = high_pressure_compressor.outputs
    
    #flow through the high pressure turbine

    high_pressure_turbine.inputs.bypass_ratio   = 0.0
    high_pressure_turbine.inputs.fan            = Data()
    high_pressure_turbine.inputs.fan.work_done  = 0.0
    high_pressure_turbine(conditions)
                    
     #link the low pressure turbine to the high pressure turbine
    low_pressure_turbine.inputs.stagnation_temperature     = high_pressure_turbine.outputs.stagnation_temperature
    low_pressure_turbine.inputs.stagnation_pressure        = high_pressure_turbine.outputs.stagnation_pressure
    
    #link the low pressure turbine to the low_pressure_compresor
    low_pressure_turbine.inputs.compressor                 = low_pressure_compressor.outputs
    
    #link the low pressure turbine to the combustor
    low_pressure_turbine.inputs.fuel_to_air_ratio          = combustor.outputs.fuel_to_air_ratio
    	
    #get the bypass ratio from the thrust component
    low_pressure_turbine.inputs.bypass_ratio               = 0.0
            
    #flow through the low pressure turbine
    low_pressure_turbine.inputs.bypass_ratio    = 0.0
    low_pressure_turbine.inputs.fan             = Data()
    low_pressure_turbine.inputs.fan.work_done   = 0.0    
    low_pressure_turbine(conditions)
    

    #----------------------
    # SEPARATION OF OPERATIONS
    #----------------------
     
    
    #-- Turbojet operation
    #---- Ignore the second combustor
    if i_tj :
        # Link the nozzle to the low pressure turbine
        core_nozzle.inputs.stagnation_temperature       = low_pressure_turbine.outputs.stagnation_temperature
        core_nozzle.inputs.stagnation_pressure          = low_pressure_turbine.outputs.stagnation_pressure
    
    
    else :
        #-- Ramjet operation
        if i_rj:
            # Link the second combustor to the network
            combustor_2.inputs.stagnation_temperature       = inlet_nozzle.outputs.stagnation_temperature
            combustor_2.inputs.stagnation_pressure          = inlet_nozzle.outputs.stagnation_pressure
        
        #-- Dual mode operation
        if i_mx:
            # mixing properties, assuming same pressure
            
            mixer.inputs.stagnation_temperature     = inlet_nozzle.outputs.stagnation_temperature
            mixer.inputs.stagnation_pressure        = inlet_nozzle.outputs.stagnation_pressure
            mixer.inputs.mach                       = inlet_nozzle.outputs.mach_number
            mixer.inputs.stagnation_temperature_2   = low_pressure_turbine.outputs.stagnation_temperature
            mixer.inputs.stagnation_pressure_2      = low_pressure_turbine.outputs.stagnation_pressure

            mixer(conditions)
            
#            
#            mixed_temperature = ram_bypass*inlet_nozzle.outputs.stagnation_temperature + (1-ram_bypass)*low_pressure_turbine.outputs.stagnation_temperature
#            mixed_pressure    = low_pressure_turbine.outputs.stagnation_pressure
            
            # link the second combustor to the network
            combustor_2.inputs.stagnation_temperature       = mixer.outputs.stagnation_temperature
            combustor_2.inputs.stagnation_pressure          = mixer.outputs.stagnation_pressure

        # flow through combustor
        combustor_2(conditions)
        
        # Link the nozzle to the combustor
        core_nozzle.inputs.stagnation_temperature       = combustor_2.outputs.stagnation_temperature
        core_nozzle.inputs.stagnation_pressure          = combustor_2.outputs.stagnation_pressure
    
    # flow through nozzle
    core_nozzle(conditions)  
    
    #link the thrust component to the core nozzle
    thrust.inputs.core_exit_velocity                       = core_nozzle.outputs.velocity
    thrust.inputs.core_area_ratio                          = core_nozzle.outputs.area_ratio
    thrust.inputs.core_nozzle                              = core_nozzle.outputs
        
    
    if i_tj :
        # if turbojet mode only, disregard fuel calculations for second combustor
        # apply correct reference parameters
        combustor_2.outputs.fuel_to_air_ratio              = 0.0
        thrust.inputs.total_temperature_reference          = low_pressure_compressor.outputs.stagnation_temperature
        thrust.inputs.total_pressure_reference             = low_pressure_compressor.outputs.stagnation_pressure

    else :
        thrust.inputs.total_temperature_reference          = inlet_nozzle.outputs.stagnation_temperature
        thrust.inputs.total_pressure_reference             = inlet_nozzle.outputs.stagnation_pressure
        
        if i_rj :
            # if ramjet mode only, disregard fuel calculations for first combustor
            combustor.outputs.fuel_to_air_ratio                = 0.0
            



    thrust.inputs.fuel_to_air_ratio                        = combustor.outputs.fuel_to_air_ratio + combustor_2.outputs.fuel_to_air_ratio


    thrust.inputs.number_of_engines                        = number_of_engines
    thrust.inputs.fan_nozzle                               = Data()
    thrust.inputs.fan_nozzle.velocity                      = 0.0
    thrust.inputs.fan_nozzle.area_ratio                    = 0.0
    thrust.inputs.fan_nozzle.static_pressure               = 0.0
    thrust.inputs.bypass_ratio                             = 0.0
    thrust.inputs.flow_through_core                        = 1.0 #scaled constant to turn on core thrust computation
    thrust.inputs.flow_through_fan                         = 0.0 #scaled constant to turn on fan thrust computation        

    #compute the thrust
    thrust.size(conditions)
    
    #update the design thrust value
    turboramjet.design_thrust = thrust.total_design

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
    results_sls = turboramjet.size(state_sls)
    turboramjet.sealevel_static_thrust = results_sls.thrust_force_vector[0,0] / number_of_engines 