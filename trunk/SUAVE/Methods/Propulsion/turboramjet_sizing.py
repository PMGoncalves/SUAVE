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

        
    # Shared components
    ram                       = turboramjet.ram
    inlet_nozzle              = turboramjet.inlet_nozzle
    combustor_2               = turboramjet.combustor_2
    core_nozzle               = turboramjet.core_nozzle

        
    #-- turboramjet components
    low_pressure_compressor   = turboramjet.low_pressure_compressor
    high_pressure_compressor  = turboramjet.high_pressure_compressor
    combustor                 = turboramjet.combustor
    high_pressure_turbine     = turboramjet.high_pressure_turbine
    low_pressure_turbine      = turboramjet.low_pressure_turbine

        
    thrust                    = turboramjet.thrust
    bypass_ratio              = turboramjet.bypass_ratio
    number_of_engines         = turboramjet.number_of_engines   
    mach_separation           = turboramjet.mach_separation
        
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
    
    print '-------------------------------------------'
    print 'INLET '
    print 'Pressure    :', inlet_nozzle.outputs.stagnation_pressure
    print 'Temperature :', inlet_nozzle.outputs.stagnation_temperature
        
        
    # Bypass system
    Mo       = conditions.freestream.mach_number
    i_tj = Mo < mach_separation
    i_rj = Mo > mach_separation
    
    print 'TURBOJET ', i_tj

        
    #----------------------
    # TURBOJET OPERATION
    #----------------------
        
    #--link low pressure compressor to the inlet nozzle
    low_pressure_compressor.inputs.stagnation_temperature  = inlet_nozzle.outputs.stagnation_temperature
    low_pressure_compressor.inputs.stagnation_pressure     = inlet_nozzle.outputs.stagnation_pressure
        
    #Flow through the low pressure compressor
    low_pressure_compressor(conditions)
    
    print '-------------------------------------------'
    print 'LOW PRESSURE COMPRESSOR '
    print 'Pressure    :', low_pressure_compressor.outputs.stagnation_pressure
    print 'Temperature :', low_pressure_compressor.outputs.stagnation_temperature
            
    #link the high pressure compressor to the low pressure compressor
    high_pressure_compressor.inputs.stagnation_temperature = low_pressure_compressor.outputs.stagnation_temperature
    high_pressure_compressor.inputs.stagnation_pressure    = low_pressure_compressor.outputs.stagnation_pressure
            
    #Flow through the high pressure compressor
    high_pressure_compressor(conditions)
    
    print '-------------------------------------------'
    print 'HIGH PRESSURE COMPRESSOR '
    print 'Pressure    :', high_pressure_compressor.outputs.stagnation_pressure
    print 'Temperature :', high_pressure_compressor.outputs.stagnation_temperature
            
    #link the combustor to the high pressure compressor
    combustor.inputs.stagnation_temperature                = high_pressure_compressor.outputs.stagnation_temperature
    combustor.inputs.stagnation_pressure                   = high_pressure_compressor.outputs.stagnation_pressure
        
    #flow through the high pressor comprresor
    combustor(conditions)
    
    print '-------------------------------------------'
    print 'COMBUSTOR 1 '
    print 'f 1 : ', combustor.outputs.fuel_to_air_ratio
            
    #link the high pressure turbine to the combustor
    high_pressure_turbine.inputs.stagnation_temperature    = combustor.outputs.stagnation_temperature
    high_pressure_turbine.inputs.stagnation_pressure       = combustor.outputs.stagnation_pressure
    high_pressure_turbine.inputs.fuel_to_air_ratio         = combustor.outputs.fuel_to_air_ratio
    
    #link the high pressuer turbine to the high pressure compressor
    high_pressure_turbine.inputs.compressor                = high_pressure_compressor.outputs
    
    #flow through the high pressure turbine

    high_pressure_turbine.inputs.bypass_ratio = 0.0
    high_pressure_turbine.inputs.fan = Data()
    high_pressure_turbine.inputs.fan.work_done = 0.0
    high_pressure_turbine(conditions)
    
    print '-------------------------------------------'
    print 'HIGH PRESSURE TURBINE '
    print 'Pressure    :', high_pressure_turbine.outputs.stagnation_pressure
    print 'Temperature :', high_pressure_turbine.outputs.stagnation_temperature
                    
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
    low_pressure_turbine.inputs.bypass_ratio = 0.0
    low_pressure_turbine.inputs.fan = Data()
    low_pressure_turbine.inputs.fan.work_done = 0.0    
    low_pressure_turbine(conditions)
    
    #link the combustor to the network   
    #-- Turbojet operation
    if i_tj :
        core_nozzle.inputs.stagnation_temperature       = low_pressure_turbine.outputs.stagnation_temperature
        core_nozzle.inputs.stagnation_pressure          = low_pressure_turbine.outputs.stagnation_pressure
    else :
        #-- Ramjet operation
        combustor_2.inputs.stagnation_temperature       = inlet_nozzle.outputs.stagnation_temperature
        combustor_2.inputs.stagnation_pressure          = inlet_nozzle.outputs.stagnation_pressure
        combustor_2(conditions)
        
        core_nozzle.inputs.stagnation_temperature       = combustor_2.outputs.stagnation_temperature
        core_nozzle.inputs.stagnation_pressure          = combustor_2.outputs.stagnation_pressure
       
    core_nozzle(conditions)  



    print '-------------------------------------------'
    print 'LOW PRESSURE TURBINE '
    print 'Pressure    :', low_pressure_turbine.outputs.stagnation_pressure
    print 'Temperature :', low_pressure_turbine.outputs.stagnation_temperature
    
    print '-------------------------------------------'
    print 'SECOND COMBUSTOR '
    print 'f 2 : ', combustor_2.outputs.fuel_to_air_ratio
    

  
    print '-------------------------------------------'
    print 'NOZZLE '
    print 'Pressure    :', core_nozzle.outputs.stagnation_pressure
    print 'Temperature :', core_nozzle.outputs.stagnation_temperature
    
    #link the thrust component to the core nozzle
    thrust.inputs.core_exit_velocity                       = core_nozzle.outputs.velocity
    thrust.inputs.core_area_ratio                          = core_nozzle.outputs.area_ratio
    thrust.inputs.core_nozzle                              = core_nozzle.outputs
        
    #link the thrust component to the combustor
    print 'F 1', combustor.outputs.fuel_to_air_ratio
    print 'F 2', combustor_2.outputs.fuel_to_air_ratio
    thrust.inputs.fuel_to_air_ratio                        = combustor.outputs.fuel_to_air_ratio + combustor_2.outputs.fuel_to_air_ratio
     
    #link the thrust component to the low pressure compressor 
    #-- Turbojet mode
    thrust.inputs.total_temperature_reference             = inlet_nozzle.outputs.stagnation_temperature
    thrust.inputs.total_pressure_reference                = inlet_nozzle.outputs.stagnation_pressure



    thrust.inputs.number_of_engines                        = number_of_engines
    thrust.inputs.fan_nozzle = Data()
    thrust.inputs.fan_nozzle.velocity = 0.0
    thrust.inputs.fan_nozzle.area_ratio = 0.0
    thrust.inputs.fan_nozzle.static_pressure = 0.0
    thrust.inputs.bypass_ratio = 0.0
    thrust.inputs.flow_through_core                        = 1./(1.+bypass_ratio) #scaled constant to turn on core thrust computation
    thrust.inputs.flow_through_fan                         = bypass_ratio/(1.+bypass_ratio) #scaled constant to turn on fan thrust computation        

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