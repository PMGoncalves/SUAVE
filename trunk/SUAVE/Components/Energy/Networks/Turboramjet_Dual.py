## @ingroup Components-Energy-Networks
#Turbofan.py
# 
# Created:  Oct 2014, A. Variyar, 
# Modified: Feb 2016, M. Vegh
#           Jul 2017, M. Clarke
#           Aug 2017, E. Botero

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# suave imports
import SUAVE

# package imports
import numpy as np

from SUAVE.Core import Data
from SUAVE.Components.Propulsors.Propulsor import Propulsor

# ----------------------------------------------------------------------
#  Turbofan Network
# ----------------------------------------------------------------------

## @ingroup Components-Energy-Networks
class Turboramjet(Propulsor):
    """ This is a turbofan. 
    
        Assumptions:
        None
        
        Source:
        Most of the componentes come from this book:
        https://web.stanford.edu/~cantwell/AA283_Course_Material/AA283_Course_Notes/
    """      
    
    def __defaults__(self):
        """ This sets the default values for the network to function.
    
            Assumptions:
            None
    
            Source:
            N/A
    
            Inputs:
            None
    
            Outputs:
            None
    
            Properties Used:
            N/A
        """           
        
        #setting the default values
        self.tag = 'Turboramjet'
        self.number_of_engines = 1.0
        self.nacelle_diameter  = 1.0
        self.engine_length     = 1.0
        self.bypass_ratio      = 0.0
        
        #areas needed for drag; not in there yet
        self.areas             = Data()
        self.areas.wetted      = 0.0
        self.areas.maximum     = 0.0
        self.areas.exit        = 0.0
        self.areas.inflow      = 0.0
    _component_root_map = None
        
    # linking the different network components
    def evaluate_thrust(self,state):
        """ Calculate thrust given the current state of the vehicle
    
            Assumptions:
            None
    
            Source:
            N/A
    
            Inputs:
            state [state()]
    
            Outputs:
            results.thrust_force_vector [newtons]
            results.vehicle_mass_rate   [kg/s]
            conditions.propulsion.acoustic_outputs:
                core:
                    exit_static_temperature      
                    exit_static_pressure       
                    exit_stagnation_temperature 
                    exit_stagnation_pressure
                    exit_velocity 
                fan:
                    exit_static_temperature      
                    exit_static_pressure       
                    exit_stagnation_temperature 
                    exit_stagnation_pressure
                    exit_velocity 
    
            Properties Used:
            Defaulted values
        """           

        #Unpack
        conditions = state.conditions
        
        # Shared components
        ram                       = self.ram
        inlet_nozzle              = self.inlet_nozzle
        combustor_2               = self.combustor_2
        core_nozzle               = self.core_nozzle

        
        #-- Turbojet components
        low_pressure_compressor   = self.low_pressure_compressor
        high_pressure_compressor  = self.high_pressure_compressor
        combustor                 = self.combustor
        high_pressure_turbine     = self.high_pressure_turbine
        low_pressure_turbine      = self.low_pressure_turbine

        
        thrust                    = self.thrust
        bypass_ratio              = self.bypass_ratio
        number_of_engines         = self.number_of_engines  
        mach_separation           = self.mach_separation
        
        #Creating the network by manually linking the different components
        
        #set the working fluid to determine the fluid properties
        ram.inputs.working_fluid                               = self.working_fluid
        
        #Flow through the ram , this computes the necessary flow quantities and stores it into conditions
        ram(conditions) 
        
        #link inlet nozzle to ram 
        inlet_nozzle.inputs.stagnation_temperature             = ram.outputs.stagnation_temperature 
        inlet_nozzle.inputs.stagnation_pressure                = ram.outputs.stagnation_pressure
        
        #Flow through the inlet nozzle
        inlet_nozzle(conditions)
        
        # Bypass system
        Mo      = conditions.freestream.mach_number
        i_tj    = Mo < mach_separation
        i_rj    = Mo > mach_separation
    
            
        #----------------------
        # TURBOJET OPERATION
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
        low_pressure_turbine(conditions) 
        
        #----------------------
        # SEPARATION OF OPERATIONS
        #----------------------
        
        # Initalize arrays
        mixed_temperature   = 0.0*Mo/Mo
        mixed_pressure      = 0.0*Mo/Mo
        mach_number         = 0.1*Mo/Mo
        
        # Combustor
        #-- Turbojet operation
        mixed_temperature[i_tj] = low_pressure_turbine.outputs.stagnation_temperature[i_tj] 
        mixed_pressure[i_tj]    = low_pressure_turbine.outputs.stagnation_pressure[i_tj] 
        
        
        #-- Ramjet operation
        mixed_temperature[i_rj] = inlet_nozzle.outputs.stagnation_temperature[i_rj]
        mixed_pressure[i_rj]    = inlet_nozzle.outputs.stagnation_pressure[i_rj]
        mach_number[i_rj]       = inlet_nozzle.outputs.mach_number[i_rj]
        
        #-- link the combustor to the correct stagnation properties
        combustor_2.inputs.stagnation_temperature   = mixed_temperature
        combustor_2.inputs.stagnation_pressure      = mixed_pressure
        #combustor_2.inputs.mach_number              = mach_number
        
        # flow through combustor
        combustor_2(conditions)
        
        # Fuel-to-air ratio adjustments
        final_f1        = combustor.outputs.fuel_to_air_ratio
        final_f2        = combustor_2.outputs.fuel_to_air_ratio
        
        #-- If ramjet mode, no fuel used for turbojet mode
        final_f1[i_rj]  = 0.0*Mo[i_rj]/Mo[i_rj]
        
        #-- If turbojet mode, no fuel used for ramjet mode
        final_f2[i_tj]  = 0.0*Mo[i_tj]/Mo[i_tj]
        
        #-- Corrected fuel-to-air ratios of both combustors
        combustor.outputs.fuel_to_air_ratio     = final_f1
        combustor_2.outputs.fuel_to_air_ratio   = final_f2
            
        # Nozzle
        #-- Turbojet operation
        mixed_temperature[i_tj] = low_pressure_turbine.outputs.stagnation_temperature[i_tj] 
        mixed_pressure[i_tj]    = low_pressure_turbine.outputs.stagnation_pressure[i_tj] 

        #-- Ramjet operation
        mixed_pressure[i_rj]    = combustor_2.outputs.stagnation_pressure[i_rj]
        mixed_temperature[i_rj] = combustor_2.outputs.stagnation_temperature[i_rj]

        #-- link the combustor to the correct stagnation temperatures
        core_nozzle.inputs.stagnation_temperature   = mixed_temperature
        core_nozzle.inputs.stagnation_pressure      = mixed_pressure
        
        #flow through nozzle
        core_nozzle(conditions)
        
#        print 'TURBOJET ', i_tj
#        print '-------------------------------------------'
#        print 'INLET '
#        print 'Pressure    :', inlet_nozzle.outputs.stagnation_pressure
#        print 'Temperature :', inlet_nozzle.outputs.stagnation_temperature
#        print '-------------------------------------------'
#        print 'LOW PRESSURE COMPRESSOR '
#        print 'Pressure    :', low_pressure_compressor.outputs.stagnation_pressure
#        print 'Temperature :', low_pressure_compressor.outputs.stagnation_temperature         
#        print '-------------------------------------------'
#        print 'HIGH PRESSURE COMPRESSOR '
#        print 'Pressure    :', high_pressure_compressor.outputs.stagnation_pressure
#        print 'Temperature :', high_pressure_compressor.outputs.stagnation_temperature       
#        print '-------------------------------------------'
#        print 'COMBUSTOR 1 '
#        print 'f 1 : ', combustor.outputs.fuel_to_air_ratio
#        print '-------------------------------------------'
#        print 'HIGH PRESSURE TURBINE '
#        print 'Pressure    :', high_pressure_turbine.outputs.stagnation_pressure
#        print 'Temperature :', high_pressure_turbine.outputs.stagnation_temperature
#                      
#        print '-------------------------------------------'
#        print 'LOW PRESSURE TURBINE '
#        print 'Pressure    :', low_pressure_turbine.outputs.stagnation_pressure
#        print 'Temperature :', low_pressure_turbine.outputs.stagnation_temperature
#    
#        print '-------------------------------------------'
#        print 'SECOND COMBUSTOR '
#        print 'f 2 : ', combustor_2.outputs.fuel_to_air_ratio
#        print 'Pressure IN :', combustor_2.inputs.stagnation_pressure
#        print 'TemperatureIN:', combustor_2.inputs.stagnation_temperature
#    
#
#        print '-------------------------------------------'
#        print 'NOZZLE '
#        print 'Pressure    :', core_nozzle.outputs.stagnation_pressure
#        print 'Temperature :', core_nozzle.outputs.stagnation_temperature
#        print '-------------------------------------------'
  
        #link the thrust component to the combustor 
        thrust.inputs.fuel_to_air_ratio                        = combustor.outputs.fuel_to_air_ratio + combustor_2.outputs.fuel_to_air_ratio
      
        # link the thrust componet to the correct referece pressures
        # Turbojet link
        mixed_temperature[i_tj]                                = low_pressure_compressor.outputs.stagnation_temperature[i_tj]
        mixed_pressure[i_tj]                                   = low_pressure_compressor.outputs.stagnation_pressure[i_tj]
        
        # Ramjet link
        mixed_temperature[i_rj]                                = combustor_2.outputs.stagnation_temperature[i_rj]
        mixed_pressure[i_rj]                                   = combustor_2.outputs.stagnation_pressure[i_rj]
        
        thrust.inputs.total_temperature_reference              = mixed_temperature
        thrust.inputs.total_pressure_reference                 = mixed_pressure

        thrust.inputs.fuel_to_air_ratio                        = combustor.outputs.fuel_to_air_ratio + combustor_2.outputs.fuel_to_air_ratio


        #link the thrust component to the core nozzle
        thrust.inputs.core_exit_velocity                       = core_nozzle.outputs.velocity
        thrust.inputs.core_area_ratio                          = core_nozzle.outputs.area_ratio
        thrust.inputs.core_nozzle                              = core_nozzle.outputs
             
        #link the thrust component to the low pressure compressor 
        #-- Turbojet mode
        thrust.inputs.total_temperature_reference             = inlet_nozzle.outputs.stagnation_temperature
        thrust.inputs.total_pressure_reference                 = inlet_nozzle.outputs.stagnation_pressure
        thrust.inputs.number_of_engines                        = number_of_engines
        thrust.inputs.bypass_ratio                             = bypass_ratio
        thrust.inputs.flow_through_core                        = 1.0 #scaled constant to turn on core thrust computation
        thrust.inputs.flow_through_fan                         = 0.0  #scaled constant to turn on fan thrust computation        

        #compute the thrust
        thrust(conditions)

        #getting the network outputs from the thrust outputs
        F            = thrust.outputs.thrust*[1,0,0]
        mdot         = thrust.outputs.fuel_flow_rate
        output_power = thrust.outputs.power
        F_vec        = conditions.ones_row(3) * 0.0
        F_vec[:,0]   = F[:,0]
        F            = F_vec
        results = Data()
        results.thrust_force_vector = F
        results.vehicle_mass_rate   = mdot
        
        # store data
        results_conditions = Data
        conditions.propulsion.acoustic_outputs.core = results_conditions(
        exit_static_temperature             = core_nozzle.outputs.static_temperature,
        exit_static_pressure                = core_nozzle.outputs.static_pressure,
        exit_stagnation_temperature         = core_nozzle.outputs.stagnation_temperature,
        exit_stagnation_pressure            = core_nozzle.outputs.static_pressure,
        exit_velocity                       = core_nozzle.outputs.velocity
        )
        
        
        return results
    
    def size(self,state):  
        """ Size the turbofan
    
            Assumptions:
            None
    
            Source:
            N/A
    
            Inputs:
            State [state()]
    
            Outputs:
            None
    
            Properties Used:
            N/A
        """             
        
        #Unpack
        conditions = state.conditions
        
        # Shared components
        ram                       = self.ram
        inlet_nozzle              = self.inlet_nozzle
        combustor_2               = self.combustor
        core_nozzle               = self.core_nozzle

        
        #-- Turbojet components
        low_pressure_compressor   = self.low_pressure_compressor
        high_pressure_compressor  = self.high_pressure_compressor
        combustor                 = self.combustor
        high_pressure_turbine     = self.high_pressure_turbine
        low_pressure_turbine      = self.low_pressure_turbine

        
        thrust                    = self.thrust
        bypass_ratio              = self.bypass_ratio
        number_of_engines         = self.number_of_engines   
        mach_separation           = self.mach_separation
        
        #Creating the network by manually linking the different components
        
        #set the working fluid to determine the fluid properties
        ram.inputs.working_fluid                               = self.working_fluid
        
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
        i_tj    = Mo < mach_separation
    
        #-- Defines ramjet operation
        i_rj    = Mo > mach_separation
    
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
    
        #-- Ramjet operation
        else :
            # Link the second combustor to the network
            combustor_2.inputs.stagnation_temperature       = inlet_nozzle.outputs.stagnation_temperature
            combustor_2.inputs.stagnation_pressure          = inlet_nozzle.outputs.stagnation_pressure
        
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
            thrust.inputs.total_temperature_reference          = inlet_nozzle.outputs.stagnation_temperature
            thrust.inputs.total_pressure_reference             = inlet_nozzle.outputs.stagnation_pressure

        else:
            # if ramjet mode only, disregard fuel calculations for first combustor
            combustor.outputs.fuel_to_air_ratio                = 0.0
            thrust.inputs.total_temperature_reference          = combustor_2.outputs.stagnation_temperature
            thrust.inputs.total_pressure_reference             = combustor_2.outputs.stagnation_pressure
    
        #link the thrust component to the combustor
        thrust.inputs.fuel_to_air_ratio                        = combustor.outputs.fuel_to_air_ratio + combustor_2.outputs.fuel_to_air_ratio

        thrust.inputs.number_of_engines                        = number_of_engines
        thrust.inputs.bypass_ratio                             = bypass_ratio
        thrust.inputs.flow_through_core                        = 1.0 #scaled constant to turn on core thrust computation
        thrust.inputs.flow_through_fan                         = 0.0 #scaled constant to turn on fan thrust computation        

        #compute the thrust
        thrust.size(conditions)
        
        #getting the network outputs from the thrust outputs
        F            = thrust.outputs.thrust*[1,0,0]
        mdot         = thrust.outputs.fuel_flow_rate
        F_vec        = conditions.ones_row(3) * 0.0
        F_vec[:,0]   = F[:,0]
        F            = F_vec

        results = Data()
        results.thrust_force_vector = F
        results.vehicle_mass_rate   = mdot
        return results
        
        
    __call__ = evaluate_thrust