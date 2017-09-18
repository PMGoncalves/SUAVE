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
        self.bypass_ratio      = 1.0
        
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
        Mo       = conditions.freestream.mach_number
        i_tj = Mo < mach_separation
        i_rj = Mo > mach_separation
    
            
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
        print 'TURBOJET ', i_tj
        print '-------------------------------------------'
        print 'INLET '
        print 'Pressure    :', inlet_nozzle.outputs.stagnation_pressure
        print 'Temperature :', inlet_nozzle.outputs.stagnation_temperature
        print '-------------------------------------------'    
        
        mixed_temperature = 0.0*Mo/Mo
        mixed_pressure    = 0.0*Mo/Mo
        
        # Combustor
        #-- Turbojet operation
        mixed_temperature[i_tj] = low_pressure_turbine.outputs.stagnation_temperature[i_tj] 
        mixed_pressure[i_tj] = low_pressure_turbine.outputs.stagnation_pressure[i_tj] 
        
        #-- Ramjet operation
        mixed_temperature[i_rj] = inlet_nozzle.outputs.stagnation_temperature[i_rj]
        mixed_pressure[i_rj] = inlet_nozzle.outputs.stagnation_pressure[i_rj]
        
        combustor_2.inputs.stagnation_temperature = mixed_temperature
        combustor_2.inputs.stagnation_pressure = mixed_pressure
        
        combustor_2(conditions)
        
        final_f1  = combustor.outputs.fuel_to_air_ratio
        final_f2  = combustor_2.outputs.fuel_to_air_ratio
        final_f1[i_rj] = 0.0*Mo[i_rj]/Mo[i_rj]
        final_f2[i_tj] = 0.0*Mo[i_tj]/Mo[i_tj]
        
        combustor.outputs.fuel_to_air_ratio = final_f1
        combustor_2.outputs.fuel_to_air_ratio = final_f2
        
        mixed_temperature = 0.0*Mo/Mo
        mixed_pressure    = 0.0*Mo/Mo
        
        # Nozzle
        #-- Turbojet operation
        mixed_temperature[i_tj] = low_pressure_turbine.outputs.stagnation_temperature[i_tj] 
        mixed_pressure[i_tj]    = low_pressure_turbine.outputs.stagnation_pressure[i_tj] 

        #-- Ramjet operation
        mixed_pressure[i_rj]    = combustor_2.outputs.stagnation_pressure[i_rj]
        mixed_temperature[i_rj] = combustor_2.outputs.stagnation_temperature[i_rj]

        
        core_nozzle.inputs.stagnation_temperature = mixed_temperature
        core_nozzle.inputs.stagnation_pressure = mixed_pressure
        core_nozzle(conditions)
        
        print 'TURBOJET ', i_tj
        print '-------------------------------------------'
        print 'INLET '
        print 'Pressure    :', inlet_nozzle.outputs.stagnation_pressure
        print 'Temperature :', inlet_nozzle.outputs.stagnation_temperature
        print '-------------------------------------------'
        print 'LOW PRESSURE COMPRESSOR '
        print 'Pressure    :', low_pressure_compressor.outputs.stagnation_pressure
        print 'Temperature :', low_pressure_compressor.outputs.stagnation_temperature         
        print '-------------------------------------------'
        print 'HIGH PRESSURE COMPRESSOR '
        print 'Pressure    :', high_pressure_compressor.outputs.stagnation_pressure
        print 'Temperature :', high_pressure_compressor.outputs.stagnation_temperature       
        print '-------------------------------------------'
        print 'COMBUSTOR 1 '
        print 'f 1 : ', combustor.outputs.fuel_to_air_ratio
        print '-------------------------------------------'
        print 'HIGH PRESSURE TURBINE '
        print 'Pressure    :', high_pressure_turbine.outputs.stagnation_pressure
        print 'Temperature :', high_pressure_turbine.outputs.stagnation_temperature
                      
        print '-------------------------------------------'
        print 'LOW PRESSURE TURBINE '
        print 'Pressure    :', low_pressure_turbine.outputs.stagnation_pressure
        print 'Temperature :', low_pressure_turbine.outputs.stagnation_temperature
    
        print '-------------------------------------------'
        print 'SECOND COMBUSTOR '
        print 'f 2 : ', combustor_2.outputs.fuel_to_air_ratio
        print 'Pressure IN :', combustor_2.inputs.stagnation_pressure
        print 'TemperatureIN:', combustor_2.inputs.stagnation_temperature
    

        print '-------------------------------------------'
        print 'NOZZLE '
        print 'Pressure    :', core_nozzle.outputs.stagnation_pressure
        print 'Temperature :', core_nozzle.outputs.stagnation_temperature
        print '-------------------------------------------'
      
        #link the thrust component to the core nozzle
        thrust.inputs.core_exit_velocity                       = core_nozzle.outputs.velocity
        thrust.inputs.core_area_ratio                          = core_nozzle.outputs.area_ratio
        thrust.inputs.core_nozzle                              = core_nozzle.outputs
        
        #link the thrust component to the combustor
        final_f1  = combustor.outputs.fuel_to_air_ratio
        final_f2  = combustor_2.outputs.fuel_to_air_ratio
        final_f1[i_rj] = 0.0*Mo[i_rj]/Mo[i_rj]
        final_f2[i_tj] = 0.0*Mo[i_tj]/Mo[i_tj]
        
        combustor.outputs.fuel_to_air_ratio = final_f1
        combustor_2.outputs.fuel_to_air_ratio = final_f2
        print 'FINAL F', (combustor.outputs.fuel_to_air_ratio + combustor_2.outputs.fuel_to_air_ratio)
        
        thrust.inputs.fuel_to_air_ratio                        = combustor.outputs.fuel_to_air_ratio + combustor_2.outputs.fuel_to_air_ratio
     
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
        print 'ola estou aqui'
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
        print 'TURBOJET ', i_tj
        print '-------------------------------------------'
        print 'INLET '
        print 'Pressure    :', inlet_nozzle.outputs.stagnation_pressure
        print 'Temperature :', inlet_nozzle.outputs.stagnation_temperature
        print '-------------------------------------------'
        print 'LOW PRESSURE COMPRESSOR '
        print 'Pressure    :', low_pressure_compressor.outputs.stagnation_pressure
        print 'Temperature :', low_pressure_compressor.outputs.stagnation_temperature         
        print '-------------------------------------------'
        print 'HIGH PRESSURE COMPRESSOR '
        print 'Pressure    :', high_pressure_compressor.outputs.stagnation_pressure
        print 'Temperature :', high_pressure_compressor.outputs.stagnation_temperature       
        print '-------------------------------------------'
        print 'COMBUSTOR 1 '
        print 'f 1 : ', combustor.outputs.fuel_to_air_ratio
        print '-------------------------------------------'
        print 'HIGH PRESSURE TURBINE '
        print 'Pressure    :', high_pressure_turbine.outputs.stagnation_pressure
        print 'Temperature :', high_pressure_turbine.outputs.stagnation_temperature
                      
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
        print '-------------------------------------------'
    
        
        #link the thrust component to the core nozzle
        thrust.inputs.core_exit_velocity                       = core_nozzle.outputs.velocity
        thrust.inputs.core_area_ratio                          = core_nozzle.outputs.area_ratio
        thrust.inputs.core_nozzle                              = core_nozzle.outputs
        
        #link the thrust component to the combustor
        thrust.inputs.fuel_to_air_ratio                        = combustor.outputs.fuel_to_air_ratio + combustor_2.outputs.fuel_to_air_ratio
        
        #link the thrust component to the low pressure compressor 
        thrust.inputs.total_temperature_reference              = low_pressure_compressor.outputs.stagnation_temperature
        thrust.inputs.total_pressure_reference                 = low_pressure_compressor.outputs.stagnation_pressure
        thrust.inputs.number_of_engines                        = number_of_engines
        thrust.inputs.bypass_ratio                             = bypass_ratio
        thrust.inputs.flow_through_core                        = 1./(1.+bypass_ratio) #scaled constant to turn on core thrust computation
        thrust.inputs.flow_through_fan                         = bypass_ratio/(1.+bypass_ratio) #scaled constant to turn on fan thrust computation        

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