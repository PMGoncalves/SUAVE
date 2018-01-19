## @ingroup Components-Energy-Networks
# Ramjet.py
# 
# Created:  Jun 2017, P. Goncalves

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# suave imports
import SUAVE
import numpy as np


from SUAVE.Core import Data, Units
from SUAVE.Components.Propulsors.Propulsor import Propulsor

# ----------------------------------------------------------------------
#  Ramjet Network
# ----------------------------------------------------------------------

## @ingroup Components-Energy-Networks
class Turbine_Based_Combined_Cycle(Propulsor):
    """ This is a ramjet for supersonic flight.

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
        self.tag = 'TBCC'
        self.number_of_engines = 1.0
        self.nacelle_diameter  = 1.0
        self.engine_length     = 1.0
        self.prop_mode         = 0.0 *conditions.propulsion.throttle 
    
    _component_root_map = None
        

#    def evaluate_thrust(self,state):
#	""" Calculate thrust given the current state of the vehicle
#    
#		Assumptions:
#		None
#    
#		Source:
#		N/A
#    
#		Inputs:
#		state [state()]
#    
#		Outputs:
#		results.thrust_force_vector                   [newtons]
#		results.vehicle_mass_rate                     [kg/s]
#		conditions.propulsion.acoustic_outputs:
#		    core:
#			exit_static_temperature                  [K] 
#			exit_static_pressure                     [K] 
#			exit_stagnation_temperature              [K] 
#			exit_stagnation_pressure                 [Pa] 
#			exit_velocity                            [m/s] 
#		    fan:
#			exit_static_temperature                  [K]  
#			exit_static_pressure                     [K] 
#			exit_stagnation_temperature              [K] 
#			exit_stagnation_pressure                 [Pa] 
#			exit_velocity                            [m/s] 
#    
#		Properties Used:
#		Defaulted values
#	    """  	
#
#        #Unpack
#        conditions = state.conditions
#        ###################################################
#        # Trbie
#        
#        
#        ###################################################
#        # Ramjet & Scramjet
#        ram                       = self.ram
#        inlet_nozzle              = self.inlet_nozzle
#        combustor                 = self.combustor
#        core_nozzle               = self.core_nozzle
#        thrust                    = self.thrust
#        
#        
#        number_of_engines         = self.number_of_engines
#        throttle                  = conditions.propulsion.throttle  
#        
#        # Vector init
#        stag_temperature_vector   = 0.0  * conditions.propulsion.throttle  
#        stag_pressure_vector      = 0.0  * conditions.propulsion.throttle  
#        stat_temperature_vector   = 0.0  * conditions.propulsion.throttle    
#        stat_pressure_vector      = 0.0  * conditions.propulsion.throttle  
#        velocity_vector           = 0.0  * conditions.propulsion.throttle 
#        mach_number_vector        = 0.0  * conditions.propulsion.throttle 
#        fuel_to_air_vector        = 0.0  * conditions.propulsion.throttle  
#        area_ratio_vector         = 0.0  * conditions.propulsion.throttle          
#        
#        
#        # Activate right propulsion
#        prop_mode                 = 0.0 
#            # turbojet  = 0.0
#            # ramjet    = 1.0
#            # scramet   = 2.0
#                
#        i_max_throttle            = throttle >= 1.0
#        i_max_mode                = prop_mode < 2.0
#        
#        prop_switch               = np.logical_and(i_max_mode,i_max_throttle)                        
#        
#        prop_mode[prop_switch]    = prop_mode[prop_switch] + 1
#          
#        #i_turbojet                = prop_mode == 0.0
#        i_ramjet                  = prop_mode == 1.0
#        i_scramjet                = prop_mode == 2.0
#        
#        #Creating the network by manually linking the different components
#        
#        #set the working fluid to determine the fluid properties
#        ram.inputs.working_fluid                               = self.working_fluid
#        
#        #Flow through the ram , this computes the necessary flow quantities and stores it into conditions
#        ram(conditions)
#        
#        #link inlet nozzle to ram 
#        inlet_nozzle.inputs.stagnation_temperature             = ram.outputs.stagnation_temperature 
#        inlet_nozzle.inputs.stagnation_pressure                = ram.outputs.stagnation_pressure
#        
#        
#        ###########################################
#        ################## inlet nozzle
#        ######
#        
#        #Flow through the inlet nozzle
#        
#            #   Ramjet 
#            #   run the conditions and retrieve flow properties under ramjet mode
#            
#        inlet_nozzle(conditions)
#        stag_temperature_vector[i_ramjet]                      = inlet_nozzle.outputs.stagnation_temperature[i_ramjet]
#        stag_pressure_vector[i_ramjet]                         = inlet_nozzle.outputs.stagnation_pressure[i_ramjet]
#        mach_number_vector[i_ramjet]                           = inlet_nozzle.outputs.mach_number[i_ramjet]
#        
#            ###########################################
#            
#            #   Scramjet 
#            #   run the conditions and retrieve flow properties under scramjet mode
#        inlet_nozzle.compute_scramjet(conditions)
#        stag_temperature_vector[i_scramjet]                    = inlet_nozzle.outputs.stagnation_temperature[i_scramjet]
#        stag_pressure_vector[i_scramjet]                       = inlet_nozzle.outputs.stagnation_pressure[i_scramjet]
#        
#        ###########################################
#        ################## combustor
#        ######
#        
#        #link the combustor to the inlet nozzle
#        combustor.inputs.stagnation_temperature                = stag_temperature_vector
#        combustor.inputs.stagnation_pressure                   = stag_pressure_vector
#        combustor.inputs.mach_number                           = mach_number_vector
#        
#       
#        #Flow through the combustor
#        
#            #   Ramjet 
#            #   run the conditions and retrieve flow properties under ramjet mode
#        combustor.compute_rayleigh(conditions)
#        stag_temperature_vector[i_ramjet]                      = inlet_nozzle.outputs.stagnation_temperature[i_ramjet]
#        stag_pressure_vector[i_ramjet]                         = inlet_nozzle.outputs.stagnation_pressure[i_ramjet]
#        fuel_to_air_vector[i_ramjet]                           = combustor.outputs.fuel_to_air_ratio[i_ramjet]
#        
#            ###########################################
#            
#            #   Scramjet 
#            #   run the conditions and retrieve flow properties under scramjet mode
#        combustor.compute_scramjet(conditions)
#        stag_temperature_vector[i_scramjet]                    = combustor.outputs.stagnation_temperature[i_scramjet]
#        stag_pressure_vector[i_scramjet]                       = combustor.outputs.stagnation_pressure[i_scramjet]
#        stat_temperature_vector[i_scramjet]                    = combustor.outputs.static_pressure[i_scramjet]        
#        stat_pressure_vector[i_scramjet]                       = combustor.outputs.static_temperature[i_scramjet] 
#        velocity_vector[i_scramjet]                            = combustor.outputs.velocity[i_scramjet]
#        fuel_to_air_vector[i_scramjet]                         = combustor.outputs.fuel_to_air_ratio[i_scramjet]
#  
#      
#        ###########################################
#        ################## core nozzle
#        ######
#
#        #link the combustor to the core nozzle
#        core_nozzle.inputs.stagnation_temperature              = stag_temperature_vector
#        core_nozzle.inputs.stagnation_pressure                 = stag_pressure_vector
#        core_nozzle.inputs.static_temperature                  = stat_temperature_vector
#        core_nozzle.inputs.static_pressure                     = stat_pressure_vector
#        core_nozzle.inputs.velocity                            = velocity_vector
#        core_nozzle.inputs.fuel_to_air_ratio                   = fuel_to_air_vector
#
#        #flow through the core nozzle
#        
#            #   Ramjet 
#            #   run the conditions and retrieve flow properties under ramjet mode
#        core_nozzle.compute_limited_geometry(conditions)        
#        velocity_vector[i_ramjet]                              = core_nozzle.outputs.velocity[i_ramjet]
#        area_ratio_vector[i_ramjet]                            = core_nozzle.outputs.area_ratio[i_ramjet]
#        
#            ###########################################
#            
#        # compute the thrust using the thrust component
#        #link the thrust component to the core nozzle
#        thrust.inputs.core_exit_velocity                       = core_nozzle.outputs.velocity
#        thrust.inputs.core_area_ratio                          = core_nozzle.outputs.area_ratio
#        thrust.inputs.core_nozzle                              = core_nozzle.outputs
#	
#        #link the thrust component to the combustor
#        thrust.inputs.fuel_to_air_ratio                        = combustor.outputs.fuel_to_air_ratio
#	
#        #link the thrust component to the low pressure compressor 
#        thrust.inputs.stag_temp_lpt_exit                       = core_nozzle.outputs.stagnation_temperature
#        thrust.inputs.stag_press_lpt_exit                      = core_nozzle.outputs.stagnation_pressure
#        thrust.inputs.number_of_engines                        = number_of_engines
#	thrust.inputs.flow_through_core                        =  1.0 #scaled constant to turn on core thrust computation
#	thrust.inputs.flow_through_fan                         =  0.0 #scaled constant to turn on fan thrust computation        
#
#        #compute the thrust
#        thrust(conditions)
#
#        
#        #flow through the core nozzle
#        
#            #   Scramjet 
#            #   run the conditions and retrieve flow properties under scramjet mode
#        core_nozzle.compute_scramjet(conditions)
#        stat_temperature_vector[i_scramjet]                    = core_nozzle.outputs.temperature [i_scramjet]
#        stat_pressure_vector[i_scramjet]                       = core_nozzle.outputs.pressure[i_scramjet]        
#        velocity_vector[i_scramjet]                            = core_nozzle.outputs.velocity[i_scramjet]
#        
#        
#        
#        #link the thrust component to the core nozzle
#        thrust.inputs.core_exit_pressure                       = core_nozzle.outputs.pressure
#        thrust.inputs.core_exit_temperature                    = core_nozzle.outputs.temperature 
#        thrust.inputs.core_exit_velocity                       = core_nozzle.outputs.velocity
#        thrust.inputs.core_area_ratio                          = core_nozzle.outputs.area_ratio
#        thrust.inputs.core_nozzle                              = core_nozzle.outputs
#
#        
#        
#
#        core_nozzle.compute_limited_geometry(conditions)
#
#        # compute the thrust using the thrust component
#        #link the thrust component to the core nozzle
#        thrust.inputs.core_exit_velocity                       = core_nozzle.outputs.velocity
#        thrust.inputs.core_area_ratio                          = core_nozzle.outputs.area_ratio
#        thrust.inputs.core_nozzle                              = core_nozzle.outputs
#	
#        #link the thrust component to the combustor
#        thrust.inputs.fuel_to_air_ratio                        = combustor.outputs.fuel_to_air_ratio
#	
#        #link the thrust component to the low pressure compressor 
#        thrust.inputs.stag_temp_lpt_exit                       = core_nozzle.outputs.stagnation_temperature
#        thrust.inputs.stag_press_lpt_exit                      = core_nozzle.outputs.stagnation_pressure
#        thrust.inputs.number_of_engines                        = number_of_engines
#	thrust.inputs.flow_through_core                        =  1.0 #scaled constant to turn on core thrust computation
#	thrust.inputs.flow_through_fan                         =  0.0 #scaled constant to turn on fan thrust computation        
#
#        #compute the thrust
#        thrust(conditions)
#        
#        #getting the network outputs from the thrust outputs
#        F            = thrust.outputs.thrust*[1,0,0]
#        mdot         = thrust.outputs.fuel_flow_rate
#        Isp          = thrust.outputs.specific_impulse
#        output_power = thrust.outputs.power
#        F_vec        = conditions.ones_row(3) * 0.0
#        F_vec[:,0]   = F[:,0]
#        F            = F_vec
#
#        results = Data()
#        results.thrust_force_vector = F
#        results.vehicle_mass_rate   = mdot
#        
#        return results



    def evaluate_thrust(self,state):
	""" Calculate thrust given the current state of the vehicle
    
		Assumptions:
		None
    
		Source:
		N/A
    
		Inputs:
		state [state()]
    
		Outputs:
		results.thrust_force_vector                   [newtons]
		results.vehicle_mass_rate                     [kg/s]
		conditions.propulsion.acoustic_outputs:
		    core:
			exit_static_temperature                  [K] 
			exit_static_pressure                     [K] 
			exit_stagnation_temperature              [K] 
			exit_stagnation_pressure                 [Pa] 
			exit_velocity                            [m/s] 
		    fan:
			exit_static_temperature                  [K]  
			exit_static_pressure                     [K] 
			exit_stagnation_temperature              [K] 
			exit_stagnation_pressure                 [Pa] 
			exit_velocity                            [m/s] 
    
		Properties Used:
		Defaulted values
	    """  	

        #Unpack
        conditions = state.conditions

           
        ###################################################
        # Ramjet & Scramjet
        ram                       = self.ram
        inlet_nozzle              = self.inlet_nozzle
        combustor                 = self.combustor
        core_nozzle               = self.core_nozzle
        thrust_rj                 = self.thrust1
        thrust_sj                 = self.thrust2
        
        number_of_engines         = self.number_of_engines
        throttle                  = conditions.propulsion.throttle  
        
        # Vector init
        fuel_to_air_vector        = 0.0  * conditions.propulsion.throttle  
        
        
        # Activate right propulsion
        prop_mode                 = 0.0 
            # turbojet  = 0.0
            # ramjet    = 1.0
            # scramet   = 2.0
                
        i_max_throttle            = throttle >= 1.0
        i_max_mode                = prop_mode < 2.0
        
        prop_switch               = np.logical_and(i_max_mode,i_max_throttle)                        
        
        prop_mode[prop_switch]    = prop_mode[prop_switch] + 1
          
        #i_turbojet                = prop_mode == 0.0
        i_ramjet                  = prop_mode == 1.0
        i_scramjet                = prop_mode == 2.0
        
        
        ###########################################
        #####################################
        ################# RAMJET NETWORK
        #####################################
        ###########################################
        
        #Creating the network by manually linking the different components
        
        #set the working fluid to determine the fluid properties
        ram.inputs.working_fluid                               = self.working_fluid
        
        #Flow through the ram , this computes the necessary flow quantities and stores it into conditions
        ram(conditions)
        
        #link inlet nozzle to ram 
        inlet_nozzle.inputs.stagnation_temperature             = ram.outputs.stagnation_temperature 
        inlet_nozzle.inputs.stagnation_pressure                = ram.outputs.stagnation_pressure
        
        #Flow through the inlet nozzle
        inlet_nozzle.compute(conditions)
        
        #link the combustor to the inlet nozzle
        combustor.inputs.stagnation_temperature                = inlet_nozzle.outputs.stagnation_temperature
        combustor.inputs.stagnation_pressure                   = inlet_nozzle.outputs.stagnation_pressure
        combustor.inputs.mach_number                           = inlet_nozzle.outputs.mach_number
        
        #flow through the high pressor comprresor
        combustor.compute_rayleigh(conditions)
        
        
        #link the core nozzle to the low pressure turbine
        core_nozzle.inputs.stagnation_temperature              = combustor.outputs.stagnation_temperature
        core_nozzle.inputs.stagnation_pressure                 = combustor.outputs.stagnation_pressure
        
        #flow through the core nozzle
        core_nozzle.compute_limited_geometry(conditions)

        # compute the thrust using the thrust component
        #link the thrust component to the core nozzle
        thrust.inputs.core_exit_velocity                       = core_nozzle.outputs.velocity
        thrust.inputs.core_area_ratio                          = core_nozzle.outputs.area_ratio
        thrust.inputs.core_nozzle                              = core_nozzle.outputs
	
        #link the thrust component to the combustor
        thrust.inputs.fuel_to_air_ratio                        = combustor.outputs.fuel_to_air_ratio
	
        #link the thrust component to the low pressure compressor 
        thrust.inputs.stag_temp_lpt_exit                       = core_nozzle.outputs.stagnation_temperature
        thrust.inputs.stag_press_lpt_exit                      = core_nozzle.outputs.stagnation_pressure
        thrust.inputs.number_of_engines                        = number_of_engines
        thrust.inputs.flow_through_core                        =  1.0 #scaled constant to turn on core thrust computation
        thrust.inputs.flow_through_fan                         =  0.0 #scaled constant to turn on fan thrust computation        

        #compute the thrust
        thrust(conditions)
        

        ###########################################
        #####################################
        ################# SCRAMJET NETWORK
        #####################################
        ###########################################
        
        #Creating the network by manually linking the different components
        
        #set the working fluid to determine the fluid properties
        ram.inputs.working_fluid                               = self.working_fluid
        
        #Flow through the ram , this computes the necessary flow quantities and stores it into conditions
        ram(conditions)
        
        #link inlet nozzle to ram 
        inlet_nozzle.inputs.stagnation_temperature             = ram.outputs.stagnation_temperature #conditions.freestream.stagnation_temperature
        inlet_nozzle.inputs.stagnation_pressure                = ram.outputs.stagnation_pressure #conditions.freestream.stagnation_pressure
        
        #Flow through the inlet nozzle
        inlet_nozzle(conditions)
    
        #link the combustor to the high pressure compressor
        combustor.inputs.stagnation_temperature                = inlet_nozzle.outputs.stagnation_temperature
        combustor.inputs.stagnation_pressure                   = inlet_nozzle.outputs.stagnation_pressure
        combustor.inputs.inlet_nozzle                          = inlet_nozzle.outputs
        
        #flow through the high pressor comprresor
        combustor.compute_scramjet(conditions)
        
        #link the core nozzle to the low pressure turbine
        core_nozzle.inputs.stagnation_temperature              = combustor.outputs.stagnation_temperature
        core_nozzle.inputs.stagnation_pressure                 = combustor.outputs.stagnation_pressure
        core_nozzle.inputs.static_temperature                  = combustor.outputs.static_temperature
        core_nozzle.inputs.static_pressure                     = combustor.outputs.static_pressure
        core_nozzle.inputs.velocity                            = combustor.outputs.velocity
        core_nozzle.inputs.fuel_to_air_ratio                   = combustor.outputs.fuel_to_air_ratio
        
        #flow through the core nozzle
        core_nozzle(conditions)
        
        #link the thrust component to the core nozzle
        thrust.inputs.core_exit_pressure                       = core_nozzle.outputs.pressure
        thrust.inputs.core_exit_temperature                    = core_nozzle.outputs.temperature 
        thrust.inputs.core_exit_velocity                       = core_nozzle.outputs.velocity
        thrust.inputs.core_area_ratio                          = core_nozzle.outputs.area_ratio
        thrust.inputs.core_nozzle                              = core_nozzle.outputs
	
        #link the thrust component to the combustor
        thrust.inputs.fuel_to_air_ratio                        = combustor.outputs.fuel_to_air_ratio
	
        #link the thrust component to the low pressure compressor 
        thrust.inputs.stag_temp_lpt_exit                       = core_nozzle.outputs.stagnation_temperature
        thrust.inputs.stag_press_lpt_exit                      = core_nozzle.outputs.stagnation_pressure
        thrust.inputs.number_of_engines                        = number_of_engines
        thrust.inputs.flow_through_core                        =  1.0 #scaled constant to turn on core thrust computation
        thrust.inputs.flow_through_fan                         =  0.0 #scaled constant to turn on fan thrust computation        

        #compute the thrust
        thrust.compute_stream_thrust(conditions)


    def size(self,state):  
        
        """ Size the ramjet
    
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
        
        #Unpack components
        conditions = state.conditions      
        ram                       = self.ram
        inlet_nozzle              = self.inlet_nozzle
        combustor                 = self.combustor
        core_nozzle               = self.core_nozzle
        thrust                    = self.thrust
        
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
          
        
        #link the combustor to the high pressure compressor
        combustor.inputs.stagnation_temperature                = inlet_nozzle.outputs.stagnation_temperature
        combustor.inputs.stagnation_pressure                   = inlet_nozzle.outputs.stagnation_pressure
        combustor.inputs.mach_number                           = inlet_nozzle.outputs.mach_number
        
        #flow through the high pressure compressor
        combustor.compute_rayleigh(conditions)
        
        #link the core nozzle to the low pressure turbine
        core_nozzle.inputs.stagnation_temperature              = combustor.outputs.stagnation_temperature
        core_nozzle.inputs.stagnation_pressure                 = combustor.outputs.stagnation_pressure
        
        #flow through the core nozzle
        core_nozzle(conditions)
        
        # compute the thrust using the thrust component
        #link the thrust component to the core nozzle
        thrust.inputs.core_exit_velocity                       = core_nozzle.outputs.velocity
        thrust.inputs.core_area_ratio                          = core_nozzle.outputs.area_ratio
        thrust.inputs.core_nozzle                              = core_nozzle.outputs
	
        #link the thrust component to the combustor
        thrust.inputs.fuel_to_air_ratio                        = combustor.outputs.fuel_to_air_ratio
	
        #link the thrust component to the low pressure compressor 
        thrust.inputs.stag_temp_lpt_exit                       = core_nozzle.outputs.stagnation_temperature
        thrust.inputs.stag_press_lpt_exit                      = core_nozzle.outputs.stagnation_pressure

        #compute the thrust
        thrust.size(conditions)

    __call__ = evaluate_thrust