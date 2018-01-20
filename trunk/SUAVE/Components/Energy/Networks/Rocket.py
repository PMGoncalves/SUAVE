## @ingroup Components-Energy-Networks
# Turbojet_Super.py
# 
# Created:  May 2015, T. MacDonald
# Modified: Aug 2017, E. Botero

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# suave imports
import SUAVE

from SUAVE.Core import Data, Units
from SUAVE.Components.Propulsors.Propulsor import Propulsor

# ----------------------------------------------------------------------
#  Turbojet Network
# ----------------------------------------------------------------------

## @ingroup Components-Energy-Networks
class Rocket(Propulsor):
    """ This is a rocket engine

        Assumptions:
        None

        Source:
        Most of the componentes come from this book:
        
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
        self.tag = 'Rocket'
        self.number_of_engines = 1.0
        self.nacelle_diameter  = 1.0
        self.engine_length     = 1.0
    
    _component_root_map = None
        

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
        
        combustor                 = self.combustor
        nozzle                    = self.nozzle
        thrust                    = self.thrust        
        number_of_engines         = self.number_of_engines        
        
        #Creating the network by manually linking the different components
        combustor.oxidizer_data                               = self.oxidizer
        combustor.fuel_data                                   = self.fuel
          
        #flow through combustor
        combustor.compute_rocket(conditions)

        #link the nozzle to the combustor
        nozzle.inputs.stagnation_temperature                = combustor.outputs.stagnation_temperature
        nozzle.inputs.stagnation_pressure                   = combustor.outputs.stagnation_pressure
        nozzle.inputs.isentropic_expansion_factor           = combustor.outputs.isentropic_expansion_factor
        nozzle.inputs.specific_gas_constant                 = combustor.outputs.specific_gas_constant 
        nozzle.outputs.specific_heat_constant_pressure      = combustor.outputs.specific_heat_constant_pressure    

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
        thrust.compute_rocket(conditions)
        
        #getting the network outputs from the thrust outputs
        F            = thrust.outputs.thrust*[1,0,0]
        mdot         = thrust.outputs.fuel_flow_rate
        Isp          = thrust.outputs.specific_impulse
        output_power = thrust.outputs.power
        F_vec        = conditions.ones_row(3) * 0.0
        F_vec[:,0]   = F[:,0]
        F            = F_vec

        results = Data()
        results.thrust_force_vector = F
        results.vehicle_mass_rate   = mdot
        
        return results
    
    def size(self,state):  
        
        #Unpack
        conditions = state.conditions
        
        combustor                 = self.combustor
        nozzle                    = self.nozzle
        thrust                    = self.thrust        
        number_of_engines         = self.number_of_engines        
        
        #Creating the network by manually linking the different components
        combustor.oxidizer_data                               = self.oxidizer
        combustor.fuel_data                                   = self.fuel
          
        #flow through combustor
        combustor.compute_rocket(conditions)

        #link the nozzle to the combustor
        nozzle.inputs.stagnation_temperature                = combustor.outputs.stagnation_temperature
        nozzle.inputs.stagnation_pressure                   = combustor.outputs.stagnation_pressure
        nozzle.inputs.isentropic_expansion_factor           = combustor.outputs.isentropic_expansion_factor
        nozzle.inputs.specific_gas_constant                 = combustor.outputs.specific_gas_constant 
        nozzle.outputs.specific_heat_constant_pressure      = combustor.outputs.specific_heat_constant_pressure    

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

    __call__ = evaluate_thrust