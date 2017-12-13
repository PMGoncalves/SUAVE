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
    """ This is a turbojet for supersonic flight.

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
        self.tag = 'Turbojet'
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
        
        thrust                    = self.thrust        
        number_of_engines         = self.number_of_engines        
        #compute the thrust
        thrust.inputs.number_of_engines = number_of_engines
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
        
        #Unpack components
        conditions = state.conditions      
        thrust                    = self.thrust
        
        #compute the thrust
        thrust.size_rocket(conditions)

    __call__ = evaluate_thrust