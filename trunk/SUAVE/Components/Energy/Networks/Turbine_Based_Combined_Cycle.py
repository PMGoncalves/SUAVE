## @ingroup Components-Energy-Networks
# Scramjet.py
#
# Created:  Jun 2017, P. Goncalves
# Modified: Jan 2018, W. Maier

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# suave imports
import SUAVE

from SUAVE.Core import Data, Units
from SUAVE.Components.Propulsors.Propulsor import Propulsor
from SUAVE.Components.Energy.Networks.Dual_Mode_Ramjet import Dual_Mode_Ramjet
from SUAVE.Components.Energy.Networks.Turbojet_Super import Turbojet_Super
from SUAVE.Components.Energy.Networks.Ramjet import Ramjet
from SUAVE.Components.Energy.Networks.Scramjet import Scramjet

import numpy as np


# ----------------------------------------------------------------------
#  Scramjet Network
# ----------------------------------------------------------------------
## @ingroup Components-Energy-Networks
class Turbine_Based_Combined_Cycle(Propulsor):
    """ This is a Dual_Mode_Ramjet for supersonic flight.

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
        self.tag = 'Dual Mode Ramjet'
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
		results.thrust_force_vector                      [newtons]
		results.vehicle_mass_rate                        [kg/s]
		results.specific_impulse                         [s]
		conditions.propulsion.acoustic_outputs:
		    core:
			exit_static_temperature                  [K]
			exit_static_pressure                     [K]
			exit_stagnation_temperature              [K]
			exit_stagnation_pressure                 [Pa]
			exit_velocity                            [m/s]
		    
		Properties Used:
		Defaulted values
	    """

        # unpack
        conditions                = state.conditions
        number_of_engines         = self.number_of_engines   
        
        # Common
        ram                       = self.ram
        inlet_nozzle              = self.inlet_nozzle
             
        # Turbojet
        low_pressure_compressor   = self.low_pressure_compressor
        high_pressure_compressor  = self.high_pressure_compressor
        low_pressure_turbine      = self.low_pressure_turbine
        high_pressure_turbine     = self.high_pressure_turbine
        combustor_2               = self.combustor_2
        core_nozzle_2             = self.core_nozzle_2    
        thrust_2                  = self.thrust_2
        
        # Dual mode ramjet
        combustor                 = self.combustor
        core_nozzle               = self.core_nozzle


        thrust                    = self.thrust

        # creating switch for operation seletion
        
        Mo = conditions.freestream.mach_number
        M_transition = 2.5
        M_transition2 = 5.0
        
        tj_mode = Mo < M_transition
        rj_mode = np.logical_and(Mo >= M_transition, Mo <= M_transition2)
        sj_mode = Mo > M_transition2
        
        # Build turbojet 
                
        turbojet =  SUAVE.Components.Energy.Networks.Turbojet_Super()
        turbojet.number_of_engines = number_of_engines
        turbojet.working_fluid = self.working_fluid
        turbojet.ram = ram
        turbojet.inlet_nozzle = inlet_nozzle
        turbojet.inlet_nozzle.compressibility_effects = False
        turbojet.low_pressure_compressor = low_pressure_compressor
        turbojet.high_pressure_compressor = high_pressure_compressor
        turbojet.combustor = combustor_2
        turbojet.high_pressure_turbine = high_pressure_turbine
        turbojet.low_pressure_turbine = low_pressure_turbine
        turbojet.core_nozzle = core_nozzle_2
        turbojet.thrust = thrust_2
        
        
        turbojet_results = Data()
        turbojet_results = turbojet.evaluate_thrust(state)
        
        # Build ramjet 
                
        ramjet =  SUAVE.Components.Energy.Networks.Ramjet()
        ramjet.number_of_engines = number_of_engines
        ramjet.working_fluid = self.working_fluid
        ramjet.ram = ram
        ramjet.inlet_nozzle = inlet_nozzle
        ramjet.inlet_nozzle.compressibility_effects = True
        ramjet.combustor = combustor
        ramjet.core_nozzle = core_nozzle
        ramjet.thrust = thrust
        
        ramjet_results = Data()
        ramjet_results = ramjet.evaluate_thrust(state)
        
        # Build scramjet
         
        scramjet =  SUAVE.Components.Energy.Networks.Scramjet()
        scramjet.number_of_engines = number_of_engines
        scramjet.working_fluid = self.working_fluid
        scramjet.ram = ram
        scramjet.inlet_nozzle = inlet_nozzle
        scramjet.combustor = combustor
        scramjet.core_nozzle = core_nozzle
        scramjet.thrust = thrust
        
        scramjet_results = Data()
        scramjet_results = scramjet.evaluate_thrust(state)

        #------------------------------
        
        # Initialize arrays
        F_sum        = np.ones_like(Mo)
        mdot_f_sum   = np.ones_like(Mo)
        Isp_sum      = np.ones_like(Mo)
        f_sum        = np.ones_like(Mo)
        TSFC_sum     = np.ones_like(Mo)
        Fsp_sum      = np.ones_like(Mo)
        mdot_c_sum   = np.ones_like(Mo)
        
        #-- Combine different modes
        
        F_sum[tj_mode] = turbojet_results.thrust_force_vector[tj_mode]
        F_sum[rj_mode] = ramjet_results.thrust_force_vector[rj_mode]
        F_sum[sj_mode] = scramjet_results.thrust_force_vector[sj_mode]
        
        F_vec        = conditions.ones_row(3) * [1,0,0]
        F_sum        = F_sum*F_vec
       
        mdot_f_sum[tj_mode] = turbojet_results.vehicle_mass_rate[tj_mode]
        mdot_f_sum[rj_mode] = ramjet_results.vehicle_mass_rate[rj_mode]
        mdot_f_sum[sj_mode] = scramjet_results.vehicle_mass_rate[sj_mode]

        Isp_sum[tj_mode] = turbojet_results.specific_impulse[tj_mode]
        Isp_sum[rj_mode] = ramjet_results.specific_impulse[rj_mode]
        Isp_sum[sj_mode] = scramjet_results.specific_impulse[sj_mode]

        #-- additional debug 

        f_sum[tj_mode] = turbojet_results.f[tj_mode]        
        f_sum[rj_mode] = ramjet_results.f[rj_mode]
        f_sum[sj_mode] = scramjet_results.f[sj_mode]

        TSFC_sum[tj_mode] = turbojet_results.tsfc[tj_mode]         
        TSFC_sum[rj_mode] = ramjet_results.tsfc[rj_mode]
        TSFC_sum[sj_mode] = scramjet_results.tsfc[sj_mode]

        Fsp_sum[tj_mode] = turbojet_results.fsp[tj_mode]        
        Fsp_sum[rj_mode] = ramjet_results.fsp[rj_mode]
        Fsp_sum[sj_mode] = scramjet_results.fsp[sj_mode]
 
#        mdot_c_sum[tj_mode] = turbojet_results.mdot_core[tj_mode]       
#        mdot_c_sum[rj_mode] = ramjet_results.mdot_core[rj_mode]
#        mdot_c_sum[sj_mode] = scramjet_results.mdot_core[sj_mode]


        # Final export
        results_sum                     = Data()
        results_sum.thrust_force_vector = F_sum
        results_sum.vehicle_mass_rate   = mdot_f_sum
        results_sum.specific_impulse    = Isp_sum
        
        #-- new results
        results_sum.fsp                 = Fsp_sum
        results_sum.mdot_core           = mdot_c_sum
        results_sum.tsfc                = TSFC_sum
        results_sum.f                   = f_sum
        
        return results_sum
        #return scramjet_results
        #return ramjet_results


#    def size(self,state):

#        """ Size the dual_mode_ramjet
#
#            Assumptions:
#            None
#
#            Source:
#            N/A
#
#            Inputs:
#            State [state()]
#
#            Outputs:
#            None
#
#            Properties Used:
#            N/A
#        """
#
#        # unpack
#        conditions                = state.conditions
#        number_of_engines         = self.number_of_engines   
#        
#        # Common
#        ram                       = self.ram
#        inlet_nozzle              = self.inlet_nozzle
#             
#        # Turbojet
#        low_pressure_compressor   = self.low_pressure_compressor
#        high_pressure_compressor  = self.high_pressure_compressor
#        low_pressure_turbine      = self.low_pressure_turbine
#        high_pressure_turbine     = self.high_pressure_turbine
#        combustor_2               = self.combustor_2
#        core_nozzle_2             = self.core_nozzle_2    
#
#        # Dual mode ramjet
#        combustor                 = self.combustor
#        core_nozzle               = self.core_nozzle
#
#
#        thrust                    = self.thrust
#        thrust                    = self.thrust_2
#        
#        # Build turbojet 
#                
#        turbojet =  SUAVE.Components.Energy.Networks.Turbojet_Super()
#        turbojet.number_of_engines = number_of_engines
#        turbojet.working_fluid = self.working_fluid
#        turbojet.ram = ram
#        turbojet.inlet_nozzle = inlet_nozzle
#        turbojet.inlet_nozzle.compressibility_effects = False
#        turbojet.low_pressure_compressor = low_pressure_compressor
#        turbojet.high_pressure_compressor = high_pressure_compressor
#        turbojet.combustor = combustor_2
#        turbojet.high_pressure_turbine = high_pressure_turbine
#        turbojet.low_pressure_turbine = low_pressure_turbine
#        turbojet.core_nozzle = core_nozzle_2
#        turbojet.thrust = thrust
#        
#        turbojet.thrust.size(conditions)
#        
#    
#        # Build ramjet 
#                
#        ramjet =  SUAVE.Components.Energy.Networks.Ramjet()
#        ramjet.number_of_engines = number_of_engines
#        ramjet.working_fluid = self.working_fluid
#        ramjet.ram = ram
#        ramjet.inlet_nozzle = inlet_nozzle
#        ramjet.inlet_nozzle.compressibility_effects = True
#        ramjet.combustor = combustor
#        ramjet.core_nozzle = core_nozzle
#        ramjet.thrust = thrust
#
#        # compute the thrust
#        thrust.size(conditions)
        
    __call__ = evaluate_thrust
