## @ingroup Methods-Missions-Segments-Common
# Energy.py
# 
# Created:  Jul 2014, SUAVE Team
# Modified: Jan 2016, E. Botero
#           Jul 2017, E. Botero

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------
#  Initialize Battery
# ----------------------------------------------------------------------

## @ingroup Methods-Missions-Segments-Common
def initialize_battery(segment,state):
    """ Sets the initial battery energy at the start of the mission
    
        Assumptions:
        N/A
        
        Inputs:
            state.initials.conditions:
                propulsion.battery_energy    [Joules]
            segment.battery_energy           [Joules]
            
        Outputs:
            state.conditions:
                propulsion.battery_energy    [Joules]

        Properties Used:
        N/A
                                
    """
    
    
    if state.initials:
        energy_initial  = state.initials.conditions.propulsion.battery_energy[-1,0]
    elif segment.has_key('battery_energy'):
        energy_initial  = segment.battery_energy
    else:
        energy_initial = 0.0
    
    state.conditions.propulsion.battery_energy[:,0] = energy_initial

    return

# ----------------------------------------------------------------------
#  Update Thrust
# ----------------------------------------------------------------------

## @ingroup Methods-Missions-Segments-Common
def update_thrust(segment,state):
    """ Evaluates the energy network to find the thrust force and mass rate

        Inputs -
            segment.analyses.energy_network    [Function]
            state                              [Data]

        Outputs -
            state.conditions:
               frames.body.thrust_force_vector [Newtons]
               weights.vehicle_mass_rate       [kg/s]


        Assumptions -


    """    
    
    # unpack
    energy_model = segment.analyses.energy

    # evaluate
    results   = energy_model.evaluate_thrust(state)

    # pack conditions
    conditions = state.conditions
    conditions.frames.body.thrust_force_vector = results.thrust_force_vector
    conditions.weights.vehicle_mass_rate       = results.vehicle_mass_rate
    
    conditions.propulsion.pto = results.pto
    conditions.propulsion.pt1 = results.pt1
    conditions.propulsion.pt2 = results.pt2
    conditions.propulsion.pt3 = results.pt3
    conditions.propulsion.pt4 = results.pt4
    conditions.propulsion.pt5 = results.pt5
    conditions.propulsion.pt6 = results.pt6
    conditions.propulsion.pt7 = results.pt7
    conditions.propulsion.pt8 = results.pt8
    
#    conditions.propulsion.po = results.po
#    conditions.propulsion.p1 = results.p1
#    conditions.propulsion.p2 = results.p2
#    conditions.propulsion.p3 = results.p3
    
    conditions.propulsion.tto = results.tto
    conditions.propulsion.tt1 = results.tt1
    conditions.propulsion.tt2 = results.tt2
    conditions.propulsion.tt3 = results.tt3
    conditions.propulsion.tt4 = results.tt4
    conditions.propulsion.tt5 = results.tt5
    conditions.propulsion.tt6 = results.tt6
    conditions.propulsion.tt7 = results.tt7
    conditions.propulsion.tt8 = results.tt8

    
#    conditions.propulsion.mo = results.mo
#    conditions.propulsion.m1 = results.m1
#    conditions.propulsion.m11 = results.m11
#    conditions.propulsion.m2 = results.m2
#    conditions.propulsion.m3 = results.m3
    
    conditions.propulsion.fsp = results.fsp
    conditions.propulsion.isp = results.isp
    conditions.propulsion.f   = results.f
    
    conditions.propulsion.no  = results.no

    