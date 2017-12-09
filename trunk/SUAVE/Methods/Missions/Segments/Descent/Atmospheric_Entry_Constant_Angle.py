## @ingroup Methods-Missions-Segments-Descent
# Constant_Speed_Constant_Angle.py
# 
# Created:  Jul 2014, SUAVE Team
# Modified: Jan 2016, E. Botero

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------
#  Initialize Conditions
# ----------------------------------------------------------------------

## @ingroup Methods-Missions-Segments-Descent
def initialize_conditions(segment,state):
    """Sets the specified conditions which are given for the segment type.

    Assumptions:
    Constant speed and constant descent angle

    Source:
    N/A

    Inputs:
    segment.descent_angle                       [radians]
    segment.altitude_start                      [meters]
    segment.altitude_end                        [meters]
    segment.air_speed                           [meters/second]
    state.numerics.dimensionless.control_points [array]

    Outputs:
    conditions.frames.inertial.velocity_vector  [meters/second]
    conditions.frames.inertial.position_vector  [meters]
    conditions.freestream.altitude              [meters]
    conditions.frames.inertial.time             [seconds]

    Properties Used:
    N/A
    """        
    
    # unpack
    descent_angle= segment.descent_angle
    entry_speed  = segment.air_speed   
    alt0         = segment.altitude_start 
    altf         = segment.altitude_end
    t_nondim     = state.numerics.dimensionless.control_points
    conditions   = state.conditions  

    ballistic_nd = 0.3
    
    # check for initial altitude
    if alt0 is None:
        if not state.initials: raise AttributeError('initial altitude not set')
        alt0 = -1.0 * state.initials.conditions.frames.inertial.position_vector[-1,2]

    # discretize on altitude
    alt = t_nondim * (altf-alt0) + alt0
                     
    print alt
    
    # process velocity vector
    #v_mag = state.ones_row(1)
    v_mag = entry_speed * alt
    
    print v_mag
    v_mag = entry_speed*np.exp((1/(2*ballistic_nd*np.sin(-descent_angle)))*np.exp(-(alt[:,0]/alt0)*(alt0/7640)))

    print v_mag

    v_x   = v_mag * np.cos(-descent_angle)
    v_z   = -v_mag * np.sin(-descent_angle)
    
    # pack conditions    
    conditions.frames.inertial.velocity_vector[:,0] = v_x#[:,0]
    conditions.frames.inertial.velocity_vector[:,2] = v_z#[:,0]
    conditions.frames.inertial.position_vector[:,2] = -alt[:,0] # z points down
    conditions.freestream.altitude[:,0]             =  alt[:,0] # positive altitude in this context