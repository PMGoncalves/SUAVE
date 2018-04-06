## @ingroup Methods-Aerodynamics-Common-Fidelity_Zero-Lift
# fuselage_correction.py
# 
# Created:  Dec 2013, A. Variyar 
# Modified: Feb 2014, A. Variyar, T. Lukaczyk, T. Orra 
#           Apr 2014, A. Variyar
#           Jan 2015, E. Botero

# ----------------------------------------------------------------------
#  Fuselage Correction
# ----------------------------------------------------------------------
import numpy as np
## @ingroup Methods-Aerodynamics-Common-Fidelity_Zero-Lift
def fuselage_correction(state,settings,geometry):  
    """Corrects aircraft lift based on fuselage effects

    Assumptions:
    None

    Source:
    adg.stanford.edu (Stanford AA241 A/B Course Notes)

    Inputs:
    settings.fuselage_lift_correction  [Unitless]
    state.conditions.
      freestream.mach_number           [Unitless]
      aerodynamics.angle_of_attack     [radians]
      aerodynamics.lift_coefficient    [Unitless]

    Outputs:
    aircraft_lift_total                [Unitless]

    Properties Used:
    N/A
    """         
   
    # unpack
    fus_correction  = settings.fuselage_lift_correction
    Mc              = state.conditions.freestream.mach_number
    AoA             = state.conditions.aerodynamics.angle_of_attack
    wings_lift_comp = state.conditions.aerodynamics.lift_coefficient
    
    # total lift, accounting one fuselage
    aircraft_lift_total = wings_lift_comp * fus_correction 
    
    aircraft_lift_total[Mc >= 5.698] = 0.032
    aircraft_lift_total[Mc > 9.712] = 0.01
    aircraft_lift_total[Mc > 13.712] = 0.0045
    aircraft_lift_total[Mc > 17.712] = 0.0035
    
    state.conditions.aerodynamics.lift_coefficient= aircraft_lift_total

    return aircraft_lift_total