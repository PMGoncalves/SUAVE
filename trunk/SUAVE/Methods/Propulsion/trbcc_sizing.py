# turbojet_sizing.py
# 
# Created:  May 2015, T. MacDonald 
# Modified: Jan 2016, E. Botero        

# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import SUAVE
import numpy as np
from SUAVE.Core import Data, Units

from SUAVE.Components.Energy.Networks.Turbojet_Super import Turbojet_Super
from SUAVE.Components.Energy.Networks.Ramjet import Ramjet
from SUAVE.Components.Energy.Networks.Scramjet import Scramjet

from SUAVE.Methods.Propulsion.turbojet_sizing import turbojet_sizing
from SUAVE.Methods.Propulsion.ramjet_sizing import ramjet_sizing
from SUAVE.Methods.Propulsion.scramjet_sizing import scramjet_sizing
from SUAVE.Methods.Propulsion.rocket_sizing import rocket_sizing

# ----------------------------------------------------------------------
#   Sizing
# ----------------------------------------------------------------------

def trbcc_sizing(tbcc,mach_number = None, altitude = None, delta_isa = 0, conditions = None):  
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
            conditions.freestream.altitude                    = np.atleast_1d(altitude)
            conditions.freestream.mach_number                 = np.atleast_1d(mach_number)
            conditions.freestream.pressure                    = np.atleast_1d(p)
            conditions.freestream.temperature                 = np.atleast_1d(T)
            conditions.freestream.density                     = np.atleast_1d(rho)
            conditions.freestream.dynamic_viscosity           = np.atleast_1d(mu)
            conditions.freestream.gravity                     = np.atleast_1d(9.81)
            conditions.freestream.isentropic_expansion_factor = np.atleast_1d(1.4)
            conditions.freestream.Cp                          = 1.4*(p/(rho*T))/(1.4-1)
            conditions.freestream.R                           = p/(rho*T)
            conditions.freestream.speed_of_sound              = np.atleast_1d(a)
            conditions.freestream.velocity                    = np.atleast_1d(a*mach_number)
            
            # propulsion conditions
            conditions.propulsion.throttle                    = np.atleast_1d(1.0)
    

    number_of_engines         = tbcc.number_of_engines   

    ram                       = tbcc.ram
    inlet_nozzle              = tbcc.inlet_nozzle

    low_pressure_compressor   = tbcc.low_pressure_compressor
    high_pressure_compressor  = tbcc.high_pressure_compressor
    low_pressure_turbine      = tbcc.low_pressure_turbine
    high_pressure_turbine     = tbcc.high_pressure_turbine
    combustor_2               = tbcc.combustor_2
    core_nozzle_2             = tbcc.core_nozzle_2    

    combustor                 = tbcc.combustor
    core_nozzle               = tbcc.core_nozzle
    thrust                    = tbcc.thrust
    thrust_2                  = tbcc.thrust_2
    
    thrust_3                  = tbcc.thrust_3
    
    # Turbojet
    turbojet = SUAVE.Components.Energy.Networks.Turbojet_Super()
    turbojet.number_of_engines = number_of_engines
    turbojet.working_fluid = tbcc.working_fluid
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
 
    altitude      = 0.0
    mach_number   = 0.01
    
    turbojet_sizing(turbojet,mach_number,altitude)  
    
#    print 'deu bem'
#    print turbojet.thrust.mass_flow_rate_design
    
    ramjet =  SUAVE.Components.Energy.Networks.Ramjet()
    ramjet.number_of_engines = number_of_engines
    ramjet.working_fluid = tbcc.working_fluid
    ramjet.ram = ram
    ramjet.inlet_nozzle = inlet_nozzle
    ramjet.inlet_nozzle.compressibility_effects = True
    ramjet.combustor = combustor
    ramjet.core_nozzle = core_nozzle
    ramjet.thrust = thrust
    
    altitude      = 76060*Units.ft
    mach_number   = 2.0

    ramjet_sizing(ramjet,mach_number,altitude)    

#    print 'deu bem'
#    print turbojet.thrust.mass_flow_rate_design


    rocket =  SUAVE.Components.Energy.Networks.Rocket()
    rocket.number_of_engines = number_of_engines
    rocket.working_fluid = tbcc.working_fluid
    rocket.combustor = combustor
    rocket.nozzle = core_nozzle
    rocket.thrust = thrust_3
    
    altitude      = 0.0*Units.ft
    mach_number   = 0.01
    
    rocket_sizing(rocket,mach_number,altitude)    
#    print 'deu bem 2'
#    print ramjet.thrust.mass_flow_rate_design
    return