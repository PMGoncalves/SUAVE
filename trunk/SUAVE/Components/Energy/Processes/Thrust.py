## @ingroup Components-Energy-Processes
# Thrust.py
#
# Created:  Jul 2014, A. Variyar
# Modified: Feb 2016, T. MacDonald, A. Variyar, M. Vegh


# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# SUAVE imports

import SUAVE

from SUAVE.Core import Units

# package imports
import numpy as np
import scipy as sp
from scipy.optimize import fsolve

from SUAVE.Core import Data
from SUAVE.Components import Component, Physical_Component, Lofted_Body
from SUAVE.Components.Energy.Energy_Component import Energy_Component
from SUAVE.Components.Propulsors.Propulsor import Propulsor


# ----------------------------------------------------------------------
#  Thrust Process
# ----------------------------------------------------------------------
## @ingroup Components-Energy-Processes
class Thrust(Energy_Component):
    """A class that handles computation of thrust and other outputs for a gas turbine engine.
    
    Assumptions:
    Perfect gas
    
    Source:
    https://web.stanford.edu/~cantwell/AA283_Course_Material/AA283_Course_Notes/
    """         
    
    def __defaults__(self):
        """This sets the default value.

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
        self.tag ='Thrust'
        self.bypass_ratio                             = 0.0
        self.compressor_nondimensional_massflow       = 0.0
        self.reference_temperature                    = 288.15
        self.reference_pressure                       = 1.01325*10**5
        self.inputs.number_of_engines                        = 0.0
        self.inputs.fuel_to_air_ratio                 = 0.0  #changed
        self.outputs.thrust                           = 0.0 
        self.outputs.thrust_specific_fuel_consumption = 0.0
        self.outputs.specific_impulse                 = 0.0
        self.outputs.non_dimensional_thrust           = 0.0
        self.outputs.core_mass_flow_rate              = 0.0
        self.outputs.fuel_flow_rate                   = 0.0
        self.outputs.fuel_mass                        = 0.0 #changed
        self.outputs.power                            = 0.0
        self.design_thrust                            = 0.0
        self.mass_flow_rate_design                    = 0.0
        self.isp_design                               = 0.0
        self.outputs.specific_thrust                  = 0.0
        self.inputs.stagnation_pressure               = 70*101325.
        self.inputs.stagnation_temperature            = 2220.
        self.inputs.expansion_ratio                   = 20.
        self.inputs.specific_gas_constant             = 280.


	
 
    def compute(self,conditions):
        """Computes thrust and other properties as below.

        Assumptions:
        Perfect gas

        Source:
        https://web.stanford.edu/~cantwell/AA283_Course_Material/AA283_Course_Notes/

        Inputs:
        conditions.freestream.
          isentropic_expansion_factor        [-] (gamma)
          specific_heat_at_constant_pressure [J/(kg K)]
          velocity                           [m/s]
          speed_of_sound                     [m/s]
          mach_number                        [-]
          pressure                           [Pa]
          gravity                            [m/s^2]
        conditions.throttle                  [-] (.1 is 10%)
        self.inputs.
          fuel_to_air_ratio                  [-]
          total_temperature_reference        [K]
          total_pressure_reference           [Pa]
          core_nozzle.
            velocity                         [m/s]
            static_pressure                  [Pa]
            area_ratio                       [-]
          fan_nozzle.
            velocity                         [m/s]
            static_pressure                  [Pa]
            area_ratio                       [-]
          number_of_engines                  [-]
          bypass_ratio                       [-]
          flow_through_core                  [-] percentage of total flow (.1 is 10%)
          flow_through_fan                   [-] percentage of total flow (.1 is 10%)

        Outputs:
        self.outputs.
          thrust                             [N]
          thrust_specific_fuel_consumption   [N/N-s]
          non_dimensional_thrust             [-]
          core_mass_flow_rate                [kg/s]
          fuel_flow_rate                     [kg/s]
          power                              [W]

        Properties Used:
        self.
          reference_temperature              [K]
          reference_pressure                 [Pa]
          compressor_nondimensional_massflow [-]
        """           
        #unpack the values
        
        #unpacking from conditions
        gamma                = conditions.freestream.isentropic_expansion_factor
        Cp                   = conditions.freestream.specific_heat_at_constant_pressure
        u0                   = conditions.freestream.velocity
        a0                   = conditions.freestream.speed_of_sound
        M0                   = conditions.freestream.mach_number
        p0                   = conditions.freestream.pressure  
        g                    = conditions.freestream.gravity
        throttle             = conditions.propulsion.throttle        
        
        #unpacking from inputs
        f                           = self.inputs.fuel_to_air_ratio
        total_temperature_reference = self.inputs.total_temperature_reference
        total_pressure_reference    = self.inputs.total_pressure_reference
        core_nozzle                 = self.inputs.core_nozzle
        fan_nozzle                  = self.inputs.fan_nozzle
        fan_exit_velocity           = self.inputs.fan_nozzle.velocity
        core_exit_velocity          = self.inputs.core_nozzle.velocity
        fan_area_ratio              = self.inputs.fan_nozzle.area_ratio
        core_area_ratio             = self.inputs.core_nozzle.area_ratio
        no_eng                      = self.inputs.number_of_engines                      
        bypass_ratio                = self.inputs.bypass_ratio  
        flow_through_core           = self.inputs.flow_through_core #scaled constant to turn on core thrust computation
        flow_through_fan            = self.inputs.flow_through_fan #scaled constant to turn on fan thrust computation
        
        #unpacking from self
        Tref                 = self.reference_temperature
        Pref                 = self.reference_pressure
        mdhc                 = self.compressor_nondimensional_massflow

        
    
        ##--------Cantwell method---------------------------------
        

        #computing the non dimensional thrust
        core_thrust_nondimensional  = flow_through_core*(gamma*M0*M0*(core_nozzle.velocity/u0-1) + core_area_ratio*(core_nozzle.static_pressure/p0-1))
        fan_thrust_nondimensional   = flow_through_fan*(gamma*M0*M0*(fan_nozzle.velocity/u0-1) + fan_area_ratio*(fan_nozzle.static_pressure/p0-1))
        
        Thrust_nd                   = core_thrust_nondimensional + fan_thrust_nondimensional
      
     
        Fsp              = 1./(gamma*M0)*Thrust_nd

        #Computing the specific impulse
        #Isp              = Fsp*a0*(1+bypass_ratio)/(f*g)
        
        #Computing the TSFC
        TSFC             = 3600.*f*g/(Fsp*a0*(1+bypass_ratio))  
        #computing the core mass flow
        mdot_core        = mdhc*np.sqrt(Tref/total_temperature_reference)*(total_pressure_reference/Pref)

        #computing the dimensional thrust
        FD2              = Fsp*a0*(1+bypass_ratio)*mdot_core*no_eng*throttle

     
        
        #fuel flow rate
        a = np.array([0.])        
        fuel_flow_rate   = np.fmax(0.1019715*FD2*TSFC/3600,a) #use units package for the constants
        
        #computing the power 
        power            = FD2*u0
        
        #pack outputs
        
        self.outputs.thrust                            = FD2 
        self.outputs.thrust_specific_fuel_consumption  = TSFC
        self.outputs.non_dimensional_thrust            = Fsp 
        self.outputs.core_mass_flow_rate               = mdot_core
        self.outputs.fuel_flow_rate                    = fuel_flow_rate    
        self.outputs.power                             = power  
    
        
    
    def size(self,conditions):
        """Sizes the core flow for the design condition.

        Assumptions:
        Perfect gas

        Source:
        https://web.stanford.edu/~cantwell/AA283_Course_Material/AA283_Course_Notes/

        Inputs:
        conditions.freestream.speed_of_sound [m/s] (conditions is also passed to self.compute(..))
        self.inputs.
          bypass_ratio                       [-]
          total_temperature_reference        [K]
          total_pressure_reference           [Pa]
          number_of_engines                  [-]

        Outputs:
        self.outputs.non_dimensional_thrust  [-]

        Properties Used:
        self.
          reference_temperature              [K]
          reference_pressure                 [Pa]
          total_design                       [N] - Design thrust
        """             
        #unpack inputs
        a0                   = conditions.freestream.speed_of_sound
        throttle             = 1.0
        
        #unpack from self
        bypass_ratio                = self.inputs.bypass_ratio
        Tref                        = self.reference_temperature
        Pref                        = self.reference_pressure
        design_thrust               = self.total_design
        
        total_temperature_reference = self.inputs.total_temperature_reference  # low pressure turbine output for turbofan
        total_pressure_reference    = self.inputs.total_pressure_reference
        no_eng                      = self.inputs.number_of_engines
        
        #compute nondimensional thrust
        self.compute(conditions)
        
        #unpack results 
        Fsp                         = self.outputs.non_dimensional_thrust

                
        #compute dimensional mass flow rates
        mdot_core                   = design_thrust/(Fsp*a0*(1+bypass_ratio)*no_eng*throttle)  
        mdhc                        = mdot_core/ (np.sqrt(Tref/total_temperature_reference)*(total_pressure_reference/Pref))
    
        #pack outputs
        self.mass_flow_rate_design               = mdot_core
        self.compressor_nondimensional_massflow  = mdhc
 
         
        
        return
    
    
    
    def size_rocket(self,conditions):         
        #unpack inputs
        throttle             = 1.0
        
        #unpack from self
        design_thrust               = self.total_design
        design_isp                  = self.isp_design
        no_eng                      = self.number_of_engines
        #compute nondimensional thrust
        
        
        self.compute_rocket(conditions)
        
        #unpack results 
        Fsp                         = self.outputs.specific_thrust
                
        #compute dimensional mass flow rates
        mdot_core                   = design_thrust/(Fsp*no_eng*throttle)  
    
        #pack outputs
        self.mass_flow_rate_design               = mdot_core
        self.isp_design                          = design_isp
 
         
        
        return
    
    
    def compute_rocket(self,conditions):
        #unpack the values
        
        #unpacking from conditions
        g                    = conditions.freestream.gravity
        throttle             = conditions.propulsion.throttle 
        u0                   = conditions.freestream.velocity
        density              = conditions.freestream.density
        Po                   = conditions.freestream.pressure

        no_eng               = self.inputs.number_of_engines                      
        Pto                  = self.inputs.stagnation_pressure
        Tto                  = self.inputs.stagnation_temperature
        exp_ratio            = self.inputs.expansion_ratio    
        Rm                   = self.inputs.specific_gas_constant


        P_out = 1.0*Po/Po
        gamma = 1.2
        #-- Method 1
                                               
        # Vandenkerckhove function
        a = gamma
        b = ((1+gamma)/2)**((1+gamma)/(1-gamma))        
        Gf = np.sqrt(a*b)
        
        # P_out calculation
        a = Gf
        b = 2*gamma/(gamma-1)
        

        func = lambda paux : exp_ratio - (  a / ( np.sqrt(b * (paux/Pto)**(2/gamma) * (1-(paux/Pto)**((gamma-1)/gamma)))))
        i = 0
        for x in Po :    
            P_out[i] = fsolve(func,0.4*Po[i],factor = 0.01) 
            i = i+1
                      
        
#        print P_out
        i = P_out < 0.4*Po
        P_out[i] = 0.4*Po[i]
             
        # CF
        a = (2*gamma**2)/(gamma-1)
        b = (2/(gamma+1))**((gamma+1)/(gamma-1))
        c = (1 - (P_out/Pto)**((gamma-1)/gamma))
        d = ((P_out - Po)/(Pto))*exp_ratio                
        CF = np.sqrt(a*b*c)+d

        # CD
        a = 2/(gamma + 1)
        b = (gamma+1)/(2*(gamma-1))
        c = np.sqrt((gamma)/(Rm*Tto))
        CD = (a**b)*c
             
        # Properties
        Isp = CF/(CD*g)*Po/Po
        
        u_out  = np.sqrt((2/(gamma-1))*(Rm)*Tto * (1 - (P_out/Pto)**((gamma-1)/gamma)))
        
        Fsp = Isp*g
#        Fsp = u_out + (1/CD)*(P_out - Po)*exp_ratio / Po
#        print 'u_out', u_out
#        print 'CD', CD
#        print 'P_out', P_out
#        print 'Po', Po
#        print 'exp_ratio', exp_ratio
#        print '---------------------------------------------'
#        print 'term 1', u_out
#        print 'term 2', (1/CD)*(P_out - Po)*exp_ratio / Po

        
                            
        mdot_core        = self.mass_flow_rate_design
        
        FD2              = Fsp*mdot_core*no_eng*throttle
#        print 'FD2', np.shape(FD2)
        TSFC             = 1/ Fsp
        
        #fuel flow rate
        fuel_flow_rate   = mdot_core*density/density #use units package for the constants
        
        print 'mdot', fuel_flow_rate
        print 'ISP', Isp
        #computing the power 
        power            = FD2*u0
        
        #pack outputs
        self.outputs.thrust                            = FD2 
        self.outputs.specific_thrust                   = Fsp
        self.outputs.thrust_specific_fuel_consumption  = TSFC
        self.outputs.core_mass_flow_rate               = mdot_core
        self.outputs.fuel_flow_rate                    = fuel_flow_rate    
        self.outputs.power                             = power  
    
        
    
    
#    def compute_rocket(self,conditions):
#   
#        #unpack the values
#        
#        #unpacking from conditions
#        g                    = conditions.freestream.gravity
#        throttle             = conditions.propulsion.throttle 
#        u0                   = conditions.freestream.velocity
#        density              = conditions.freestream.density
#        
#        #unpacking from inputs
#        no_eng               = self.inputs.number_of_engines                      
#        
#        #unpacking from self
#        Isp                  = self.isp_design
#      
#        Fsp                 = Isp*g*density/density
#        
#        
#        mdot_core           = self.mass_flow_rate_design
#
#        #computing the dimensional thrust
#        FD2              = Fsp*mdot_core*no_eng*throttle
#
#        TSFC             = 1/ Fsp
#        
#        #fuel flow rate
#        fuel_flow_rate   = mdot_core*density/density #use units package for the constants
#        
#        #computing the power 
#        power            = FD2*u0
#        
#        #pack outputs
#        self.outputs.thrust                            = FD2 
#        self.outputs.specific_thrust                   = Fsp
#        self.outputs.thrust_specific_fuel_consumption  = TSFC
#        self.outputs.core_mass_flow_rate               = mdot_core
#        self.outputs.fuel_flow_rate                    = fuel_flow_rate    
#        self.outputs.power                             = power  
#    
#############################

 
         
        
        return

    
    __call__ = compute         

############################################3




#if __name__ == '__main__':
#    
#    import pylab as plt
#    gamma = 1.2
#    Po = 101325.0
#    Pto = 10*Po
#    max_exp = 100.
#    P_out = np.linspace(0, 100*Po, 200)
#    Rm = 280.2
#    Tto = 2000.
#    g = 9.81
#    
#    # Vandenkerckhove function
#    a = gamma
#    b = ((1+gamma)/2)**((1+gamma)/(1-gamma))
#        
#    Gf = np.sqrt(a*b)
#    exp_ratio2 = np.linspace(1, max_exp, 200)
#    
#    P_out2 = exp_ratio2/exp_ratio2*1
#
#    a = Gf
#    b = 2*gamma/(gamma-1)
#    i=0
#    for x in exp_ratio2 :
#        func = lambda paux : exp_ratio2[i] - (  a / ( np.sqrt(b * (paux/Pto)**(2/gamma) * (1-(paux/Pto)**((gamma-1)/gamma)))))
#        P_out2[i] = fsolve(func,0.4*Po,factor = 0.01)
#        i=i+1
#        
#    i = P_out2 >= 0.4*Po
#
#    # CF
#    a = (2*gamma**2)/(gamma-1)
#    b = (2/(gamma+1))**((gamma+1)/(gamma-1))
#    c = (1 - (P_out2/Pto)**((gamma-1)/gamma))
#    d = ((P_out2 - Po)/(Pto))*exp_ratio2
#            
#    CF = np.sqrt(a*b*c)+d
#    plt.figure(1)
#    plt.plot(exp_ratio2[i],CF[i], 'g*')
#                
#    Pto = 100*Po
#
#    # Vandenkerckhove function
#    a = gamma
#    b = ((1+gamma)/2)**((1+gamma)/(1-gamma))
#        
#    Gf = np.sqrt(a*b)
#    exp_ratio2 = np.linspace(1, max_exp, 200)
#    
#    Pto = 6.89e6
#    P_out2 = exp_ratio2/exp_ratio2*1
#
#    a = Gf
#    b = 2*gamma/(gamma-1)
#    i=0
#    for x in exp_ratio2 :
#        func = lambda paux : exp_ratio2[i] - (  a / ( np.sqrt(b * (paux/Pto)**(2/gamma) * (1-(paux/Pto)**((gamma-1)/gamma)))))
#        P_out2[i] = fsolve(func,0.4*Po,factor = 0.01)
#        i=i+1
#        
#    i = P_out2 >= 0.4*Po
#
#    # CF
#    a = (2*gamma**2)/(gamma-1)
#    b = (2/(gamma+1))**((gamma+1)/(gamma-1))
#    c = (1 - (P_out2/Pto)**((gamma-1)/gamma))
#    d = ((P_out2 - Po)/(Pto))*exp_ratio2
#            
#    CF2 = np.sqrt(a*b*c)+d
#    plt.figure(1)
#    plt.plot(exp_ratio2[i],CF2[i], 'k*')
#    
#    
#    Pto = 1000*Po
#
#    # Vandenkerckhove function
#    a = gamma
#    b = ((1+gamma)/2)**((1+gamma)/(1-gamma))
#        
#    Gf = np.sqrt(a*b)
#    exp_ratio2 = np.linspace(1, max_exp, 200)
#    
#    Pto = 6.89e6
#    P_out2 = exp_ratio2/exp_ratio2*1
#
#    a = Gf
#    b = 2*gamma/(gamma-1)
#    i=0
#    for x in exp_ratio2 :
#        func = lambda paux : exp_ratio2[i] - (  a / ( np.sqrt(b * (paux/Pto)**(2/gamma) * (1-(paux/Pto)**((gamma-1)/gamma)))))
#        P_out2[i] = fsolve(func,0.4*Po,factor = 0.01)
#        i=i+1
#        
#    i = P_out2 >= 0.4*Po
#
#    # CF
#    a = (2*gamma**2)/(gamma-1)
#    b = (2/(gamma+1))**((gamma+1)/(gamma-1))
#    c = (1 - (P_out2/Pto)**((gamma-1)/gamma))
#    d = ((P_out2 - Po)/(Pto))*exp_ratio2
#            
#    CF3 = np.sqrt(a*b*c)+d
##    plt.plot(exp_ratio2[i],CF3[i], 'r*')     
#    plt.ylim([0.5,2.3])
#    plt.xlabel('exp')
#    plt.ylabel('CF')
#    plt.show()        

if __name__ == '__main__':

    #        vander_f    = np.sqrt(gamma * ( (1 + gamma) / 2)**((1+gamma)/(1-gamma)))
#        
#        
#        func = lambda pratio : exp_ratio - (vander_f / np.sqrt(((2*gamma/(gamma-1)) * pratio**(2/gamma) * (1 - (pratio)**((gamma-1)/gamma))   )))
#        pratio = fsolve(func,0.01,factor = 0.1)     
#        P_out  = pratio*Pto
    import pylab as plt
    
    
    gamma = 1.2
    Po = 101325.0
    Pto = 100*Po
    max_exp = 100.
    P_out = np.linspace(0, 100*Po, 200)
    Rm = 280.2
    Tto = 2500.
    g = 9.8066
    
    minp = 0.05
    
    # Vandenkerckhove function
    a = gamma
    b = ((1+gamma)/2)**((1+gamma)/(1-gamma))
        
    Gf = np.sqrt(a*b)
        
    # Expansion ratio
    a = Gf
    b = 2*gamma/(gamma-1)
    c = (P_out/Pto)**(2/gamma)
    d = 1 - (P_out/Pto)**((gamma-1)/gamma)
    exp_ratio = a/(np.sqrt(b*c*d))
    
    # CF
    a = (2*gamma**2)/(gamma-1)
    b = (2/(gamma+1))**((gamma+1)/(gamma-1))
    c = (1 - (P_out/Pto)**((gamma-1)/gamma))
    d = ((P_out - Po)/(Pto))*exp_ratio
      
        
    x1 = (P_out >= minp*Po) & (exp_ratio <= max_exp) & (exp_ratio >= 1)
    
    CF = np.sqrt(a*b*c)+d         
    plt.figure(1)
    plt.plot(Pto/P_out[x1],exp_ratio[x1])
    plt.xlabel('PC/POUT')
    plt.ylabel('EXP')

    
    # ---------------
    exp_ratio2 = np.linspace(1, max_exp, 200)
    P_out2 = exp_ratio2/exp_ratio2*1

    a = Gf
    b = 2*gamma/(gamma-1)
    i=0
    for x in exp_ratio2 :
        func = lambda paux : exp_ratio2[i] - (  a / ( np.sqrt(b * (paux/Pto)**(2/gamma) * (1-(paux/Pto)**((gamma-1)/gamma)))))
        P_out2[i] = fsolve(func,Po,factor = 0.01)
        i=i+1
        
    i = P_out2 >= minp*Po

    # CF
    a = (2*gamma**2)/(gamma-1)
    b = (2/(gamma+1))**((gamma+1)/(gamma-1))
    c = (1 - (P_out2/Pto)**((gamma-1)/gamma))
    d = ((P_out2 - Po)/(Pto))*exp_ratio2
            
    CF2 = np.sqrt(a*b*c)+d

    # CD
    a = 2/(gamma + 1)
    b = (gamma+1)/(2*(gamma-1))
    c = np.sqrt((gamma)/(Rm*Tto))
    CD = (a**b)*c
             
    # Properties
    Isp = CF2/(CD*g)
        
    u_out  = np.sqrt((2/(gamma-1))*(Rm)*Tto * (1 - (P_out2/Pto)**((gamma-1)/gamma)))
        
        
    Fsp = Isp*g
    Isp = CF2/(CD*g)*Po/Po
    
    plt.figure(1)
    plt.plot(Pto/P_out2[i],exp_ratio2[i])
    plt.plot(Pto/P_out[x1],exp_ratio[x1])
    plt.xlabel('PC/POUT')
    plt.ylabel('EXP')


    plt.figure(2)
    plt.plot(exp_ratio[x1],CF[x1], 'k^')
    plt.plot(exp_ratio2[i],CF2[i], 'g*')
    plt.ylim([0.5,2.3])
    plt.xlabel('exp')
    plt.ylabel('CF')
    
    plt.figure(3)
    plt.plot(exp_ratio2[i],Fsp[i], 'g*')
    plt.xlabel('exp')
    plt.ylabel('Fsp')

    plt.figure(4)
    plt.plot(exp_ratio2[i],Isp[i], 'g*')
    plt.xlabel('exp')
    plt.ylabel('Isp')
    plt.show()