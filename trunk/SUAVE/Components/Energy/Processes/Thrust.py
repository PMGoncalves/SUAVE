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
        self.number_of_engines                        = 0.0
        self.inputs.fuel_to_air_ratio                 = 0.0
        self.outputs.thrust                           = 0.0 
        self.outputs.thrust_specific_fuel_consumption = 0.0
        self.outputs.specific_impulse                 = 0.0
        self.outputs.non_dimensional_thrust           = 0.0
        self.outputs.core_mass_flow_rate              = 0.0
        self.outputs.fuel_flow_rate                   = 0.0
        self.outputs.fuel_mass                        = 0.0
        self.outputs.power                            = 0.0
        self.design_thrust                            = 0.0
        self.mass_flow_rate_design                    = 0.0
       
        
        #Silly stuff
        self.outputs.heat_flux_and                    = 0.0

	
    def Cpw(T):
    
        T = T - 273.15
        
        if T < 507 :
            Cpw     = 526.727 - 0.0476148 * T - 8.92138 * 10e6 / (273.15 + T)**2.
        elif T < 787 :
            Cpw     = 2.71*(-133.249 + T)*(399965. - 711.241*T + T**2) / (273.15 + T)**2.
        else:
            Cpw     = 725.
            
        return Cpw


    
    def thermo_anderson(self,conditions):
        
        Vo = conditions.freestream.velocity
        rho = conditions.freestream.density
        Rle = 1
        
        q_and  = (1.83e-8*(rho/Rle)**0.5*Vo**3.) #W/m2
        
        self.outputs.heat_flux_and = q_and
                
        return
#        from SUAVE.Methods.Propulsion.oblique_shock import theta_beta_mach
#        from SUAVE.Attributes.Gases import Air
#        # unpacking from conditions
#        gamma                = conditions.freestream.isentropic_expansion_factor
#        Cp                   = conditions.freestream.specific_heat_at_constant_pressure
#        R                    = conditions.freestream.specific_gas_constant
#        Vo                   = conditions.freestream.velocity
#        To                   = conditions.freestream.temperature
#        a0                   = conditions.freestream.speed_of_sound
#        M                 = conditions.freestream.mach_number
#        Po                   = conditions.freestream.pressure  
#        g                    = conditions.freestream.gravity
#        a                    = conditions.freestream.speed_of_sound
#        rho                  = conditions.freestream.density
#        mew                  = conditions.freestream.dynamic_viscosity
#        Tto                   = conditions.freestream.stagnation_temperature
#        Pto                  = conditions.freestream.stagnation_pressure
#        gas = Air()
#
#        Pr = 0.715
#        Tto = conditions.freestream.stagnation_temperature   
#        Pto = conditions.freestream.stagnation_pressure
#        
#        #---- Geometry
#        Rle     = 3.0e-3                  # Leading edge radius
#        L_flat  = 82.3e-3                # Lenght flat plate   
#        theta_f = 84.*np.pi/180.         # radians
#        phi     = np.pi/2. - theta_f    # radians
#        
#        #---- Material
#        eps     = 0.95      # emissivity
#        
#        #---- Constants
#        sig     = 5.67e-8   # Steffan-Boltzman constant
#    
#    
#        #-- Stagnation
#        rhoo    = Pto/(R*Tto)
#        
#        #-- Shock
#        Tt1     = Tto
#        Pt1     = Pto * ((((gamma+1.)*M*M)/(2.+(gamma-1.)*M*M))**(gamma/(gamma-1.)))*((gamma+1.)/(2.*gamma*M*M-(gamma-1.)))**(1./(gamma-1.))
#        rho1    = Pt1/(R*Tt1)
#        mew1    = gas.compute_absolute_viscosity(Tt1)
#        
#        #-- Walls
#        haw     = Cp*To + Vo*Vo/2.       #adiabatic wall enthalpy
#        duds    = (1/Rle)*np.sqrt(2*(Pt1-Po)/rho1)
#    
#    
#    
#        #------------------------
#    
#        def integrand(s, Tiso):
#            
#            #-- Wall Cp    
#            Cpw1 = Cpw(Tiso)
#            
#                
#            #-- Position
#            theta = s/Rle
#            #-- Stagnation point heat flux
#            q_st_gw = 0.57 * Pr**(-0.6) * (rho1*mew1)**0.5 * (haw - Cpw1*Tiso) * duds**0.5
#            
#            #-- Around the tip curvature heat flux
#            G           = (1. - 1./(gamma*M*M))*(theta*theta - theta/2.*np.sin(4.*theta) + (1.-np.cos(4.*theta))/8.) + 4./(gamma*M*M)*((theta*theta - theta*np.sin(2.*theta) + (1.-np.cos(2.*theta))/2.))
#            
#            qt_qst  = 2.*theta*np.sin(theta)*((1.-(1./(gamma*M*M)))*(np.cos(theta))**2. + 1./(gamma*M*M))*G**-0.5
#    
#            qr      = sig*eps*Tiso**4.
#            
#            
#            return q_st_gw*qt_qst - qr
#    
#        
#            
#        def integrand2(y,Tiso) :
#            #-- Flat plate
#            
#            #-- Wall Cp    
#            Cpw1 = Thrust.Cpw(Tiso)
#            
#            #---- After shock conditions
#            beta    = theta_beta_mach(M, phi, gamma,1 )
#            Msin    = M*np.sin(beta)
#            M2      = ((1/(np.sin(beta-phi))**2.)*((Msin)**2.+(2./(gamma-1.)))/((2.*gamma*Msin**2.)/(gamma-1.)-1.))**0.5
#            P2      = p*(1.+((2.*gamma)/(1.+gamma))*(Msin**2.-1.))
#            rho2    = rhoo*(gamma+1.)*Msin*Msin/(2.+(gamma-1.)*Msin*Msin)
#            V2      = Vo*(M2/M)*((P2/p)*(rhoo/rho2))**(0.5)
#            T2      = ((P2*rhoo)/(p*rho2))*T          
#            mew2    = gas.compute_absolute_viscosity(T2)
#                
#            #---- Wall enthalpy
#            h2aw    = Cp*T2 + r*V2*V2/2.
#            
#            #---- Find wall temperature
#            
#            Cpw2      = Cpw(3000)
#            T2aw    = h2aw / Cpw2
#                 
#            Pr_star = (((T2aw/T2)-1.)*(2./((gamma-1.)*M2*M2)))**2.
#            Re_y    = rho2*V2*y/mew2
#            
#            Tstar   = T2*(0.5 + 0.039*M2*M2 + 0.5*Tiso/T2)
#            Cfstar  = 0.664*(Tstar/T2)**(-1./6.)*Re_y**-0.5
#            CH      = Cfstar/(2.*Pr_star**(2./3.))
#            q_flat  = CH*rho2*V2*(h2aw-Cpw1*Tiso)
#            
#            qr      = sig*eps*Tiso**4.
#            
#            return q_flat - qr
#
#
        
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
          Specific Impulse                   [s]

        Properties Used:
        self.
          reference_temperature              [K]
          reference_pressure                 [Pa]
          compressor_nondimensional_massflow [-]
        """           
        # unpack the values
        
        # unpacking from conditions
        gamma                = conditions.freestream.isentropic_expansion_factor
        Cp                   = conditions.freestream.specific_heat_at_constant_pressure
        u0                   = conditions.freestream.velocity
        a0                   = conditions.freestream.speed_of_sound
        M0                   = conditions.freestream.mach_number
        p0                   = conditions.freestream.pressure  
        g                    = conditions.freestream.gravity
        throttle             = conditions.propulsion.throttle        
        
        # unpacking from inputs
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
        
        # unpacking from self
        Tref                 = self.reference_temperature
        Pref                 = self.reference_pressure
        mdhc                 = self.compressor_nondimensional_massflow

        
    
        ##--------Cantwell method---------------------------------
        

        # computing the non dimensional thrust
        core_thrust_nondimensional  = flow_through_core*(gamma*M0*M0*(core_nozzle.velocity/u0-1.) + core_area_ratio*(core_nozzle.static_pressure/p0-1.))
        fan_thrust_nondimensional   = flow_through_fan*(gamma*M0*M0*(fan_nozzle.velocity/u0-1.) + fan_area_ratio*(fan_nozzle.static_pressure/p0-1.))
        
        Thrust_nd                   = core_thrust_nondimensional + fan_thrust_nondimensional
      
     
        Fsp              = 1./(gamma*M0)*Thrust_nd

        # computing the specific impulse
        Isp              = Fsp*a0*(1.+bypass_ratio)/(f*g)
        
        # computing the TSFC
        TSFC             = 3600.*f*g/(Fsp*a0*(1.+bypass_ratio))  
        
        # computing the core mass flow
        mdot_core        = mdhc*np.sqrt(Tref/total_temperature_reference)*(total_pressure_reference/Pref)

        # computing the dimensional thrust
        FD2              = Fsp*a0*(1+bypass_ratio)*mdot_core*no_eng*throttle

     
        
        # fuel flow rate
        a = np.array([0.])        
        fuel_flow_rate   = np.fmax(0.1019715*FD2*TSFC/3600.,a) #use units package for the constants
        
        # computing the power 
        power            = FD2*u0
        
        # pack outputs
        self.outputs.thrust                            = FD2 
        self.outputs.thrust_specific_fuel_consumption  = TSFC
        self.outputs.non_dimensional_thrust            = Fsp 
        self.outputs.core_mass_flow_rate               = mdot_core
        self.outputs.fuel_flow_rate                    = fuel_flow_rate    
        self.outputs.power                             = power  
        self.outputs.specific_impulse                  = Isp
        
                     
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
        # unpack inputs
        a0                   = conditions.freestream.speed_of_sound
        throttle             = 1.0
        
        # unpack from self
        bypass_ratio                = self.inputs.bypass_ratio
        Tref                        = self.reference_temperature
        Pref                        = self.reference_pressure
        design_thrust               = self.total_design
        
        total_temperature_reference = self.inputs.total_temperature_reference  # low pressure turbine output for turbofan
        total_pressure_reference    = self.inputs.total_pressure_reference
        no_eng                      = self.inputs.number_of_engines
        
        # compute nondimensional thrust
        self.compute(conditions)
        
        # unpack results 
        Fsp                         = self.outputs.non_dimensional_thrust
     
        # compute dimensional mass flow rates
        mdot_core                   = design_thrust/(Fsp*a0*(1.+bypass_ratio)*no_eng*throttle)  
        mdhc                        = mdot_core/ (np.sqrt(Tref/total_temperature_reference)*(total_pressure_reference/Pref))
    
        # pack outputs
        self.mass_flow_rate_design               = mdot_core
        self.compressor_nondimensional_massflow  = mdhc   
        
        return
    
    __call__ = compute         

