## @ingroup Analyses-Atmospheric
# US_Standard_1976.py
#
# Created: 
# Modified: Feb 2016, Andrew Wendorff

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
from warnings import warn

import SUAVE

from SUAVE.Analyses.Atmospheric import Atmospheric

from SUAVE.Attributes.Gases import Air
from SUAVE.Attributes.Planets import Earth

from SUAVE.Analyses.Mission.Segments.Conditions import Conditions

from SUAVE.Core import Units
from SUAVE.Core.Arrays import atleast_2d_col


# ----------------------------------------------------------------------
#  Classes
# ----------------------------------------------------------------------

## @ingroup Analyses-Atmospheric
class US_Standard_1976(Atmospheric):

    """ Implements the U.S. Standard Atmosphere (1976 version)
        
    Assumptions:
    None
    
    Source:
    U.S. Standard Atmosphere, 1976, U.S. Government Printing Office, Washington, D.C., 1976
    """
    
    def __defaults__(self):
        """This sets the default values for the analysis to function.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Output:
        None

        Properties Used:
        None
        """     
        
        atmo_data = SUAVE.Attributes.Atmospheres.Earth.US_Standard_1976()
        self.update(atmo_data)        
    
    def compute_values(self,altitude,temperature_deviation=0.0):

        """Computes atmospheric values.

        Assumptions:
        US 1976 Standard Atmosphere

        Source:
        U.S. Standard Atmosphere, 1976, U.S. Government Printing Office, Washington, D.C., 1976

        Inputs:
        altitude                                 [m]
        temperature_deviation                    [K]

        Output:
        atmo_data.
          pressure                               [Pa]
          temperature                            [K]
          speed_of_sound                         [m/s]
          dynamic_viscosity                      [kg/(m*s)]

        Properties Used:
        self.
          fluid_properties.gas_specific_constant [J/(kg*K)]
          planet.sea_level_gravity               [m/s^2]
          planet.mean_radius                     [m]
          breaks.
            altitude                             [m]
            temperature                          [K]
            pressure                             [Pa]
        """

        # unpack
        zs        = altitude
        gas       = self.fluid_properties
        planet    = self.planet
        grav      = self.planet.sea_level_gravity        
        Rad       = self.planet.mean_radius
        gamma     = gas.gas_specific_constant
        delta_isa = temperature_deviation
        
        # check properties
        if not gas == Air():
            warn('US Standard Atmosphere not using Air fluid properties')
        if not planet == Earth():
            warn('US Standard Atmosphere not using Earth planet properties')          
        
        # convert input if necessary
        zs = atleast_2d_col(zs)

        # get model altitude bounds
        zmin = self.breaks.altitude[0]
        zmax = self.breaks.altitude[-1]   
        
        # convert geometric to geopotential altitude
        zs = zs/(1 + zs/Rad)
        
        # check ranges
        if np.amin(zs) < zmin:
            print "Warning: altitude requested below minimum for this atmospheric model; returning values for h = -2.0 km"
            zs[zs < zmin] = zmin
        if np.amax(zs) > zmax:
            print "Warning: altitude requested above maximum for this atmospheric model; returning values for h = 86.0 km"   
            zs[zs > zmax] = zmax        

        # initialize return data
        zeros = np.zeros_like(zs)
        p     = zeros * 0.0
        T     = zeros * 0.0
        rho   = zeros * 0.0
        a     = zeros * 0.0
        mew   = zeros * 0.0
        z0    = zeros * 0.0
        T0    = zeros * 0.0
        p0    = zeros * 0.0
        alpha = zeros * 0.0
        
        # populate the altitude breaks
        # this uses >= and <= to capture both edges and because values should be the same at the edges
        for i in xrange( len(self.breaks.altitude)-1 ): 
            i_inside = (zs >= self.breaks.altitude[i]) & (zs <= self.breaks.altitude[i+1])
            z0[ i_inside ]    = self.breaks.altitude[i]
            T0[ i_inside ]    = self.breaks.temperature[i]
            p0[ i_inside ]    = self.breaks.pressure[i]
            alpha[ i_inside ] = -(self.breaks.temperature[i+1] - self.breaks.temperature[i])/ \
                                 (self.breaks.altitude[i+1]    - self.breaks.altitude[i])
        
        # interpolate the breaks
        dz = zs-z0
        i_isoth = (alpha == 0.)
        i_adiab = (alpha != 0.)
        p[i_isoth] = p0[i_isoth] * np.exp(-1.*dz[i_isoth]*grav/(gamma*T0[i_isoth]))
        p[i_adiab] = p0[i_adiab] * ( (1.-alpha[i_adiab]*dz[i_adiab]/T0[i_adiab]) **(1.*grav/(alpha[i_adiab]*gamma)) )
        
        T   = T0 - dz*alpha + delta_isa
        rho = gas.compute_density(T,p)
        a   = gas.compute_speed_of_sound(T)
        mew = gas.compute_absolute_viscosity(T)
                
        atmo_data = Conditions()
        atmo_data.expand_rows(zs.shape[0])
        atmo_data.pressure          = p
        atmo_data.temperature       = T
        atmo_data.density           = rho
        atmo_data.speed_of_sound    = a
        atmo_data.dynamic_viscosity = mew
        
        return atmo_data


# ----------------------------------------------------------------------
#   Module Tests
# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    import pylab as plt
    import scipy as sp
    from scipy.optimize import fsolve

    
    h = np.linspace(-1.,60.,200) * Units.km
    delta_isa = 0.
    atmosphere = US_Standard_1976()
    
    data = atmosphere.compute_values(h,delta_isa)
    p   = data.pressure
    T   = data.temperature
    rho = data.density
    a   = data.speed_of_sound
    mew = data.dynamic_viscosity
    
    # ---------------
    exp_ratio = 69.
    gamma = 1.196
    Pto = 206.4*101325.
    Tto  = 2800.
    minp = 0.2
    Rm = 8134/13.6
    g = 9.8066
    P_out = np.ones((200,1))
    

    # Vandenkerckhove function
    a = gamma
    b = ((1+gamma)/2)**((1+gamma)/(1-gamma))
        
    Gf = np.sqrt(a*b)

    a = Gf
    b = 2*gamma/(gamma-1)
    print 'p', p
    print np.shape(p)

    func = lambda paux : exp_ratio - (  a / ( np.sqrt(b * (paux/Pto)**(2/gamma) * (1-(paux/Pto)**((gamma-1)/gamma)))))
    P_aux = fsolve(func,p,factor = 0.01)
    
    P_out[:,0] = P_aux
    
    print 'P_out', P_out
    print np.shape(P_out)
        
    i = P_out >= minp*p
    # CF
    a = (2*gamma**2)/(gamma-1)
    b = (2/(gamma+1))**((gamma+1)/(gamma-1))
    c = (1 - (P_out/Pto)**((gamma-1)/gamma))
    d = ((P_out - p)/(Pto))*exp_ratio
            
    CF2 = np.sqrt(a*b*c)+d
    print 'CF2', CF2
    print np.shape(CF2)  
    # CD
    a = 2/(gamma + 1)
    b = (gamma+1)/(2*(gamma-1))
    c = np.sqrt((gamma)/(Rm*Tto))
    CD = (a**b)*c
    print np.shape(CD)  
    u_out  = np.sqrt((2/(gamma-1))*(Rm)*Tto * (1 - (P_out/Pto)**((gamma-1)/gamma)))
        
    Isp = CF2/(CD*g)
    Fsp = Isp*g
    print np.shape(Isp)  
              
    plt.figure(1)
    plt.plot(h/Units.km,Isp, 'k')
    plt.xlabel('Altitude (km)')
    plt.ylabel('Isp (s)')
    plt.xlim([1,60])
    plt.ylim([330,450])
    
    plt.figure(2)
    plt.plot(h/Units.km, p/101325.)
    plt.xlabel('Altitude (km)')
    plt.ylabel('Pressure (atm)')
    
    plt.figure(3)
    plt.plot(h/Units.km, CF2)
    plt.xlabel('Altitude (km)')
    plt.ylabel('CF')
    print data
    