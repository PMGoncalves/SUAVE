## @ingroup Attributes-Atmospheres-Earth
#US_Standard_1976.py

# Created:  Mar, 2014, SUAVE Team
# Modified: Feb, 2015, T. MacDonald

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
from SUAVE.Attributes.Gases import Air
from SUAVE.Attributes.Atmospheres import Atmosphere
from SUAVE.Attributes.Planets import Earth
from SUAVE.Core import Data
from SUAVE.Core import Units

# ----------------------------------------------------------------------
#  US_Standard_1976 Atmosphere Class
# ----------------------------------------------------------------------
## @ingroup Attributes-Atmospheres-Earth
class US_Standard_1976(Atmosphere):
    """Contains US Standard 1976 values.
    
    Assumptions:
    None
    
    Source:
    U.S. Standard Atmosphere (1976 version)
    """
    
    def __defaults__(self):
        """This sets the default values at breaks in the atmosphere.

        Assumptions:
        None

        Source:
        U.S. Standard Atmosphere (1976 version)

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        None
        """          
        self.tag = ' U.S. Standard Atmosphere (1976)'

        # break point data: 
        self.fluid_properties = Air()
        self.planet = Earth()
        self.breaks = Data()
        self.breaks.altitude    = np.array( [-2.00    , 0.00,     11.00,      20.00,      32.00,      47.00,      51.00,      71.00,      84.852]) * Units.km     # m, geopotential altitude
        self.breaks.temperature = np.array( [301.15   , 288.15,   216.65,     216.65,     228.65,     270.65,     270.65,     214.65,     186.95])      # K
        self.breaks.pressure    = np.array( [127774.0 , 101325.0, 22632.1,    5474.89,    868.019,    110.906,    66.9389,    3.95642,    0.3734])      # Pa
        self.breaks.density     = np.array( [1.47808e0, 1.2250e0, 3.63918e-1, 8.80349e-2, 1.32250e-2, 1.42753e-3, 8.61606e-4, 6.42099e-5, 6.95792e-6])  # kg/m^3
    
        #86- 120 km
        altitude_append     = np.array([86.,    89.,    91.,  93.,    95.,    97., 99., 101., 103., 105., 107., 110., 112., 114., 116., 118., 120.]) *Units.km
        temperature_append   = np.array([186.87, 186.87, 186.87, 187.25, 188.42, 190.40, 193.27, 197.16, 202.23, 208.84, 217.63, 240., 264., 288., 312., 336., 360.])
        density_append      = np.array([6.92789, 4.09588, 2.86182, 1.99463, 1.38963, 9.68809e-1, 6.76231e-1, 4.7264e-1, 3.30816e-1, 2.31952e-1, 1.63057e-1, 9.69560e-2, 6.9276e-2, 5.01920e-2, 3.71263e-2, 2.82884e-2, 2.2456e-2])*10e-6
        pressure_append     = np.array([0.3741, 0.2195, 0.1540, 0.1081, 0.0760, 0.0536, 0.038, 0.027, 0.019, 0.014, 0.0107, 7.133e-3, 5.5708e-3, 4.4433e-3, 3.619e-3, 3.0065e-3, 2.541e-3])

        self.breaks.altitude        = np.append(self.breaks.altitude, altitude_append)
        self.breaks.temperature     = np.append(self.breaks.temperature,temperature_append)
        self.breaks.pressure        = np.append(self.breaks.pressure,pressure_append)
        self.breaks.density         = np.append(self.breaks.temperature,density_append)
        
        np.set_printoptions(threshold=np.nan)
      
# ----------------------------------------------------------------------
#   Module Tests
# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    import pylab as plt
    
    h = np.linspace(-1.,120,200) * Units.km
    
    atmosphere = US_Standard_1976()
    
    atmo_data = atmosphere.compute_values(h)

    p   = atmo_data.pressure          
    T   = atmo_data.temperature       
    rho = atmo_data.density          
    a   = atmo_data.speed_of_sound    
    mu  = atmo_data.dynamic_viscosity   
    
    plt.figure(1)
    plt.plot(p,h)
    plt.xlabel('Pressure (Pa)')
    plt.ylabel('Altitude (km)')
    
    plt.figure(2)
    plt.plot(T,h)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Altitude (km)')    
    
    plt.figure(3)
    plt.plot(rho,h)
    plt.xlabel('Density (kg/m^3)')
    plt.ylabel('Altitude (km)')       
    
    plt.figure(4)
    plt.plot(a,h)
    plt.xlabel('Speed of Sound (m/s)')
    plt.ylabel('Altitude (km)') 
    
    plt.figure(6)
    plt.plot(mu)
    plt.xlabel('Viscosity (kg/m-s)')
    plt.ylabel('Altitude (km)')   

    plt.show()