## @ingroup Attributes-Propellants
#Liquid H2
#
# Created:  Unk 2013, SUAVE TEAM
# Modified: Apr 2015, SUAVE TEAM
#           Feb 2016, M.Vegh
# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from Propellant import Propellant

# ----------------------------------------------------------------------
#  Jet_A1 Propellant Class
# ----------------------------------------------------------------------
## @ingroup Attributes-Propellants
class Liquid_H2(Propellant):
    """Holds values for this propellant
    
    Assumptions:
    None
    
    Source:
    None
    """

    def __defaults__(self):
        """This sets the default values.

        Assumptions:
        None

        Source:
        Values commonly available

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        None
        """    
        self.tag                       = 'Liquid H2'
        self.reactant                  = 'O2'
        self.density                   = 59.9                             # kg/m^3 (15 C, 1 atm)
        self.specific_energy           = 141.86e6                          # J/kg
        self.energy_density            = 8491.0e6                        # J/m^3
        self.stoichiometric_air        = 0.0291        