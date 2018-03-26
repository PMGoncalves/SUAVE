# compute_isothermal_wall.py
# 
# Created:  May 2015, T. MacDonald 
# Modified: Jan 2016, E. Botero        

# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import SUAVE
import numpy as np
from SUAVE.Core import Data

# ----------------------------------------------------------------------
#   Sizing
# ----------------------------------------------------------------------

def compute_isothermal_wall(vehicle, conditions):  
    """ create and evaluate a gas turbine network
    """    
    
    
    #-- Atmospheric model
    Po = conditions.freestream.pressure
    To = conditions.freestream.temperature
    rho = conditions.freestream.density
    a  = conditions.freestream.speed_of_sound
    mew = conditions.freestream.dynamic_viscosity
    Mach = conditions.freestream.mach_number
    gamma = 1.4
    R  = 287.
    Pr = 0.715
    Cp = 1006.
    Vo = conditions.freestream.velocity
    Tto = conditions.freestream.stagnation_temperature    
    # -- Geometry
    R = vehicle.fuselages['fuselage'].leading_edge_radius
    L_flat = vehicle.fuselages['fuselage'].nose_length
    theta_f = 84.*np.pi/180.         # radians
    phi     = np.pi/2. - theta_f    # radians
    
    # -- Material
    eps     = 0.95      # emissivity
    
    #---- Constants
    sig     = 5.67e-8   # Steffan-Boltzman constant
    
    
    print R
    
    return Mach
    
    # Freestream conditions
    Tto = conditions.freestream.stagnation_temperature    
#    
#    #-- Stagnation
#    Tto     = T*(1.+((gamma-1.)/2. *M*M))
#    Pto     = p* ((1.+(gamma-1.)/2. *M*M )**(gamma/(gamma-1.)))
#    rhoo    = Pto/(R*Tto)
#    
#    #-- Shock
#    Tt1     = Tto
#    Pt1     = Pto * ((((gamma+1.)*M*M)/(2.+(gamma-1.)*M*M))**(gamma/(gamma-1.)))*((gamma+1.)/(2.*gamma*M*M-(gamma-1.)))**(1./(gamma-1.))
#    rho1    = Pt1/(R*Tt1)
#    mew1    = gas.compute_absolute_viscosity(Tt1)
#    
#    #-- Walls
#    haw     = Cp*T + Vo*Vo/2.       #adiabatic wall enthalpy
#    duds    = (1/Rle)*np.sqrt(2*(Pt1-p)/rho1)
#
#
#    #------------------------
#
#    def integrand(s, Tiso):
#        
#        #-- Wall Cp    
#        Cpw1 = Cpw(Tiso)
#        
#            
#        #-- Position
#        theta = s/Rle
#        #-- Stagnation point heat flux
#        q_st_gw = 0.57 * Pr**(-0.6) * (rho1*mew1)**0.5 * (haw - Cpw1*Tiso) * duds**0.5
#        
#        #-- Around the tip curvature heat flux
#        G           = (1. - 1./(gamma*M*M))*(theta*theta - theta/2.*np.sin(4.*theta) + (1.-np.cos(4.*theta))/8.) + 4./(gamma*M*M)*((theta*theta - theta*np.sin(2.*theta) + (1.-np.cos(2.*theta))/2.))
#        
#        qt_qst  = 2.*theta*np.sin(theta)*((1.-(1./(gamma*M*M)))*(np.cos(theta))**2. + 1./(gamma*M*M))*G**-0.5
#
#        qr      = sig*eps*Tiso**4.
#        
#        
#        return q_st_gw*qt_qst - qr
#
#    
#        
#    def integrand2(y,Tiso) :
#        #-- Flat plate
#        
#        #-- Wall Cp    
#        Cpw1 = Cpw(Tiso)
#        
#        #---- After shock conditions
#        beta    = theta_beta_mach(M, phi, gamma,1 )
#        Msin    = M*np.sin(beta)
#        M2      = ((1/(np.sin(beta-phi))**2.)*((Msin)**2.+(2./(gamma-1.)))/((2.*gamma*Msin**2.)/(gamma-1.)-1.))**0.5
#        P2      = p*(1.+((2.*gamma)/(1.+gamma))*(Msin**2.-1.))
#        rho2    = rhoo*(gamma+1.)*Msin*Msin/(2.+(gamma-1.)*Msin*Msin)
#        V2      = Vo*(M2/M)*((P2/p)*(rhoo/rho2))**(0.5)
#        T2      = ((P2*rhoo)/(p*rho2))*T          
#        mew2    = gas.compute_absolute_viscosity(T2)
#            
#        #---- Wall enthalpy
#        h2aw    = Cp*T2 + r*V2*V2/2.
#        
#        #---- Find wall temperature
#        
#        Cpw2      = Cpw(3000)
#        T2aw    = h2aw / Cpw2
#             
#        Pr_star = (((T2aw/T2)-1.)*(2./((gamma-1.)*M2*M2)))**2.
#        Re_y    = rho2*V2*y/mew2
#        
#        Tstar   = T2*(0.5 + 0.039*M2*M2 + 0.5*Tiso/T2)
#        Cfstar  = 0.664*(Tstar/T2)**(-1./6.)*Re_y**-0.5
#        CH      = Cfstar/(2.*Pr_star**(2./3.))
#        q_flat  = CH*rho2*V2*(h2aw-Cpw1*Tiso)
#        
#        qr      = sig*eps*Tiso**4.
#        
#        return q_flat - qr
#
#
#    def func(Tiso):
#        integral, err   = quad(integrand,0,Rle*theta_f, args=(Tiso,))
#        integral2, err2 = quad(integrand2,Rle*np.tan(theta_f),L_flat + Rle*np.tan(theta_f), args=(Tiso,))
#        
#        
#        print '--- ITERATION'
#        print integral
#        print '---'
#        print integral2
#        print '---'
#        print 'Tiso :', Tiso
#        print '------------------'
#        
#        return integral + integral2 
#    
#    T_iso_final = fsolve(func,300, factor=1.)
#    #--------------------    
##        theta_show = np.linspace(0,theta_f,20)
#    Cpw1 = Cpw(T_iso_final)
#    q_st_gw = 0.57 * Pr**(-0.6) * (rho1*mew1)**0.5 * (haw - Cpw1*T_iso_final) * duds**0.5
#    
#    hwc = 1 - Cpw1*T_iso_final/haw
#    q_and  = (1.83e-8*(rho/Rle)**0.5*Vo**3.)*hwc 
#    print '===> Isothermal wall'
#    print q_st_gw *1e-6
#    print q_and * 1e-2
#
#    print '-------------------------'
#    T_iso_final = 0.01
#    Cpw1 = Cpw(T_iso_final)
#    q_st_gw = 0.57 * Pr**(-0.6) * (rho1*mew1)**0.5 * (haw ) * duds**0.5
#    
#    hwc = 1
#    q_and  = (1.83e-8*(rho/Rle)**0.5*Vo**3.)*hwc 
#    print '===> Cold wall'
#    print q_st_gw *1e-6
#    print q_and * 1e-2
#        

