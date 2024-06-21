"""Pressure drop and heat transfer calculation.
"""

from math import pi, sin, log, log10, sqrt, tan, tanh     
from .std_conditions import ureg, Q_, P_NTP
from .functions import Re, Ra, Nu_vcyl, Nu_hcyl     
#from .functions import Material, Property     
#from .functions import nist_property, conduction_cyl     
#from scipy.optimize import root_scalar, minimize  
#from .functions import AIR
from .functions import heat_trans_coef, Ra, Nu_vcyl, Nu_hcyl, Pr
#from .cp_wrapper import ThermState
from .piping import Mach, Mach_total, K_lim, ChokedFlow, HydraulicError, velocity, dP_Darcy, dP_adiab, Pipe, Tube, CopperTube

def laminar_flow(Re_, Pr_, L_ID):
    # Non dimentional calcul of the Nusselt and fiction factor in pipe in laminar flow following Section 5.2.4 of Nellis and Klein (2020)
    # Verify the input conditions     
    if Pr_ < 0.1:
        raise ValueError(f'Prandtl number (Pr) must be > 0.1. The value is {Pr}')
        
    #Inverse Graetz numbers verification
    GZ = L_ID / (Re_ * Pr_)
    SGZ = L_ID / Re_
    if GZ < 1e-6:
        raise ValueError(f'Inverse Graetz number (GZ) must be > 1e-6. The value is {GZ}')
         
    # Calculate friction factor (f)
    f = 4 * (3.44 / sqrt(SGZ) + (1.25 / (4 * SGZ) + 16 - 3.44 / sqrt(SGZ)) / (1 + 0.00021 * SGZ**(-2))) / Re_
    
    # Calculate Nusselt numbers, temperature constant and flux constant
    Nu_T = ((5.001 / GZ**1.119 + 136.0)**0.2978 - 0.6628) / tanh(2.444 * SGZ**(1 / 6) * (1 + 0.565 * SGZ**(1 / 3)))
    Nu_Q = ((6.562 / GZ**1.137 + 220.4)**0.2932 - 0.5003) / tanh(2.530 * SGZ**(1 / 6) * (1 + 0.639 * SGZ**(1 / 3)))
    
    return Nu_T, Nu_Q, f

def turbulent_flow(Re_, Pr_, L_ID, eps):
    # Non dimentional calcul of the Nusselt and fiction factor in pipe in turbulent flow following Section 5.2.3 of Nellis and Klein (2020)
    # Verify the input conditions
    if Pr_ < 0.004 or Pr_ > 2000:
        raise ValueError(f'Prandtl number (Pr) must be between 0.004 and 2000. The value is {Pr}')
    if L_ID <= 1:  
        if L_ID < 0:  ###not inferior to zero - make no sense
            raise ValueError('L/ID ratio < 0. Not possible')
        print('L/ID ratio should be > 1. The value is {L_ID}')
        L_ID = 1
        
    # Friction Factor 
    if eps > 1e-5:
        #Offor & Alabi, Advances in Chemical Engineering and Science, 2016
        friction = (-2 * log10(eps / 3.71 - 1.975 / Re_ * log((eps / 3.93)**1.092 + 7.627 / (Re_ + 395.9))))**(-2)
    else:
        #Li & Seem correlation, A New Explicity Equation for Accurate Friction Factor Calculation for Smooth Tubes, 2011
        friction = (-0.001570232 / log(Re_) + 0.394203137 / log(Re_)**2 + 2.534153311 / log(Re_)**3) * 4
        
    # Nusselt, Gnielinski, Int. Chem. Eng., 1976
    Nusselt = ((friction / 8) * (Re_ - 1000) * Pr_) / (1 + 12.7 * sqrt(friction / 8) * (Pr_**(2 / 3) - 1))
    
    # Correct Nusselt number for low Prandtl numbers, Notter & Sleicher, Chem. Eng. Sci., 1972
    if Pr_ < 0.5:
        Nusselt_lp = 4.8 + 0.0156 * Re_ ** 0.85 * Pr_ ** 0.93
        if Pr_ < 0.1:
            Nusselt = Nusselt_lp
        else:
            Nusselt = Nusselt_lp + (Pr_ - 0.1) * (Nusselt - Nusselt_lp) / 0.4   
            
    # Correct f and Nusselt for developing flow
    f = friction * (1 + (1 / L_ID)**0.7)
    Nu = Nusselt * (1 + (1 / L_ID)**0.7)
    
    return Nu, f

def dP_Pipe(m_dot, fluid, pipe):
    """Calculate pressure drop for flow with heat transfer at the surface of the pipe.
    Section 5.2.3 and 5.2.4 of Nellis and Klein (2020)
    
    Parameters
    ----------
    m_dot : Quantity {mass: 1, time: -1}
        mass flow rate
    fluid : ThermState
        Inlet fluid conditions
    pipe : Pipe

    Returns
    -------
    Quantity {length: -1, mass: 1, time: -2}
        Pressure drop
 
        Heat_transfer_coefficients
    """

    #Parameters
    Re_ = Re(fluid, m_dot, pipe.ID, pipe.area)
    L_ID = pipe.L.m_as(ureg.m)/pipe.ID.m_as(ureg.m)
    eps = (pipe.eps/pipe.ID)
    w = velocity(fluid, m_dot, pipe.area)
    Pr = fluid.Prandtl
    
    #Verify two phase flow and Chockedflow
    if (fluid.phase == 0 or fluid.phase == 6) and fluid.Q < 0.9:
        Phase = 'liquid or two-phase'
    else:
        if Mach(fluid, w) > 1/(fluid.gamma):
            raise ChokedFlow(' Reduce hydraulic resistance or mass flow.')

    # Check Flow conditions
    if Re_ < 0.001 or Re_ > 5e7: 
        raise ValueError(f'Reynolds number (Re) must be > 0.001. The value is {Re_}')
    if eps < 0 or eps > 0.05:
        raise ValueError(f'Relative roughness (eps) should be between 0 and 0.05. The value is {eps}')

    if Re_ > 3000:  # Turbulent flow
        Nu_T, f = turbulent_flow(Re_, Pr, L_ID, eps)
        Nu_Q = Nu_T
        
    elif Re_ < 2300:  # Laminar flow
        Nu_T, Nu_Q, f = laminar_flow(Re_, Pr, L_ID)
        
    else:  # Transitional flow (Re between 2300 and 3000)
        Nu_T_turbulent, f_turbulent = turbulent_flow(3000, Pr, L_ID, eps)
        Nu_lam_T, Nu_lam_Q, f_lam = laminar_flow(2300, Pr, L_ID)
        
        # Interpolate between laminar and turbulent values
        alpha = (Re_ - 2300) / (3000 - 2300)
        Nu_T = Nu_lam_T + alpha * (Nu_T_turbulent - Nu_lam_T)
        Nu_Q = Nu_lam_Q + alpha * (Nu_T_turbulent - Nu_lam_Q)
        f = f_lam + alpha * (f_turbulent - f_lam)
     
    #Pressure drop
    dP = dP_Darcy(f*L_ID, fluid.Dmass, w)    
    
    #Heat transfer
    h_T = heat_trans_coef(fluid, Nu_T, pipe.ID)
    h_Q = heat_trans_coef(fluid, Nu_Q, pipe.ID)
    
    return dP.to(ureg.pascal), h_T, h_Q
