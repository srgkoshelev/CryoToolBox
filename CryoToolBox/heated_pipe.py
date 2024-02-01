# -*- coding: utf-8 -*-
"""Utilities for hydraulics calculations on heated pipe.

Contains functions for thermo hydrodynamic calculations. 
"""
from math import pi, sin, log, log10, sqrt, tan, tanh
from .std_conditions import ureg, Q_, P_NTP
from .functions import Re, Ra, Pr, Nu_vcyl, Nu_hcyl_2
from .functions import Material, Property
from .functions import nist_property, conduction_cyl
from scipy.optimize import root_scalar, minimize
from .functions import AIR
from .cp_wrapper import ThermState
from .piping import Mach_total, K_lim, ChokedFlow, HydraulicError, velocity, dP_Darcy, heat_trans_coef, dP_adiab, Pipe, Tube, CopperTube

### New class extension for pipe, copper tube and tube
class Tube_HP(Tube):
    def __init__(self, OD, wall=0*ureg.m, L=0*ureg.m, c=Q_('0 mm'),
                 eps=0.0018*ureg.inch, material = Material.SS304, 
                 T_p_out = 0 * ureg.K, T_p_out_Q = 0 * ureg.K, 
                 T_p_in = 0 * ureg.K, Q = 0 * ureg.W/ureg.m ** 2):
        # Pipe class constructor
        super().__init__(OD, wall, L, c, eps)

        # new attribute 
        self.material = material
        self.T_p_out = T_p_out
        self.T_p_out_Q = T_p_out_Q
        self.T_p_in = T_p_in
        self.Q = Q

class Pipe_HP(Pipe):
    def __init__(self, D_nom, SCH=40, L=0*ureg.m, c=Q_('0 mm'),
                 eps=0.0018*ureg.inch, material = Material.SS304, 
                 T_p_out = 0 * ureg.K, T_p_out_Q = 0 * ureg.K, 
                 T_p_in = 0 * ureg.K, Q = 0 * ureg.W/ureg.m ** 2):
        # Pipe class constructor

        super().__init__(D_nom, SCH, L, c, eps)

        # new attribute 
        self.material = material
        self.T_p_out = T_p_out
        self.T_p_out_Q = T_p_out_Q
        self.T_p_in = T_p_in
        self.Q = Q
        
class CopperTube_HP(CopperTube):
    def __init__(self, D_nom, type_='K', L=0*ureg.m,
                 eps=0.0018*ureg.inch, material = Material.OFHC, 
                 RRR_OFHC=None, T_p_out = 0 * ureg.K, 
                 T_p_out_Q = 0 * ureg.K, T_p_in = 0 * ureg.K, 
                 Q = 0 * ureg.W/ureg.m ** 2):
        # Pipe class constructor
        super().__init__(D_nom, type_, L, eps)

        # new attribute 
        self.material = material
        self.RRR_OFHC = RRR_OFHC
        self.T_p_out = T_p_out
        self.T_p_out_Q = T_p_out_Q
        self.T_p_in = T_p_in
        self.Q = Q

###additional fct from EES for pipes with heat transfer from EES

def pipeflow_turbulent(Re, Pr, L_D, relRough):
    # Check input conditions
    if Re > 5e7:
        raise ValueError(f'Reynolds number (Re) must be < 5E7 for Nusselt number and < 1E8 for the friction factor. The value is {Re}')
    if Re < 2300:
        raise ValueError(f'Reynolds number (Re) must be > 2300. The value is {Re}')
    if Pr < 0.004 or Pr > 2000:
        raise ValueError(f'Prandtl number (Pr) must be between 0.004 and 2000. The value is {Pr}')
    if L_D <= 1:
        if L_D < 0:
            raise ValueError('L/D ratio should be > 1. The value is {L_D}')
        L_D = 1
    if relRough < 0 or relRough > 0.05:
        raise ValueError(f'Relative roughness (relRough) should be between 0 and 0.05. The value is {relRough}')

    # Calculate f_fd using the equation for turbulent flow
    f_fd = (-0.001570232 / log(Re) + 0.394203137 / log(Re)**2 + 2.534153311 / log(Re)**3) * 4

    # Calculate f using Zigrang & Sylvester correlation
    if relRough > 1e-5:
        f_fd = (-2 * log10(relRough / 3.71 - 1.975 / Re * log((relRough / 3.93)**1.092 + 7.627 / (Re + 395.9))))**(-2)

    # Calculate Nusselt number for laminar flow
    Nusselt_L = ((f_fd / 8) * (Re - 1000) * Pr) / (1 + 12.7 * sqrt(f_fd / 8) * (Pr**(2 / 3) - 1))

    # Correct Nusselt number for low Prandtl numbers
    if Pr < 0.5:
        Nusselt_L_lp = 4.8 + 0.0156 * Re**0.85 * Pr**0.93
        if Pr < 0.1:
            Nusselt_L = Nusselt_L_lp
        else:
            Nusselt_L = Nusselt_L_lp + (Pr - 0.1) * (Nusselt_L - Nusselt_L_lp) / 0.4

    # Correct f for developing flow
    f = f_fd * (1 + (1 / L_D)**0.7)

    # Correct Nusselt number for developing flow
    Nusselt = Nusselt_L * (1 + (1 / L_D)**0.7)

    return Nusselt, f

def pipeflow_laminar(Re, Pr, L_D):
    # Check input conditions
    if Re > 2300:
        raise ValueError(f'Reynolds number (Re) must be < 2300. The value is {Re}')
    if Pr < 0.1:
        raise ValueError(f'Prandtl number (Pr) must be > 0.1. The value is {Pr}')
    
    Z_H = L_D / (Re * Pr)
    if Z_H < 1e-6:
        raise ValueError(f'Inverse Graetz number (Z_H) must be > 1e-6. The value is {Z_H}')
    
    Z_M = L_D / Re

    # Calculate friction factor (f)
    f = 4 * (3.44 / sqrt(Z_M) + (1.25 / (4 * Z_M) + 16 - 3.44 / sqrt(Z_M)) / (1 + 0.00021 * Z_M**(-2))) / Re

    # Calculate Nusselt numbers
    Nusselt_T = ((5.001 / Z_H**1.119 + 136.0)**0.2978 - 0.6628) / tanh(2.444 * Z_M**(1 / 6) * (1 + 0.565 * Z_M**(1 / 3)))
    Nusselt_H = ((6.562 / Z_H**1.137 + 220.4)**0.2932 - 0.5003) / tanh(2.530 * Z_M**(1 / 6) * (1 + 0.639 * Z_M**(1 / 3)))

    return Nusselt_T, Nusselt_H, f


def pipeflow_nd(Re, Pr, L_D, relRough):
    # Check input conditions
    if Re < 0.001:
        raise ValueError(f'Reynolds number (Re) must be > 0.001. The value is {Re}')
    if relRough < 0 or relRough > 0.05:
        raise ValueError(f'Relative roughness (relRough) should be between 0 and 0.05. The value is {relRough}')

    if Re > 3000:  # Turbulent flow
        Nusselt_T, f = pipeflow_turbulent(Re, Pr, L_D, relRough)
        Nusselt_H = Nusselt_T
    elif Re < 2300:  # Laminar flow
        Nusselt_T, Nusselt_H, f = pipeflow_laminar(Re, Pr, L_D)
    else:  # Transitional flow (Re between 2300 and 3000)
        Nusselt_T_turbulent, f_turbulent = pipeflow_turbulent(3000, Pr, L_D, relRough)
        Nusselt_lam_T, Nusselt_lam_H, f_lam = pipeflow_laminar(2300, Pr, L_D)
        
        # Interpolate between laminar and turbulent values
        alpha = (Re - 2300) / (3000 - 2300)
        Nusselt_T = Nusselt_lam_T + alpha * (Nusselt_T_turbulent - Nusselt_lam_T)
        Nusselt_H = Nusselt_lam_H + alpha * (Nusselt_T_turbulent - Nusselt_lam_H)
        f = f_lam + alpha * (f_turbulent - f_lam)

    return Nusselt_T, Nusselt_H, f

def dP_HT(m_dot, fluid, pipe):
    """Calculate pressure drop for flow with heat transfer at the surface.

    """

    Re_ = Re(fluid, m_dot, pipe.ID, pipe.area)
    K_pipe = pipe.K(Re_)
    if (fluid.phase == 0 or fluid.phase == 6) and fluid.Q < 0.9:
        Phase = 'liquid or two-phase'
    else:
        M = Mach_total(fluid, m_dot, pipe.area)
        K_limit = K_lim(M, fluid.gamma)
        K_left = K_limit - K_pipe
        if K_left < 0:
            raise ChokedFlow(f'Flow is choked at K={float(K_limit):.3g} with given '
                             f'K={float(K_pipe):.3g}. Reduce hydraulic resistance or'
                             ' mass flow.')
            
    L_D = pipe.L.m_as(ureg.m)
    I_D = pipe.ID.m_as(ureg.m)
    relRough = pipe.eps.m_as(ureg.m)
    w = velocity(fluid, m_dot, pipe.area)
    Nusselt_T, Nusselt_H, f = pipeflow_nd(Re_, fluid.Prandtl, L_D, relRough)
    dP = dP_Darcy(f*L_D/I_D, fluid.Dmass, w)

    h_T = heat_trans_coef(fluid, Nusselt_T, pipe.ID)
    h_H = heat_trans_coef(fluid, Nusselt_H, pipe.ID)

    if pipe.T_p_out.m_as(ureg.K) != 0:
        try:    
            def find_Tp(x):
                property = Property.TC  # Thermal conductivity
                k = nist_property(pipe.material, property, T1 = pipe.T_p_out, T2 = x * ureg.K)
                dT1 = (x - pipe.T_p_out.m_as(ureg.K)) * ureg.K
                dQ = - conduction_cyl(pipe.ID.to(ureg.m), pipe.OD.to(ureg.m), pipe.L.to(ureg.m), k, dT1)
                #conduction_cyl(D_i, D_o, L, k, dT)
                dT2 = (fluid.T.m_as(ureg.K) - x) * ureg.K
                dH = - h_T * dT2 * pipe.ID.to(ureg.m) * pipe.L.to(ureg.m) * 3.14
                return dH - dQ
            T = root_scalar(find_Tp, x0 = fluid.T.m_as(ureg.K), x1 = (pipe.T_p_out.m_as(ureg.K) + fluid.T.m_as(ureg.K))/2)
            pipe.T_p_in = T.root.magnitude * ureg.K
            dT = (fluid.T.m_as(ureg.K) - T.root.magnitude) * ureg.K
        except:
            print('material do not exist')
            dT = (fluid.T.m_as(ureg.K) - pipe.T_p_out.m_as(ureg.K)) * ureg.K
        dH = - h_T * dT * pipe.ID.to(ureg.m) * pipe.L.to(ureg.m) * 3.14
    elif pipe.Q.m_as(ureg.W/ureg.m ** 2) != 0:  
        dH = pipe.Q * pipe.ID.to(ureg.m) * pipe.L.to(ureg.m) * 3.14
        pipe.T_p_in = fluid.T + pipe.Q/h_H
        try:    
            def find_Tp(x):
                property = Property.TC  # Thermal conductivity
                k = nist_property(pipe.material, property, T1 = pipe.T_p_in, T2 = x * ureg.K)
                dT1 = (x - pipe.T_p_in.m_as(ureg.K)) * ureg.K
                dQ = conduction_cyl(pipe.ID.to(ureg.m), pipe.OD.to(ureg.m), pipe.L.to(ureg.m), k, dT1)
                return dH - dQ
            T = root_scalar(find_Tp, x0 = pipe.T_p_in.m_as(ureg.K)+1, x1 = pipe.T_p_in.m_as(ureg.K)+10)
            pipe.T_p_out_Q = T.root.magnitude * ureg.K
            dT = (fluid.T.m_as(ureg.K) - T.root.magnitude) * ureg.K
        except:
            print('material do not exist')
    else:
        dH = 0 * ureg.W

    return dP.to(ureg.pascal), dH

def dP_HT_Foam(m_dot, fluid, pipe, k_Foam, OD_Foam, T_ambient = 300 * ureg.K):
    """Calculate pressure drop for flow with heat transfer at the surface.

    """
    Re_ = Re(fluid, m_dot, pipe.ID, pipe.area) 
    K_pipe = pipe.K(Re_)
    if (fluid.phase == 0 or fluid.phase == 6) and fluid.Q < 0.9:
        Phase = 'liquid or two-phase'
    else:
        M = Mach_total(fluid, m_dot, pipe.area)
        K_limit = K_lim(M, fluid.gamma)
        K_left = K_limit - K_pipe
        if K_left < 0:
            raise ChokedFlow(f'Flow is choked at K={float(K_limit):.3g} with given '
                             f'K={float(K_pipe):.3g}. Reduce hydraulic resistance or'
                             ' mass flow.')

    L_D = pipe.L.m_as(ureg.m)
    I_D = pipe.ID.m_as(ureg.m)
    relRough = pipe.eps.m_as(ureg.m)
    w = velocity(fluid, m_dot, pipe.area)
    Nusselt_T, Nusselt_H, f = pipeflow_nd(Re_, fluid.Prandtl, L_D, relRough)
    dP = dP_Darcy(f*L_D/I_D, fluid.Dmass, w)

    h_T = heat_trans_coef(fluid, Nusselt_T, pipe.ID)
    #h_H = heat_trans_coef(fluid, Nusselt_H, pipe.ID)

    try:    
        def find_Tp(x):
            T1 = ( k_Foam * (T_ambient - x * ureg.K) / h_T / pipe.ID / log(OD_Foam/pipe.OD) ) + fluid.T
            if T1 > T_ambient:
                T1 = T_ambient
            if T1 < fluid.T:
                T1 = fluid.T
            #print(T1.to(ureg.K))
            property = Property.TC  # Thermal conductivity
            k = nist_property(pipe.material, property, T1 = T1, T2 = x * ureg.K)
            #print(x * ureg.K)
            dT1 = x * ureg.K - T1
            dQ = conduction_cyl(pipe.ID.to(ureg.m), pipe.OD.to(ureg.m), pipe.L.to(ureg.m), k, dT1)
            #print(dQ)
            #conduction_cyl(D_i, D_o, L, k, dT)
            dT2 = fluid.T - T1
            dH = - h_T * dT2 * pipe.ID.to(ureg.m) * pipe.L.to(ureg.m) * 3.14
            #print(dH)
            return (dH - dQ)**2
        fact = (OD_Foam.m_as(ureg.m) - pipe.OD.m_as(ureg.m)) / k_Foam.m_as(ureg.W/ureg.m/ureg.K)
        fact = min(1/2, fact)  # Use min for the upper limit      
        T_v2 = minimize(find_Tp, x0=T_ambient.m_as(ureg.K) - (T_ambient.m_as(ureg.K) - fluid.T.m_as(ureg.K)) * fact, bounds=[(fluid.T.m_as(ureg.K), T_ambient.m_as(ureg.K))])
        pipe.T_p_out_Q = T_v2.x[0] * ureg.K      
        pipe.T_p_in = (( k_Foam * (T_ambient - pipe.T_p_out_Q) / h_T / pipe.ID / log(OD_Foam/pipe.OD) ) + fluid.T).to(ureg.K)
        dT = (fluid.T - pipe.T_p_in) 
    except:
        print('material do not exist')
        dT = (fluid.T.m_as(ureg.K) - pipe.T_p_out.m_as(ureg.K)) * ureg.K
    dH = - h_T * dT * pipe.ID.to(ureg.m) * pipe.L.to(ureg.m) * 3.14

    return dP.to(ureg.pascal), dH

def dP_HT_ambient(m_dot, fluid, pipe, safety_fact = 1, T_ambient = 300 * ureg.K):
    """Calculate pressure drop for flow with heat transfer at the surface.

    """
    Re_ = Re(fluid, m_dot, pipe.ID, pipe.area)
    K_pipe = pipe.K(Re_)
    if (fluid.phase == 0 or fluid.phase == 6) and fluid.Q < 0.9:
        Phase = 'liquid or two-phase'
    else:
        M = Mach_total(fluid, m_dot, pipe.area)
        K_limit = K_lim(M, fluid.gamma)
        K_left = K_limit - K_pipe
        if K_left < 0:
            raise ChokedFlow(f'Flow is choked at K={float(K_limit):.3g} with given '
                             f'K={float(K_pipe):.3g}. Reduce hydraulic resistance or'
                             ' mass flow.')

    L_D = pipe.L.m_as(ureg.m)
    I_D = pipe.ID.m_as(ureg.m)
    relRough = pipe.eps.m_as(ureg.m)
    w = velocity(fluid, m_dot, pipe.area)
    Nusselt_T, Nusselt_H, f = pipeflow_nd(Re_, fluid.Prandtl, L_D, relRough)
    dP = dP_Darcy(f*L_D/I_D, fluid.Dmass, w)

    h_T = heat_trans_coef(fluid, Nusselt_T, pipe.ID)
    #h_H = heat_trans_coef(fluid, Nusselt_H, pipe.ID)
    
    fluid2 = ThermState('air', T= T_ambient, P=1. * ureg.bar)

    try:    
        def find_Tp(x):
            Ra1 = Ra(fluid2, x * ureg.K, pipe.OD)
            Pr1 = Pr(fluid2)
            Nu1 = max (Nu_vcyl(Pr1, Ra1, pipe.OD, pipe.L), Nu_hcyl_2(Pr1, Ra1))
            h_ext = heat_trans_coef(fluid2, Nu1 * safety_fact, pipe.OD)
            T1 = ( h_ext * pipe.OD * (T_ambient - x * ureg.K) / h_T / pipe.ID  ) + fluid.T
            if T1 > T_ambient:
                T1 = T_ambient
            if T1 < fluid.T:
                T1 = fluid.T
            #print(T1.to(ureg.K))
            property = Property.TC  # Thermal conductivity
            k = nist_property(pipe.material, property, T1 = T1, T2 = x * ureg.K)
            #print(x * ureg.K)
            dT1 = x * ureg.K - T1
            dQ = conduction_cyl(pipe.ID.to(ureg.m), pipe.OD.to(ureg.m), pipe.L.to(ureg.m), k, dT1)
            #print(dQ)
            #conduction_cyl(D_i, D_o, L, k, dT)
            dT2 = fluid.T - T1
            dH = - h_T * dT2 * pipe.ID.to(ureg.m) * pipe.L.to(ureg.m) * 3.14
            return (dH - dQ)**2

        #fact = min(1/2, h_ext.m_as(ureg.W/(ureg.m**2*ureg.K)))  # Use min for the upper limit      
        T_v2 = minimize(find_Tp, x0=(fluid2.T + fluid.T) / 2, bounds=[(fluid.T.m_as(ureg.K), T_ambient.m_as(ureg.K))])
        #print(T_v2)
        #print(AIR.T)
        pipe.T_p_out_Q = T_v2.x[0] * ureg.K
        #print(pipe.T_p_out_Q)
        Ra1 = Ra(fluid2, pipe.T_p_out_Q, pipe.OD)
        Pr1 = Pr(fluid2)
        Nu1 = max (Nu_vcyl(Pr1, Ra1, pipe.OD, pipe.L), Nu_hcyl_2(Pr1, Ra1))
        h_ext = heat_trans_coef(fluid2, Nu1 * safety_fact, pipe.OD)
        pipe.T_p_in = ( h_ext * pipe.OD * (T_ambient - pipe.T_p_out_Q) / h_T / pipe.ID  ) + fluid.T
        #print(pipe.T_p_in)        
        dT = (fluid.T - pipe.T_p_in) 
    except:
        print('material do not exist')
        dT = (fluid.T.m_as(ureg.K) - pipe.T_p_out.m_as(ureg.K)) * ureg.K
    dH = - h_T * dT * pipe.ID.to(ureg.m) * pipe.L.to(ureg.m) * 3.14

    return dP.to(ureg.pascal), dH
