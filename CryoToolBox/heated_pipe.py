"""Pressure drop and heat transfer calculation.
"""

from math import pi, sin, log, log10, sqrt, tan, tanh     
from .std_conditions import ureg, Q_, P_NTP
from .functions import Re, Ra, Nu_vcyl, Nu_hcyl     
from .functions import Material, Property     
from .functions import nist_property, conduction_cyl     
from scipy.optimize import root_scalar, minimize  
from .functions import AIR
from .functions import heat_trans_coef, Ra, Nu_vcyl, Nu_hcyl, Pr
from .cp_wrapper import ThermState
from .piping import Mach, Mach_total, K_lim, ChokedFlow, HydraulicError, velocity, dP_Darcy, dP_adiab, Pipe, Tube, CopperTube

class pipe_isolation:
    ###  class to define the necessary isolation imputs
    def __init__(self, k, OD, T_ext = 293 * ureg.K):
        self.k = k
        self.OD = OD
        self.T_ext = T_ext

def laminar_flow(Re_, Pr_, L_ID):
    # Non dimentional calcul of the Nusselt and fiction factor in pipe in laminar flow following Section 5.2.4 of Nellis and Klein (2020)
    # Verify the input conditions     
    if Pr_ < 0.1:
        raise ValueError(f'Prandtl number (Pr) must be > 0.1. The value is {Pr_}')
        
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
        raise ValueError(f'Prandtl number (Pr) must be between 0.004 and 2000. The value is {Pr_}')
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

def ht_def(pipe):
    """Determine the heated status of the piping component
    """
    ###heat flux defined on the wall of the pipe
    if hasattr(pipe, 'Q_def') and pipe.Q_def != None :
        try: 
            pipe.Q_def.m_as(ureg.W / ureg.m ** 2)
        except:
            raise ValueError(f"the Q_def is not properly defined in component {pipe}" )
        pipe.ht_status = 1
        
    ###Temperature defined on the wall of the pipe
    elif hasattr(pipe, 'Tw_def') and pipe.Tw_def != None :
        try: 
            pipe.Tw_def.m_as(ureg.K)
        except:
            raise ValueError(f"the Tw_def is not properly defined in component {pipe}" )
        pipe.ht_status = 2
        
    ###Heat transfer defined on the wall of the pipe
    elif hasattr(pipe, 'h_ext') and pipe.h_ext != None :
        try: 
            pipe.h_ext.m_as(ureg.W / ureg.m ** 2 / ureg.K)
        except:
            raise ValueError(f"the h_ext is not properly defined in component {pipe}" )
        try:
            pipe.T_ext.m_as(ureg.K)
        except:
            pipe.T_ext = 293 * ureg.K
        pipe.ht_status = 3
        
    ###Ambiant/external temperature defined 
    elif hasattr(pipe, 'T_ext') and pipe.T_ext != None :
        try: 
            pipe.T_ext.m_as(ureg.K)
        except:
            raise ValueError(f"the T_ext is not properly defined in component {pipe}" )
        try:
            pipe.safety_fact 
        except:
            pipe.safety_fact = 1
        pipe.ht_status = 3
        
    ###Isolation on the external of the pipe 
    elif hasattr(pipe, 'isolation') and pipe.isolation != None :
        try: 
            pipe.isolation.k.m_as(ureg.W / ureg.m / ureg.K)
            pipe.isolation.OD.m_as(ureg.m)
        except:
            raise ValueError(f"the isolation (k, OD) is not properly defined in component {pipe}" )
        pipe.ht_status = 4
        
    ###Other
    else: 
        pipe.ht_status = 0

def pipe_Q_def(fluid, pipe, m_dot, dP, h_Q):
    """Calculate the inlet and outlet average temperature of the wall of the component.

    """
    ### Calculate the average temperature of the fluid inside the component
    fluid_temp = fluid.copy()
    dH = (pipe.Q_def * pipe.ID.to(ureg.m) * pipe.L.to(ureg.m) * 3.14) / m_dot
    fluid_temp.update('P', fluid.P - dP,'Hmass', fluid.Hmass + dH.to(ureg.J/ureg.kg))
    T_avg = (fluid.T + fluid_temp.T)/2
    
    ### internal wall temperature
    Tw_i = T_avg + pipe.Q_def/h_Q
    
    ### external wall temperature (root_scalar method, can be changed)
    def find_Tw_o(x):
        k = k_pipe(pipe, Tw_i, x * ureg.K)
        dT1 = (x - Tw_i.m_as(ureg.K)) * ureg.K
        dQ = conduction_cyl(pipe.ID.to(ureg.m), pipe.OD.to(ureg.m), pipe.L.to(ureg.m), k, dT1)
        return dH*m_dot - dQ
    T = root_scalar(find_Tw_o, x0 = Tw_i.m_as(ureg.K)+1, x1 = Tw_i.m_as(ureg.K)+10)
    Tw_o = T.root.magnitude * ureg.K
    
    return Tw_i, Tw_o

def pipe_Tw_def(fluid, pipe, m_dot, dP, h_T):
    """Calculate the inlet and outlet average temperature of the wall of the component,
        as weel as the heat load reaching the fluid.

    """
    H = fluid.Hmass
    fluid_temp = fluid.copy()
    T_avg = fluid.T
    Tw_o = pipe.Tw_def
    res = 1
    j = 0
    
    while res>0.0001:                 
        def find_Tw_i(x):
            k = k_pipe(pipe, Tw_o, x * ureg.K)
            dT1 = x * ureg.K - Tw_o
            dQ = (conduction_cyl(pipe.ID.to(ureg.m), pipe.OD.to(ureg.m), pipe.L.to(ureg.m), k, dT1)).m_as(ureg.watt)
            dT2 = T_avg - x * ureg.K
            dH = (h_T * dT2 * pipe.ID.to(ureg.m) * pipe.L.to(ureg.m) * 3.14).m_as(ureg.watt)
            return (dH ** 2 - dQ ** 2)
        if fluid.T < Tw_o:
            bracket = [T_avg.m_as(ureg.K), Tw_o.m_as(ureg.K)-0.00001]  # Limits search range
        else:
            bracket = [Tw_o.m_as(ureg.K)+0.00001, T_avg.m_as(ureg.K)]
        solution = root_scalar(find_Tw_i, x0 = T_avg.m_as(ureg.K), x1 = (Tw_o.m_as(ureg.K) + T_avg.m_as(ureg.K))/2, bracket=bracket)
        Tw_i = solution.root * ureg.K   
        
        ####### Second part identical
        dT = T_avg - Tw_i
        dH = - h_T * dT * pipe.ID * pipe.L * 3.14 / m_dot
        fluid_temp.update('P', fluid.P - dP,'Hmass', H + dH.to(ureg.J/ureg.kg))
        T_out = fluid_temp.T                
        T_avg_new = (fluid.T + T_out)/2
        res = ((T_avg_new - T_avg)**2 / (T_out - fluid.T)**2)
        T_avg = T_avg_new
        if fluid.T < Tw_o and T_out > Tw_o:
            if j>0:
                raise Exception('the pipe is too long')
            j=j+1
            T_avg = (fluid.T + Tw_o) / 2
        if fluid.T > Tw_o and T_out < Tw_o:
            if j>0:
                raise Exception('the pipe is too long')
            j=j+1
            T_avg = (fluid.T + Tw_o) / 2
            
    Q = (- h_T * dT ).to(ureg.W/ureg.m ** 2)
    
    return Tw_i, Tw_o, Q

def pipe_h_ext(fluid, pipe, m_dot, dP, h_T):
    """Calculate the inlet and outlet average temperature of the wall of the component,
        as weel as the heat load reaching the fluid.

    """
    H = fluid.Hmass
    fluid_temp = fluid.copy()
    T_avg = fluid.T
    res = 1
    j = 0
    
    while res>0.0001:       
        fluid2 = ThermState('air', T= pipe.T_ext, P=1. * ureg.bar)
        
        def find_Tw_o(x):
            h_ext = h_ext_(fluid2, pipe, x * ureg.K)
            T1 = ( h_ext * pipe.OD * (pipe.T_ext - x * ureg.K) / h_T / pipe.ID  ) + T_avg
            if T1 > pipe.T_ext:
                T1 = pipe.T_ext
            if T1 < T_avg:
                T1 = T_avg
            k = k_pipe(pipe, T1, x * ureg.K)
            dT1 = x * ureg.K - T1
            dQ = conduction_cyl(pipe.ID.to(ureg.m), pipe.OD.to(ureg.m), pipe.L.to(ureg.m), k, dT1).m_as(ureg.watt)
            dT2 = T_avg - T1
            dH = - (h_T * dT2 * pipe.ID.to(ureg.m) * pipe.L.to(ureg.m) * 3.14).m_as(ureg.watt)
            return (dH - dQ) ** 2
        T_v2 = minimize(find_Tw_o, x0=(fluid2.T + T_avg) / 2, bounds=[(T_avg.m_as(ureg.K), pipe.T_ext.m_as(ureg.K))])   
        Tw_o = T_v2.x[0] * ureg.K    

        h_ext = h_ext_(fluid2, pipe, Tw_o)
        Tw_i = ( h_ext * pipe.OD * (pipe.T_ext - Tw_o) / h_T / pipe.ID  ) + T_avg
        
        ####### Second part identical
        dT = T_avg - Tw_i
        dH = - h_T * dT * pipe.ID * pipe.L * 3.14 / m_dot
        fluid_temp.update('P', fluid.P - dP,'Hmass', H + dH.to(ureg.J/ureg.kg))
        T_out = fluid_temp.T                
        T_avg_new = (fluid.T + T_out)/2
        res = ((T_avg_new - T_avg)**2 / (T_out - fluid.T)**2)
        T_avg = T_avg_new
        if fluid.T < Tw_o and T_out > Tw_o:
            if j>0:
                raise Exception('the pipe is too long')
            j=j+1
            T_avg = (fluid.T + Tw_o) / 2
        if fluid.T > Tw_o and T_out < Tw_o:
            if j>0:
                raise Exception('the pipe is too long')
            j=j+1
            T_avg = (fluid.T + Tw_o) / 2
            
    Q = (- h_T * dT ).to(ureg.W/ureg.m ** 2)
    
    return Tw_i, Tw_o, Q

def pipe_insulated(fluid, pipe, m_dot, dP, h_T):
    """Calculate the inlet and outlet average temperature of the wall of the component,
        as weel as the heat load reaching the fluid.

    """
    H = fluid.Hmass
    fluid_temp = fluid.copy()
    T_avg = fluid.T
    res = 1
    j = 0
    
    while res>0.0001:   
        
        def find_Tw_o(x):
            T1 = ( pipe.isolation.k * (pipe.isolation.T_ext - x * ureg.K) / h_T / pipe.ID / log(pipe.isolation.OD/pipe.OD) ) + T_avg
            if T1 > pipe.isolation.T_ext:
                T1 = pipe.isolation.T_ext
            if T1 < T_avg:
                T1 = T_avg
            k = k_pipe(pipe, T1, x * ureg.K)
            dT1 = x * ureg.K - T1
            dQ = conduction_cyl(pipe.ID.to(ureg.m), pipe.OD.to(ureg.m), pipe.L.to(ureg.m), k, dT1)
            dT2 = T_avg - T1
            dH = - h_T * dT2 * pipe.ID.to(ureg.m) * pipe.L.to(ureg.m) * 3.14
            return (dH - dQ)**2
        
        fact = (pipe.isolation.OD.m_as(ureg.m) - pipe.OD.m_as(ureg.m)) / pipe.isolation.k.m_as(ureg.W/ureg.m/ureg.K)
        fact = min(1/2, fact)  # Use min for the upper limit      
        T_v2 = minimize(find_Tw_o, x0=pipe.isolation.T_ext.m_as(ureg.K) - (pipe.isolation.T_ext.m_as(ureg.K) - T_avg.m_as(ureg.K)) * fact, bounds=[(T_avg.m_as(ureg.K), pipe.isolation.T_ext.m_as(ureg.K))])
        Tw_o = T_v2.x[0] * ureg.K      
        Tw_i = (( pipe.isolation.k * (pipe.isolation.T_ext - Tw_o) / h_T / pipe.ID / log(pipe.isolation.OD/pipe.OD) ) + T_avg)
        
        ####### Second part identical
        dT = T_avg - Tw_i
        dH = - h_T * dT * pipe.ID * pipe.L * 3.14 / m_dot
        fluid_temp.update('P', fluid.P - dP,'Hmass', H + dH.to(ureg.J/ureg.kg))
        T_out = fluid_temp.T                
        T_avg_new = (fluid.T + T_out)/2
        res = ((T_avg_new - T_avg)**2 / (T_out - fluid.T)**2)
        T_avg = T_avg_new
        if fluid.T < Tw_o and T_out > Tw_o:
            if j>0:
                raise Exception('the pipe is too long')
            j=j+1
            T_avg = (fluid.T + Tw_o) / 2
        if fluid.T > Tw_o and T_out < Tw_o:
            if j>0:
                raise Exception('the pipe is too long')
            j=j+1
            T_avg = (fluid.T + Tw_o) / 2
            
    Q = (- h_T * dT ).to(ureg.W/ureg.m ** 2)
    
    return Tw_i, Tw_o, Q

def pipe_heat(pipe, fluid, m_dot):
    # try:
    #     pipe.ht_status
    # except:
    ht_def(pipe)

    ### Calculate pressure drop and heat transfer coefficient
    dP, h_T, h_Q = dP_Pipe(m_dot, fluid, pipe)  
    
    ### fix heat flux on the wall
    if pipe.ht_status == 1:
        Tw_i, Tw_o = pipe_Q_def(fluid, pipe, m_dot, dP, h_Q)
        Q = pipe.Q_def
        
    elif pipe.ht_status == 2:
        Tw_i, Tw_o, Q = pipe_Tw_def(fluid, pipe, m_dot, dP, h_Q)
        
    elif pipe.ht_status == 3:
        Tw_i, Tw_o, Q = pipe_h_ext(fluid, pipe, m_dot, dP, h_T)
        
    elif pipe.ht_status == 4:
        Tw_i, Tw_o, Q = pipe_insulated(fluid, pipe, m_dot, dP, h_T)
                
    else :
        pipe.Q_def = Q = 0 * ureg.W/ ureg.m ** 2
        Tw_i, Tw_o = pipe_Q_def(fluid, pipe, m_dot, dP, h_Q)

    return Tw_i.to(ureg.K), Tw_o.to(ureg.K), dP.to(ureg.bar), Q.to(ureg.W/ ureg.m ** 2)

def k_pipe(pipe, T_wall, T_ext = 293 * ureg.K):    ### you should add the possibility to define a table of values to do that
    try:
        k = pipe.k.m_as(ureg.W / ureg.K / ureg.m) * ureg.W / ureg.K / ureg.m
    except:
        try:
            mat = pipe.material
        except:
            if isinstance(pipe, CopperTube):
                mat = Material.OFHC
            else:
                mat = Material.SS304
        if T_wall < T_ext:
            k = nist_property(mat, Property.TC, T1 = T_wall, T2 = T_ext)
        else:
            k = nist_property(mat, Property.TC, T1 = T_ext, T2 = T_wall)
    return k.to(ureg.W / ureg.K / ureg.m)

def h_ext_(fluid, pipe, T_wall):  ### you should add the possibility to define a table of values to do that
    try:
        h = pipe.h_ext.m_as(ureg.W / ureg.K / ureg.m ** 2) * ureg.W / ureg.K / ureg.m ** 2
    except:
        try:
            type = pipe.h_type
        except:   
            raise Exception('You should define an external heat transfer type pipe.h_type') #or define type = 1
        if type >= 1:
            try:
                orientation = pipe.orientation
            except:
                orientation = None
            Ra_ = Ra(fluid, T_wall, pipe.OD)
            Pr_ = Pr(fluid)
            if orientation == 'vertical':        
                Nu_ = Nu_vcyl(Pr_, Ra_, pipe.OD, pipe.L)
            elif orientation == 'horizontal':
                Nu_ = Nu_hcyl(Pr_, Ra_)
            else:
                Nu_ = max (Nu_vcyl(Pr_, Ra_, pipe.OD, pipe.L), Nu_hcyl(Pr_, Ra_))
            h = heat_trans_coef(fluid, Nu_, pipe.OD)
        if type == 2:
            try:
                epsilon = pipe.epsilon
            except:
                epsilon = 0.075 # polished stainless steel emissivity
            sigma = 5.670373e-8 * ureg.W/ureg.m ** 2 / ureg.K ** 4
            h_rad = epsilon * sigma * (fluid.T ** 4 - T_wall ** 4) / (fluid.T - T_wall) 
            h = h + h_rad

        if type == 3: ###including Rad ice and h_ice specific in the problem to do 
            print('to do')
    return h.to(ureg.W / ureg.K / ureg.m ** 2)