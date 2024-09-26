
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
import numpy as np

class pipe_isolation:
    ### class to define the necessary isolation imputs
    # testing
    def __init__(self, k, OD, T_ext = 293 * ureg.K):
        self.k = k
        self.OD = OD
        self.T_ext = T_ext

def laminar_flow(Re_, Pr_, L_ID):
    # Non dimentional calculation of the Nusselt and fiction factor in pipe in laminar flow  
    # Section 5.2.4 of Nellis and Klein (2020)
    
    # Verify the input conditions     
    if Pr_ < 0.1:
        raise ValueError(f'Prandtl number (Pr) must be > 0.1. The value is {Pr}')
        
    # Calculate Graetz number and Inverse Graetz number verification
    SGZ = L_ID / Re_
    GZ = L_ID / (Re_ * Pr_)
    if GZ < 1e-6:
        raise ValueError(f'Inverse Graetz number (GZ) must be > 1e-6. The value is {GZ}')
         
    # Calculate friction factor (f)
    f = 4 * (3.44 / sqrt(SGZ) + (1.25 / (4 * SGZ) + 16 - 3.44 / sqrt(SGZ)) / (1 + 0.00021 * SGZ**(-2))) / Re_
    
    # Calculate Nusselt numbers: temperature constant and flux constant
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
        Defining the pipe characteristics 
    Returns
    -------
    dP: Quantity {length: -1, mass: 1, time: -2}
        Pressure drop
    h_T, h_Q: Quantity: {mass : 1, temperature : -1, time : -3}
        Heat transfer coefficients
    """

    # Calculate fluid pameters
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
        raise ValueError(f'Reynolds number (Re) must be between 0.001 and 5E7. The value is {Re_}')
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
     
    # Calculate pressure drop
    dP = dP_Darcy(f*L_ID, fluid.Dmass, w)    
    
    # Calculate heat transfer coefficient: temperature and heat flux constant 
    h_T = heat_trans_coef(fluid, Nu_T, pipe.ID)
    h_Q = heat_trans_coef(fluid, Nu_Q, pipe.ID)
    
    return dP.to(ureg.pascal), h_T, h_Q


def find_Tw(x, T_avg, pipe, h_coeff, m_dot):
  
    """Calculate the average temperature of the inner or outer wall of the component.
    
    Parameters
    ----------
    dH : Quantity {length: 2, mass: 1, time: -2}
        heat rate calculated from convection equation 
    dQ : Quantity {length: 2, mass: 1, time: -2}
        heat rate calculated from thermal resistance equation
    dT : Quantity {temperature: 1}
        temperature difference
    fluid : ThermState
        Inlet fluid conditions
    fluid_external : Thermstate
        External fluid conditions
    h_coeff: Quantity {mass : 1, temperature : -1, time : -3}
            heat transfer coefficient: chosen to be either h_T or h_Q
    h_ext: Quantity {mass : 1, temperature : -1, time : -3}
            heat transfer coefficient for external fluid
    k : Quantity {length: 1, mass: 1, temperature: -1, time: -3}
        material thermal conductivity
    m_dot : Quantity { mass: 1, time: -1}
            mass flow rate    
    pipe : Pipe
        defining the pipe characteristics      
    Q_def : Quantity {mass: 1, time: -3}
            Heat load reaching the fluid  
    Tw_i: Quantity {temperature: 1}
        inside temperature of wall  
    Tw_o: Quantity {temperature: 1}
        outside temperature of wall         
    Returns
    -------
    (dH - dQ) ** 2 : Equation
            Quadratic Expression that computes wall temperature when the minimum is solved
    """
    
    if pipe.Q_def != None:
        #For a system with defined heat load: pipe_Q_def
        Tw_i = T_avg + pipe.Q_def/h_coeff
        Tw_o =  x * ureg.K
        
    elif pipe.Tw_def != None:
        #For a system with defined external temperature: pipe_Tw_def
        Tw_i = x * ureg.K
        Tw_o = pipe.Tw_def
        
    elif pipe.T_ext != None:
        #For a system with defined external heat transfer coeff: pipe_h_ext
        fluid_external = ThermState('air', T= pipe.T_ext, P=1 * ureg.bar) #to do: improve structure 
        h_ext = h_ext_(fluid_external, pipe, x * ureg.K)     
        
        Tw_i = (h_ext * pipe.OD * (pipe.T_ext - x * ureg.K) / h_coeff / pipe.ID) + T_avg
        Tw_i = max(min(Tw_i, pipe.T_ext), T_avg)
        Tw_o =  x * ureg.K
        
    elif pipe.isolation.k != None: 
        #For a defined insulated system: pipe_insulated
        Tw_i = (pipe.isolation.k * (pipe.isolation.T_ext - x * ureg.K) / h_coeff / pipe.ID / log(pipe.isolation.OD / pipe.OD)) + T_avg
        Tw_i = max(min(Tw_i, pipe.isolation.T_ext), T_avg)
        Tw_o =  x * ureg.K
    else:
        raise ValueError("Insufficient or invalid parameters provided.")                
    
    k = k_pipe(pipe, Tw_o, Tw_i)
    dQ = conduction_cyl(pipe.ID.to(ureg.m), pipe.OD.to(ureg.m), pipe.L.to(ureg.m), k, (Tw_o - Tw_i))
    dH = -h_coeff * (T_avg - Tw_i) * pipe.ID.to(ureg.m) * pipe.L.to(ureg.m) * 3.14  
    return (dH - dQ).m_as(ureg.W) ** 2  
    


def pipe_Q_def(fluid, pipe, m_dot, dP, h_Q):
    
    """Calculate the inner and outer wall temperatures of a system with a defined heat flux.

    Parameters
    ----------
    dH : Quantity {length: 2, time: -2}
        specific enthalpy of fluid    

    fluid : ThermState
        Inlet fluid conditions
    fluid_downstream : Thermstate
        Outlet fluid conditions
    m_dot : Quantity { mass: 1, time: -1}
            mass flow rate       
    pipe : Pipe
        Defining the pipe characteristics 
    Q_def : Quantity {mass: 1, time: -3}
            Heat load reaching the fluid
    T_avg : Quantity {temperature: 1}
            average fluid temperature 
    Returns
    -------
    Tw_i, Tw_o : Quantity {temperature: 1}
        Inside temperature of the wall, Outside temperature of the wall
    """
    #Calculate downstream conditions
    fluid_downstream = fluid.copy()
    dH = (pipe.Q_def * pipe.ID.to(ureg.m) * pipe.L.to(ureg.m) * 3.14) / m_dot 
    fluid_downstream.update('P', fluid.P - dP, 'Hmass' , fluid.Hmass + dH.to(ureg.J/ureg.kg)) 
    
    ##Calculate the average temperature of the fluid inside the component
    T_avg = (fluid.T + fluid_downstream.T)/2
          
    #Calculate Tw_i and Tw_o: minimum of the quadratic find_Tw
    Tw_i = T_avg + pipe.Q_def/h_Q
    dT = log(pipe.OD/pipe.ID) * (pipe.Q_def * pipe.OD ) / (2 * k_pipe(pipe, Tw_i))
    # Q * pi * D_o * L = 2 * pi * L * k * dT / log(D_o/D_i)
    
    if pipe.Q_def > 0 * ureg.W / ureg.m ** 2 :
        bounds = [(Tw_i, Tw_i + dT * 2)]
    else:
        bounds =[(Tw_i - dT * 2, Tw_i)]
    Tw_o = minimize(find_Tw, x0=T_avg.m_as(ureg.K) + 1, args=(T_avg, pipe, h_Q, m_dot), bounds=bounds).x[0] * ureg.K
    
    return Tw_i, Tw_o  
 

# Define main functions
def pipe_Tw_def(fluid, pipe, m_dot, dP, h_T):
    """Calculate the heat flux and the inner wall temperature for a pipe with a defined external wall temperature. 
    Parameters
    ----------
    dH : Quantity {length: 2, time: -2}
        specific enthalpy of fluid
    dT : Quantity {temperature: 1}
        temperature difference
    fluid : ThermState
       Inlet fluid conditions
    fluid_downstream : Thermstate
       Outlet fluid conditions
    h_T : Quantity : {mass : 1, temperature : -1, time : -3}
            heat transfer coefficient 
    m_dot : Quantity { mass: 1, time: -1}
            mass flow rate            
    pipe : Pipe
        Defining the pipe characteristics 
    T_ds : Quantity {temperature: 1}
            temperature of fluid downstream     
    Returns
    -------
    Tw_i and Tw_o : Quantity {temperature: 1}
        Inside temperature of the wall, Outside temperature of the wall
    Q : Quantity { mass: 1, time: -3}
        Heat load reaching the fluid
    """
    #### Initial conditions and parameters
    H = fluid.Hmass
    fluid_downstream = fluid.copy()
    T_avg = fluid.T
    Tw_o = pipe.Tw_def
    res = 1
    j = 0
    dT = 'none'
    
    while res>0.0001:  

        # Limits search range for Tw_o
        if T_avg < Tw_o:
            bracket = (T_avg.m_as(ureg.K), Tw_o.m_as(ureg.K)) 
        else:
            bracket = (Tw_o.m_as(ureg.K), T_avg.m_as(ureg.K))

        # Calculate Tw_i:  minimum of the quadratic find_Tw
        Tw_i = minimize(find_Tw,  x0 = T_avg.m_as(ureg.K), args = (T_avg, pipe, h_T, m_dot), bounds=[bracket]).x[0] *ureg.K        

        ### Calculate downstream fluid conditions
        dT = T_avg - Tw_i
        dH = - h_T * dT * pipe.ID * pipe.L * 3.14 / m_dot
        fluid_downstream.update('P', fluid.P - dP,'Hmass', H + dH.to(ureg.J/ureg.kg))
        T_ds = fluid_downstream.T 

        ###Check convergence of T_average               
        T_avg_new = (fluid.T + T_ds)/2
        res = ((T_avg_new - T_avg)**2 / (T_ds - fluid.T)**2)
        
        ##Update average temperature
        T_avg = T_avg_new
        
        ### Eliminate nonphysical solutions         
        if (fluid.T < Tw_o and T_ds > Tw_o) or (fluid.T > Tw_o and T_ds < Tw_o):
            if j > 0:
                raise Exception('the pipe is too long')
            j += 1
            T_avg = (fluid.T + Tw_o) / 2

        ### Calculate heat flux 
        Q = (- h_T * dT ).to(ureg.W/ureg.m ** 2)       
 
        return Tw_i, Tw_o, Q


def pipe_h_ext(fluid, pipe, m_dot, dP, h_T): 
    """Calculate the heat flux and the inner and outer wall temperatures for a pipe with a defined
        external heat transfer coefficient.
     Parameters
     ----------
     dH : Quantity {length: 2, time: -2}
         specific enthalpy of fluid
     dT : Quantity {temperature: 1}
         temperature difference
     fluid : ThermState
        Inlet fluid conditions
     fluid_downstream : Thermstate
         outlet fluid conditions
     fluid_external : Thermstate
         external fluid conditions
     h_ext : Quantity : {mass : 1, temperature : -1, time : -3}
        transfer coefficient of external fluid
     m_dot : Quantity { mass: 1, time: -1}
        mass flow rate        
     pipe : Pipe
        defining the pipe characteristics 
    T_avg: Quantity {temperature: 1}
        average temperature of the fluid
    T_ds : Quantity {temperature: 1}
        temperature of downstream fluid
     Returns
     -------
     Tw_i and Tw_o : Quantity {temperature: 1}
         Inside temperature of the wall, Outside temperature of the wall
     Q : Quantity { mass: 1, time: -3}
         Heat load reaching the fluid
    """
    ### Initial conditions and parameters
    H = fluid.Hmass
    fluid_downstream = fluid.copy()
    T_avg = fluid.T
    res = 1
    j = 0
    
    while res>0.0001:  
        
        ##Define external fluid
        fluid_external = ThermState('air', T= pipe.T_ext, P=1 * ureg.bar) 
        
        # Limits search range for Tw_o
        if T_avg < pipe.T_ext:   
            bracket = [(T_avg.m_as(ureg.K)+0.0001, pipe.T_ext.m_as(ureg.K)-0.0001)]   
        else:
            bracket = [(pipe.T_ext.m_as(ureg.K)+0.0001, T_avg.m_as(ureg.K)-0.0001)]       
        
        #Calculate Tw_o: minimum of the quadratic find_Tw_o
        Tw_o = minimize(find_Tw, x0=(fluid_external.T + T_avg).m_as(ureg.K)/ 2, args = (T_avg, pipe, h_T, m_dot), bounds=bracket).x[0] * ureg.K                                                  
        
        # Caclulate external heat transfer coefficient for system: calculates h_ext for system with h_type defined otherwise uses defined h_ext
        h_ext = h_ext_(fluid_external, pipe, Tw_o)
        
        # Calculatue inner wall temperature
        Tw_i = ( h_ext * pipe.OD * (pipe.T_ext - Tw_o) / h_T / pipe.ID  ) + T_avg
        
        ###Calculate downstream flow conditions
        dT = T_avg - Tw_i
        dH = - h_T * dT * pipe.ID * pipe.L * 3.14 / m_dot
        fluid_downstream.update('P', fluid.P - dP,'Hmass', H + dH.to(ureg.J/ureg.kg))
        T_ds = fluid_downstream.T 

        ###Check convergence of T_average and update T_avg value            
        T_avg_new = (fluid.T + T_ds)/2
        res = ((T_avg_new - T_avg)**2 / (T_ds - fluid.T)**2)
        T_avg = T_avg_new
    
        ###Eliminate nonphysical solutions         
        if (fluid.T < Tw_o and T_ds > Tw_o) or (fluid.T > Tw_o and T_ds < Tw_o):
            if j > 0:
                raise Exception('the pipe is too long')
            j += 1
            T_avg = (fluid.T + Tw_o) / 2
       
        ###Calculate heat flux  
        Q = (- h_T * dT ).to(ureg.W/ureg.m ** 2)
    
        return Tw_i, Tw_o, Q

def pipe_insulated(fluid, pipe, m_dot, dP, h_T): 
    """Calculate the heat flux and the inner and outer wall temperatures for a pipe with defined insulation.
     Parameters
     ----------
     dH : Quantity {length: 2, time: -2}
         specific enthalpy of fluid
     dT : Quantity {temperature: 1}
         temperature difference
     fluid : ThermState
        Inlet fluid conditions
     fluid_downstream : Thermstate
         outlet fluid conditions
     h_ext : Quantity : {mass : 1, temperature : -1, time : -3}
        transfer coefficient of external fluid
     m_dot : Quantity { mass: 1, time: -1}
        mass flow rate        
     pipe : Pipe
        defining the pipe characteristics 
     T_avg: Quantity {temperature: 1}
        average temperature of the fluid
     T_ds : Quantity {temperature: 1}
        temperature of downstream fluid
     Returns
     -------
     Tw_i and Tw_o : Quantity {temperature: 1}
         Inside temperature of the wall, Outside temperature of the wall
     Q : Quantity { mass: 1, time: -3}
         Heat load reaching the fluid
    """
    ###Initial conditions and parameters
    H = fluid.Hmass
    fluid_downstream = fluid.copy()
    T_avg = fluid.T
    res = 1
    j = 0
    
    while res>0.0001:   
        # Limits search range for Tw_o       
        if T_avg < pipe.isolation.T_ext:  
            bracket = [(T_avg.m_as(ureg.K)+0.0001, pipe.isolation.T_ext.m_as(ureg.K)-0.0001)]  
        else:
            bracket = [(pipe.isolation.T_ext.m_as(ureg.K)+0.0001, T_avg.m_as(ureg.K)-0.0001)]
         
        #Calculate Tw_i and Tw_o: minimum of the quadratic find_Tw
        Tw_o = minimize(find_Tw, x0=pipe.isolation.T_ext.m_as(ureg.K) - (pipe.isolation.T_ext.m_as(ureg.K) - T_avg.m_as(ureg.K)), args = (T_avg, pipe, h_T, m_dot), bounds=bracket).x[0] * ureg.K 
        Tw_i = (( pipe.isolation.k * (pipe.isolation.T_ext - Tw_o) / h_T / pipe.ID / log(pipe.isolation.OD/pipe.OD) ) + T_avg)
        
        ###Calculate downstream flow conditions
        dT = T_avg - Tw_i
        dH = - h_T * dT * pipe.ID * pipe.L * 3.14 / m_dot
        fluid_downstream.update('P', fluid.P - dP,'Hmass', H + dH.to(ureg.J/ureg.kg))
        T_ds = fluid_downstream.T 

        ###Check convergence of T_average and update T_avg value       
        T_avg_new = (fluid.T + T_ds)/2
        res = ((T_avg_new - T_avg)**2 / (T_ds - fluid.T)**2)
        T_avg = T_avg_new
        
        ###Eliminate nonphysical solutions
        if (fluid.T < Tw_o and T_ds > Tw_o) or (fluid.T > Tw_o and T_ds < Tw_o):
            if j > 0:
                raise Exception('the pipe is too long')
            j += 1
            T_avg = (fluid.T + Tw_o) / 2
        
        ### Calculate heat flux 
        Q = (- h_T * dT ).to(ureg.W/ureg.m ** 2)

        return Tw_i, Tw_o, Q


def pipe_heat(pipe, fluid, m_dot):
    """Determine the heated status of the piping component
    """
    ### Calculate pressure drop and heat transfer coefficient
    dP, h_T, h_Q = dP_Pipe(m_dot, fluid, pipe)  
    
    ###heat flux defined on the wall of the pipe
    if hasattr(pipe, 'Q_def') and pipe.Q_def != None :
        try: 
            pipe.Q_def.m_as(ureg.W / ureg.m ** 2)
        except:
            raise ValueError(f"the Q_def is not properly defined in component {pipe}" )
        Tw_i, Tw_o = pipe_Q_def(fluid, pipe, m_dot, dP, h_Q)
        Q = pipe.Q_def
        
    ###Temperature defined on the wall of the pipe
    elif hasattr(pipe, 'Tw_def') and pipe.Tw_def != None :
        try: 
            pipe.Tw_def.m_as(ureg.K)
        except:
            raise ValueError(f"the Tw_def is not properly defined in component {pipe}" )
        Tw_i, Tw_o, Q = pipe_Tw_def(fluid, pipe, m_dot, dP, h_Q)

    ### Heat transfer or ambient/external temperature defined on the wall of the pipe
    elif (hasattr(pipe, 'h_ext') and pipe.h_ext is not None) or (hasattr(pipe, 'h_type') and pipe.h_type is not None):
        if hasattr(pipe, 'h_ext') and pipe.h_ext is not None:
            try:
                pipe.h_ext.m_as(ureg.W / ureg.m ** 2 / ureg.K)
            except:
                raise ValueError(f"the h_ext is not properly defined in component {pipe}")
    
            try:
                pipe.T_ext.m_as(ureg.K)
            except:
                pipe.T_ext = 293 * ureg.K
    
        if hasattr(pipe, 'h_type') and pipe.h_type is not None:
            try:
                pipe.safety_fact
            except:
                pipe.safety_fact = 1
    
        Tw_i, Tw_o, Q = pipe_h_ext(fluid, pipe, m_dot, dP, h_T)  

    ###Isolation on the external of the pipe 
    elif hasattr(pipe, 'isolation') and pipe.isolation != None :
        try: 
            pipe.isolation.k.m_as(ureg.W / ureg.m / ureg.K)
            pipe.isolation.OD.m_as(ureg.m)
        except:
            raise ValueError(f"the isolation (k, OD) is not properly defined in component {pipe}" )
        Tw_i, Tw_o, Q = pipe_insulated(fluid, pipe, m_dot, dP, h_T)
                
    ###Other
    else: 
        pipe.Q_def = Q = 0 * ureg.W/ ureg.m ** 2
        Tw_i, Tw_o = pipe_Q_def(fluid, pipe, m_dot, dP, h_Q)
    
    return Tw_i.to(ureg.K), Tw_o.to(ureg.K), dP.to(ureg.bar), Q.to(ureg.W/ ureg.m ** 2)
        

def k_pipe(pipe, T_wall, T_ext=293 * ureg.K):     ### you should add the possibility to define a table of values to do that (for materials not in the database)
    try:
        # Check if the 'pipe' object has the attribute 'k'
        if hasattr(pipe, 'k'):
            k = pipe.k.m_as(ureg.W / ureg.K / ureg.m) * ureg.W / ureg.K / ureg.m
        else:
            raise AttributeError("Pipe object has no attribute 'k'")
    except AttributeError:
        # Handle the case where 'k' is not an attribute of 'pipe'
        try:
            mat = pipe.material
        except AttributeError:
            # If 'pipe' has no 'material' attribute, default to a material
            if isinstance(pipe, CopperTube):
                mat = Material.OFHC
            else:
                mat = Material.SS304
       ### Determine the thermal conductivity using the 'nist_property' function
        if np.abs((T_wall - T_ext).m_as(ureg.K)) < 0.001:
            k = nist_property(mat, Property.TC, T1=T_wall)
        else:
            if T_wall < T_ext:
                k = nist_property(mat, Property.TC, T1=T_wall, T2=T_ext)
            else:
                k = nist_property(mat, Property.TC, T1=T_ext, T2=T_wall)
    return k.to(ureg.W / ureg.K / ureg.m)

def h_ext_(fluid, pipe, T_wall):  ### you should add the possibility to define a table of values to do that
    try:
        h = pipe.h_ext.m_as(ureg.W / ureg.K / ureg.m ** 2) * ureg.W / ureg.K / ureg.m ** 2
    except:
        try:
            type = pipe.h_type
        except:   
            raise Exception('You should define an external heat transfer type pipe.h_type') #or define type = 1
        if type >= 1: #Convection considered
            
            #Calculate Rayleigh and Prandtl numbers
            Ra_ = Ra(fluid, T_wall, pipe.OD)
            Pr_ = Pr(fluid)
            
            #Determine orientation of the pipe      
            try:
                orientation = pipe.orientation
            except:
                orientation = None
            #Calculatue Nusselt numbers  
            if orientation == 'vertical':        
                Nu_ = Nu_vcyl(Pr_, Ra_, pipe.OD, pipe.L)
            elif orientation == 'horizontal':
                Nu_ = Nu_hcyl(Pr_, Ra_)
            else:
                #Use maximum Nusselt if orientation is not defined
                Nu_ = max (Nu_vcyl(Pr_, Ra_, pipe.OD, pipe.L), Nu_hcyl(Pr_, Ra_))
            
            #Calculate heat transfer coefficient of external fluid
            h = heat_trans_coef(fluid, Nu_, pipe.OD)
        
        if type == 2: #Convection and radiation considered
            
            #Determine emissivity of material
            try:
                epsilon = pipe.epsilon
            except:
                #Assume material is polished steel if not defined
                epsilon = 0.075 
            
            #Calculate radiation heat transfer coefficient
            sigma = 5.670373e-8 * ureg.W/ureg.m ** 2 / ureg.K ** 4
            h_rad = epsilon * sigma * (fluid.T ** 4 - T_wall ** 4) / (fluid.T - T_wall) 
            
            #Calculate total heat transfer coefficient 
            h = h + h_rad

        if type == 3: ###including Rad ice and h_ice specific in the problem: to do 
            print('to do')
    return h.to(ureg.W / ureg.K / ureg.m ** 2)
