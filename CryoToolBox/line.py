# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:10:00 2023

@author: rbruce
"""

from . import piping
from . import heated_pipe
import numpy as np
import heat_transfer as ht
u = ht.ureg
from scipy.optimize import fsolve

class Line:
    '''Class defining a single line of cryo components.
    
    :Parameters:

    :param compList:    list of components (types Pipe, Valve, etc.)
    :param fluid:       fluid inlet
    :param mFlow_in:    inlet mass flow imposed in g/s if no upstream
    '''
    
    
    def __init__(self, compList):
        self.compList   = compList      # list of components
        self.P_out      = None          # outlet pressure fo the line
        self.P_out_comp = None          # outlet pressure for each components
        self.T_out      = None          # outlet temperature for the line
        self.T_out_comp = None          # outlet temperature for each components
        self.T_p_in_comp = None
        self.T_p_out_comp = None

            
    def dP_isot_line(self, fluid, mFlow_in):    #without temperature variation
            
        fluid_temp = fluid.copy()
        self.P_out_comp = np.zeros(np.size(self.compList)) * u.bar
        i = 0  
        D = 0
        for comp in self.compList:
            if i > 0 and comp.ID > D:
                enl = piping.Enlargement(D, comp.ID)
                dP = piping.dP_isot(mFlow_in, fluid_temp, enl)
                fluid_temp.update('P', fluid_temp.P - dP,'T', fluid_temp.T)
            if i > 0 and comp.ID < D:
                con = piping.Contraction(D, comp.ID)
                dP = piping.dP_isot(mFlow_in, fluid_temp, con)
                fluid_temp.update('P', fluid_temp.P - dP,'T', fluid_temp.T)
            if 'Xt' in comp.__dict__ and comp.Xt != 0:
                dP = piping.dP_control_valve(mFlow_in, fluid_temp, comp)
            else:
                dP = piping.dP_isot(mFlow_in, fluid_temp, comp)
            fluid_temp.update('P', fluid_temp.P - dP,'T', fluid_temp.T)
            self.P_out_comp[i] = fluid_temp.P
            i = i+1
            D = comp.ID
        self.P_out = fluid_temp.P
        
    # def dP_isot_wt_line(self, fluid, mFlow_in):    #with temperature variation (isenthalpic) no really usefull
            
    #     fluid_temp = fluid.copy()
    #     H = fluid_temp.Hmass  #enthalpy of the fluid
    #     self.P_out_comp = np.zeros(np.size(self.compList)) * u.bar
    #     self.T_out_comp = np.zeros(np.size(self.compList)) * u.K
    #     i = 0  
    #     D = 0
    #     for comp in self.compList:
    #         if i > 0 and comp.ID > D:
    #             enl = piping.Enlargement(D, comp.ID)
    #             dP = piping.dP_isot(mFlow_in, fluid_temp, enl)
    #             fluid_temp.update('P', fluid_temp.P - dP,'Hmass', H)
    #         if i > 0 and comp.ID < D:
    #             con = piping.Contraction(D, comp.ID)
    #             dP = piping.dP_isot(mFlow_in, fluid_temp, con)
    #             fluid_temp.update('P', fluid_temp.P - dP,'Hmass', H)
    #         if 'Xt' in comp.__dict__ and comp.Xt != 0:
    #             dP = piping.dP_control_valve(mFlow_in, fluid_temp, comp)
    #         else:
    #             dP = piping.dP_isot(mFlow_in, fluid_temp, comp)
    #         fluid_temp.update('P', fluid_temp.P - dP,'Hmass', H)
    #         self.P_out_comp[i] = fluid_temp.P
    #         self.T_out_comp[i] = fluid_temp.T
    #         i = i+1
    #         D = comp.ID
    #     self.P_out = fluid_temp.P
    #     self.T_out = fluid_temp.T
        
    def dP_adiab_line(self, fluid, mFlow_in):    #without temperature variation
            
        fluid_temp = fluid.copy()
        self.P_out_comp = np.zeros(np.size(self.compList)) * u.bar
        i = 0  
        D = 0
        for comp in self.compList:
            if i > 0 and comp.ID > D:
                enl = piping.Enlargement(D, comp.ID)
                dP = piping.dP_adiab(mFlow_in, fluid_temp, enl)
                fluid_temp.update('P', fluid_temp.P - dP,'T', fluid_temp.T)
            if i > 0 and comp.ID < D:
                con = piping.Contraction(D, comp.ID)
                dP = piping.dP_adiab(mFlow_in, fluid_temp, con)
                fluid_temp.update('P', fluid_temp.P - dP,'T', fluid_temp.T)
            if 'Xt' in comp.__dict__ and comp.Xt != 0:
                dP = piping.dP_control_valve(mFlow_in, fluid_temp, comp)
            else:
                dP = piping.dP_adiab(mFlow_in, fluid_temp, comp)
            fluid_temp.update('P', fluid_temp.P - dP,'T', fluid_temp.T)
            self.P_out_comp[i] = fluid_temp.P
            i = i+1
            D = comp.ID
        self.P_out = fluid_temp.P
        
    # def dP_adiab_wt_line(self, fluid, mFlow_in):    #with temperature variation (isenthalpic)
            
    #     fluid_temp = fluid.copy()
    #     H = fluid_temp.Hmass  #enthalpy of the fluid
    #     self.P_out_comp = np.zeros(np.size(self.compList)) * u.bar
    #     self.T_out_comp = np.zeros(np.size(self.compList)) * u.K
    #     i = 0  
    #     D = 0
    #     for comp in self.compList:
    #         if i > 0 and comp.ID > D:
    #             enl = piping.Enlargement(D, comp.ID)
    #             dP = piping.dP_adiab(mFlow_in, fluid_temp, enl)
    #             fluid_temp.update('P', fluid_temp.P - dP,'Hmass', H)
    #         if i > 0 and comp.ID < D:
    #             con = piping.Contraction(D, comp.ID)
    #             dP = piping.dP_adiab(mFlow_in, fluid_temp, con)
    #             fluid_temp.update('P', fluid_temp.P - dP,'Hmass', H)
    #         if 'Xt' in comp.__dict__ and comp.Xt != 0:
    #             dP = piping.dP_control_valve(mFlow_in, fluid_temp, comp)
    #         else:
    #             dP = piping.dP_adiab(mFlow_in, fluid_temp, comp)
    #         fluid_temp.update('P', fluid_temp.P - dP,'Hmass', H)
    #         self.P_out_comp[i] = fluid_temp.P
    #         self.T_out_comp[i] = fluid_temp.T
    #         i = i+1
    #         D = comp.ID
    #     self.P_out = fluid_temp.P
    #     self.T_out = fluid_temp.T
        
    def mf_isot_line(self, fluid, dP, m_guess):    #max mass flow without temperature variation
        def Delta_P(x):
            m_dot = x * u.g/u.s
            self.dP_isot_line(fluid, m_dot)
            return abs((dP-(fluid.P - self.P_out)).magnitude)
        
        return float(fsolve(Delta_P, m_guess)) * u.g/u.s

    def mf_isot_it_line(self, fluid, dP, step = 0.1):    #max mass flow without temperature variation with scpecifc iterations to help
        def Delta_P(x):
            m_dot = x * u.g/u.s
            self.dP_isot_line(fluid, m_dot)
            return abs((dP-(fluid.P - self.P_out)).magnitude)
        
        dP_i = 0 * u.bar
        m_guess = 0
        while dP_i < dP:
            m_guess = m_guess + step
            m_dot = m_guess * u.g/u.s
            self.dP_isot_line(fluid, m_dot)
            dP_i = (fluid.P - self.P_out)
            print(dP_i)
        
        return float(fsolve(Delta_P, m_guess)) * u.g/u.s
    
    # def mf_isot_wt_line(self, fluid, dP, m_guess):    #max mass flow without temperature variation
    #     def Delta_P(x):
    #         m_dot = x * u.g/u.s
    #         self.dP_isot_wt_line(fluid, m_dot)
    #         return abs((dP-(fluid.P - self.P_out)).magnitude)
        
    #     return float(fsolve(Delta_P, m_guess)) * u.g/u.s

    def mf_adiab_line(self, fluid, dP, m_guess):    #max mass flow without temperature variation
        def Delta_P(x):
            m_dot = x * u.g/u.s
            self.dP_adiab_line(fluid, m_dot)
            return abs((dP-(fluid.P - self.P_out)).magnitude)
        
        return float(fsolve(Delta_P, m_guess)) * u.g/u.s
    
    def mf_adiab_it_line(self, fluid, dP, step = 0.1):    #max mass flow without temperature variation with scpecifc iterations to help
        def Delta_P(x):
            m_dot = x * u.g/u.s
            self.dP_adiab_line(fluid, m_dot)
            return abs((dP-(fluid.P - self.P_out)).magnitude)
        
        dP_i = 0 * u.bar
        m_guess = 0
        while dP_i < dP:
            m_guess = m_guess + step
            m_dot = m_guess * u.g/u.s
            self.dP_adiab_line(fluid, m_dot)
            dP_i = (fluid.P - self.P_out)
            print(dP_i)
        
        return float(fsolve(Delta_P, m_guess)) * u.g/u.s

    # def mf_adiab_wt_line(self, fluid, dP, m_guess):    #max mass flow without temperature variation
    #     def Delta_P(x):
    #         m_dot = x * u.g/u.s
    #         self.dP_adiab_wt_line(fluid, m_dot)
    #         return abs((dP-(fluid.P - self.P_out)).magnitude)
        
    #     return float(fsolve(Delta_P, m_guess)) * u.g/u.s
    
    def dP_ht_line(self, fluid, mFlow_in):    
            
        fluid_temp = fluid.copy()
        #H = fluid_temp.Hmass  #enthalpy of the fluid
        self.P_out_comp = np.zeros(np.size(self.compList)) * u.bar
        self.T_out_comp = np.zeros(np.size(self.compList)) * u.K
        self.T_p_in_comp = np.zeros(np.size(self.compList)) * u.K
        self.T_p_out_comp = np.zeros(np.size(self.compList)) * u.K
        self.Q_out_comp = np.zeros(np.size(self.compList)) * u.dimensionless
        i = 0  
        D = 0
        for comp in self.compList:
            H = fluid_temp.Hmass  #enthalpy of the fluid
            if i > 0 and comp.ID > D:
                enl = piping.Enlargement(D, comp.ID)
                dP = piping.dP_adiab(mFlow_in, fluid_temp, enl)
                fluid_temp.update('P', fluid_temp.P - dP,'Hmass', H)
            if i > 0 and comp.ID < D:
                con = piping.Contraction(D, comp.ID)
                dP = piping.dP_adiab(mFlow_in, fluid_temp, con)
                fluid_temp.update('P', fluid_temp.P - dP,'Hmass', H)
            #dP = piping.dP_isot(mFlow_in, fluid_temp, comp)
            if 'T_p_out' in comp.__dict__ and comp.T_p_out.m_as(u.K) != 0:
                dP, dH = heated_pipe.dP_HT(mFlow_in, fluid_temp, comp)
                dH1 = dH.to(u.J/u.s)/mFlow_in
                self.T_p_in_comp[i] = comp.T_p_in
            elif 'Q' in comp.__dict__ and comp.Q.m_as(u.W/u.m ** 2) != 0:
                dP, dH = heated_pipe.dP_HT(mFlow_in, fluid_temp, comp)
                dH1 = dH.to(u.J/u.s)/mFlow_in
                self.T_p_in_comp[i] = comp.T_p_in
                self.T_p_out_comp[i] = comp.T_p_out_Q
            elif 'T_ambient' in comp.__dict__ and comp.T_ambient.m_as(u.K) != 0:
                if 'k_Foam' in comp.__dict__ and 'OD_Foam' in comp.__dict__:
                    dP, dH = heated_pipe.dP_HT_Foam(mFlow_in, fluid_temp, comp, comp.k_Foam, comp.OD_Foam, comp.T_ambient)
                else:
                    if 'safety_fact' in comp.__dict__:
                        safety_fact = comp.safety_fact
                    else:
                        safety_fact = 1
                    dP, dH = heated_pipe.dP_HT_ambient(mFlow_in, fluid_temp, comp, safety_fact, comp.T_ambient)
                dH1 = dH.to(u.J/u.s)/mFlow_in
                self.T_p_in_comp[i] = comp.T_p_in
                self.T_p_out_comp[i] = comp.T_p_out_Q
            else:
                if 'Xt' in comp.__dict__ and comp.Xt != 0:
                    dP = piping.dP_control_valve(mFlow_in, fluid_temp, comp)
                else:
                    dP = piping.dP_adiab(mFlow_in, fluid_temp, comp)
                dH1 = 0*u.J/u.kg
            fluid_temp.update('P', fluid_temp.P - dP,'Hmass', H + dH1.to(u.J/u.kg))
            self.P_out_comp[i] = fluid_temp.P
            self.T_out_comp[i] = fluid_temp.T

            if fluid_temp.Q > 0 and fluid_temp.Q < 1:
                self.Q_out_comp[i] = fluid_temp.Q
            i = i+1
            D = comp.ID
        self.P_out = fluid_temp.P
        self.T_out = fluid_temp.T
        if fluid_temp.Q > 0 and fluid_temp.Q < 1:
            self.Q_out = fluid_temp.Q
        
    def dP_incomp_line(self, fluid, mFlow_in):    #for liquid incompressible no temp variation
            
        fluid_temp = fluid.copy()
        self.P_out_comp = np.zeros(np.size(self.compList)) * u.bar
        i = 0  
        D = 0
        for comp in self.compList:
            if i > 0 and comp.ID > D:
                enl = piping.Enlargement(D, comp.ID)
                dP = piping.dP_incomp_1(mFlow_in, fluid_temp, enl)
                fluid_temp.update('P', fluid_temp.P - dP,'T', fluid_temp.T)
            if i > 0 and comp.ID < D:
                con = piping.Contraction(D, comp.ID)
                dP = piping.dP_incomp_1(mFlow_in, fluid_temp, con)
                fluid_temp.update('P', fluid_temp.P - dP,'T', fluid_temp.T)
            if 'Xt' in comp.__dict__ and comp.Xt != 0:
                dP = piping.dP_control_valve(mFlow_in, fluid_temp, comp)
            else:
                dP = piping.dP_incomp_1(mFlow_in, fluid_temp, comp)
            fluid_temp.update('P', fluid_temp.P - dP,'T', fluid_temp.T)
            self.P_out_comp[i] = fluid_temp.P
            i = i+1
            D = comp.ID
        self.P_out = fluid_temp.P