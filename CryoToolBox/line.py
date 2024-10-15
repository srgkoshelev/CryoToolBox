# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:10:00 2023

@author: rbruce
"""

from . import piping
from . import heated_pipe
import numpy as np
from .std_conditions import ureg
#from scipy.optimize import fsolve
from math import pi

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
        self.Tw_i_comp = None
        self.Tw_o_comp = None

            
    def dP_isot_line(self, fluid, mFlow_in):    #without temperature variation
            
        fluid_temp = fluid.copy()
        self.P_out_comp = np.zeros(np.size(self.compList)) * ureg.bar
        i = 0  
        D = 0
        for comp in self.compList:
            #if the next component diameter is bigger
            if i > 0 and comp.ID > D:
                enl = piping.Enlargement(D, comp.ID)
                dP = piping.dP_isot(mFlow_in, fluid_temp, enl)
                fluid_temp.update('P', fluid_temp.P - dP,'T', fluid_temp.T)
            #if the next component diameter is smaller
            if i > 0 and comp.ID < D:
                con = piping.Contraction(D, comp.ID)
                dP = piping.dP_isot(mFlow_in, fluid_temp, con)
                fluid_temp.update('P', fluid_temp.P - dP,'T', fluid_temp.T)
            #if the component is a control valve 
            if 'Xt' in comp.__dict__ and comp.Xt != 0:
                dP = piping.dP_control_valve(mFlow_in, fluid_temp, comp)
            else:
                dP = piping.dP_isot(mFlow_in, fluid_temp, comp)
            fluid_temp.update('P', fluid_temp.P - dP,'T', fluid_temp.T)
            self.P_out_comp[i] = fluid_temp.P
            i = i+1
            D = comp.ID
        self.P_out = fluid_temp.P
    
        
    def dP_incomp_line(self, fluid, mFlow_in):    #for liquid incompressible no temp variation
            
        fluid_temp = fluid.copy()
        self.P_out_comp = np.zeros(np.size(self.compList)) * ureg.bar
        i = 0  
        D = 0
        for comp in self.compList:
            #if the next component diameter is bigger
            if i > 0 and comp.ID > D:
                enl = piping.Enlargement(D, comp.ID)
                dP = piping.dP_incomp(mFlow_in, fluid_temp, enl)
                fluid_temp.update('P', fluid_temp.P - dP,'T', fluid_temp.T)
            #if the next component diameter is smaller
            if i > 0 and comp.ID < D:
                con = piping.Contraction(D, comp.ID)
                dP = piping.dP_incomp(mFlow_in, fluid_temp, con)
                fluid_temp.update('P', fluid_temp.P - dP,'T', fluid_temp.T)
            #if the component is a control valve 
            if 'Xt' in comp.__dict__ and comp.Xt != 0:
                dP = piping.dP_control_valve(mFlow_in, fluid_temp, comp)
            else:
                dP = piping.dP_incomp(mFlow_in, fluid_temp, comp)
            fluid_temp.update('P', fluid_temp.P - dP,'T', fluid_temp.T)
            self.P_out_comp[i] = fluid_temp.P
            i = i+1
            D = comp.ID
        self.P_out = fluid_temp.P
        
    def dP_ht_line(self, fluid, mFlow_in):    #with temp variation
            
        fluid_temp = fluid.copy()
        self.P_out_comp = np.zeros(np.size(self.compList)) * ureg.bar
        self.T_out_comp = np.zeros(np.size(self.compList)) * ureg.K
        self.Tw_i_comp = np.zeros(np.size(self.compList)) * ureg.K
        self.Tw_o_comp = np.zeros(np.size(self.compList)) * ureg.K
        self.Q_out_comp = np.zeros(np.size(self.compList)) * ureg.W/ureg.m**2
        i = 0  
        D = 0
        for comp in self.compList:
            H = fluid_temp.Hmass  #enthalpy of the fluid
            if i > 0 and comp.ID > D:
                enl = piping.Enlargement(D, comp.ID)
                dP = piping.dP_isot(mFlow_in, fluid_temp, enl)
                fluid_temp.update('P', fluid_temp.P - dP,'Hmass', H)
            if i > 0 and comp.ID < D:
                con = piping.Contraction(D, comp.ID)
                dP = piping.dP_isot(mFlow_in, fluid_temp, con)
                fluid_temp.update('P', fluid_temp.P - dP,'Hmass', H)
            if 'Xt' in comp.__dict__ and comp.Xt != 0:
                dP = piping.dP_control_valve(mFlow_in, fluid_temp, comp)
                dH = 0*ureg.J/ureg.kg
            else:
                try:
                    Tw_i, Tw_o, dP, Q = heated_pipe.pipe_heat(comp, fluid_temp, mFlow_in) #change to dP, and dH
                    self.Tw_i_comp[i] = Tw_i.m_as(ureg.K) * ureg.K
                    self.Tw_o_comp[i] = Tw_o.m_as(ureg.K) * ureg.K
                    self.Q_out_comp[i] = Q.m_as(ureg.W/ureg.m**2) * ureg.W/ureg.m**2       #to modify
                    dH = Q * comp.ID * comp.L * pi / mFlow_in
                except:
                    dP = piping.dP_isot(mFlow_in, fluid_temp, comp)
                    dH = 0*ureg.J/ureg.kg
            fluid_temp.update('P', fluid_temp.P - dP,'Hmass', H + dH.to(ureg.J/ureg.kg))
            self.P_out_comp[i] = fluid_temp.P
            self.T_out_comp[i] = fluid_temp.T
            
            i = i+1
            D = comp.ID
        self.P_out = fluid_temp.P
        self.T_out = fluid_temp.T