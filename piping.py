#' % Utilities for hydraulics calculations
#' % by Sergey Koshelev

#' This document describes the functions that are used for hydraulic analysis. The code is written in python and utilizes refprop for fluid properties. A separate wrapper for refprop functions allowing transparent usage of units is used. The main source of the equations is Crane TP-410.
#'
#' Importing math and thermodynamic functions, and NPS tabulated data. Setting up physical quantities for operations with units and logging.
from math import pi, log10, sin, log, sqrt
from . import logger
from . import ureg, Q_
from .cp_wrapper import *
from .functions import *
from . import T_NTP, P_NTP
from . import os, __location__
from pint import set_application_registry
import pickle
from scipy.optimize import root_scalar
from abc import ABC, abstractmethod

set_application_registry(ureg) #Should be used for both pickling and unpickling
NPS_table = pickle.load(open(os.path.join(__location__, "NPS.pkl"), "rb")) 

#' Pipe class is used to create Pipe objects, representing actual parts of the pipeline. The Pipe object will contain information such as OD, ID, wall thickness, and can be used to calculate flow coefficient K that is used for flow calculations.
class Pipe:
    """
    NPS pipe class.
    """
    def __init__ (self, D_nom, SCH=40, L=0*ureg.m, c=Q_('0 mm')):
        """
        Initiate instance of the Pipe class.

        :D_nom: nominal diameter of piping; can be dimensionless or having a unit
        :SCH: pip schedule
        :L: pipe length
        :c: sum of the mechanical allowances plus corrosion and erosion allowances
        :returns: None
        """
        try:
            self.D = D_nom.magnitude #Nominal diameter
        except AttributeError:
            self.D = D_nom #Nominal diameter
        self._OD = NPS_table[self.D]['OD']
        self.SCH = SCH
        self._wall = NPS_table[self.D].get(self.SCH)
        self._ID = self.OD - 2*self.wall
        self.L = L
        self._K = self.f_T()*self.L/self.ID
        self.c = c
        #c = Q_('0.5 mm') for unspecified machined surfaces
        # TODO add calculation for c based on thread depth c = h of B1.20.1
        self._type = 'NPS Pipe'

    @property
    def OD(self):
        """
        Return OD of the Pipe based on NPS table.

        :returns: pipe OD in unit of length
        """
        return self._OD

    @property
    def wall(self):
        """
        Return wall thickness of Pipe based on NPS table.

        :returns: pipe wall in unit of length
        """
        return self._wall

    @property
    def ID(self):
        """
        Return ID of the Pipe based on NPS table.

        :returns: pipe ID in unit of length
        """
        return self._ID

    @property
    def area(self):
        """
        Calculate cross sectional area of pipe.

        :returns: pipe cross section area
        """
        return pi * self.ID**2 / 4

    @property
    def volume(self):
        return self.area * self.L

    def f_T(self):
        """
        Calculate Darcy friction factor for complete turbulence for clean steel pipe.
        Fitted logarithmic function to data from A-25.

        :returns: Darcy friction factor
        """
        if self.ID<0.2*ureg.inch or self.ID>48*ureg.inch:
            logger.debug('''Tabulated friction data is given for 
                           ID = 0.2..48 inch, given {:.2~}'''.format(self.ID))
        ln_ID = log(self.ID.to(ureg.inch).magnitude)
        return 0.0236-6.36e-3*ln_ID+8.12e-4*ln_ID**2 #Fitting by S. Koshelev

    @property
    def K(self):
        """
        Calculate resistance coefficient for the Pipe element.

        :returns: resistance coefficient
        """
        return self._K

    def pressure_design_thick(self, P_int, P_ext=Q_('0 psig')):
        """
        Calculate pressure design thickness for given pressure and pipe material.
        Based on B31.3 304.1.
        :P: internal design pressure, gauge
        """
        if self.check_material_defined(): # Check whether S, E, W, and Y are defined
            pass
        D = self.OD
        d = self.ID
        S, E, W, Y = self.S, self.E, self.W, self.Y
        P = P_int - P_ext # Differential pressure across the wall
        t = P * D / (2*(S*E*W + P*Y))
        # TODO add 3b equation handling:
        # t = P * (d+2*c) / (2*(S*E*W-P*(1-Y)))
        if (t >= D/6) or (P/(S*E)) > 0.385:
            logger.error('Calculate design thickness in accordance \
            with B31.3 304.1.2 (b)')
            return None
        tm = t + self.c
        return tm

    def update(self, **kwargs):
        """
        Add attributes, e.g. material stress or weld joint strength to the pipe.
        """
        self.__dict__.update(kwargs)

    def check_material_defined(self):
        """Check whether instance has following material properties defined:
        :S: Stress value for pipe material
        :E: Quality factor from Table A-1A or A-1B
        :W: Weld joint strength reduction factor in accordance with 302.3.5 (e)
        :Y: coefficient from Table 304.1.1
        """
        try:
            self.S; self.E; self.W; self.Y
        except AttributeError:
            raise AttributeError('S, E, W, Y properties of the piping need to be \
            defined')

    def branch_reinforcement(self, BranchPipe, P, beta=Q_('90 deg'), d_1=None,
                             T_r=Q_('0 in')):
        """
        Calculate branch reinforcement status for given BranchPipe and reinforcing ring thickness, Tr.
        :BranchPipe: Branch pipe/tube instance with S, E, W, Y properties defined
        :P: Design pressure
        :beta: Smaller angle between axes of branch and run
        :d_1: Effective length removed from pipe at branch (opening for branch)
        :T_r: Minimum thickness of reinforcing ring
        """
        t_h = self.pressure_design_thick(P)
        if d_1 is None:
            d_1 = BranchPipe.OD
        D_h = self.OD
        T_h = self.wall
        T_b = BranchPipe.wall
        t_b = BranchPipe.pressure_design_thick(P)
        c = max(self.c, BranchPipe.c) # Max allowance is used for calculation
        # B31.3 has no specification which allowance to use
        d_2 = half_width(d_1, T_b, T_h, c, D_h)
        A_1 = t_h * d_1 * (2-sin(beta))
        A_2 = (2*d_2-d_1) * (T_h - t_h - c)
        L_4 = min(2.5*(T_h-c), 2.5*(T_b-c)) # height of reinforcement zone outside of run pipe
        A_3 = 2 * L_4 * (T_b-t_b-c)/sin(beta) * BranchPipe.S / self.S
        print(f'Required Reinforcement Area A_1: {A_1.to(ureg.inch**2):.3g~}')
        A_avail = A_2 + A_3 # Ignoring welding reinforcement
        print(f'Available Area A_3+A_3: {A_avail.to(ureg.inch**2):.3g~}')
        print(f'Weld branch connection is safe: {A_avail>A_1}')

    def __str__(self):
        return f'NPS {self.D}" SCH {self.SCH}'

#' Other Pipe elements are based on Pipe class and only specify the difference in regards to parent class.
class VJ_Pipe(Pipe):
    """
    Vacuum jacketed pipe.
    """
    def __init__ (self, D_nom, SCH, L, VJ_D, VJ_SCH=5):
        super().__init__(D_nom, SCH, L)
        self.VJ = Pipe(VJ_D, VJ_SCH, L)
        self._type = 'Vacuum jacketed pipe'

class Corrugated_Pipe(Pipe):
    '''
    Corrugated pipe class.
    '''
    def __init__ (self, D, L=0*ureg.m):
        super().__init__(D, None, L)
        self._K = 4*super().K #Multiplier 4 is used for corrugated pipe
        self._type = 'Corrugated pipe'

    @property
    def K(self):
        return self._K

    @property
    def OD(self):
        logger.debug('For corrugated piping assumed OD = D')
        return Q_(self.D*ureg.inch)

    @property
    def wall(self):
        logger.debug('For corrugated piping assumed wall = 0')
        return 0*ureg.m

class Entrance (Pipe):
    """
    Pipe entrance, flush, sharp edged.
    """
    def __init__ (self, ID):
        self._ID = ID
        self._type = 'Entrance'
        self._K = 0.5 #Crane TP-410, A-29

    @property
    def volume(self):
        return 0 * ureg.ft**3

class Exit (Pipe):
    """
    Pipe exit, projecting or sharp-edged, or rounded.
    """
    def __init__ (self, ID):
        self._ID = ID
        self._type = 'Exit'
        self._K = 1 #Crane TP-410, A-29

    @property
    def volume(self):
        return 0 * ureg.ft**3

class Orifice(Pipe):
    """
    Square-edged orifice plate
    """
    def __init__(self, ID):
        self.Cd = 0.61 #Thin sharp edged orifice plate
        self._ID = ID
        self._type = 'Orifice'

    @property
    def K(self):
        return 1/self.Cd**2

    @property
    def volume(self):
        return 0 * ureg.ft**3

class Conic_Orifice(Orifice):
    """
    Conic orifice
    """
    def __init__(self, D, ID):
        super().__init__(ID)
        if NPS_table[D]['OD'] >= 1*ureg.inch: 
            #For a smaller diameter using value for 
            #square-edged plate (unfounded assumption)
            self.Cd = 0.73 #Flow Measurements Engineering Handbook, Table 9.1, p. 9.16
        self._type = 'Conic orifice'

    @property
    def volume(self):
        return 0 * ureg.ft**3

class Tube(Pipe):
    """
    Tube, requires OD and wall thickness specified
    """
    def __init__(self, OD, wall=0*ureg.m, L=0*ureg.m, c=0*ureg.m):
        self._OD = OD
        self.D = OD.to(ureg.inch).magnitude
        self._wall = wall
        self._ID = self.OD - 2*self.wall
        self.L = L
        self._K = self.f_T()*self.L/self.ID
        self.c = c
        self._type = 'Tube'

class AbstractElbow(ABC):
    """
    Abstract Elbow class. __init__ method is abstract to avoid instantiation of this class.
    method K defines flow resistance calculation.
    R_D: elbow radius/diameter ratio
    N: number of elbows in the pipeline (to be used with lumped Darcy equation)
    """
    def __init__(self, R_D, N, angle):
        self.R_D = R_D
        self.N = N
        self.angle = angle
        self.L = R_D*self.ID*angle
        self._type = 'Elbow'

    @property
    def K(self):
        """
        Pressure drop in an elbow fitting. 
        Based on Handbook of Hydraulic Resistance by I.E. Idelchik.
        """
        if self.angle <= 70*ureg.deg:
            A1 = 0.9*sin(self.angle)
        elif self.angle == 90*ureg.deg:
            A1 = 1
        elif self.angle >= 100*ureg.deg:
            A1 = 0.7+0.35*self.angle/(90*ureg.deg)
        else:
            logger.error('''Improper bend angle for elbow. 
            90 degrees used instead: {}'''.format(self.angle))
            A1 = 1

        if self.R_D < 1:
                B1 = 0.21*(self.R_D)**(-0.25)
        else:
                B1 = 0.21*(self.R_D)**(-0.5)

        C1 = 1 #use different value for non-axis symmetric
        return (A1*B1*C1+super().K)*self.N

class PipeElbow(AbstractElbow, Pipe): #MRO makes method K from Elbow class to override method from Pipe class
    """
    NPS Tee fitting.
    """
    def __init__(self, D_nom, SCH=40, R_D=1.5, N=1, angle=90*ureg.deg):
        Pipe.__init__(self, D_nom, SCH)
        super().__init__(R_D, N, angle)

class TubeElbow(AbstractElbow, Tube): #MRO makes method K from Elbow class to override method from Pipe class
    """
    NPS Tee fitting.
    """
    def __init__(self, OD, wall, R_D=1.5, N=1, angle=90*ureg.deg):
        Tube.__init__(self, OD, wall)
        super().__init__(R_D, N, angle)

class AbstractTee(ABC):
    """
    Abstract Tee class. __init__ method is abstract to avoid instantiation of this class.
    method K defines flow resistance calculation.
    """
    @abstractmethod
    def __init__(self, direction):
        if direction in ['thru', 'through', 'run']:
            self.direction = 'run'
        elif direction in ['branch', 'side']:
            self.direction = 'branch'
        else:
            logger.error('''Tee direction is not recognized, 
                         try "thru" or "branch": {}'''.format(direction))
        self._type = 'Tee'

    @property
    def K(self):
        if self.direction == 'run':
            return 20*self.f_T() #Crane TP-410 p. A-29
        elif self.direction == 'branch':
            return 60*self.f_T() #Crane TP-410 p. A-29

class PipeTee(AbstractTee, Pipe): #MRO makes method K from Tee class to override method from Pipe class
    """
    NPS Tee fitting.
    """
    def __init__(self, D_nom, SCH=40, direction='thru'):
        Pipe.__init__(self, D_nom, SCH)
        super().__init__(direction)

class TubeTee(AbstractTee, Tube):
    """
    Tee fitting based.
    """
    def __init__(self, OD, wall, direction='thru'):
        Tube.__init__(self, OD, wall)
        super().__init__(direction)

class Valve(Pipe):
    """
    Generic valve with known Cv.
    """
    def __init__(self, ID, Cv):
        self._Cv = Cv
        self._ID = ID
        self._type = 'Valve'
        self._K = Cv_to_K(self._Cv, self.ID) 

    @property
    def volume(self):
        return 0 * ureg.ft**3

class Globe_valve(Pipe):
    """
    Globe valve.
    """
    def __init__(self, D):
        super().__init__(D, None, None)
        #ID for the valve is assumed equal to SCH40 ID:
        self._ID = self.OD - 2*NPS_table[D].get(40) 
        self._type = 'Globe valve'
        self._K = 340*self.f_T() #Horizontal ball valve with beta = 1

    @property
    def volume(self):
        return 0 * ureg.ft**3

class V_Cone(Pipe):
    """
    McCrometer V-Cone flowmeter.
    """
    def __init__(self, D, beta, Cf, SCH=40):
        super().__init__(D, SCH, None)
        self._beta = beta
        self._Cf = Cf
        self._type = 'V-cone flow meter'
        #Equation is reverse-engineered from McCrometer V-Cone equations
        self._K = 1/(self._beta**2/(1-self._beta**4)**0.5*self._Cf)**2 

class Contraction(Pipe):
    """
    Sudden and gradual contraction based on Crane TP-410.
    """
    def __init__(self, Pipe1, Pipe2, theta=ureg('180 deg')):
        """
        Pipe1: upstream pipe
        Pipe2: downstream pipe
        theta: contraction angle
        """
        ID1 = Pipe1.ID
        ID2 = Pipe2.ID
        self._beta = beta(ID1, ID2)
        self._theta = theta
        self._type = 'Contraction'
        self._ID = min(ID1, ID2)
        if theta <= 45*ureg.deg:
            self._K = 0.8 * sin(theta/2) * (1-self._beta**2) #Crane TP-410 A-26, Formula 1 for K1 (smaller dia)
        elif theta <= 180*ureg.deg:
            self._K = 0.5 * (1-self._beta**2) * sqrt(sin(theta/2)) #Crane TP-410 A-26, Formula 2 for K1 (smaller dia)
        else:
            logger.error(f'Theta cannot be greater than {180*ureg.deg} (sudden contraction): {theta}')

    @property
    def volume(self):
        return 0 * ureg.ft**3

class Enlargement(Pipe):
    """
    Sudden and gradual enlargement based on Crane TP-410.
    """
    def __init__(self, Pipe1, Pipe2, theta=ureg('180 deg')):
        """
        Pipe1: upstream pipe
        Pipe2: downstream pipe
        theta: contraction angle
        """
        ID1 = Pipe1.ID
        ID2 = Pipe2.ID
        self._beta = beta(ID1, ID2)
        self._theta = theta
        self._type = 'Enlargement'
        self._ID = min(ID1, ID2)
        if theta <= 45*ureg.deg:
            self._K = 2.6 * sin(theta/2) * (1-self._beta**2)**2 #Crane TP-410 A-26, Formula 3 for K1 (smaller dia)
        elif theta <= 180*ureg.deg:
            self._K = (1-self._beta**2)**2 #Crane TP-410 A-26, Formula 4 for K1 (smaller dia)
        else:
            logger.error(f'Theta cannot be greater than {180*ureg.deg} (sudden contraction): {theta}')

    @property
    def volume(self):
        return 0 * ureg.ft**3

#' Piping is modeled as a list of Pipe objects with given conditions at the beginning. Implemented methods allow to calculate pressure drop for given mass flow rate or mass flow rate for given pressure drop using lumped Darcy equation. All flow coefficients K are converted to the same base and added together to calculate single K value for the whole piping. This K value is used with Darcy equation to calculate pressure drop or mass flow. 
class Piping (list):
    '''
    Piping system defined by initial conditions and structure of 
    pipe elements.
    '''
    def __init__ (self, Fluid, Pipes=[]):
        self.Fluid = Fluid
        self.extend(Pipes)

    def add(self, *Pipes):
        self.extend(Pipes)

    def K(self):
        """
        Converting and adding flow coefficients K to same area.
        """
        K0 = 0*ureg.dimensionless
        try:
            A0 = self[0].area #using area of the first element as base
        except IndexError:
            raise IndexError('''Piping has no elements! 
                                Use Piping.add to add sections to piping.''')
        for section in self:
            K0 += section.K*(A0/section.area)**2
        return (K0, A0)

    @property
    def volume(self):
        self._volume = 0 * ureg.ft**3
        for pipe in self:
            self._volume += pipe.volume
        return self._volume

    def dP(self, m_dot):
        '''
        Calculate pressure drop through piping. 
        Lumped method using Darcy equation is used.
        The pressure dropped is checked for choked condition.
        '''
        P_0 = self.Fluid.P
        T_0 = self.Fluid.T
        rho_0 = self.Fluid.Dmass
        K, area = self.K()
        w = m_dot / (rho_0*area)
        dP = dP_darcy (K, rho_0, w) #first iteration
        P_out = P_0 - dP
        k = self.Fluid.gamma #adiabatic coefficient
        #Critical pressure drop;
        #Note: according to Crane TP-410 should be dependent on
        #the hydraulic resistance of the flow path
        rc = (2/(k+1))**(k/(k-1))
        if self.Fluid.Q < 0 or dP/P_0 <= 0.1: #if q<0 then fluid is a liquid
            return dP
        elif dP/P_0 <= 0.4:
            TempState = ThermState(Fluid.name, backend=Fluid.backend) #Only working for pure fluids and pre-defined mixtures
            TempState.update('T', T_0, 'P', P_out)
            rho_out = TempState.Dmass
            rho_ave = (rho_0+rho_out) / 2
            w = m_dot/(rho_ave*area)
            return dP_darcy (K, rho_ave, w)
        elif 0.4<dP/P_0<(1-rc): #Subsonic flow
            logger.warning('Pressure drop too high for Darcy equation!')
            #Complete isothermal equation, Crane TP-410, p. 1-8, eq. 1-6:
            w = (1/rho_0*(K+2*log(P_0/P_out))*(P_0**2-P_out**2)/P_0)**0.5
            return dP_darcy (K, rho_0, w)
        else:
            logger.warning('''Sonic flow developed. Calculated value ignores
                           density changes. Consider reducing mass flow:
                           {:.3~}'''.format(m_dot))
            return dP_darcy (K, rho_0, w)

    def m_dot(self, P_out=0*ureg.psig):
        '''
        Calculate mass flow through the piping using initial conditions 
        at the beginning of piping.
        Simple solution using Darcy equation is used.
        '''
        P_0 = self.Fluid.P
        if P_0 <= P_out:
            logger.warning('Input pressure less or equal to output: {P_0:.3g}, {P_out:.3g}')
            return Q_('0 g/s')
        rho = self.Fluid.Dmass
        K, area = self.K()
        k = self.Fluid.gamma #adiabatic coefficient
        #Critical pressure drop
        #Note: according to Crane TP-410 should be dependent on 
        #the hydraulic resistance of the flow path
        rc = (2/(k+1))**(k/(k-1)) 
        if P_out/P_0 > rc: #Subsonic flow
            delta_P = P_0-P_out
        else: #Sonic flow
            logger.warning('''End pressure creates sonic flow. 
                           Max possible dP will be used''')
            delta_P = P_0*(1-rc) #Crane TP-410, p 2-15
        #Net expansion factor for discharge is assumed to be 1 
        #(conservative value):
        m_dot_ = area * (2*delta_P*rho/K)**0.5 
        return m_dot_.to(ureg.g/ureg.s)

    def _solver_func(self, P_in_Pa, m_dot, P_out_act):
        """
        Solver function for calculating upstream pressure given flow and downstream pressure.

        :P_in_Pa: input pressure in Pa, float
        :args: calculation parameters:
            :m_dot: mass flow
            :P_out_act: actual downstream pressure
        """
        P_in = Q_(P_in_Pa, ureg.Pa)
        self.Fluid.update('P', P_in, 'Smass', self.Fluid.Smass)
        P_out_calc = P_in - self.dP(m_dot)
        P_out_calc_Pa = P_out_calc.to(ureg.Pa).magnitude
        P_out_act_Pa = P_out_act.to(ureg.Pa).magnitude
        return P_out_calc_Pa - P_out_act_Pa

    def P_in(self, m_dot, P_out=ureg('0 psig'), P_min=ureg('0 psig'), P_max=ureg('200 bar')):
        """
        Calculate upstream pressure given mass flow and downstream pressure.

        :m_dot: mass flow
        :P_out: downstream pressure
        :P_min: min expected pressure; should be lower than actual value
        :P_max: max expected pressure; should be higher than actual value
        """
        args = (m_dot, P_out) #arguments passed to _solver_func
        P_min_Pa = P_min.to(ureg.Pa).magnitude #Convert pressure to dimensionless form
        P_max_Pa = P_max.to(ureg.Pa).magnitude
        bracket = [P_min_Pa, P_max_Pa]
        logger_level = logger.getEffectiveLevel()
        logger.setLevel(40) #ERROR and CRITICAL only will be shown; WARNING is suppressed
        solution = root_scalar(self._solver_func, args, bracket=bracket, method='brentq')
        logger.setLevel(logger_level)
        P_in = Q_(solution.root, ureg.Pa)
        self.Fluid.update('P', P_in, 'Smass', self.Fluid.Smass)
        logger.debug(f'Comparing pressure drop:\n    dP method:\n        {self.dP(m_dot).to(ureg.psi)}\n    P_in - P_out:\n        {(P_in-P_out).to(ureg.psi)}')
        logger.info(f'Calculated initial pressure: {P_in.to(ureg.psi):.1~f}')

#' Supporting functions used for flow rate and pressure drop calculations.
def dP_darcy (K, rho, w):
    '''
    Darcy equation for pressure drop.
    K - resistance coefficient
    rho - density of flow at entrance
    w - flow speed
    '''
    d_P = K*rho*w**2/2
    return d_P.to(ureg.psi)

def K_to_Cv(K, ID):
    """
    Calculate flow coefficient Cv based on resistance coefficient value K.
    Based on definition:
    Cv = Q*sqrt(rho/(d_P*rho_w))
    where Q - volumetric flow, rho - flow density
    rho_w - water density at 60 F
    d_P - pressure drop through the valve.
    [Cv] = gal/(min*(psi)**0.5)
    """
    A = pi*ID**2/4
    rho_w = 999*ureg('kg/m**3') #Water density at 60 F
    Cv = A*(2/(K*rho_w))**0.5 #Based on Crane TP-410 p. 2-10 and Darcy equation
    Cv.ito(ureg('gal/(min*(psi)**0.5)')) #Convention accepted in the US
    return Cv

def Cv_to_K(Cv, ID):
    """
    Calculate resistance coefficient K based on flow coefficient value Cv.
    Based on definition:
    Cv = Q*sqrt(rho/(d_P*rho_w))
    where Q - volumetric flow, rho - flow density
    rho_w - water density at 60 F
    d_P - pressure drop through the valve.
    [Cv] = gal/(min*(psi)**0.5)
    """
    if isinstance(Cv, (int, float)): # Check if dimensionless input
        Cv = Cv*ureg('gal/(min*(psi)**0.5)') #Convention accepted in the US
    A = pi*ID**2/4
    rho_w = 999*ureg('kg/m**3') #Water density at 60 F
    K = 2*A**2/(Cv**2*rho_w) #Based on Crane TP-410 p. 2-10 and Darcy equation
    return K.to(ureg.dimensionless)

def beta(d1, d2):
    """
    Calculate beta = d/D for contraction or enlargement.
    """
    return min(d1, d2)/max(d1, d2)

def to_standard_flow(flow_rate, Fluid):
    '''
    Converting volumetric flow at certain conditions or mass flow to 
    flow at NTP.
    '''
    Fluid_NTP = ThermState(Fluid.name)
    Fluid_NTP.update('T', T_NTP, 'P', P_NTP)
    if flow_rate.dimensionality == ureg('kg/s').dimensionality: 
        #mass flow, flow conditions are unnecessary
        q_std = flow_rate / Fluid_NTP.Dmass
    elif flow_rate.dimensionality == ureg('m^3/s').dimensionality: 
        #volumetric flow given, converting to standard pressure and temperature
        if Fluid.Dmass != -float('Inf'): #By default ThermState is initialized with all fields == -inf
            q_std = flow_rate * Fluid.Dmass / Fluid_NTP.Dmass
        else:
            logger.warning('''Flow conditions for volumetric flow {:.3~} 
                           are not set. Assuming standard flow at NTP.
                           '''.format(flow_rate))
            q_std = flow_rate
    else:
        logger.warning('''Flow dimensionality is not supported: {:.3~}.
                       '''.format(flow_rate.dimensionality))
    q_std.ito(ureg.ft**3/ureg.min)
    return q_std

def equivalent_orifice(m_dot, dP, Fluid=Air):
    """
    Calculate ID for the equivalent square edge orifice (Cd = 0.61) for given flow and pressure drop.
    """
    Cd = 0.61
    rho = Fluid.Dmass
    ID = 2 * (m_dot/(pi*Cd*(2*dP*rho)**0.5))**0.5
    return ID.to(ureg.inch)

def to_mass_flow(Q_std, Fluid=Air):
    """
    Calculate mass flow for given volumetric flow at standard conditions.
    """
    Fluid_NTP = ThermState(Fluid.name)
    Fluid_NTP.update('T', T_NTP, 'P', P_NTP)
    m_dot = Q_std * Fluid_NTP.Dmass
    return m_dot.to(ureg.g/ureg.s)

#Parallel Plate relief valve designed by Fermilab
class ParallelPlateRelief:
    def __init__(self, Springs, Plate, Supply_pipe):
        """
        Initiate Parallel Plate instance.

        Springs: dictionary containing 'N' - number of springs, 'k' - spring rate, and 'dx_precomp' - pre-compression length
        Plate: dictionary containing 'OD_plate' - OD of the plate, 'OD_O_ring' - OD of the O-Ring installed, and 'W_plate' - plate weight
        Supply_pipe: Pipe/Tube object of upstream pipe
        """
        self.Springs = Springs
        self.Plate = Plate
        self.Supply_pipe = Supply_pipe
    def P_set(self):
        """"
        Calculate set (lift) pressure of the Parallel Plate relief
        """
        A_lift = pi * self.Plate['OD_O_ring']**2 / 4 #Before lifting pressure affects area within O-Ring OD
        F_precomp = self.Springs['N'] * self.Springs['k'] * self.Springs['dx_precomp']
        F_lift = F_precomp + self.Plate['W_plate'] #Force required to lift the plate
        self.F_lift = F_lift
        P_set = F_lift / A_lift
        return P_set.to(ureg.psi)
    def P_open(self):
        """
        Calculate pressure required to fully open Parallel Plate relief
        """
        dx_open = self.Supply_pipe.area / (pi*self.Plate['OD_plate']) #compression required to provide vent area equal to supply pipe area
        F_open = self.F_lift = self.Springs['N']*self.Springs['k']*dx_open #Force at fully open
        #At fully open pressure is distributed as:
        #Full pressure for up to supply pipe diameter
        #Linear fall off up to plate OD
        A_open = self.Supply_pipe.area + pi/8*(self.Plate['OD_plate']**2 - self.Supply_pipe.ID**2)
        P_open = F_open / A_open
        return P_open.to(ureg.psi)

def half_width(d_1, T_b, T_h, c, D_h):
    """
    Calculate 'half width' of reinforcement zone per B31.3 304.3.3.
    """
    d_2_a = (T_b-c) + (T_h-c) + d_1/2
    return min(max(d_1, d_2_a), D_h)

def PRV_flow(ID, Kd, Fluid):
    """
    Calculate mass flow through the relief valve based on BPVC VIII div. 1 UG-131 (e) (2).
    """
    A = pi * ID**2 /4
    C = Fluid.C_gas_constant
    P = Fluid.P
    W_T = C * A * P * Fluid.MZT # Theoretical flow
    W_a = W_T * Kd # Actual flow
    return W_a.to(ureg.g/ureg.s)
