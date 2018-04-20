#python3
from math import pi, log10, sin, log
import logging
from .functions import *
from .NPS_data import NPS_table
Q_ = ureg.Quantity
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.WARNING)
logger = logging.getLogger(__name__)

class Pipe:
    '''
    General NPS pipe class. All basic methods implemented.
    '''
    def __init__ (self, D_nom, SCH=10, L=0*ureg.m):
        self.D = D_nom #Nominal diameter
        self.SCH = SCH
        self.L = L

    @property
    def OD(self):
            """
            Return OD of the Pipe element based on NPS table
            """
            try:
                return self._OD
            except AttributeError:
                self._OD = NPS_table[self.D]['OD']
                return self._OD

    @property
    def wall(self):
            """
            Return wall thickness of Pipe element based on NPS table
            """
            try:
                return self._wall
            except AttributeError:
                self._wall = NPS_table[self.D].get(self.SCH)
                return self._wall

    @property
    def ID(self):
            """
            Return ID of the Pipe element based on NPS table
            """
            try:
                return self._ID
            except AttributeError:
                self._ID = self.OD - 2*self.wall
                return self._ID

    @property
    def Area(self):
        """
        Calculate cross sectional area of pipe
        """
        try:
            return self._Area
        except AttributeError:
            self._Area = pi*self.ID**2/4
        return self._Area

    def f_T(self):
        '''
        Friction factor for complete turbulence for clean steel pipe.
        Fitted logarithmic function to data from A-25.
        '''
        if self.ID<0.2*ureg.inch or self.ID>48*ureg.inch:
            logger.warning('Tabulated friction data is given for ID = 0.2..48 inch, given {:.2~}'.format(self.ID))
        ln_ID = log(self.ID.to(ureg.inch).magnitude)
        return 0.0236-6.36e-3*ln_ID+8.12e-4*ln_ID**2 #Fitting by S. Koshelev

    @property
    def K(self):
        try:
            return self._K
        except AttributeError:
            self._K = self.f_T()*self.L/self.ID
            return self._K
    #TODO Implement more accurate method of friction factor estimation
class VJ_Pipe(Pipe):
    '''
    Vacuum jacketed pipe
    '''
    def __init__ (self, D_nom, SCH, L, VJ_D, VJ_SCH = 5):
        super().__init__(D_nom, SCH, L) 
        self.VJ = Pipe(VJ_D, VJ_SCH, L)

class Corrugated_Pipe(Pipe):
    '''
    Corrugated pipe class.
    '''
    def __init__ (self, D, L=0*ureg.m):
        super().__init__(D, None, L) 
        self.corrugated = True #deprecated

    @property
    def K(self):
        try:
            return self._K
        except AttributeError:
            self._K = 4*super().K #Multiplier 4 is used for corrugated pipe
        return self._K

    @property
    def OD(self):
        logger.warning('For corrugated piping assumed OD = D')
        return self.D
    @property
    def wall(self):
        logger.warning('For corrugated piping assumed wall = 0')
        return 0*ureg.m

class Openning (Pipe):
    def __init__ (self, ID):
        self._ID = ID

    @property
    def K(self):
        return 1 #For piping end

class Tube (Pipe):
    def __init__(self, OD, wall, L=0*ureg.m):
        self._OD = OD
        self.D = OD.to(ureg.inch).magnitude
        self._wall = wall
        self.L = L

class Elbow(Pipe):
    def __init__(self, D_nom, SCH=10, R_D=1.5, N=1, angle=90*ureg.deg):
        super().__init__(D_nom, SCH)
        self.R_D = R_D
        self.N = N
        self.angle = angle
        self.L = R_D*self.ID*angle

    @property
    def K(self):
        """
        Pressure drop in an elbow fitting. Based on Handbook of Hydraulic Resistance by I.E. Idelchik.
        """
        if self.angle <= 70*ureg.deg:
            A1 = 0.9*sin(self.angle)
        elif self.angle == 90*ureg.deg:
            A1 = 1
        elif self.angle >= 100*ureg.deg:
            A1 = 0.7+0.35*self.angle/(90*ureg.deg)
        else:
            logger.error('Improper bend angle for elbow. 90 degrees used instead: {}'.format(self.angle))
            A1 = 1

        if self.R_D < 1:
                B1 = 0.21*(self.R_D)**(-0.25)
        else:
                B1 = 0.21*(self.R_D)**(-0.5)

        C1 = 1 #use different value for non-axis symmetric
        return (A1*B1*C1+super().K)*self.N
 

class Piping (list):
    '''
    Piping system defined by nintial conditions and structure of pipe elements.
    '''
    def __init__ (self, Init_fdata, *Pipes):
        self.init_cond = Init_fdata 
        self.extend(Pipes)

    def add(self, *Pipes):
        self.extend(Pipes)

    def K(self):
        if len(self) > 0:
            K0 = 0*ureg.dimensionless
            A0 = self[0].Area
            for section in self:
                K0 += section.K*(A0/section.Area)**2
            return (K0, A0)
        else:
            logger.error('Piping has no elements! Use Piping.add to add sections to piping.')


    def dP(self, m_dot):
        '''
        Calculate pressure drop through piping. Lumped method using Darcy equation is used.
        '''
        (_, M, D_in) = rp_init(self.init_cond)
        P_0 = self.init_cond['P']
        rho = D_in*M
        K, Area = self.K()
        w = m_dot/(rho*Area)
        dP = dP_darcy (K, rho, w)
        P_out = P_0 - dP
        k = gamma(self.init_cond) #adiabatic coefficient
        rc = (2/(k+1))**(k/(k-1)) #Critical pressure drop; Note: according to Crane TP-410 should be depndent on the hydraulic resistance of the flow path
        if dP/P_0 <= 0.1:
            return dP
        elif dP/P_0 <= 0.4:
            (x, _, D_out) = rp_init(self.init_cond)
            T_0 = self.init_cond['T']
            D_out = flsh ("TP", T_0, P_out, x)['D']
            rho = (D_in+D_out)/2*M
            w = m_dot/(rho*Area)
            return dP_darcy (K, rho, w)
        elif 0.4<dP/P_0<(1-rc): #Subsonic flow
            logger.warning('Pressure drop too high for Darcy equation!')
            w = A*(rho/(K+2*log(P_0/P_out))*(P_0**2-P_out**2)/P_0)**0.5 #Complete isothermal equation, Crane TP-410, p. 1-8, eq. 1-6
            return dP_darcy (K, rho, w)
        else:
            logger.warning('Sonic flow developed. Consider reducing massflow: {:.3~}'.format(m_dot))

    def m_dot(self, P_out=0*ureg.psig):
        '''
        Calculate mass flow through the piping using initial conditions at the beginning of piping.
        Simple solution using Darcy equation is used.
        '''
        (x, M, D) = rp_init(self.init_cond)
        P_0 = self.init_cond['P']
        rho = D*M
        K, Area = self.K()
        k = gamma(self.init_cond) #adiabatic coefficient
        rc = (2/(k+1))**(k/(k-1)) #Critical pressure drop; Note: according to Crane TP-410 should be depndent on the hydraulic resistance of the flow path
        if P_out/P_0 > rc: #Subsonic flow
            delta_P = P_0-P_out
        else: #Sonic flow
            logger.warning('End pressure creates sonic flow. Max possible dP will be used')
            delta_P = P_0*(1-rc) #Crane TP-410, p 2-15
        return Area*(2*delta_P*rho/K)**0.5 #Net expansion factor for discharge is assumed to be 1 (conservative value)






#Hydraulic functions
def Re (M_dot = 0.01*ureg('kg/s'), Fluid_data = {'fluid':'air', 'P':101325*ureg.Pa, 'T':Q_(38, ureg.degC)}, Dim = 1.097*ureg.inch):
        """
        Reynolds number
        """
        fluid, T_fluid, P_fluid = unpack_fluid(Fluid_data)
        (x, M, D_fluid) = rp_init(Fluid_data)
        fluid_trans_prop = trnprp(T_fluid, D_fluid, x)
        mu_fluid = fluid_trans_prop['eta']*ureg('uPa*s') #dynamic viscosity

        d = Dim
        A = pi*d**2/4
        rho_fluid = D_fluid*M
        w_flow = M_dot/(rho_fluid*A)
        Re_number = w_flow*d*rho_fluid/mu_fluid
        return Re_number.to(ureg.dimensionless)
        #TODO Make Re() a simple function; move more complex function to pipe class or create a class containing also Fluid_data


#Pressure drops for different configurations of piping

def f_friction(M_dot, pipe, Fluid_data):
        """
        Calculate friction coefficient for pressure drop calculation. Based on Handbook of Hydraulic Resistance by I.E. Idelchik.
        """
        Re_num = Re(M_dot, Fluid_data, pipe.ID())
        mult = 1 #Default value for multiplier
        try:
            if pipe.corrugated:
                mult = 4 #Using 4x multiplicator compared to straight pipe
        except AttributeError:
            pass
        if Re_num < 2000:
            return 64/Re_num*mult
        elif Re_num > 4000:
            return 1/(1.8*log10(Re_num)-1.64)**2*mult
        else:
            return max(64/Re_num*mult, 1/(1.8*log10(Re_num)-1.64)**2*mult)
        
def dP_darcy (K, rho, w):
    '''
    Darcy equation for pressure drop.
    K - resistance coefficient
    rho - density of flow at entrance
    w - flow speed
    '''
    d_P = K*rho*w**2/2
    return d_P.to(ureg.psi)


#def dp_elbow (M_dot = 0.01*ureg('kg/s'), Elbow = {'R/D':1, 'D_nom':1, 'SCH':10, 'L':10*ureg.ft}, Fluid_data = {'fluid':'air', 'P':101325*ureg.Pa, 'T':Q_(38, ureg.degC)}):
#        """
#        Pressure drop in an elbow fitting. Based on Handbook of Hydraulic Resistance by I.E. Idelchik.
#        """
#                
#        delta = Elbow.get('Angle', 90*ureg.deg)
#        if delta <= 70*ureg.deg:
#                A1 = 0.9*sin(delta/rad)
#        elif delta == 90*ureg.deg:
#                A1 = 1
#        elif delta >= 100*ureg.deg:
#                A1 = 0.7+0.35*delta/(90*ureg.deg)
#        else:
#                raise BaseException ('Non standard angle is used: {:g}.'.format(delta))
#
#        R_frac_D = Elbow.get('R/D', 1) #Most of the elbows we use have R/D = 1.5 but it is not covered by this calculation. Using conservative value instead
#        if (R_frac_D >= 0.5) and (R_frac_D <= 1):
#                B1 = 0.21*(R_frac_D)**(-2.5)
#        elif R_frac_D >= 10:
#                B1 = 0.21*(R_frac_D)**(-0.5)
#        else:
#                raise BaseException ('Non standard elbow R/D is used: {:g}.'.format(R_frac_D))
#
#        a0_frac_b0 = Elbow.get('a0/b0', 1)
#        if a0_frac_b0 == 1:
#                C1 = 1
#        elif a0_frac_b0 <=4:
#                C1 = 0.85+0.125/a0_frac_b0
#        else:
#                C1 = 1.115-0.84/a0_frac_b0
#
#        zeta = A1*B1*C1
#
#        fluid = Fluid_data['fluid']
#        T_fluid = Fluid_data['T']
#        P_fluid = Fluid_data['P']
#
#        ID_pipe = ID(Elbow)
#        L_pipe = pi*R_frac_D*ID_pipe/2
#
#
#        (x, M, D_fluid) = rp_init(Fluid_data)
#        rho_fluid = D_fluid*M
#
#        A = Area(Elbow)
#        w_flow = M_dot/(rho_fluid*A)
#
#        delta_P_local = zeta*rho_fluid*w_flow**2/2
#
#        delta_P_frict = dP_pipe (M_dot, Elbow, Fluid_data)
#
#        return delta_P_frict+delta_P_local

