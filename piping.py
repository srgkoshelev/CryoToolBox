#python3
#from natu.units import *
#from natu.units import kureg.Pa, uureg.Pa, kJ
from math import pi, log10, sin, log
from pyrefprop import refprop as rp
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
    def __init__ (self, D_nom, SCH, L=0*ureg.m):
        self.D = D_nom #Nominal diameter
        self.SCH = SCH
        self.L = L
        self.corrugated = False

    def OD(self):
            """
            Return OD of the Pipe element based on NPS table
            """
            try:
                return self._OD_
            except AttributeError:
                self._OD_ = NPS_table[self.D]['OD']
                return self._OD_
#            elif piping_type == 'tube':
#                    OD_pipe = D_nom*ureg.inch
#                    Pipe.update({'OD':OD_pipe})
#                    return OD_pipe

    def wall(self):
            """
            Return wall thickness of Pipe element based on NPS table
            """
            try:
                return self._wall_
            except AttributeError:
                self._wall_ = NPS_table[self.D].get(self.SCH)
                return self._wall_

    def ID(self):
            """
            Return ID of the Pipe element based on NPS table
            """
            try:
                return self._ID_
            except AttributeError:
                self._ID_ = self.OD() - 2*self.wall()
                return self._ID_

    def Area(self):
        """
        Calculate cross sectional area of pipe
        """
        return pi*self.ID()**2/4

    def f_T(self):
        '''
        Friction factor for complete turbulence for clean steel pipe.
        Fitted logarithmic function to data from A-25.
        '''
        if self.ID()<0.2*ureg.inch or self.ID()>48*ureg.inch:
            input('WARNING: Tabulated friction data is given for ID = 1/2..24 inch, given {.2~}'.format(self.ID()))
        ln_ID = log(self.ID().to(ureg.inch).magnitude)
        return 0.0236-6.36e-3*ln_ID+8.12e-4*ln_ID**2 #Fitting by S. Koshelev

    def K(self):
        return self.f_T()*self.L/self.ID()
    #TODO Implement more accurate method of friction factor estimation






#def make_surface (Pipe, method = 'OD'):
#        """
#        Make surface element for convection heat load calculation.
#        Method determines which surface is considered. Orientation changes which dimension should be used for Nu  calculation. 
#        """
#        T = Pipe['fluid']['T']
#        if method == 'OD':
#                Diam = OD(Pipe)
#        elif method == 'VJ':
#                Diam = VJOD(Pipe)
#        elif method == 'average':
#                Diam = (OD(Pipe) + VJOD(Pipe))/2
#
#        if Pipe['Orientation'] == 'Horizontal':
#                Dim = Diam
#                Dim_sec = Pipe['L']
#        elif Pipe['Orientation'] == 'Vertical':
#                Dim = Pipe['L']
#                Dim_sec = Diam
#        return {'T':T, 'Dim':Dim, 'Dim_sec':Dim_sec}

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
    def __init__ (self, D_nom, SCH, L=0*ureg.m):
        super().__init__(D_nom, SCH, L) 
        self.corrugated = True

    def K(self):
        return 4*super().K() #Multiplier 4 is used for corrugated pipe

class Openning (Pipe):
    def __init__ (self, ID):
        self._ID_ = ID
        self.corrugated = False

    def K(self):
        return 1*ureg.dimensionless #For piping end

class Tube (Pipe):
    def __init__(self, OD, wall, L=0*ureg.m):
        self._OD_ = OD
        self.D = OD.to(ureg.inch).magnitude
        self._wall_ = wall
        self.L = L
        self.corrugated = False


class Piping (list):
    '''
    Piping system defined by nintial conditions and structure of pipe elements.
    '''
    def __init__ (self, Pipe, Fluid_data):
        self.append(Pipe)
        self[0].fdata = Fluid_data #Initial conditions

    def add(self, *Pipes):
        self.extend(Pipes)

    def K(self):
        K0 = 0*ureg.dimensionless
        A0 = self[0].Area()
        for section in self:
            K0 += section.K()*(A0/section.Area())**2
        return (K0, A0)

    def dP(self, m_dot):
        '''
        Calculate pressure drop through piping. Lumped method using Darcy equation is used.
        '''
        (_, M, D_in) = rp_init(self[0].fdata)
        P_0 = self[0].fdata['P']
        rho = D_in*M
        K, Area = self.K()
        w = m_dot/(rho*Area)
        dP = dP_darcy (K, rho, w)
        P_out = P_0 - dP
        k = gamma(self[0].fdata) #adiabatic coefficient
        rc = (2/(k+1))**(k/(k-1)) #Critical pressure drop; Note: according to Crane TP-410 should be depndent on the hydraulic resistance of the flow path
        if dP/P_0 <= 0.1:
            return dP
        elif dP/P_0 <= 0.4:
            (_, _, D_out) = rp_init(self[0].fdata)
            D_out = rp.flsh ("TP", T_fluid, P_fluid_out, x)['D']
            rho = (D_in+D_out)/2*M
            w = m_dot/(rho*Area)
            return dP_darcy (K, rho, w)
        elif 0.4<dP/P_0<(1-rc): #Subsonic flow
            logger.warning('Pressure drop too high for Darcy equation!')
            w = A*(rho/(K+2*log(P_0/P_out))*(P_0**2-P_out**2)/P_0)**0.5 #Complete isothermal equation, Crane TP-410, p. 1-8, eq. 1-6
        else:
            logger.warning('Sonic flow developed. Consider reducing massflow: {:.3~}'.format(m_dot))




    def m_dot(self, P_out=0*ureg.psig):
        '''
        Calculate mass flow through the piping using initial conditions at the beginning of piping.
        Simple solution using Darcy equation is used.
        '''
        (x, M, D) = rp_init(self[0].fdata)
        P_0 = self[0].fdata['P']
        rho = D*M
        K, Area = self.K()
        k = gamma(self[0].fdata) #adiabatic coefficient
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
        if pipe.corrugated:
                mult = 4 #Using 4x multiplicator compared to straight pipe
        else:
                mult = 1
        if Re_num < 2000:
                return 64/Re_num*mult
        elif Re_num > 4000:
                return 1/(1.8*log10(Re_num)-1.64)**2*mult
        else:
                # print ("Warning: Re = {:g}, transition flow. Maximum value between pure laminar or pure turbulent flow will be used".format(Reynolds))
                return max(64/Re_num*mult, 1/(1.8*log10(Re_num)-1.64)**2*mult)
        

def dP_pipe (M_dot, Pipe, Fluid_data={'fluid':'air', 'P':101325*ureg.Pa, 'T':Q_(38, ureg.degC)}):
        """
        Calculate pressure drop for a straight Pipe element. Works with hard pipe and corrugated hose (4x coefficient used). Calculation is based on Crane TP-410.
        """
        fluid, T_fluid, P_fluid = unpack_fluid(Fluid_data)
        x, M, D_fluid = rp_init(Fluid_data)
        rho_fluid = D_fluid*M

        w_flow = M_dot/(rho_fluid*Pipe.Area())
        f = f_friction(M_dot, Pipe, Fluid_data)
        K = f*Pipe.L/Pipe.ID()

        delta_P = dP_darcy(K, rho_fluid, w_flow)

        if delta_P/P_fluid < 0.1:
                return delta_P
        elif delta_P/P_fluid < 0.4: #If pressure drop above 10% of input pressure, recalculating for average density
                D_in = D_fluid
                P_fluid_out = P_fluid - delta_P
                D_out = rp.flsh ("TP", T_fluid, P_fluid_out, x)['D']
                rho_fluid = (D_in+D_out)/2*M
                w_flow = M_dot/(rho_fluid*Pipe.Area())
                delta_P = dP_darcy(K, rho_fluid, w_flow)
                return delta_P
        else:
                raise BaseException ('Pressure drop is {:.0%} which is greater than 40% recommended by Crane TP-410. Consider separating pipeline into sections.'.format(delta_P/P_fluid))

def dP_darcy (K, rho, w):
    '''
    Darcy equation for pressure drop.
    K - resistance coefficient
    rho - density of flow at entrance
    w - flow speed
    '''
    d_P = K*rho*w**2/2
    return d_P.to(ureg.psi)

def dP_openning (M_dot, openning, Fluid_data):
    '''
    Calculate pressure drop through the openning
    '''
    x, M, D_fluid = rp_init(Fluid_data)
    rho = D_fluid*M
    w = M_dot/(rho*openning.Area())
    K = openning.K
    return dP (K, rho, w)


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


#Simple solver for calculating flow
def derivative(f, x, h):
                return (f(x+h) - f(x-h)) / (2*h)  # might want to return a small non-zero if ==0


def solve(f, x0, h):
        lastX = x0
        nextX = lastX + 10* h  # different than lastX so loop starts OK
        i = 0
        i_max = 100
        while (abs(lastX - nextX) > h): 
            print (nextX)
            try:
                newY = f(nextX)                     
                print ("f(", nextX, ") = ", newY)     # print out progress... again just debug
                lastX = nextX
                nextX = lastX - newY / derivative(f, lastX, h)  # update estimate using N-R
                #i += 1
                #if i > i_max:
                #        raise BaseException ("To many cycles! Possible instability.")
            except:
                nextX *= 0.9
        return nextX

def flow_solve (f, x0, step):
        """
        Simpe iteration-based solver, backup if Newton solver fails.
        """
        X = x0
        while f(X) < -0.5*ureg.psi:
                X += step
                print (X, f(X))
        return X


def calculate_flow(piping, P_in, P_out=0*ureg.psig, M_dot_0 = 1*ureg('g/s'), M_step = 0.1*ureg('g/s')):
        def to_solve (M_dot):
                return dP_piping(M_dot, piping)-(P_in-P_out)
        return solve(to_solve, M_dot_0, M_step)

        #try:
        #except BaseException:
        #return flow_solve(to_solve, M_dot_0, M_step)


def dP_piping(M_dot, piping):
    """
    Calculate pressure drop for whole piping.
    Enthalpy is assumed constant (case of adiabatic flow with no work).
    """
    (x, M, D_fluid) = rp_init(piping[0].fdata)
    fluid, T0, P0 = unpack_fluid(piping[0].fdata)
    piping[0].fdata['h'] = flsh('TP',T0, P0, x)['h']
    for sec_num, Section in enumerate (piping):
        try:
            delta_P = dP_pipe(M_dot, Section, Section.fdata)
        except AttributeError:
            delta_P = dP_orifice(M_dot, Section, Section.fdata)
        P_next = Section.fdata['P'] - delta_P
#        if P_next < 0*ureg.psi:
#            raise BaseException ('Pressure drop {:.3} is greater than input pressure {:.3} for flow {:.3}! '.format(delta_P, Section.fdata['P'], M_dot))
        h_next = Section.fdata['h']
        T_next = flsh('PH',P_next, h_next, x)['t']
        try:
            piping[sec_num+1].fdata = {'fluid':fluid, 'T':T_next, 'P':P_next, 'h':h_next}
        except IndexError:
            Delta_P = piping[0].fdata['P'] - P_next
            return Delta_P.to(ureg.psi)
        #TODO return auto-elbow (if needed) calculation


#if __name__ == "__main__":
#        print("Testing piping module...")
##        print ("D_nom    OD    SCH    wall")
##        for D in [.125, .25, .375, .5, .75, 1, 1.125, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9]:
##                for sch in [5, 10, 30, 40, 80]:
##                        Test_pipe = Pipe(D, sch)
##                        print (D, Test_pipe.OD(), sch, Test_pipe.wall())
#
#        Test_pipe = Pipe(1, 10, 1*ureg.ft)
#        Test_pipe_corr = Corrugated_Pipe(1,10, 1*ureg.ft)
#        M_test = Q_(10, ureg.g/ureg.s)
#        Fluid = {'fluid':'air', 'P':101325*ureg.Pa, 'T':Q_(38, ureg.degC)}
#
#        print("Here is a dimensionless Reynolds number for 0.01 kg/s of standard air in 1 in SCH 10 pipe: {:g}.".format(Re()))
#        print("Pressure drop for 1 ft pipe would be {:.4~}, while for corrugated hose would 4 times bigger: {:.3~}".format (dP_pipe(M_test, Test_pipe, Fluid), dP_pipe(M_test, Test_pipe_corr, Fluid)))
#        #print ("And for a 90 deg elbow, pressure drop is {:g}.".format(dp_elbow()))
#
#        Test_piping = Piping(Test_pipe, Fluid)
#        print ("Flow for pressure drop for the straight pipe above: {:.3~}".format(calculate_flow(Test_piping, P_in=0.004828*ureg.psig)))
#        Test_orifice = Orifice(100*ureg.mm)
#        print("Pressure drop through 10 mm diameter orifice for the same flow will be {:.3~}".format (dP_orifice(M_test, Test_orifice, Fluid)))
#        Test_piping.add(Test_orifice)
#        print ('Pressure drop through pipe + orifice fo 10 g/s: {:.3~}'.format(dP_piping(M_test, Test_piping)))

