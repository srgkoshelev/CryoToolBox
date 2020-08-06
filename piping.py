"""Utilities for hydraulics calculations.

Contains functions for hydrodynamic calculations. The main source of the
equations is Crane TP-410.
"""
from math import pi, sin, log, sqrt, tan
from . import logger
from . import ureg, Q_
from .cp_wrapper import ThermState
from .functions import Air
from . import T_NTP, P_NTP
from . import os, __location__
from pint import set_application_registry
from serialize import load
from scipy.optimize import root_scalar

set_application_registry(ureg)


def _load_table(table_name):
    yaml_table = load((os.path.join(__location__, table_name)))
    result = {}
    for D, sub_table in yaml_table.items():
        result.update({D: {k: v*ureg.inch for k, v in sub_table.items()}})
    return result


NPS_table = _load_table('NPS_table.yaml')
COPPER_TABLE = _load_table('copper_table.yaml')


class Tube:
    """
    Tube, requires OD and wall thickness specified
    """
    def __init__(self, OD, wall=0*ureg.m, L=0*ureg.m, c=0*ureg.m):
        """Generate tube object.

        Parameters
        ----------
        OD : ureg.Quantity {length: 1}
            Outer diameter of the tube.
        wall : ureg.Quantity {length: 1}
            Wall thickness of the tube.
        L : ureg.Quantity {length: 1}
            Length of the tube.
        c : ureg.Quantity {length: 1}
            Sum of the mechanical allowances plus corrosion and erosion
            allowances.
        """
        self.OD = OD
        self.D = OD.to(ureg.inch).magnitude
        self.wall = wall
        self.ID = self.OD - 2*self.wall
        self.L = L
        self.area = self.calculate_area()
        self.volume = self.calculate_volume()
        self.K = self.calculate_K()
        # c = Q_('0.5 mm') for unspecified machined surfaces
        # TODO add calculation for c based on thread depth c = h of B1.20.1
        # Wall thickness under tolerance is 12.5% as per ASTM A999
        wall_tol = 0.125 * self.wall
        self.c = c + wall_tol
        self.type = 'Tube'

    def calculate_ID(self):
        """ureg.Quantity {length: 1} : Wall thickness of Pipe based on NPS table.
        """
        return self.OD - 2*self.wall

    def calculate_area(self):
        """ureg.Quantity {length: 2} : Cross sectional area of pipe.
        """
        return pi * self.ID**2 / 4

    def calculate_volume(self):
        """ureg.Quantity {length: 3} : Pipe inner volume.
        """
        return self.area * self.L

    def calculate_K(self):
        """ureg.Quantity {length: 1}: Resistance coefficient.
        """
        return self.f_T()*self.L/self.ID

    def f_T(self):
        return Pipe.f_T(self)

    def pressure_design_thick(self, P_int, P_ext=Q_('0 psig')):
        """Calculate pressure design thickness for given pressure and pipe material.

        Based on B31.3 304.1.

        Parameters
        ----------
        P_int : ureg.Quantity {length: -1, mass: 1, time: -2}
            Internal pressure, absolute
        P_ext : ureg.Quantity {length: -1, mass: 1, time: -2}
            External pressure, absolute

        Returns
        -------
        ureg.Quantity {length: 1}
            Minimum required wall thickness.
        """
        if self.check_material_defined():
            # Check whether S, E, W, and Y are defined
            pass
        D = self.OD
        # d = self.ID
        S, E, W, Y = self.S, self.E, self.W, self.Y
        P = P_int - P_ext  # Differential pressure across the wall
        t = P * D / (2*(S*E*W + P*Y))
        # TODO add 3b equation handling:
        # t = P * (d+2*c) / (2*(S*E*W-P*(1-Y)))
        if (t >= D/6) or (P/(S*E)) > 0.385:
            logger.error('Calculate design thickness in accordance \
            with B31.3 304.1.2 (b)')
            return None
        tm = t + self.c
        return tm

    def pressure_rating(self):
        """Calculate internal pressure rating.

        Based on B31.3 304.1.

        Returns
        -------
        ureg.Quantity {length: 1}
            Minimum required wall thickness.
        """
        if self.check_material_defined():
            # Check whether S, E, W, and Y are defined
            pass
        D = self.OD
        S, E, W, Y = self.S, self.E, self.W, self.Y
        t = self.wall - self.c
        P = 2 * t * S * E * W / (D-2*Y*t)
        return P

    def update(self, **kwargs):
        """ Add attributes, e.g. material stress or weld joint strength to the pipe.

        Parameters
        ----------
        **kwargs
            Parameters of the pipe, e.g. S=Q_('16700 psi'), Y=0.4.

        Returns
        -------
        None
        """
        # TODO Change to set material for predefined materials; check_material_defined will no longer be needed
        self.__dict__.update(kwargs)

    def check_material_defined(self):
        """Check whether following properties S, E, W, Y for stress
        calculations are defined.

        Parameters
        ----------
        * S: Stress value for pipe material
        * E: Quality factor from Table A-1A or A-1B
        * W: Weld joint strength reduction factor in accordance with 302.3.5 (e)
        * Y: coefficient from Table 304.1.1

        Returns
        -------
        None
        """
        try:
            self.S; self.E; self.W; self.Y
        except AttributeError:
            raise AttributeError('S, E, W, Y properties of the piping need to be \
            defined')

    def branch_reinforcement(self, BranchPipe, P, beta=Q_('90 deg'), d_1=None,
                             T_r=Q_('0 in')):
        """Calculate branch reinforcement status for given BranchPipe and
        reinforcing ring thickness, Tr.

        Parameters
        ----------
        BranchPipe : Pipe
            Branch pipe/tube instance with S, E, W, Y properties defined
        P : ureg.Quantity {length: -1, mass: 1, time: -2}
            Design pressure
        beta : ureg.Quantity {dimensionless}
            Smaller angle between axes of branch and run
        d_1 : ureg.Quantity {length: 1}
            Effective length removed from pipe at branch (opening for branch)
        T_r : ureg.Quantity {length: 1}
            Minimum thickness of reinforcing ring

        Returns
        -------
        None

        """
        # TODO Rename to reinforcement area?
        t_h = self.pressure_design_thick(P)
        if d_1 is None:
            d_1 = BranchPipe.OD
        D_h = self.OD
        T_h = self.wall
        T_b = BranchPipe.wall
        t_b = BranchPipe.pressure_design_thick(P)
        c = max(self.c, BranchPipe.c)  # Max allowance is used for calculation
        # B31.3 has no specification which allowance to use
        d_2 = half_width(d_1, T_b, T_h, c, D_h)
        A_1 = t_h * d_1 * (2-sin(beta))
        A_2 = (2*d_2-d_1) * (T_h - t_h - c)
        # height of reinforcement zone outside of run pipe
        L_4 = min(2.5*(T_h-c), 2.5*(T_b-c))
        A_3 = 2 * L_4 * (T_b-t_b-c)/sin(beta) * BranchPipe.S / self.S
        print(f'Required Reinforcement Area A_1: {A_1.to(ureg.inch**2):.3g~}')
        A_avail = A_2 + A_3  # Ignoring welding reinforcement
        print(f'Available Area A_2+A_3: {A_avail.to(ureg.inch**2):.3g~}')
        print(f'Weld branch connection is safe: {A_avail>A_1}')

    def info(self):
        return f'{self.type}, {self.OD:.3g~}x{self.wall:.3g~}, ' + \
            f'L={self.L:.3g~}'

    def __str__(self):
        return f'{self.D:.3g}" {self.type}'


class Pipe(Tube):
    """NPS pipe class.

    Pipe objects represent actual parts of the pipeline. The Pipe object
    will contain information such as OD, ID, wall thickness, and can be
    used to calculate flow coefficient K that is used for flow calculations.
    """
    def __init__(self, D_nom, SCH=40, L=0*ureg.m, c=Q_('0 mm')):
        """Generate `Pipe` object.

        Parameters
        ----------
        D_nom : float or ureg.Quantity {length: 1}
            Nominal diameter of piping; can be dimensionless or having a unit
            of length.
        SCH : int
            Pipe schedule. Default value is SCH 40 (STD).
        L : ureg.Quantity {length: 1}
            Pipe length
        c : ureg.Quantity {length: 1}
            Sum of the mechanical allowances plus corrosion and erosion
            allowances.
        """
        if isinstance(D_nom, ureg.Quantity):
            D = D_nom.magnitude
        else:
            D = D_nom
        OD = NPS_table[D]['OD']
        self.SCH = SCH
        # TODO Should I be using get here?
        wall = NPS_table[D].get(SCH)
        super().__init__(OD, wall, L, c)
        self.D = D
        self.type = 'NPS pipe'

        # """ureg.Quantity {length: 1} : Pipe OD based on NPS table.
        # """


        """ureg.Quantity {length: 1} : ID of the Pipe based on NPS table.
        """

    def f_T(self):
        """Calculate Darcy friction factor for complete turbulence for clean
        steel pipe.

        Fitted logarithmic function to data from A-25.

        Returns
        -------
        ureg.Quantity {dimensionless}
            Darcy friction factor.
        """
        if self.ID < 0.2*ureg.inch or self.ID > 48*ureg.inch:
            logger.debug('''Tabulated friction data is given for
                           ID = 0.2..48 inch, given {:.2~}'''.format(self.ID))
        ln_ID = log(self.ID.to(ureg.inch).magnitude)
        return 0.0236-6.36e-3*ln_ID+8.12e-4*ln_ID**2  # Fitting by S. Koshelev

    def info(self):
        return f'{self.type} {self.D}" SCH {self.SCH}, L={self.L:.3~g}'


class CopperTube(Tube):
    """Copper tube.

    Parameters
    ----------
    OD : ureg.Quantity {length: 1}
        Outer diameter of the tube.
    wall : ureg.Quantity {length: 1}
        Wall thickness of the tube.
    L : ureg.Quantity {length: 1}
        Length of the tube.
    """
    def __init__(self, D_nom, type_='Type K', L=0*ureg.m):
        if isinstance(D_nom, ureg.Quantity):
            D = D_nom.magnitude
        else:
            D = D_nom
        OD = COPPER_TABLE[D]['OD']
        wall = COPPER_TABLE[D][type_]
        c = 0 * ureg.inch  # Not affected by corrosion
        super().__init__(OD, wall, L, c)
        self.D = D
        self.type = 'Copper tube ' + type_


class VJPipe(Pipe):
    """Vacuum jacketed pipe.
    """
    def __init__(self, D_nom, *, SCH=5, L=0*ureg.m,
                 VJ_D, VJ_SCH=5, c=0*ureg.inch):
        """Generate Vacuum jacketed pipe object.

        Parameters
        ----------
        D_nom : float or ureg.Quantity {length: 1}
            Nominal diameter of the inner pipe.
        SCH : int
            Inner pipe schedule. Default value is SCH 40 (STD).
        L : ureg.Quantity {length: 1}
            Length of the inner pipe.
        VJ_D : float or ureg.Quantity {length: 1}
            Nominal diameter of the vacuum jacket.
        VJ_SCH : int
            Vacuum jacket pipe schedule. Default value is SCH 40 (STD).
        c : ureg.Quantity {length: 1}
            Sum of the mechanical allowances plus corrosion and erosion
            allowances of the inner pipe.
        """
        super().__init__(D_nom, SCH, L)
        self.VJ = Pipe(VJ_D, VJ_SCH, L)
        self.type = 'VJ pipe'

    def info(self):
        return f'NPS {self.D}" SCH {self.SCH} with VJ {self.VJ.D}", ' + \
            f'SCH {self.VJ.SCH}, L={self.L:.3~g}'


class CorrugatedPipe(Tube):
    """Corrugated pipe class.
    """
    def __init__(self, D_nom, L=0*ureg.m):
        """Generate corrugated pipe object.

        Parameters
        ----------
        D_nom : ureg.Quantity {length: 1}
            Nominal diameter of the inner pipe.
        L : ureg.Quantity {length: 1}
            Length of the pipe.
        """
        # TODO DRY
        OD = D_nom
        logger.debug('For corrugated piping assumed OD = D')
        wall = 0 * ureg.m
        c = 0 * ureg.inch
        super().__init__(OD, wall, L, c)
        self.K = 4*self.K  # Multiplier 4 is used for corrugated pipe
        self.type = 'Corrugated pipe'
        logger.debug('For corrugated piping assumed wall = 0')

    def branch_reinforcement(self):
        raise NotImplementedError('Branch reinforcement not implemented for'
                                  ' corrugated pipe')

    def pressure_design_thick(self):
        raise NotImplementedError('Pressure design thickness not implemented'
                                  ' for corrugated pipe')

    def info(self):
        return f'Corrugated pipe D={self.OD:.3g~}, L={self.L:.3g~}'


class Entrance ():
    """Pipe entrance, flush, sharp edged.
    """
    def __init__(self, ID):
        """Generate an pipe entrance.

        Parameters
        ----------
        ID : ureg.Quantity {length: 1}
            Inside diameter of the entrance.
        """
        self.ID = ID
        self.area = Tube.calculate_area(self)
        self.type = 'Entrance'
        self.K = 0.5  # Crane TP-410, A-29
        self.volume = 0 * ureg.ft**3

    def info(self):
        return f'{self.type}, {self.ID:.3g~}'

    def __str__(self):
        return self.type


class Exit (Entrance):
    """Pipe exit, projecting or sharp-edged, or rounded.
    """
    def __init__(self, ID):
        """Generate pipe exit.

        Parameters
        ----------
        ID : ureg.Quantity {length: 1}
            Inside diameter of the exit.
        """
        self.ID = ID
        self.area = Pipe.calculate_area(self)
        self.type = 'Exit'
        self.K = 1  # Crane TP-410, A-29
        self.volume = 0 * ureg.ft**3

    def info(self):
        return f'Exit opening, {self.ID:.3g~}'


class Orifice():
    """Square-edged orifice plate
    """
    def __init__(self, ID):
        """Generate orifice.

        Parameters
        ----------
        ID : ureg.Quantity {length: 1}
            Inside diameter of the orifice.
        """
        self.Cd = 0.61  # Thin sharp edged orifice plate
        self.ID = ID
        self.area = Pipe.calculate_area(self)
        self.type = 'Orifice'
        self.K = 1/self.Cd**2
        self.volume = 0 * ureg.ft**3

    def info(self):
        return f'Orifice, {self.ID:.3g~}'

    def __str__(self):
        return f'{self.ID:.3g~} orifice'


class ConicOrifice(Orifice):
    """Conic orifice
    """
    def __init__(self, D, ID):
        """Generate conic orifice.

        Parameters
        ----------
        D :
        ID : ureg.Quantity {length: 1}
            Inside diameter of the orifice.
        """
        super().__init__(ID)
        if NPS_table[D]['OD'] >= 1*ureg.inch:
            # For a smaller diameter using value for
            # square-edged plate (unfounded assumption)
            self.Cd = 0.73
            # Flow Measurements Engineering Handbook, Table 9.1, p. 9.16
        self.type = 'Conic orifice'

    def info(self):
        return f'Conic orifice, {self.ID:.3g~}'

    def __str__(self):
        return f'{self.ID:.3g~} conic orifice'


class Annulus():
    """"Annulus or tube in tube.

    Parameters
    ----------
    D1 : ureg.Quantity {length: 1}
        ID of the outer tube.
    D2 : ureg.Quantity {length: 1}
        OD of the inner tube.
    L : ureg.Quantity {length: 1}
        Length of the annulus tube.
    """
    def __init__(self, D1, D2, L=Q_('0m')):
        self.D1 = D1
        self.D2 = D2
        assert D1 > D2, 'D1 should be larger than D2'
        self.L = L
        assert D1 > D2
        self.area = pi / 4 * (D1**2 - D2**2)
        self.volume = Pipe.calculate_volume(self)
        self.ID = D1 - D2  # Hydraulic diameter
        self.f_T = lambda: Tube.f_T(self)
        self.K = Tube.calculate_K(self)

    def info(self):
        return f'Annulus D1={self.D1:.3g~}, D2={self.D2:.3g~}, L={self.L:.3g~}'

    def __str__(self):
        return f'{self.D1:.3g~} x {self.D2:.3g~} annulus'


class Elbow(Tube):
    """
    NPS Tee fitting.
    MRO makes method K from PipeElbow class to override method from Pipe class.
    """
    def __init__(self, OD, wall=0*ureg.inch, c=0*ureg.inch, R_D=1.5, N=1,
                 angle=90*ureg.deg):
        """Generate a tube elbow object.

        Parameters
        ----------
        OD : ureg.Quantity {length: 1}
            Outer diameter of the tube.
        wall : ureg.Quantity {length: 1}
            Wall thickness of the tube.
        R_D : ureg.Quantity {length: 1}
            Elbow radius/diameter ratio
        N : int
            Number of elbows in the pipeline
        angle : ureg.Quantity {dimensionless}
            Number of elbows in the pipeline
        """
        self.R_D = R_D
        self.N = N
        self.angle = angle
        super().__init__(OD, wall, L=0*ureg.m, c=c)
        self.L = R_D*self.ID*angle
        self.type = 'Tube elbow'

    def calculate_K(self):
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

        C1 = 1  # use different value for non-axis symmetric
        # Friction losses in the elbow
        K_frict = super().calculate_K()
        return (A1*B1*C1+K_frict)*self.N

    def info(self):
        return f'{self.N}x {self.type}, {self.OD}"x{self.wall}", ' + \
            f'{self.angle.to(ureg.deg)}, R_D = {self.R_D}'


class PipeElbow(Elbow, Pipe):
    """
    NPS Tee fitting.
    MRO makes method K from Elbow class to override method from Pipe class.
    """
    def __init__(self, D_nom, SCH=40, c=0*ureg.inch, R_D=1.5, N=1, angle=90*ureg.deg):
        """Generate a pipe elbow object.

        Parameters
        ----------
        D_nom : float or ureg.Quantity {length: 1}
            Nominal diameter of piping; can be dimensionless or having a unit
            of length.
        SCH : int
            Pipe schedule. Default value is SCH 40 (STD).
        R_D : ureg.Quantity {length: 1}
            Elbow radius/diameter ratio
        N : int
            Number of elbows in the pipeline
        angle : ureg.Quantity {dimensionless}
            Number of elbows in the pipeline
        """
        # D_nom and SCH go as positional arguments to Pipe __init__
        super().__init__(D_nom, SCH, c=c, R_D=R_D, N=N, angle=angle)
        self.type = 'NPS elbow'

    def info(self):
        return f'{self.N}x {self.type}, {self.D}" SCH {self.SCH}, ' + \
            f'{self.angle.to(ureg.deg)}, R_D = {self.R_D}'


class Tee(Tube):
    """
    Tee fitting based.
    """
    def __init__(self, OD, wall=0*ureg.inch, c=0*ureg.inch, direction='thru'):
        if direction in ['thru', 'through', 'run']:
            self.direction = 'run'
        elif direction in ['branch', 'side']:
            self.direction = 'branch'
        else:
            logger.error('''Tee direction is not recognized,
                         try "thru" or "branch": {}'''.format(direction))
        L = 0*ureg.m
        super().__init__(OD, wall, L, c)
        self.type = 'Tube tee'

    def calculate_K(self):
        if self.direction == 'run':
            return 20*self.f_T()  # Crane TP-410 p. A-29
        elif self.direction == 'branch':
            return 60*self.f_T()  # Crane TP-410 p. A-29

    def info(self):
        return f'{self.type}, {self.OD}x{self.wall}, {self.direction}'


class PipeTee(Tee, Pipe):
    """
    NPS Tee fitting.
    MRO makes method K from Tee class to override method from Pipe class.
    """
    def __init__(self, D_nom, SCH=40, c=0*ureg.inch, direction='thru'):
        # D_nom and SCH go as positional arguments to Pipe __init__
        super().__init__(D_nom, SCH, c, direction)
        self.type = 'NPS tee'

    def info(self):
        return f'{self.type}, {self.D}" SCH {self.SCH}, ' + \
            f'{self.direction}'


class Valve():
    """
    Generic valve with known Cv.
    """
    def __init__(self, D, Cv):
        # TODO DRY
        self.D = D
        self._Cv = Cv
        self.OD = None
        self.ID = self.D
        self.area = Pipe.calculate_area(self)
        self.L = None
        self.type = 'Valve'
        self.K = Cv_to_K(self._Cv, self.D)
        self.volume = 0 * ureg.ft**3

    def info(self):
        return f'{self.type}, {self.D}", Cv = {self._Cv:.3g}'

    def __str__(self):
        # TODO remove after init updated
        return f'{self.D:.3g~} valve'


# class GlobeValve(Pipe):
#     """
#     Globe valve.
#     """
#     def __init__(self, D):
#         super().__init__(D, None, None)
#         # ID for the valve is assumed equal to SCH40 ID:
#         self.ID = self.OD - 2*NPS_table[D].get(40)
#         self.type = 'Globe valve'
#         self._K = 340*self.f_T()  # Horizontal ball valve with beta = 1

#     @property
#     def volume(self):
#         return 0 * ureg.ft**3

#     def __str__(self):
#         return f'{self.type}, {self.D}"'


# class VCone(Pipe):
#     """
#     McCrometer V-Cone flowmeter.
#     """
#     def __init__(self, D, beta, Cf, SCH=40):
#         super().__init__(D, SCH, None)
#         self.beta = beta
#         self._Cf = Cf
#         self.type = 'V-cone flow meter'
#         # Equation is reverse-engineered from McCrometer V-Cone equations
#         self._K = 1/(self.beta**2/(1-self.beta**4)**0.5*self._Cf)**2


class Contraction():
    """
    Sudden and gradual contraction based on Crane TP-410.
    """
    def __init__(self, Pipe1, Pipe2, theta=ureg('180 deg')):
        """
        Pipe1: upstream pipe
        Pipe2: downstream pipe
        theta: contraction angle
        """
        self._Pipe1 = Pipe1
        self._Pipe2 = Pipe2
        ID1 = Pipe1.ID
        ID2 = Pipe2.ID
        self.beta = beta(ID1, ID2)
        self.theta = theta
        self.type = 'Contraction'
        self.L = None
        self.OD = None
        self.ID = min(ID1, ID2)
        self.area = Pipe.calculate_area(self)
        self.L = abs(ID1 - ID2) / tan(theta/2)
        self.volume = pi * self.L / 3 * (ID1**2 + ID1*ID2 + ID2**2)
        self.K = self.calculate_K()

    def calculate_K(self):
        if self.theta <= 45*ureg.deg:
            K_ = 0.8 * sin(self.theta/2) * (1-self.beta**2)
            # Crane TP-410 A-26, Formula 1 for K1 (smaller dia)
        elif self.theta <= 180*ureg.deg:
            K_ = 0.5 * (1-self.beta**2) * sqrt(sin(self.theta/2))
            # Crane TP-410 A-26, Formula 2 for K1 (smaller dia)
        else:
            logger.error(f'Theta cannot be greater than {180*ureg.deg} '
                         f'(sudden contraction): {self.theta}')
        return K_

    def info(self):
        return f'{self.type}, {self.theta.to(ureg.deg)} from {self._Pipe1} ' + \
            f'to {self._Pipe2}'

    def __str__(self):
        return f'{self.type}'


class Enlargement(Contraction):
    """
    Sudden and gradual enlargement based on Crane TP-410.
    """
    def __init__(self, Pipe1, Pipe2, theta=ureg('180 deg')):
        """
        Pipe1: upstream pipe
        Pipe2: downstream pipe
        theta: contraction angle
        """
        super().__init__(Pipe1, Pipe2, theta)
        self.type = 'Enlargement'

    def calculate_K(self):
        if self.theta <= 45*ureg.deg:
            K_ = 2.6 * sin(self.theta/2) * (1-self.beta**2)**2
            # Crane TP-410 A-26, Formula 3 for K1 (smaller dia)
        elif self.theta <= 180*ureg.deg:
            K_ = (1-self.beta**2)**2
            # Crane TP-410 A-26, Formula 4 for K1 (smaller dia)
        else:
            logger.error(f'Theta cannot be greater than {180*ureg.deg} \
            (sudden enlargement): {self.theta}')
        return K_


class Piping(list):
    '''
    Piping system defined by initial conditions and structure of
    pipe elements.
    Piping is modeled as a list of Pipe objects with given conditions at the
    beginning. Implemented methods allow to calculate pressure drop for given
    mass flow rate or mass flow rate for given pressure drop using lumped
    Darcy equation. All flow coefficients K are converted to the same base and
    added together to calculate single K value for the whole piping. This K
    value is used with Darcy equation to calculate pressure drop or mass flow.
    '''
    def __init__(self, fluid, pipes=[]):
        self.fluid = fluid
        self.extend(pipes)

    def add(self, *pipes):
        self.extend(pipes)

    def K(self):
        """
        Converting and adding flow coefficients K to same area.
        """
        K0 = 0*ureg.dimensionless
        try:
            A0 = self[0].area  # using area of the first element as base
        except IndexError:
            raise IndexError('Piping has no elements! '
                             'Use Piping.add to add sections to piping.')
        K0 = sum([section.K*(A0/section.area)**2 for section in self])
        return (K0, A0)

    def volume(self):
        result = []
        pipes = (Pipe, VJPipe, CorrugatedPipe, Tube, Annulus)
        # TODO remove elbows and tees after merge
        for pipe in self:
            if isinstance(pipe, pipes):
                result.append((str(pipe),
                               f'{pipe.L.to(ureg.ft).magnitude:.3g}',
                               f'{pipe.volume.to(ureg.ft**3).magnitude:.3g}'))
        return result

    def dP(self, m_dot):
        '''
        Calculate pressure drop through piping.
        Lumped method using Darcy equation is used.
        The pressure dropped is checked for choked condition.
        '''
        P_0 = self.fluid.P
        T_0 = self.fluid.T
        rho_0 = self.fluid.Dmass
        K, area = self.K()
        w = m_dot / (rho_0*area)
        dP = dP_darcy(K, rho_0, w)  # first iteration
        P_out = P_0 - dP
        k = self.fluid.gamma  # adiabatic coefficient
        # Critical pressure drop;
        # Note: according to Crane TP-410 should be dependent on
        # the hydraulic resistance of the flow path
        rc = (2/(k+1))**(k/(k-1))
        if self.fluid.Q < 0 or dP/P_0 <= 0.1:  # if q<0 then fluid is a liquid
            return dP
        elif dP/P_0 <= 0.4:
            TempState = ThermState(self.fluid.name, backend=self.fluid.backend)
            # Only working for pure fluids and pre-defined mixtures
            TempState.update('T', T_0, 'P', P_out)
            rho_out = TempState.Dmass
            rho_ave = (rho_0+rho_out) / 2
            w = m_dot/(rho_ave*area)
            return dP_darcy(K, rho_ave, w)
        elif 0.4 < dP/P_0 < (1-rc):  # Subsonic flow
            logger.warning('Pressure drop too high for Darcy equation!')
            # Complete isothermal equation, Crane TP-410, p. 1-8, eq. 1-6:
            w = (1/rho_0*(K+2*log(P_0/P_out))*(P_0**2-P_out**2)/P_0)**0.5
            return dP_darcy(K, rho_0, w)
        else:
            logger.warning('''Sonic flow developed. Calculated value ignores
            density changes. Consider reducing mass flow:
                           {:.3~}'''.format(m_dot))
            return dP_darcy(K, rho_0, w)

    def m_dot(self, P_out=0*ureg.psig):
        '''
        Calculate mass flow through the piping using initial conditions
        at the beginning of piping.
        Simple solution using Darcy equation is used.
        '''
        P_0 = self.fluid.P
        if P_0 <= P_out:
            logger.warning('Input pressure less or equal to output: \
            {P_0:.3g}, {P_out:.3g}')
            return Q_('0 g/s')
        rho = self.fluid.Dmass
        K, area = self.K()
        k = self.fluid.gamma  # adiabatic coefficient
        # Critical pressure drop
        # Note: according to Crane TP-410 should be dependent on
        # the hydraulic resistance of the flow path
        rc = (2/(k+1))**(k/(k-1))
        if P_out/P_0 > rc:  # Subsonic flow
            delta_P = P_0-P_out
        else:  # Sonic flow
            # logger.warning('''End pressure creates sonic flow.
            #                Max possible dP will be used''')
            delta_P = P_0*(1-rc)  # Crane TP-410, p 2-15
        # Net expansion factor for discharge is assumed to be 1
        # (conservative value):
        m_dot_ = area * (2*delta_P*rho/K)**0.5
        return m_dot_.to(ureg.g/ureg.s)

    def _solver_func(self, P_in_Pa, m_dot, P_out_act):
        """
        Solver function for calculating upstream pressure given flow and
        downstream pressure.

        :P_in_Pa: input pressure in Pa, float
        :args: calculation parameters:
            :m_dot: mass flow
            :P_out_act: actual downstream pressure
        """
        P_in = Q_(P_in_Pa, ureg.Pa)
        self.fluid.update('P', P_in, 'Smass', self.fluid.Smass)
        P_out_calc = P_in - self.dP(m_dot)
        P_out_calc_Pa = P_out_calc.to(ureg.Pa).magnitude
        P_out_act_Pa = P_out_act.to(ureg.Pa).magnitude
        return P_out_calc_Pa - P_out_act_Pa

    def P_in(self, m_dot, P_out=ureg('0 psig'), P_min=ureg('0 psig'),
             P_max=ureg('200 bar')):
        """
        Calculate upstream pressure given mass flow and downstream pressure.

        :m_dot: mass flow
        :P_out: downstream pressure
        :P_min: min expected pressure; should be lower than actual value
        :P_max: max expected pressure; should be higher than actual value
        """
        args = (m_dot, P_out)  # arguments passed to _solver_func
        # Convert pressure to dimensionless form
        P_min_Pa = P_min.to(ureg.Pa).magnitude
        P_max_Pa = P_max.to(ureg.Pa).magnitude
        bracket = [P_min_Pa, P_max_Pa]
        logger_level = logger.getEffectiveLevel()
        # ERROR and CRITICAL only will be shown; WARNING is suppressed
        logger.setLevel(40)
        solution = root_scalar(self._solver_func, args, bracket=bracket,
                               method='brentq')
        logger.setLevel(logger_level)
        P_in = Q_(solution.root, ureg.Pa)
        self.fluid.update('P', P_in, 'Smass', self.fluid.Smass)
        logger.debug(f'Comparing pressure drop:\n    dP method:\n        \
        {self.dP(m_dot).to(ureg.psi)}\n    \
        P_in - P_out:\n        \
        {(P_in-P_out).to(ureg.psi)}')
        logger.info(f'Calculated initial pressure: {P_in.to(ureg.psi):.1~f}')

# Supporting functions used for flow rate and pressure drop calculations.


def dP_darcy(K, rho, w):
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
    rho_w = 999*ureg('kg/m**3')  # Water density at 60 F
    # TODO move density to global const (or use ThermState)
    # Based on Crane TP-410 p. 2-10 and Darcy equation:
    Cv = A*(2/(K*rho_w))**0.5
    Cv.ito(ureg('gal/(min*(psi)**0.5)'))  # Convention accepted in the US
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
    if isinstance(Cv, (int, float)):  # Check if dimensionless input
        Cv = Cv*ureg('gal/(min*(psi)**0.5)')  # Convention accepted in the US
    A = pi*ID**2/4
    rho_w = 999*ureg('kg/m**3')  # Water density at 60 F
    # Based on Crane TP-410 p. 2-10 and Darcy equation:
    K = 2*A**2/(Cv**2*rho_w)
    return K.to(ureg.dimensionless)


def beta(d1, d2):
    """
    Calculate beta = d/D for contraction or enlargement.
    """
    return min(d1, d2)/max(d1, d2)


def to_standard_flow(flow_rate, fluid):
    '''
    Converting volumetric flow at certain conditions or mass flow to
    flow at NTP.
    '''
    fluid_NTP = ThermState(fluid.name)
    fluid_NTP.update('T', T_NTP, 'P', P_NTP)
    if flow_rate.dimensionality == ureg('kg/s').dimensionality:
        # mass flow, flow conditions are unnecessary
        q_std = flow_rate / fluid_NTP.Dmass
    elif flow_rate.dimensionality == ureg('m^3/s').dimensionality:
        # volumetric flow given, converting to standard pressure and
        # temperature
        if fluid.Dmass != -float('Inf')*ureg.kg/ureg.m**3:
            # By default ThermState is initialized with all fields == -inf
            q_std = flow_rate * fluid.Dmass / fluid_NTP.Dmass
        else:
            logger.warning('''Flow conditions for volumetric flow {:.3~}
                           are not set. Assuming standard flow at NTP.
                           '''.format(flow_rate))
            q_std = flow_rate
    else:
        logger.error('''Flow dimensionality is not supported: {:.3~}.
                       '''.format(flow_rate.dimensionality))
    q_std.ito(ureg.ft**3/ureg.min)
    return q_std


def equivalent_orifice(m_dot, dP, fluid=Air):
    """
    Calculate ID for the equivalent square edge orifice (Cd = 0.61) for given
    flow and pressure drop.
    """
    Cd = 0.61
    rho = fluid.Dmass
    ID = 2 * (m_dot/(pi*Cd*(2*dP*rho)**0.5))**0.5
    return ID.to(ureg.inch)


def to_mass_flow(Q_std, fluid=Air):
    """
    Calculate mass flow for given volumetric flow at standard conditions.
    """
    fluid_NTP = ThermState(fluid.name)
    fluid_NTP.update('T', T_NTP, 'P', P_NTP)
    m_dot = Q_std * fluid_NTP.Dmass
    return m_dot.to(ureg.g/ureg.s)


class ParallelPlateRelief:
    def __init__(self, Springs, Plate, Supply_pipe):
        """
        Initiate Parallel Plate instance.
        Parallel Plate relief valve designed by Fermilab

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
        # Before lifting pressure affects area within O-Ring OD
        A_lift = pi * self.Plate['OD_O_ring']**2 / 4
        F_precomp = self.Springs['N'] * self.Springs['k'] *\
                    self.Springs['dx_precomp']
        # Force required to lift the plate
        F_lift = F_precomp + self.Plate['W_plate']
        self.F_lift = F_lift
        P_set = F_lift / A_lift
        return P_set.to(ureg.psi)

    def P_open(self):
        """
        Calculate pressure required to fully open Parallel Plate relief
        """
        # compression required to provide vent area equal to supply pipe area
        dx_open = self.Supply_pipe.area / (pi*self.Plate['OD_plate'])
        # Force at fully open
        F_open = self.F_lift = self.Springs['N'] * self.Springs['k'] * dx_open
        # At fully open pressure is distributed as:
        # Full pressure for up to supply pipe diameter
        # Linear fall off up to plate OD
        A_open = self.Supply_pipe.area + pi/8*(self.Plate['OD_plate']**2 -
                                               self.Supply_pipe.ID**2)
        P_open = F_open / A_open
        return P_open.to(ureg.psi)


def half_width(d_1, T_b, T_h, c, D_h):
    """
    Calculate 'half width' of reinforcement zone per B31.3 304.3.3.
    """
    d_2_a = (T_b-c) + (T_h-c) + d_1/2
    return min(max(d_1, d_2_a), D_h)


def PRV_flow(ID, Kd, fluid):
    """
    Calculate mass flow through the relief valve based on
    BPVC VIII div. 1 UG-131 (e) (2).
    """
    A = pi * ID**2 / 4
    C = fluid.C_gas_const
    P = fluid.P
    W_T = C * A * P * fluid.MZT  # Theoretical flow
    W_a = W_T * Kd  # Actual flow
    return W_a.to(ureg.g/ureg.s)
