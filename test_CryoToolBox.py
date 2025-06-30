import CryoToolBox as ctb
from CryoToolBox import odh
import pprint
import unittest
import doctest
from math import pi, sin, sqrt, tan
from unittest.mock import MagicMock

pp = pprint.PrettyPrinter()

ht = ctb
u = ht.ureg
Q_ = ht.Q_

def load_tests(loader, tests, ignore):
    """
    Load unit tests and doctests with ellipsis support.

    This function is automatically called by unittest when discovering tests.
    It adds doctests from the ctb.piping module with ellipsis directive enabled.
    """
    # Add doctests with ellipsis and other useful flags
    tests.addTests(doctest.DocTestSuite(
        ctb.piping,
        optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE | doctest.IGNORE_EXCEPTION_DETAIL
    ))
    return tests


# Checking if REFPROP works to run the tests
try:
    ctb.ThermState('helium', P=1*u.bar, T=ctb.T_NTP, backend='REFPROP')
    skip_RP_test = False
except ValueError:
    skip_RP_test = True

# Checking if Heprop works to run the tests
try:
    ctb.ThermState('helium', T= 300*u.K, P=1*u.bar, backend='HEPROP')
    skip_HP_test = False
except ValueError:
    skip_HP_test = True

@unittest.skipIf(skip_RP_test, 'REFPROP failed to initialize, skipping REFPROP tests.')
class RefPropTest(unittest.TestCase):
    """Verify fluid properties are calculated accurately.

    Tests borrowed from pyrefprop package.
    """

    criteria = 0.00001  # NIST acceptance criteria

    def assertNIST(self, nist, calc, criteria=None):
        criteria = criteria or self.criteria
        assert abs(nist-calc) / nist < criteria, \
            f'Calculated value {calc} is not within {criteria:%} of NIST value ' + \
            f'{nist}'

    # def test_air(self):
    #     fluid = ht.ThermState('air', backend='REFPROP')
    #     self.assertNIST(28.958600656, fluid.M)
    # TODO Check why different values

    def test_argon(self):
        fluid = ht.ThermState('argon', backend='REFPROP')
        P = 2 * 1000 * u.kPa
        D = 15 / fluid.M * u.mol/u.L
        fluid.update_kw(P=P, Dmolar=D)
        self.assertNIST(637.377588657857, fluid.T.to(u.K).magnitude)

    def test_r134a(self):
        fluid = ht.ThermState('r134a', backend='REFPROP')
        T = 400 * u.K
        D = 50 / fluid.M * u.mol/u.L
        fluid.update_kw(T=T, Dmolar=D)
        self.assertNIST(1.45691892789737, fluid.P.to(u.kPa).magnitude/1000)

    def test_oxygen(self):
        fluid = ht.ThermState('oxygen', backend='REFPROP')
        T = 100 * u.K
        P = 1000 * u.kPa
        fluid.update_kw(T=T, P=P)
        self.assertNIST(153.886680663753,
                        fluid.viscosity.to(u.uPa*u.s).magnitude)

    def test_nitrogen(self):
        fluid = ht.ThermState('nitrogen', backend='REFPROP')
        T = 100 * u.K
        fluid.update_kw(T=T, Q=0)
        self.assertNIST(100.111748964945,
                        fluid.conductivity.to(u.W/(u.m*u.K)).magnitude*
                        1000)

    # TODO Check why different values
    # def test_air_density(self):
    #     fluid = ht.ThermState('air', backend='REFPROP')
    #     T = (((70 - 32) * 5 / 9) + 273.15) * u.K
    #     P = 14.7 / 14.50377377 * (10**5) / 1000 * u.kPa
    #     fluid.update_kw(P=P, T=T)
    #     self.assertNIST(0.0749156384666842,
    #                     fluid.Dmolar.to(u.mol/u.L).magnitude*fluid.M *
    #                     0.062427974)

    # def test_freon_mixture(self):
    #     fluid = ht.ThermState('R32&R125', backend='REFPROP')
    #     fluid.set_mole_fractions(0.3, 0.7)
    #     P = 10 * 1000 * u.kPa
    #     Smolar = 110 * u.J / (u.mol * u.K)
    #     fluid.update_kw(P=P, Smolar=Smolar)
    #     self.assertNIST(23643.993624382,
    #                     fluid.Hmolar.to(u.J/u.mol).magnitude)
    # TODO Figure out the issue with mixtures (possibly refprop_setref/ixflag thing)

    # def test_ethane_butane(self):
    #     fluid = ht.ThermState('ethane&butane', backend='REFPROP')
    #     fluid.set_mole_fractions(0.5, 0.5)
    #     Dmolar = 30 * 0.45359237 / 0.028316846592 / fluid.M * u.mol / u.L
    #     Hmolar = 283 * 1.05435026448889 / 0.45359237 * fluid.M * u.J / u.mol
    #     fluid.update_kw(Dmolar=Dmolar, Hmolar=Hmolar)
    #     self.assertNIST(298.431320311048,
    #                     fluid.T.to(u.degF).magnitude)
    # TODO Figure out the issue with mixtures (possibly refprop_setref/ixflag thing)

    # def test_ammonia_water(self):
    #     fluid = ht.ThermState('ammonia&water', backend='REFPROP')
    #     fluid.set_mole_fractions(0.4, 0.6)
    #     T = (((300 - 32) * 5 / 9) + 273.15) * u.K
    #     P = 10000 / 14.50377377 * (10**5) / 1000 * u.kPa
    #     fluid.update_kw(T=T, P=P)
    #     self.assertNIST(5536.79144924071,
    #                     fluid.speed_sound.to(u.m/u.s).magnitude *
    #                     1000 / 25.4 / 12)
    # TODO Check why the values are different

    # def test_octane(self):
    #     fluid = ht.ThermState('octane', backend='REFPROP')
    #     T = (100+273.15) * u.K
    #     fluid.update_kw(T=T, Q=0)
    #     self.assertNIST(319.167499870568,
    #                     fluid.latent_heat.to(u.J/u.g).magnitude)
    # TODO Check why the values are different

    # def test_R410A_mole_fraction(self):
    #     fluid = ht.ThermState('R410A.MIX', backend='REFPROP')
    #     self.assertNIST(0.697614699375863,
    #                     fluid.mole_fractions[0])

    # def test_R410A_mass_fraction(self):
    #     fluid = ht.ThermState('R410A.MIX', backend='REFPROP')
    #     self.assertNIST(0.5,
    #                     fluid.mass_fractions[0])


class CpWrapperTest(unittest.TestCase):
    def test_super_critical(self):
        fluid = ctb.ThermState('helium', P=1*u.bar, T=300*u.K)
        self.assertTrue(fluid.is_super_critical)
        fluid.update_kw(T=2.2*u.K, P=100*u.bar)
        self.assertTrue(fluid.is_super_critical)
        fluid.update_kw(T=300*u.K, P=100*u.bar)
        self.assertTrue(fluid.is_super_critical)

class FunctionsTest(unittest.TestCase):
    def assertApproxEqual(self, data, calc, uncertainty=0.1):
        if isinstance(data, u.Quantity):
            calc.ito(data.units)
        assert abs(data-calc) / data < uncertainty, \
            f'Calculated value {calc} is not within {uncertainty:.1%} of data ' + \
            f'value {data}'

    def test_rad_hl_1(self):
        eps = 0.02
        T1 = Q_(27, u.degC)
        T2 = Q_(-183, u.degC)
        A = 0.05 * u.m**2
        heat_flow = abs(ht.rad_hl(T1, eps, T2, eps)) * A
        self.assertAlmostEqual(0.23054, heat_flow.to(u.W).magnitude, 5)

    def test_rad_hl_2(self):
        eps = 0.55
        T1 = Q_(27, u.degC)
        T2 = Q_(-183, u.degC)
        A = 0.05 * u.m**2
        heat_flow = abs(ht.rad_hl(T1, eps, T2, eps)) * A
        self.assertAlmostEqual(8.65727, heat_flow.to(u.W).magnitude, 5)

    def test_rad_hl_3(self):
        eps = 0.8
        T1 = Q_(327, u.degC)
        T2 = Q_(127, u.degC)
        eps_b = 0.05
        heat_flow = abs(ht.rad_hl(T1, eps, T2, eps,
                                  baffles={'N': 1, 'eps': eps_b}))
        self.assertAlmostEqual(145.7373409,
                               heat_flow.to(u.W/u.m**2).magnitude, 5)

    def test_rad_hl_4(self):
        T1 = Q_(327, u.degC)
        eps_1 = 0.8
        T2 = Q_(127, u.degC)
        eps_2 = 0.4
        eps_b = 0.05
        heat_flow = abs(ht.rad_hl(T1, eps_1, T2, eps_2,
                                  baffles={'N': 1, 'eps': eps_b}))
        self.assertAlmostEqual(141.3739475,
                               heat_flow.to(u.W/u.m**2).magnitude, 5)

    @unittest.skipIf(skip_RP_test, 'REFPROP failed to initialize, skipping REFPROP tests.')
    def test_API_subsonic(self):
        """Example from API 5.6"""
        m_dot = 53500 * u.lb/u.hr
        fluid = ht.ThermState('propane&butane', backend='REFPROP')
        fluid.set_mole_fractions(0.5, 0.5)
        fluid.update_kw(T=627*u.degR, P=97.2*u.psi)
        P_back = 77.2 * u.psi
        A_expect = 6.55 * u.inch**2
        A_calc = ht.A_relief_API(m_dot, fluid, P_back=P_back)
        self.assertApproxEqual(A_expect, A_calc, uncertainty=0.05)

    @unittest.skipIf(skip_RP_test, 'REFPROP failed to initialize, skipping REFPROP tests.')
    def test_API_sonic(self):
        """Example from API 5.6"""
        m_dot = 53500 * u.lb/u.hr
        fluid = ht.ThermState('propane&butane', backend='REFPROP')
        fluid.set_mole_fractions(0.5, 0.5)
        fluid.update_kw(T=627*u.degR, P=97.2*u.psi)
        P_back = Q_(0, u.psig)
        A_expect = 5.73 * u.inch**2
        A_calc = ht.A_relief_API(m_dot, fluid, P_back=P_back)
        self.assertApproxEqual(A_expect, A_calc, uncertainty=0.05)

    def test_to_standard_flow(self):
        self.assertRaises(ValueError, ctb.to_standard_flow, 1*u.inch, ctb.AIR)

    def test_theta(self):
        # Subcritical test
        fluid = ctb.ThermState('helium', P=1*u.bar, Q=0.5)
        self.assertAlmostEqual(4.20982594*u.K, ctb.cga.theta(fluid))
        # Supercritical test
        fluid.update_kw(P=20*u.bar, T=10*u.K)
        self.assertAlmostEqual(11.9574303*u.K, ctb.cga.theta(fluid))

    def test_calculate_fluid_FR(self):
        # Subcritical test
        fluid = ctb.ThermState('helium', P=0.9090909090909091*u.bar, Q=0.5)
        fluid_FR = ctb.cga.calculate_fluid_FR(fluid)
        self.assertAlmostEqual(4.20982594*u.K, fluid_FR.T)
        # Supercritical test
        fluid.update_kw(P=18.18181818181818*u.bar, T=10*u.K)
        fluid_FR = ctb.cga.calculate_fluid_FR(fluid)
        self.assertAlmostEqual(11.9574303*u.K, fluid_FR.T)

    def test_G_i(self):
        test_fluid = ctb.ThermState('air', P=Q_(100, u.psig), Q=1)
        self.assertApproxEqual(10.2, ctb.cga.G_i(test_fluid), uncertainty=0.05)
        test_fluid = ctb.ThermState('hydrogen', P=100*u.psi, Q=1)
        self.assertApproxEqual(8.272, ctb.cga.G_i(test_fluid), uncertainty=0.01)
        test_fluid = ctb.ThermState('hydrogen', P=200*u.psi, T=ctb.T_NTP)
        T = ctb.cga.theta(test_fluid)
        test_fluid.update_kw(P=test_fluid.P, T=T)
        self.assertApproxEqual(13.52, ctb.cga.G_i(test_fluid), uncertainty=0.01)

    def test_G_u(self):
        test_fluid = ctb.ThermState('air', P=Q_(100, u.psig), Q=1)
        self.assertApproxEqual(59, ctb.cga.G_u(test_fluid), uncertainty=0.05)
        test_fluid = ctb.ThermState('air', P=Q_(200, u.psig), Q=1)
        self.assertApproxEqual(69, ctb.cga.G_u(test_fluid), uncertainty=0.05)

    def test_average_spec_heat(self):
        # CGA S-1.3 2008 Table 5
        table_5 = {'nitrogen': (2170*u.kPa, 1.2167*u.kJ/u.kg/u.K),
                   'argon': (2170*u.kPa, 0.6225*u.kJ/u.kg/u.K),
                   'hydrogen': (1067*u.kPa, 13.3701*u.kJ/u.kg/u.K),
                   'neon': (2170*u.kPa, 1.1245*u.kJ/u.kg/u.K),
                   }
        for gas, (P, Cp) in table_5.items():
            fluid = ctb.ThermState(gas, P=P, Q=1)
            Cp_calc = ctb.cga._average_spec_heat(fluid)
            self.assertApproxEqual(Cp, Cp_calc, 0.15)

    def test_calculate_inlet_temp(self):
        # Test of basic functionality
        fluid = ctb.ThermState('argon', P=2*u.bar, Q=0)
        m_dot = 100 * u.g/u.s
        pipe = ctb.piping.Pipe(1, L=10*u.m)
        Ti = ctb.cga.calculate_inlet_temp(fluid, m_dot, pipe, condition=None)

    def test_mean_free_path(self):
        """Barron, Example 9.1"""
        fluid = ctb.ThermState('helium', T=290*u.K, P=5*u.mtorr)
        mfp_exp = 28.2*u.mm
        self.assertApproxEqual(mfp_exp, ctb.mean_free_path(fluid), 0.01)

    def test_Kv_Cv(self):
        """
        Test conversion between flow coefficients Kv and Cv.

        Verifies that:
        - Cv_to_Kv conversion is accurate
        - Round-trip conversion (Cv_to_Kv and Kv_to_Cv) is accurate
        """
        Cv = 1
        self.assertApproxEqual(0.86, ctb.Cv_to_Kv(Cv), 0.86)
        self.assertApproxEqual(Cv, ctb.Kv_to_Cv(ctb.Cv_to_Kv(Cv)), Cv)

    def test_kt_Dacron(self):
        """Test the thermal conductivity of Dacron spacer and MLI"""
        tc = 80 * u.K
        ks_exp = 0.1409 * u.W/u.m/u.K
        ks = ctb.ks_Dacron(tc)
        self.assertApproxEqual(ks_exp, ks, uncertainty=0.001)
        ts = 0.06 * u.mm
        hc_exp = 1.3526 * u.W/u.m**2/u.K
        hc = ctb.hc_Dacron(tc, ts=ts)
        self.assertApproxEqual(hc_exp, hc, uncertainty=0.001)
        kt_exp = 0.465 * u.mW/u.m/u.K
        kt = ctb.k_MLI(30/u.cm, hc, 80*u.K)
        self.assertApproxEqual(kt_exp, kt, uncertainty=0.001)

class PipingTest(unittest.TestCase):
    """Piping checks, mostly taken from textbooks.
    """

    def assertApproxEqual(self, expected, calculated, uncertainty=0.05):
        error = float(abs(expected-calculated)/expected)
        if isinstance(expected, u.Quantity) and isinstance(calculated, u.Quantity):
            unit = expected.units
            calculated.ito(unit)
        assert error < uncertainty, \
            f'Calculated value {calculated:g} is not within {uncertainty:.1%} of ' + \
            f'expected value {expected:g}'

    def test_create_pipes(self):
        vj_pipe = ht.piping.VJPipe(Q_('1 inch'), VJ_D=Q_('2 inch'))
        entrance = ht.piping.Entrance(Q_('1 inch'))
        pipe_exit = ht.piping.Exit(Q_('1 inch'))
        orifice = ht.piping.Orifice(Q_('1 inch'))
        c_orifice = ht.piping.ConicOrifice(1, Q_('3/4 inch'))
        annulus = ht.piping.Annulus(Q_('1 inch'), Q_('3/4 inch'))
        pipe_tee = ht.piping.PipeTee(Q_('1 inch'), N=2)
        tee = ht.piping.Tee(Q_('1 inch'), N=2)
        valve = ht.piping.Valve(Q_('1 inch'), 1)
        # g_valve = ht.piping.GlobeValve(Q_('1 inch'))
        # print(f'Generated {g_valve}')
        # v_cone = ht.piping.VCone(Q_('1 inch'), 0.7, 1)
        # print(f'Generated {v_cone}')
        cont = ht.piping.Contraction(1*u.inch, 0.5*u.inch)
        enl = ht.piping.Enlargement(0.5*u.inch, 1*u.inch)
        test_state = ht.AIR
        piping = [
            vj_pipe, entrance, pipe_exit, orifice, c_orifice,
            annulus, pipe_tee, tee, valve,
            # g_valve, v_cone,
            cont, enl
        ]
        ht.piping.volume(piping)
        dP = ht.piping.dP_incomp(Q_('10 g/s'), test_state, piping)
        self.assertApproxEqual(21.2*u.psi, dP)

    def test_f_Darcy(self):
        eps_smooth = 0.0018 * u.inch
        Re = 1e8
        ID = 0.2 * u.inch
        eps_r = eps_smooth/ID
        self.assertApproxEqual(0.0368, ht.piping.f_Darcy(Re, eps_r))
        ID = 0.4 * u.inch
        eps_r = eps_smooth/ID
        self.assertApproxEqual(0.0296, ht.piping.f_Darcy(Re, eps_r))
        ID = 1 * u.inch
        eps_r = eps_smooth/ID
        self.assertApproxEqual(0.0228, ht.piping.f_Darcy(Re, eps_r))
        ID = 2 * u.inch
        eps_r = eps_smooth/ID
        self.assertApproxEqual(0.0191, ht.piping.f_Darcy(Re, eps_r))
        eps_r = 0.006
        Re = 1e5
        self.assertApproxEqual(0.033, ht.piping.f_Darcy(Re, eps_r))
        eps_r = 0.006
        Re = 1e3
        self.assertApproxEqual(64/Re, ht.piping.f_Darcy(Re, eps_r))

    # def test_Crane_4_22(self):
    #     test_air = ht.ThermState('air', P=2.343*u.bar, T=40*u.degC)
    #     pipe = ht.piping.Pipe(1/2, SCH=80, L=3*u.m)
    #     piping = ht.piping.Piping(
    #         test_air,
    #         [pipe,
    #          ht.piping.Exit(pipe.ID)])
    #     Y = 0.76  # Taken from Crane; temporary stub
    #     flow = ht.to_standard_flow(piping.m_dot(), test_air) * Y
    #     # TODO Check this test (might have to do with subsonic flow)
    #     self.assertAlmostEqual(1.76, flow.to(u.m**3/u.min).magnitude)

    def test_Crane_7_16(self):
        air = ht.ThermState('air', P=Q_(65, u.psig), T=Q_(110, u.degF))
        pipe = ht.piping.Pipe(1, SCH=40, L=75*u.ft)
        flow = 100 * u.ft**3/u.min
        m_dot = flow * ht.AIR.Dmass
        piping = [pipe, ht.piping.Exit(pipe.ID)]
        dP = ht.piping.dP_incomp(m_dot, air, piping)
        self.assertApproxEqual(2.61, dP.m_as(u.psi))
        dP_isot = ht.piping.dP_isot(m_dot, air, pipe)
        self.assertApproxEqual(2.61, dP_isot.m_as(u.psi))

    def test_Crane_7_22(self):
        air = ht.ThermState('air', P=Q_(19.3, u.psig), T=Q_(100, u.degF))
        pipe = ht.piping.Pipe(1/2, SCH=80, L=7.04/6.04*10*u.ft)  # Adjust for entrance K = 1, total 7.04
        dP = 19.3 * u.psi
        P_out = air.P - dP
        q_expected = 3762 * u.ft**3/u.hr  # STD flow
        m_dot_expected = ht.to_mass_flow(q_expected, air)
        Re = ht.Re(air, m_dot_expected, pipe.ID, pipe.area)
        f = 0.0275
        m_dot_isot = ht.piping.m_dot_isot(air, pipe)
        self.assertApproxEqual(m_dot_expected, m_dot_isot)
        # M = Mach(air, m_dot_expected, pipe.area)
        # # print('Mach at inlet: ', M)
        # K_lim = K_lim(M, air.gamma)
        # # print('K limit: ', K_lim)
        # # print('K pipe: ', pipe.K(Re).to_base_units())
        # check(M, M_Klim(K_lim, air.gamma))
        # K_left = K_lim - pipe.K(Re)
        # # print('K left: ', K_left)
        # M_end = M_Klim(K_left, air.gamma)
        # # print('M end: ', M_end)
        # M_out = Mach(ht.ThermState('air', P=P_out, T=air.T), m_dot_expected, pipe.area)
        # check(M_out, M_end)
        # P_static_end = P_from_M(air.P, M, M_end, air.gamma)
        # # print('P static end: ', P_static_end)
        # P_static_out = P_from_M(air.P, M, M_out, air.gamma)
        # check(P_static_end, P_static_out)
        # P_total_end = P_total(P_static_end, M, air.gamma)
        # # print('P total end: ', P_total_end)
        # print(f"Mach dP works: {check(P_out, P_total_end)}")

    def test_Rennels_4_5(self):
        pipe = ht.piping.Pipe(4, SCH=40, L=100*u.ft)
        fluid = ht.ThermState('nitrogen', P=100*u.psi, T=530*u.degR)
        piping = [pipe]
        P_out = 84.056 * u.psi
        m_dot_isot = ht.piping.m_dot_isot(fluid, pipe, P_out)
        self.assertApproxEqual(10, m_dot_isot.m_as(u.lb/u.s))

    def test_White_9_3(self):
        v = 240 * u.m/u.s
        T1 = 320 * u.K
        P1 = 170 * u.kPa
        T01 = 349 * u.K
        P01 = 230 * u.kPa
        fluid_static = ht.ThermState('air', P=P1, T=T1)
        area = 0.05 * u.m**2  # Assumed value
        M = ht.piping.Mach(fluid_static, v)
        M_exp = 0.67
        self.assertApproxEqual(M_exp, M, uncertainty=0.005)
        fluid_total = ht.ThermState('air', P=P01, T=T01)
        m_dot = v * area * fluid_static.Dmass
        M_total = ht.piping.Mach_total(fluid_total, m_dot, area)
        self.assertApproxEqual(M_exp, M_total, uncertainty=0.005)

    def test_White_9_4(self):
        fluid_static = ht.ThermState('air', P=500*u.kPa, T=470*u.K)
        area = 0.05 * u.m**2
        v = 180 * u.m/u.s
        m_dot = v * area * fluid_static.Dmass
        M = ht.piping.Mach(fluid_static, v)
        M_exp = 0.414
        self.assertApproxEqual(M_exp, M, uncertainty=0.005)
        P0 = 563 * u.kPa
        T0 = 486 * u.K
        fluid_total = ht.ThermState('air', P=P0, T=T0)
        M_total = ht.piping.Mach_total(fluid_total, m_dot, area)
        self.assertApproxEqual(M_exp, M_total, uncertainty=0.005)

    def test_White_9_10(self):
        Ma1 = 0.1
        Ma2 = 0.5
        Ma3 = 1
        k = ht.AIR.gamma
        f = 0.024
        d = 2 * u.cm
        K_lim1 = ht.piping.K_lim(Ma1, k)
        L1 = K_lim1 * d / f
        K_lim2 = ht.piping.K_lim(Ma2, k)
        L2 = K_lim2 * d / f
        K_lim3 = ht.piping.K_lim(Ma3, k)
        L3 = K_lim3 * d / f
        self.assertApproxEqual(55*u.m, L1-L2)
        self.assertApproxEqual(0.9*u.m, L2-L3)

    def test_White_9_11(self):
        f = 0.024
        d = 2 * u.cm
        Ma1 = 0.1
        Ma2 = 0.5
        P1 = 600 * u.kPa
        T1 = 450 * u.K
        fluid = ht.ThermState('air', P=P1, T=T1)
        k = fluid.gamma
        P2 = ht.piping.P_from_M(P1, Ma1, Ma2, k)
        self.assertApproxEqual(117*u.kPa, P2)
        P02 = ht.piping.P_total(P2, Ma2, k)
        self.assertApproxEqual(139*u.kPa, P02)

    def test_White_9_12(self):
        fluid_total = ht.ThermState('air', P=200*u.kPa, T=500*u.K)
        v = 100 * u.m/u.s
        pipe = ht.piping.Tube(3*u.cm, L=15*u.m)
        m_dot = 0.0961 * u.kg/u.s  # Based on static density from the source
        f = 0.02
        Ma = ht.piping.Mach_total(fluid_total, m_dot, pipe.area)
        K_lim = ht.piping.K_lim(Ma, fluid_total.gamma)
        Lstar = K_lim * pipe.ID / f
        self.assertApproxEqual(16.5*u.m, Lstar)

    def test_White_9_13(self):
        fluid_static = ht.ThermState('air', P=220*u.kPa, T=300*u.K)
        P_out = 140 * u.kPa
        pipe = ht.piping.Tube(1*u.cm, wall=0*u.m, L=1.2*u.m)
        K = 0.025 * pipe.L/pipe.ID
        pipe.K = MagicMock(return_value=K)
        m_dot_adiab = ht.piping.m_dot_adiab(fluid_static, pipe, P_out, state='static')
        self.assertApproxEqual(0.0233*u.kg/u.s, m_dot_adiab, uncertainty=0.001)

    def test_White_P_9_100(self):
        fluid = ht.ThermState('methane', T=Q_(68, u.degF), P=5*u.bar+101325*u.Pa)
        P_out = 1*u.bar + 101325*u.Pa
        pipe = ht.piping.Pipe(6, L=31*u.mile)
        K = 0.019 * pipe.L/pipe.ID
        pipe.K = MagicMock(return_value=K)
        m_dot_isot = ht.piping.m_dot_isot(fluid, pipe, P_out)
        self.assertApproxEqual(0.345*u.kg/u.s, m_dot_isot)
        dP_isot = ht.piping.dP_isot(m_dot_isot, fluid, pipe)
        self.assertApproxEqual(fluid.P-P_out, dP_isot)

    def test_White_P_9_101(self):
        air = ht.ThermState('air', T=Q_(20, u.degC), P=102*u.kPa)
        pipe = ht.piping.Tube(1*u.cm, L=3*u.m)
        P_out = 100*u.kPa
        K = 0.028 * pipe.L/pipe.ID
        pipe.K = MagicMock(return_value=K)
        m_dot_incomp = ht.piping.m_dot_incomp(air, [pipe], P_out=P_out)
        m_dot_isot = ht.piping.m_dot_isot(air, pipe, P_out)
        self.assertApproxEqual(m_dot_incomp, m_dot_isot, uncertainty=0.01)
        m_dot_adiab = ht.piping.m_dot_adiab(air, pipe, P_out, state='static')
        self.assertApproxEqual(m_dot_incomp, m_dot_adiab, uncertainty=0.01)
        dP_incomp = ht.piping.dP_incomp(m_dot_incomp, air, [pipe])
        self.assertApproxEqual(air.P-P_out, dP_incomp, uncertainty=0.01)
        dP_isot = ht.piping.dP_isot(m_dot_isot, air, pipe)
        self.assertApproxEqual(air.P-P_out, dP_isot, uncertainty=0.01)
        dP_adiab = ctb.piping.dP_adiab(m_dot_adiab, air, pipe, state='static')
        self.assertApproxEqual(air.P-P_out, dP_adiab, uncertainty=0.01)


    def test_Cengel_12_15(self):
        pipe = ctb.piping.Tube(3*u.cm)
        fluid = ctb.ThermState('air', P=150*u.kPa, T=300*u.K)
        Ma1 = 0.4
        v1 = fluid.speed_sound * Ma1
        self.assertApproxEqual(139*u.m/u.s, v1, uncertainty=0.001)
        m_dot = fluid.Dmass*v1*pipe.area
        Re1 = 2.637e5  # Differences in viscosity compared to the source
        f = ctb.piping.serghide(Re1, 0)
        self.assertApproxEqual(0.0148, f, uncertainty=0.01)
        K_lim = ctb.piping.K_lim(Ma1, fluid.gamma)
        L = K_lim * pipe.ID / f
        self.assertApproxEqual(4.68*u.m, L, uncertainty=0.01)


    def test_M_from_K(self):
        air = ht.ThermState('air', P=300*u.kPa, T=500*u.K)
        self.assertApproxEqual(0.225, ht.piping.M_Klim(
            11, air.gamma))
        self.assertApproxEqual(0.2, ht.piping.M_Klim(
            14.533, air.gamma))
        self.assertApproxEqual(0.15, ht.piping.M_Klim(
            27.932, air.gamma))
        self.assertApproxEqual(0.17, ht.piping.M_Klim(
            21.115, air.gamma))
        self.assertApproxEqual(0.175, ht.piping.M_Klim(
            19.772, air.gamma))
        self.assertApproxEqual(0.1741, ht.piping.M_Klim(
            20, air.gamma))
        self.assertApproxEqual(1, ht.piping.M_Klim(
            0, air.gamma))

    def test_piping_req_thick(self):
        pipe = ctb.piping.Tube(1*u.inch, wall=0.065*u.inch)
        S = 16700*u.psi
        E = 0.8
        W = 1
        Y = 0.4
        P_max = ctb.piping.pressure_rating(pipe, S=S, E=E, W=W, Y=Y)
        t_m = ctb.piping.pressure_req_thick(pipe, P_max, S=S, E=E, W=W, Y=Y)
        self.assertAlmostEqual(t_m, pipe.T)

    def test_direct_integration(self):
        fluid = ht.ThermState('helium', P=200*u.psi, T=7*u.K)
        self.assertApproxEqual(13920*u.kg/(u.s*u.m**2),
                               ht.piping.G_nozzle(fluid, P_out=1*u.atm),
                               uncertainty=0.01)

    # TODO Add Crane examples: 4-22 (may need Y implementation),
    # 4-20, 4-19, 4-18, 4-16, 4-12?, 4-10?

    def test_packed_bed(self):
        fluid = ht.ThermState('argon', T=Q_(80, u.degF), P=126.2*u.psi)
        filter_shell = ht.piping.Pipe(12, SCH=10)
        filter_L = 3.33 * u.ft
        eps = 0.37  # filter void fraction
        m_dot = 419.5 * u.lb/u.hr
        x_ox = 0.00336 * u.ft  # oxygen filter particle size
        x_ms = 0.00666 * u.ft  # molsieve filter particle size

        ox_filter = ht.piping.PackedBed(filter_shell.ID, filter_L, x_ox, eps)
        ms_filter = ht.piping.PackedBed(filter_shell.ID, filter_L, x_ms, eps)
        # TODO use dP_incomp once it is simplified to only incompressible
        rho = fluid.Dmass
        mu = fluid.viscosity
        U_s_ox = m_dot / (ox_filter.area*rho)
        Re_s = U_s_ox*ox_filter.ID*rho/mu
        dP_ox = ht.piping.dP_Darcy(ox_filter.K(Re_s), rho, U_s_ox)
        self.assertApproxEqual(0.261*u.psi, dP_ox, uncertainty=0.1)
        U_s_ms = m_dot / (ms_filter.area*rho)
        Re_s = U_s_ms*ms_filter.ID*rho/mu
        dP_ms = ht.piping.dP_Darcy(ms_filter.K(Re_s), rho, U_s_ms)
        self.assertApproxEqual(0.092*u.psi, dP_ms, uncertainty=0.1)

    def test_dP_incomp(self):
        m_dot = 1 * u.kg/u.s
        fluid = ctb.AIR
        piping = [ctb.piping.Pipe(1/2, L=100*u.m)]
        # Execution below should not produce an error
        # despite process not being physical
        ctb.piping.dP_incomp(m_dot, fluid, piping)

    def test_reinforcement_area_H301(self):
        header = ctb.piping.Pipe(8, c=2.5*u.mm)
        branch = ctb.piping.Pipe(4, c=2.5*u.mm)
        P_des = 2068 * u.kPa
        S = 110 * u.MPa
        E = 1
        W = 1
        Y = 0.4
        self.assertAlmostEqual(header.T, 7.16*u.mm, delta=0.005*u.mm)
        self.assertAlmostEqual(branch.T, 5.27*u.mm, delta=0.005*u.mm)
        result = ctb.piping.calculate_branch_reinforcement(header,
                                                         branch,
                                                         P_des,
                                                         S=S,
                                                         E=E,
                                                         W=W,
                                                         Y=Y)
        tmh = ctb.piping.pressure_req_thick(header, P_des, S=S, E=E, W=W, Y=Y).to(u.mm)
        self.assertAlmostEqual(tmh-header.c,
                               2.04*u.mm, delta=0.005*u.mm)
        tmb = ctb.piping.pressure_req_thick(branch, P_des, S=S, E=E, W=W, Y=Y).to(u.mm)
        self.assertAlmostEqual(tmb-branch.c,
                               1.07*u.mm, delta=0.005*u.mm)
        A1, A2, A3, A4, d2, safe = result
        self.assertAlmostEqual(A1.to(u.mm**2), 222*u.mm**2, delta=0.5*u.mm**2)
        self.assertAlmostEqual(A2.to(u.mm**2), 285*u.mm**2, delta=1*u.mm**2)
        self.assertAlmostEqual(A3.to(u.mm**2), 24*u.mm**2, delta=0.5*u.mm**2)
        self.assertAlmostEqual(d2.to(u.mm), 108.8*u.mm, delta=0.05*u.mm)
        self.assertAlmostEqual(safe, True)

    def test_reinforcement_area_H313(self):
        header = ctb.piping.Pipe(16, SCH=40, c=0.1*u.inch)
        branch = ctb.piping.Pipe(6, SCH=40, c=0.1*u.inch)
        beta = 60 * u.deg
        Lr = 12 * u.inch
        Tr = 0.5 * u.inch
        P_des = Q_(500, u.psi)
        E = 1
        W = 1
        Y = 0.4
        S = 14.4 * u.kpsi
        self.assertAlmostEqual(header.T, 0.438*u.inch, delta=0.001*u.inch)
        self.assertAlmostEqual(branch.T, 0.245*u.inch, delta=0.001*u.inch)
        tmh = ctb.piping.pressure_req_thick(header, P_des, S=S, E=E, W=W, Y=Y).to(u.inch)
        self.assertAlmostEqual(tmh-header.c, 0.274*u.inch, delta=0.001*u.inch)
        tmb = ctb.piping.pressure_req_thick(branch, P_des, S=S, E=E, W=W, Y=Y).to(u.inch)
        self.assertAlmostEqual(tmb-header.c, 0.113*u.inch, delta=0.001*u.inch)
        result = ctb.piping.calculate_branch_reinforcement(header,
                                                           branch,
                                                           P_des,
                                                           beta=beta,
                                                           Tr=Tr,
                                                           Lr=Lr,
                                                           S=S,
                                                           E=E,
                                                           W=W,
                                                           Y=Y)
        A1, A2, A3, A4, d2, safe = result
        self.assertAlmostEqual(A1.to(u.inch**2), 2.27*u.inch**2, delta=0.01*u.inch**2)
        self.assertAlmostEqual(A2.to(u.inch**2), 0.468*u.inch**2, delta=0.005*u.inch**2)
        self.assertAlmostEqual(A3.to(u.inch**2), 0.062*u.inch**2, delta=0.001*u.inch**2)
        self.assertAlmostEqual(A4.to(u.inch**2), 2.175*u.inch**2, delta=0.001*u.inch**2)

    def test_piping_stored_energy(self):
        fluid = ctb.ThermState('nitrogen', P=100*u.psi, T=ctb.T_NTP)
        short_pipe = ctb.piping.Pipe(1, L=2*u.inch)
        exp_volume = pi*short_pipe.ID**2/4 * short_pipe.L
        self.assertAlmostEqual(ctb.stored_energy(fluid, exp_volume),
                               ctb.piping.piping_stored_energy(fluid, [short_pipe]))
        long_pipe = ctb.piping.Pipe(1, L=100*u.ft)
        exp_volume = pi*long_pipe.ID**2/4 * 8*long_pipe.ID
        self.assertAlmostEqual(ctb.stored_energy(fluid, exp_volume),
                               ctb.piping.piping_stored_energy(fluid, [long_pipe]))

    def test_Kv_Cv(self):
        """
        Test the conversion functions between Cv and Kv.

        This test verifies that:

        1. Converting Cv (value 1) to Kv using Cv_to_Kv returns approximately 0.86,
            within a relative error tolerance of 1e-2.
        2. Converting Cv to Kv and then back to Cv (using Kv_to_Cv) returns the original Cv,
            with a relative difference less than 1e-6.

        Raises
        ------
        AssertionError
            If either conversion does not meet the specified tolerance levels.

        Notes
        -----
        The test relies on the correct implementation of Cv_to_Kv and Kv_to_Cv.
        """
        Cv = 1
        computed_Kv = ctb.piping.Cv_to_Kv(Cv)
        expected_Kv = 0.86
        rel_error_Kv = abs((computed_Kv - expected_Kv) / expected_Kv)
        self.assertLess(rel_error_Kv, 1e-2,
                        msg=f"Expected Kv ~ {expected_Kv}, but got {computed_Kv} with relative error {rel_error_Kv}")

        Cv_converted = ctb.piping.Kv_to_Cv(computed_Kv)
        rel_error_Cv = abs((Cv_converted - Cv) / Cv)
        self.assertLess(rel_error_Cv, 1e-6,
                        msg=f"Round-trip conversion error is {rel_error_Cv}, which exceeds 1e-6")


class CPWrapperTest(unittest.TestCase):
    """Test for additional methods of ThermState class"""
    def test_copy(self):
        """Minimal test of copy method"""
        test_state = ht.ThermState('nitrogen')
        test_state.copy()
        test_state.update_kw(T=300*u.K, P=1*u.MPa)
        test_state.copy()

    def test_standard(self):
        air_std = ht.AIR.to_standard(conditions='NTP')
        self.assertAlmostEqual(ht.P_NTP, air_std.P)
        self.assertAlmostEqual(ht.T_NTP, air_std.T)
        air_std = ht.AIR.to_standard(conditions='MSC')
        self.assertAlmostEqual(ht.P_MSC, air_std.P)
        self.assertAlmostEqual(ht.T_MSC, air_std.T)
        air_std = ht.AIR.to_standard(conditions='STD')
        self.assertAlmostEqual(ht.P_STD, air_std.P)
        self.assertAlmostEqual(ht.T_STD, air_std.T)

    def test_str(self):
        argon = ctb.ThermState('argon', P=ctb.P_NTP, T=ctb.T_NTP)
        self.assertEqual('argon at T: 293 K and P: 1.01 bar', argon.__str__())
        argon.update_kw(P=argon.P, Q=0)
        self.assertEqual('saturated argon liquid at P: 1.01 bar', argon.__str__())
        argon.update_kw(P=argon.P, Q=1)
        self.assertEqual('saturated argon vapor at P: 1.01 bar', argon.__str__())
        argon.update_kw(P=argon.P, Q=0.5)
        self.assertEqual('two-phase argon at P: 1.01 bar and Q: 0.50', argon.__str__())


class ODHTest(FunctionsTest):
    """Set of tests for basic functionality of ODH analysis
    """
    def test_hole_leak_Crane_7_23(self):
        fluid = ht.ThermState('water', P=Q_(4.4, u.inHg)+ht.P_NTP,
                              T=Q_(60, u.degF))
        d_hole = 2 * u.inch
        area = pi * d_hole**2 / 4
        self.assertApproxEqual(106*u.gal/u.min,
                               odh.hole_leak(area, fluid),
                               uncertainty=0.05)
        # Crane gives 131 gal/min for flow
        # difference in equations gives 20 % difference

    def test_hole_leak_voirin(self):
        P_150 = Q_(150, u.psig)
        P_15 = Q_(15, u.psig)
        T = Q_(60, u.degF)
        V_flux_150 = 2.856 * u.ft**3 / (u.min * u.mm**2)
        V_flux_15 = 0.515 * u.ft**3 / (u.min * u.mm**2)

        area = 1*u.mm**2  # Chosen for simplicity
        ID = (4*area/pi)**0.5  # Assuming tube ID based on area
        tube = ht.piping.Tube(ID, wall=0*u.m)
        fluid = ht.ThermState('nitrogen', P=P_150, T=T)
        m_dot = ht.odh.hole_leak(area, fluid)
        self.assertApproxEqual(V_flux_150*area, m_dot)
        fluid = ht.ThermState('nitrogen', P=P_15, T=T)
        m_dot = ht.odh.hole_leak(area, fluid)
        self.assertApproxEqual(V_flux_15*area, m_dot)

    def test_pipe_failure(self):
        fluid = ht.ThermState('nitrogen', P=100*u.psi, Q=0)
        source = odh.Source('Test source', fluid, 100*u.L)
        zero_length_tube = ht.piping.Pipe(3, SCH=10)
        self.assertRaises(odh.ODHError,
                          source.add_pipe_failure, zero_length_tube, fluid)
        N = 2
        tube = ht.piping.Pipe(3, SCH=10, L=N*1*u.m)
        D_t = tube.OD / tube.wall
        source.add_pipe_failure(tube, fluid, N_welds=N)
        leaks = [
            odh.Leak('', N*1e-9/u.hr, fluid,
                     odh.hole_leak(10*u.mm**2, fluid), 0*u.s, 1),
            odh.Leak('', N*1e-10/u.hr, fluid,
                     odh.hole_leak(1000*u.mm**2, fluid), 0*u.s, 1),
            odh.Leak('', N*3e-11/u.hr, fluid,
                     odh.hole_leak(tube.area, fluid), 0*u.s, 1),
            odh.Leak('', N*2e-11*D_t/u.hr, fluid,
                     odh.hole_leak(10*u.mm**2, fluid), 0*u.s, 1),
            odh.Leak('', N*2e-12*D_t/u.hr, fluid,
                     odh.hole_leak(1000*u.mm**2, fluid), 0*u.s, 1),
            odh.Leak('', N*6e-13*D_t/u.hr, fluid,
                     odh.hole_leak(tube.area, fluid), 0*u.s, 1),
        ]
        self.assertEqual(len(leaks), len(source.leaks))
        for leak1, leak2 in zip(leaks, source.leaks):
            self.assertAlmostEqual(leak1.failure_rate, leak2.failure_rate)
            self.assertAlmostEqual(leak1.q_std, leak2.q_std)
            self.assertAlmostEqual(source.volume/leak1.q_std, leak2.tau)

    def test_flange_failure(self):
        fluid = ht.ThermState('nitrogen', P=100*u.psi, Q=0)
        source = odh.Source('Test source', fluid, 100*u.L)
        tube = ht.piping.Pipe(3, SCH=10, L=1*u.m)
        N = 2
        source.add_flange_failure(tube, fluid, N=N)
        leaks = [
            odh.Leak('', N*4e-7/u.hr, fluid, odh.hole_leak(10*u.mm**2, fluid),
                     0*u.s, 1),
            odh.Leak('', N*1e-9/u.hr, fluid, odh.hole_leak(tube.area, fluid),
                     0*u.s, 1),
        ]
        self.assertEqual(len(leaks), len(source.leaks))
        for leak1, leak2 in zip(leaks, source.leaks):
            self.assertAlmostEqual(leak1.failure_rate, leak2.failure_rate)
            self.assertAlmostEqual(leak1.q_std, leak2.q_std)
            self.assertAlmostEqual(source.volume/leak1.q_std, leak2.tau)

    def test_valve_failure(self):
        fluid = ht.ThermState('nitrogen', P=100*u.psi, Q=0)
        source = odh.Source('Test source', fluid, 100*u.L)
        tube = ht.piping.Pipe(3, SCH=10, L=1*u.m)
        N = 2
        source.add_valve_failure(tube, fluid, N=N)
        leaks = [
            odh.Leak('', N*1e-8/u.hr, fluid,
                     odh.hole_leak(10*u.mm**2, fluid), 0*u.s, 1),
            odh.Leak('', N*5e-10/u.hr, fluid,
                     odh.hole_leak(tube.area, fluid), 0*u.s, 1),
        ]
        self.assertEqual(len(leaks), len(source.leaks))
        for leak1, leak2 in zip(leaks, source.leaks):
            self.assertAlmostEqual(leak1.failure_rate, leak2.failure_rate)
            self.assertAlmostEqual(leak1.q_std, leak2.q_std)
            self.assertAlmostEqual(source.volume/leak1.q_std, leak2.tau)

    def test_line_failure(self):
        fluid = ht.ThermState('nitrogen', P=100*u.psi, Q=0)
        source = odh.Source('Test source', fluid, 100*u.L)
        tube = ht.piping.Pipe(3, SCH=10, L=1*u.m)
        source.add_line_failure(tube, fluid, N_welds=1, N_reinforced=1, N_soft=0, N_valves=1)
        D_t = tube.OD / tube.wall
        leaks = [
            odh.Leak('', 1e-9/u.hr, fluid,
                     odh.hole_leak(10*u.mm**2, fluid), 0*u.s, 1),
            odh.Leak('', 1e-10/u.hr, fluid,
                     odh.hole_leak(1000*u.mm**2, fluid), 0*u.s, 1),
            odh.Leak('', 3e-11/u.hr, fluid,
                     odh.hole_leak(tube.area, fluid), 0*u.s, 1),
            odh.Leak('', 2e-11*D_t/u.hr, fluid,
                     odh.hole_leak(10*u.mm**2, fluid), 0*u.s, 1),
            odh.Leak('', 2e-12*D_t/u.hr, fluid,
                     odh.hole_leak(1000*u.mm**2, fluid), 0*u.s, 1),
            odh.Leak('', 6e-13*D_t/u.hr, fluid,
                     odh.hole_leak(tube.area, fluid), 0*u.s, 1),
            odh.Leak('', 4e-7/u.hr, fluid,
                     odh.hole_leak(10*u.mm**2, fluid), 0*u.s, 1),
            odh.Leak('', 1e-9/u.hr, fluid,
                     odh.hole_leak(tube.area, fluid), 0*u.s, 1),
            odh.Leak('', 1e-8/u.hr, fluid,
                     odh.hole_leak(10*u.mm**2, fluid), 0*u.s, 1),
            odh.Leak('', 5e-10/u.hr, fluid,
                     odh.hole_leak(tube.area, fluid), 0*u.s, 1),
        ]
        self.assertEqual(len(leaks), len(source.leaks))
        for leak1, leak2 in zip(leaks, source.leaks):
            self.assertAlmostEqual(leak1.failure_rate, leak2.failure_rate)
            self.assertAlmostEqual(leak1.q_std, leak2.q_std)
            self.assertAlmostEqual(source.volume/leak1.q_std, leak2.tau)

    def test_const_leak(self):
        fluid = ht.ThermState('nitrogen', P=100*u.psi, Q=0)
        source = odh.Source('Test source', fluid, 100*u.L)
        q_leak = 1*u.ft**3/u.min
        N_leaks = 5
        source.add_const_leak('Constant leak', fluid, q_leak, N=N_leaks)
        self.assertEqual(1, len(source.leaks))
        self.assertAlmostEqual(N_leaks*q_leak, source.leaks[0].q_std)


class Nu_test(FunctionsTest):
    def test_Nu_vplate(self):
        fluid = ht.AIR
        Pr = fluid.Prandtl
        # Estimated values from Figure 4.7 of Handbook of Heat Transfer
        figure_4_7_HHT = [
            (100, 2.5),
            (1e5, 10),
            (1e6, 17),
            (1e10, 200)
        ]
        for Ra, Nu in figure_4_7_HHT:
            self.assertApproxEqual(Nu, ht.Nu_vplate(Pr, Ra), 0.1)

    def test_C_t_bar_cyl(self):
        self.assertApproxEqual(0.103, ht.C_t_bar_cyl(0.71), uncertainty=1e-5)
        self.assertApproxEqual((0.103+0.109)/2, ht.C_t_bar_cyl((0.71+6)/2), uncertainty=1e-5)

    def test_Nu_hcyl(self):
        Pr = 0.71  # Taken as Pr for air
        figure_4_17_HHT = [
            (1, 1),
            (1e-6, 0.3),
            (1e-9, 0.2),
            (1e5, 8)
        ]
        for Ra, Nu in figure_4_17_HHT:
            self.assertApproxEqual(Nu, ht.Nu_hcyl(Pr, Ra), 0.1)

    def test_Nu_vcyl(self):
        Pr = 0.7
        Ra = 1e6
        thesis_data = {
            10: [
                (1e6, 23),
                (1e4, 11),
                # (1e3, 8.5),
                # (1e2, 7),
                (1e7, 38),
                (1e8, 60),
                (1e9, 100)
            ],
            5: [
                (1e4, 8.2),
                # (1e5, 10.3),
                (1e6, 20),
                (1e7, 35),
                (1e8, 60)
            ],
            1: [
                (1e4, 6),
                (1e5, 10),
                # (1e6, 10.8),
                (1e7, 30)
            ],
            0.1: [
                # (1e2, 1.3),
                # (1e3, 2.6),
                (1e4, 5)
            ]
        }
        for AR, plot_data in thesis_data.items():
            for Ra, Nu in plot_data:
                D = 1 * u.inch
                L = D * AR
                self.assertApproxEqual(Nu, ht.Nu_vcyl(Pr, Ra, D, L), uncertainty=0.4)


class Tubestest(unittest.TestCase):
    def test_tube(self):
        tubes = [
            (1*u.inch, 0*u.inch, 0*u.inch),
            (1*u.inch, 1*u.inch, 0*u.inch),
            (2*u.inch, 1*u.inch, 1*u.inch)
        ]
        for OD, wall, L in tubes:
            tube = ht.piping.Tube(OD, wall=wall, L=L)
            self.assertEqual(tube.OD, OD)
            self.assertEqual(tube.wall, wall)
            ID = OD-wall*2
            self.assertEqual(tube.ID, ID)
            A = pi * ID**2 / 4
            self.assertEqual(tube.area, A)
            V = A * L
            self.assertEqual(tube.volume, V)
            # Check for eps
            # Check for c

    def test_pipe(self):
        pipes = [
            (1, 10, 1.315*u.inch, 0.109*u.inch),
            (1, 40, 1.315*u.inch, 0.133*u.inch),
        ]
        # Check for proper error on non-existent schedules
        for D, SCH, OD, wall in pipes:
            pipe = ht.piping.Pipe(D, SCH=SCH)
            self.assertEqual(pipe.OD, OD)
            self.assertEqual(pipe.wall, wall)
            ID = OD-wall*2
            self.assertEqual(pipe.ID, ID)
            A = pi * ID**2 / 4
            self.assertEqual(pipe.area, A)

    def test_CopperTube(self):
        copper_tube = ht.piping.CopperTube(1)
        OD = 1.125 * u.inch
        wall = 0.065 * u.inch
        self.assertEqual(copper_tube.OD, OD)
        self.assertEqual(copper_tube.wall, wall)
        ID = OD-wall*2
        self.assertEqual(copper_tube.ID, ID)
        A = pi * ID**2 / 4
        self.assertEqual(copper_tube.area, A)

    def test_CorrugatedPipe(self):
        OD = 1 * u.inch
        corr_pipe = ht.piping.CorrugatedPipe(OD)
        self.assertEqual(corr_pipe.OD, OD)
        wall = 0 * u.inch
        self.assertEqual(corr_pipe.wall, wall)
        ID = OD-wall*2
        self.assertEqual(corr_pipe.ID, ID)
        A = pi * ID**2 / 4
        self.assertEqual(corr_pipe.area, A)

    def test_Elbow(self):
        OD = 1*u.inch
        wall = 0.049 * u.inch
        R_D = 1.5
        N = 1
        angle = 90 * u.degrees
        elbow = ht.piping.Elbow(OD, wall=wall, R_D=R_D, N=N, angle=angle)
        self.assertEqual(elbow.OD, OD)
        self.assertEqual(elbow.wall, wall)
        ID = OD-wall*2
        self.assertEqual(elbow.ID, ID)
        A = pi * ID**2 / 4
        self.assertEqual(elbow.area, A)
        L = R_D * ID * angle.m_as(u.degree)/180*pi
        self.assertEqual(elbow.L, L)
        V = A * L
        self.assertEqual(elbow.volume, V)

    def test_PipeElbow(self):
        D_nom = 1
        SCH = 10
        OD = 1.315 * u.inch
        wall = 0.109 * u.inch

        R_D = 1.5
        N = 1
        angle = 90 * u.degrees
        elbow = ht.piping.PipeElbow(D_nom, SCH, R_D=R_D, N=N, angle=angle)
        self.assertEqual(elbow.OD, OD)
        self.assertEqual(elbow.wall, wall)
        ID = OD-wall*2
        self.assertEqual(elbow.ID, ID)
        A = pi * ID**2 / 4
        self.assertEqual(elbow.area, A)
        L = R_D * ID * angle
        self.assertEqual(elbow.L, L)
        V = A * L
        self.assertEqual(elbow.volume, V)

    def test_Entrance(self):
        ID = 1 * u.inch
        K = 0.78
        entrance = ht.piping.Entrance(ID, K=K)
        self.assertEqual(entrance.ID, ID)
        A = pi * ID**2 / 4
        self.assertEqual(entrance.area, A)
        L = 0*u.inch
        self.assertEqual(entrance.L, L)
        V = A * L
        self.assertEqual(entrance.volume, V)
        self.assertEqual(entrance.K(), K)

    def test_Exit(self):
        ID = 1 * u.inch
        exit_ = ht.piping.Exit(ID)
        self.assertEqual(exit_.ID, ID)
        A = pi * ID**2 / 4
        self.assertEqual(exit_.area, A)
        L = 0*u.inch
        self.assertEqual(exit_.L, L)
        V = A * L
        self.assertEqual(exit_.volume, V)
        self.assertEqual(exit_.K(), 1)

    def test_Orifice(self):
        ID = 1 * u.inch
        orifice = ht.piping.Orifice(ID)
        self.assertEqual(orifice.ID, ID)
        A = pi * ID**2 / 4
        self.assertEqual(orifice.area, A)
        L = 0*u.inch
        self.assertEqual(orifice.L, L)
        V = A * L
        self.assertEqual(orifice.volume, V)
        self.assertEqual(orifice.K(), 1/.61**2)


class geometry(unittest.TestCase):
    def test_circle_area(self):
        D = 1 * u.inch
        A = pi * D**2 / 4
        self.assertEqual(A, ctb.geometry.circle_area(D))

    def test_cylinder_volume(self):
        D = 1 * u.inch
        H = 10 * u.inch
        A = pi * D**2 / 4
        V = A * H
        self.assertEqual(V, ctb.geometry.cylinder_volume(D, H))


@unittest.skipIf(skip_HP_test, 'Heprop failed to initialize, skipping Heprop tests.')
class Heprop(unittest.TestCase):
    def test_fluid_name(self):
        self.assertRaises(ValueError, ctb.ThermState, 'nitrogen', T= 300*u.K, P=1*u.bar, backend='HEPROP')

    def test_create(self):
        CP_helium_state = ctb.ThermState('helium', T= 300*u.K, P=1*u.bar)
        heprop_state = ctb.ThermState('helium', T= 300*u.K, P=1*u.bar, backend='HEPROP')
        self.assertAlmostEqual(CP_helium_state.T_critical / heprop_state.T_critical, 1, delta=0.01)
        self.assertAlmostEqual(CP_helium_state.P_critical / heprop_state.P_critical, 1, delta=0.01)
        self.assertAlmostEqual(CP_helium_state.Dmass / heprop_state.Dmass, 1, delta=0.01)
        self.assertAlmostEqual(CP_helium_state.Dmolar / heprop_state.Dmolar, 1, delta=0.01)
        CP_helium_state_2 = ctb.ThermState('helium', T= 100*u.K, P=1*u.bar)
        heprop_state_2 = ctb.ThermState('helium', T= 100*u.K, P=1*u.bar, backend='HEPROP')
        self.assertAlmostEqual((CP_helium_state.Smass - CP_helium_state_2.Smass) / (heprop_state.Smass - heprop_state_2.Smass), 1, delta=0.01)
        self.assertAlmostEqual((CP_helium_state.Smolar - CP_helium_state_2.Smolar)/ (heprop_state.Smolar - heprop_state_2.Smolar), 1, delta=0.01)
        self.assertAlmostEqual((CP_helium_state.Hmass - CP_helium_state_2.Hmass) / (heprop_state.Hmass - heprop_state_2.Hmass), 1, delta=0.01)
        self.assertAlmostEqual((CP_helium_state.Hmolar - CP_helium_state_2.Hmolar)/ (heprop_state.Hmolar - heprop_state_2.Hmolar), 1, delta=0.01)
        self.assertAlmostEqual(CP_helium_state.viscosity / heprop_state.viscosity, 1, delta=0.001)
        self.assertAlmostEqual(CP_helium_state.conductivity / heprop_state.conductivity, 1, delta=0.001)
        self.assertAlmostEqual(CP_helium_state.Prandtl / heprop_state.Prandtl, 1, delta=0.001)
        self.assertAlmostEqual(CP_helium_state.isothermal_compressibility / heprop_state.isothermal_compressibility, 1, delta=0.001)
        self.assertAlmostEqual(CP_helium_state.isobaric_expansion_coefficient / heprop_state.isobaric_expansion_coefficient, 1, delta=0.001)
        self.assertAlmostEqual(CP_helium_state.specific_heat_input / heprop_state.specific_heat_input, 1, delta=0.001)
        self.assertAlmostEqual(CP_helium_state.compressibility_factor / heprop_state.compressibility_factor, 1, delta=0.001)


class TestContraction(unittest.TestCase):
    """Unit tests for the Contraction class."""

    def test_contraction_basic_properties(self):
        """Test basic contraction creation and properties."""
        c = ctb.piping.Contraction(4*u.inch, 2*u.inch)

        self.assertEqual(c.beta, 0.5)
        self.assertEqual(c.theta, 180*u.deg)
        self.assertEqual(c.ID, 2*u.inch)
        self.assertEqual(c.ID1, 4*u.inch)
        self.assertEqual(c.ID2, 2*u.inch)
        self.assertEqual(c.type, 'Contraction')
        self.assertIsNone(c.OD)

    def test_contraction_gradual_30_degrees(self):
        """Test gradual contraction with 30° angle."""
        c = ctb.piping.Contraction(4*u.inch, 2*u.inch, 30*u.deg)
        K = c.K()

        # Expected K ≈ 0.8 * sin(15°) * (1-0.25) ≈ 0.155
        expected_K = 0.8 * sin(15*u.deg) * (1 - 0.5**2)
        self.assertAlmostEqual(K, expected_K, places=4)
        self.assertTrue(0.15 < K < 0.16)

    def test_contraction_gradual_45_degrees(self):
        """Test gradual contraction at 45° boundary."""
        c = ctb.piping.Contraction(4*u.inch, 2*u.inch, 45*u.deg)
        K = c.K()

        # Should use first formula: K = 0.8 * sin(22.5°) * (1-0.25)
        expected_K = 0.8 * sin(22.5*u.deg) * (1 - 0.5**2)
        self.assertAlmostEqual(K, expected_K, places=4)

    def test_contraction_sudden_180_degrees(self):
        """Test sudden contraction (180°)."""
        c = ctb.piping.Contraction(4*u.inch, 2*u.inch)
        K = c.K()

        # Expected K = 0.5 * (1-0.25) * sqrt(sin(90°)) = 0.375
        expected_K = 0.5 * (1 - 0.5**2) * sqrt(sin(90*u.deg))
        self.assertAlmostEqual(K, expected_K, places=6)
        self.assertAlmostEqual(K, 0.375, places=6)

    def test_contraction_90_degrees(self):
        """Test contraction with 90° angle."""
        c = ctb.piping.Contraction(4*u.inch, 2*u.inch, 90*u.deg)
        K = c.K()

        # Should use second formula: K = 0.5 * (1-0.25) * sqrt(sin(45°))
        expected_K = 0.5 * (1 - 0.5**2) * sqrt(sin(45*u.deg))
        self.assertAlmostEqual(K, expected_K, places=4)

    def test_contraction_different_diameter_ratio(self):
        """Test contraction with different diameter ratio."""
        c = ctb.piping.Contraction(6*u.inch, 3*u.inch, 60*u.deg)
        K = c.K()

        beta = 3.0/6.0  # 0.5
        expected_K = 0.5 * (1 - beta**2) * sqrt(sin(30*u.deg))
        self.assertAlmostEqual(K, expected_K, places=4)
        self.assertEqual(c.beta, 0.5)

    def test_contraction_small_diameter_ratio(self):
        """Test contraction with small diameter ratio (large contraction)."""
        c = ctb.piping.Contraction(10*u.inch, 2*u.inch, 120*u.deg)
        K = c.K()

        beta = 2.0/10.0  # 0.2
        expected_K = 0.5 * (1 - beta**2) * sqrt(sin(60*u.deg))
        self.assertAlmostEqual(K, expected_K, places=4)
        self.assertEqual(c.beta, 0.2)

    def test_contraction_validation_wrong_direction(self):
        """Test input validation - ID1 must be greater than ID2."""
        with self.assertRaises(ValueError) as context:
            ctb.piping.Contraction(2*u.inch, 4*u.inch)

        self.assertIn("ID1", str(context.exception))
        self.assertIn("greater than ID2", str(context.exception))

    def test_contraction_validation_equal_diameters(self):
        """Test input validation - equal diameters not allowed."""
        with self.assertRaises(ValueError) as context:
            ctb.piping.Contraction(4*u.inch, 4*u.inch)

        self.assertIn("greater than ID2", str(context.exception))

    def test_contraction_validation_zero_angle(self):
        """Test input validation - zero angle not allowed."""
        with self.assertRaises(ValueError) as context:
            ctb.piping.Contraction(4*u.inch, 2*u.inch, 0*u.deg)

        self.assertIn("angle must be in range", str(context.exception))

    def test_contraction_validation_negative_angle(self):
        """Test input validation - negative angle not allowed."""
        with self.assertRaises(ValueError) as context:
            ctb.piping.Contraction(4*u.inch, 2*u.inch, -30*u.deg)

        self.assertIn("angle must be in range", str(context.exception))

    def test_contraction_validation_excessive_angle(self):
        """Test input validation - angle > 180° not allowed."""
        with self.assertRaises(ValueError) as context:
            ctb.piping.Contraction(4*u.inch, 2*u.inch, 200*u.deg)

        self.assertIn("angle must be in range", str(context.exception))

    def test_contraction_K_with_reynolds_warning(self):
        """Test K calculation with Reynolds number (should warn)."""
        c = ctb.piping.Contraction(4*u.inch, 2*u.inch)

        # Should not raise exception but may warn
        K = c.K(Re=10000)
        self.assertIsInstance(K, float)
        self.assertGreater(K, 0)

    def test_contraction_length_calculation(self):
        """Test axial length calculation."""
        c = ctb.piping.Contraction(4*u.inch, 2*u.inch, 60*u.deg)

        # L = |ID1 - ID2| / (2 * tan(theta/2))
        expected_L = abs(4 - 2) / (2 * tan(30*u.deg))
        expected_L_with_units = expected_L * u.inch

        self.assertAlmostEqual(c.L.to(u.inch).magnitude, expected_L, places=4)

    def test_contraction_volume_calculation(self):
        """Test internal volume calculation."""
        c = ctb.piping.Contraction(4*u.inch, 2*u.inch, 90*u.deg)

        # Volume should be positive and reasonable for a truncated cone
        self.assertGreater(c.volume.to(u.inch**3).magnitude, 0)

        # Volume of truncated cone: V = π*L/3 * (r1² + r1*r2 + r2²)
        r1 = 2.0  # radius in inches
        r2 = 1.0  # radius in inches
        L = c.L.to(u.inch).magnitude
        expected_volume = pi * L / 3 * (r1**2 + r1*r2 + r2**2)

        self.assertAlmostEqual(c.volume.to(u.inch**3).magnitude, expected_volume, places=4)

    def test_contraction_area_property(self):
        """Test area property returns downstream area."""
        c = ctb.piping.Contraction(4*u.inch, 2*u.inch)

        # Area should be based on downstream (smaller) diameter
        expected_area = pi * (1.0)**2 * u.inch**2  # radius = 1 inch
        self.assertAlmostEqual(c.area.to(u.inch**2).magnitude,
                             expected_area.to(u.inch**2).magnitude, places=4)

    def test_contraction_str_representation(self):
        """Test string representation."""
        c = ctb.piping.Contraction(4*u.inch, 2*u.inch, 60*u.deg)
        str_repr = str(c)

        self.assertIn("Contraction", str_repr)
        self.assertIn("60", str_repr)
        self.assertIn("4.000", str_repr)
        self.assertIn("2.000", str_repr)

    def test_contraction_repr_representation(self):
        """Test detailed representation."""
        c = ctb.piping.Contraction(4*u.inch, 2*u.inch, 60*u.deg)
        repr_str = repr(c)

        self.assertIn("Contraction(", repr_str)
        self.assertIn("ID1=", repr_str)
        self.assertIn("ID2=", repr_str)
        self.assertIn("theta=", repr_str)

    def test_contraction_edge_case_very_small_angle(self):
        """Test contraction with very small angle."""
        c = ctb.piping.Contraction(4*u.inch, 2*u.inch, 1*u.deg)
        K = c.K()

        # Should use first formula and give very small K value
        self.assertGreater(K, 0)
        self.assertLess(K, 0.01)  # Very gradual contraction

    def test_contraction_edge_case_nearly_sudden(self):
        """Test contraction with angle just under 180°."""
        c = ctb.piping.Contraction(4*u.inch, 2*u.inch, 179*u.deg)
        K = c.K()

        # Should be very close to sudden contraction value
        sudden_c = ctb.piping.Contraction(4*u.inch, 2*u.inch, 180*u.deg)
        sudden_K = sudden_c.K()

        self.assertAlmostEqual(K, sudden_K, places=3)

class TestDPDarcy(unittest.TestCase):
    """Unit tests for dP_Darcy function."""

    def test_basic_calculation(self):
        """Test basic pressure drop calculation."""
        K = 0.9  # 90° elbow
        rho = 1000 * u.kg / u.m**3  # Water density
        w = 2.0 * u.m / u.s  # Flow velocity

        dp = ctb.piping.dP_Darcy(K, rho, w)

        # Expected: K * rho * w^2 / 2 = 0.9 * 1000 * 4 / 2 = 1800 Pa
        expected_pa = 1800 * u.Pa

        self.assertAlmostEqual(dp.m_as(u.Pa), expected_pa.m_as(u.Pa), places=5)

    def test_zero_velocity(self):
        """Test with zero velocity should give zero pressure drop."""
        dp = ctb.piping.dP_Darcy(0.9, 1000*u.kg/u.m**3, 0*u.m/u.s)
        self.assertEqual(dp.magnitude, 0.0)

    def test_zero_K_coefficient(self):
        """Test with zero K coefficient should give zero pressure drop."""
        dp = ctb.piping.dP_Darcy(0.0, 1000*u.kg/u.m**3, 2*u.m/u.s)
        self.assertEqual(dp.magnitude, 0.0)

    def test_liquid_nitrogen_example(self):
        """Test with liquid nitrogen properties."""
        K = 0.375  # Sudden contraction 4" to 2"
        rho = 808 * u.kg / u.m**3  # Liquid nitrogen
        w = 1.5 * u.m / u.s

        dp = ctb.piping.dP_Darcy(K, rho, w)

        self.assertAlmostEqual(dp.m_as(u.psi), 0.0494397, places=5)

    def test_high_velocity_flow(self):
        """Test with high velocity flow (quadratic relationship)."""
        K = 1.0
        rho = 1000 * u.kg / u.m**3
        w1 = 1 * u.m / u.s
        w2 = 2 * u.m / u.s  # Double velocity

        dp1 = ctb.piping.dP_Darcy(K, rho, w1)
        dp2 = ctb.piping.dP_Darcy(K, rho, w2)

        # Pressure drop should be 4x for 2x velocity (quadratic)
        ratio = float(dp2 / dp1)
        self.assertAlmostEqual(ratio, 4.0, places=2)

    def test_validation_negative_K(self):
        """Test validation of negative K coefficient."""
        with self.assertRaises(ValueError) as context:
            ctb.piping.dP_Darcy(-0.5, 1000*u.kg/u.m**3, 2*u.m/u.s)

        self.assertIn("non-negative", str(context.exception))

    def test_validation_negative_density(self):
        """Test validation of negative density."""
        with self.assertRaises(ValueError) as context:
            ctb.piping.dP_Darcy(0.9, -1000*u.kg/u.m**3, 2*u.m/u.s)

        self.assertIn("positive", str(context.exception))

    def test_validation_zero_density(self):
        """Test validation of zero density."""
        with self.assertRaises(ValueError) as context:
            ctb.piping.dP_Darcy(0.9, 0*u.kg/u.m**3, 2*u.m/u.s)

        self.assertIn("positive", str(context.exception))

    def test_validation_negative_velocity(self):
        """Test validation of negative velocity."""
        with self.assertRaises(ValueError) as context:
            ctb.piping.dP_Darcy(0.9, 1000*u.kg/u.m**3, -2*u.m/u.s)

        self.assertIn("non-negative", str(context.exception))

    def test_multiple_fittings_additive(self):
        """Test that K values are additive for multiple fittings."""
        # Individual fittings
        K1, K2, K3 = 0.5, 0.9, 0.15  # entrance, elbow, valve
        rho = 1000 * u.kg / u.m**3
        w = 2 * u.m / u.s

        dp1 = ctb.piping.dP_Darcy(K1, rho, w)
        dp2 = ctb.piping.dP_Darcy(K2, rho, w)
        dp3 = ctb.piping.dP_Darcy(K3, rho, w)

        # Combined fitting
        K_total = K1 + K2 + K3
        dp_total = ctb.piping.dP_Darcy(K_total, rho, w)

        # Total should equal sum of individuals
        dp_sum = dp1 + dp2 + dp3
        self.assertAlmostEqual(dp_total.magnitude, dp_sum.magnitude, places=6)

if __name__ == '__main__':
    unittest.main()
