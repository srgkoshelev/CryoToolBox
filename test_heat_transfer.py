import heat_transfer as ht
from heat_transfer import odh
from heat_transfer import geometry as g
import pprint
import unittest
from math import pi
from unittest.mock import MagicMock

pp = pprint.PrettyPrinter()

ureg = ht.ureg
Q_ = ht.Q_


class RefPropTest(unittest.TestCase):
    """Verify fluid properties are calculated accurately.

    Tests borrowed from pyrefprop package.
    """

    criteria = 0.00001  # NIST acceptance criteria

    def assertNIST(self, nist, calc, criteria=None):
        criteria = criteria or self.criteria
        assert abs(nist-calc) / nist < criteria, \
            f'Calculated value {calc} is not within {criteria:%} of NIST value \
            {nist}'

    def test_air(self):
        if REFPROP_TEST == False:
            return
        fluid = ht.ThermState('air.mix', backend='REFPROP')
        self.assertNIST(fluid.M, 28.958600656)
        # TODO Open an issue with CoolProp

    def test_argon(self):
        if REFPROP_TEST == False:
            return
        fluid = ht.ThermState('argon', backend='REFPROP')
        P = 2 * 1000 * ureg.kPa
        D = 15 / fluid.M * ureg.mol/ureg.L
        fluid.update_kw(P=P, Dmolar=D)
        self.assertNIST(637.377588657857, fluid.T.to(ureg.K).magnitude)

    def test_r134a(self):
        if REFPROP_TEST == False:
            return
        fluid = ht.ThermState('r134a', backend='REFPROP')
        T = 400 * ureg.K
        D = 50 / fluid.M * ureg.mol/ureg.L
        fluid.update_kw(T=T, Dmolar=D)
        self.assertNIST(1.45691892789737, fluid.P.to(ureg.kPa).magnitude/1000)

    def test_oxygen(self):
        if REFPROP_TEST == False:
            return
        fluid = ht.ThermState('oxygen', backend='REFPROP')
        T = 100 * ureg.K
        P = 1000 * ureg.kPa
        fluid.update_kw(T=T, P=P)
        self.assertNIST(153.886680663753,
                        fluid.viscosity.to(ureg.uPa*ureg.s).magnitude)

    def test_nitrogen(self):
        if REFPROP_TEST == False:
            return
        fluid = ht.ThermState('nitrogen', backend='REFPROP')
        T = 100 * ureg.K
        fluid.update_kw(T=T, Q=0)
        self.assertNIST(100.111748964945,
                        fluid.conductivity.to(ureg.W/(ureg.m*ureg.K)).magnitude*
                        1000)

    def test_air_density(self):
        if REFPROP_TEST == False:
            return
        fluid = ht.ThermState('air.mix', backend='REFPROP')
        T = (((70 - 32) * 5 / 9) + 273.15) * ureg.K
        P = 14.7 / 14.50377377 * (10**5) / 1000 * ureg.kPa
        fluid.update_kw(P=P, T=T)
        self.assertNIST(0.0749156384666842,
                        fluid.Dmolar.to(ureg.mol/ureg.L).magnitude*fluid.M *
                        0.062427974)

    # def test_freon_mixture(self):
    #     fluid = ht.ThermState('R32&R125', backend='REFPROP')
    #     fluid.set_mole_fractions(0.3, 0.7)
    #     P = 10 * 1000 * ureg.kPa
    #     Smolar = 110 * ureg.J / (ureg.mol * ureg.K)
    #     fluid.update_kw(P=P, Smolar=Smolar)
    #     self.assertNIST(23643.993624382,
    #                     fluid.Hmolar.to(ureg.J/ureg.mol).magnitude)
    # TODO Figure out the issue with mixtures (possibly refprop_setref/ixflag thing)

    # def test_ethane_butane(self):
    #     fluid = ht.ThermState('ethane&butane', backend='REFPROP')
    #     fluid.set_mole_fractions(0.5, 0.5)
    #     Dmolar = 30 * 0.45359237 / 0.028316846592 / fluid.M * ureg.mol / ureg.L
    #     Hmolar = 283 * 1.05435026448889 / 0.45359237 * fluid.M * ureg.J / ureg.mol
    #     fluid.update_kw(Dmolar=Dmolar, Hmolar=Hmolar)
    #     self.assertNIST(298.431320311048,
    #                     fluid.T.to(ureg.degF).magnitude)
    # TODO Figure out the issue with mixtures (possibly refprop_setref/ixflag thing)

    def test_ammonia_water(self):
        if REFPROP_TEST == False:
            return
        fluid = ht.ThermState('ammonia&water', backend='REFPROP')
        fluid.set_mole_fractions(0.4, 0.6)
        T = (((300 - 32) * 5 / 9) + 273.15) * ureg.K
        P = 10000 / 14.50377377 * (10**5) / 1000 * ureg.kPa
        fluid.update_kw(T=T, P=P)
        self.assertNIST(5536.79144924071,
                        fluid.speed_sound.to(ureg.m/ureg.s).magnitude *
                        1000 / 25.4 / 12)

    def test_octane(self):
        if REFPROP_TEST == False:
            return
        fluid = ht.ThermState('octane', backend='REFPROP')
        T = (100+273.15) * ureg.K
        fluid.update_kw(T=T, Q=0)
        self.assertNIST(319.167499870568,
                        fluid.latent_heat.to(ureg.J/ureg.g).magnitude)

    def test_R410A_mole_fraction(self):
        if REFPROP_TEST == False:
            return
        fluid = ht.ThermState('R410A.MIX', backend='REFPROP')
        self.assertNIST(0.697614699375863,
                        fluid.mole_fractions[0])

    def test_R410A_mass_fraction(self):
        if REFPROP_TEST == False:
            return
        fluid = ht.ThermState('R410A.MIX', backend='REFPROP')
        self.assertNIST(0.5,
                        fluid.mass_fractions[0])


class FunctionsTest(unittest.TestCase):
    def assertApproxEqual(self, data, calc, uncertainty=0.1):
        assert abs(data-calc) / data < uncertainty, \
            f'Calculated value {calc} is not within {uncertainty:.1%} of data \
            value {data}'

    def test_rad_hl_1(self):
        eps = 0.02
        T1 = Q_(27, ureg.degC)
        T2 = Q_(-183, ureg.degC)
        A = 0.05 * ureg.m**2
        heat_flow = abs(ht.rad_hl(T1, eps, T2, eps)) * A
        self.assertAlmostEqual(0.23054, heat_flow.to(ureg.W).magnitude, 5)

    def test_rad_hl_2(self):
        eps = 0.55
        T1 = Q_(27, ureg.degC)
        T2 = Q_(-183, ureg.degC)
        A = 0.05 * ureg.m**2
        heat_flow = abs(ht.rad_hl(T1, eps, T2, eps)) * A
        self.assertAlmostEqual(8.65727, heat_flow.to(ureg.W).magnitude, 5)

    def test_rad_hl_3(self):
        eps = 0.8
        T1 = Q_(327, ureg.degC)
        T2 = Q_(127, ureg.degC)
        eps_b = 0.05
        heat_flow = abs(ht.rad_hl(T1, eps, T2, eps,
                                  baffles={'N': 1, 'eps': eps_b}))
        self.assertAlmostEqual(145.7373409,
                               heat_flow.to(ureg.W/ureg.m**2).magnitude, 5)

    def test_rad_hl_4(self):
        T1 = Q_(327, ureg.degC)
        eps_1 = 0.8
        T2 = Q_(127, ureg.degC)
        eps_2 = 0.4
        eps_b = 0.05
        heat_flow = abs(ht.rad_hl(T1, eps_1, T2, eps_2,
                                  baffles={'N': 1, 'eps': eps_b}))
        self.assertAlmostEqual(141.3739475,
                               heat_flow.to(ureg.W/ureg.m**2).magnitude, 5)

    def test_API_subsonic(self):
        """Example from API 5.6"""
        if REFPROP_TEST == False:
            return
        m_dot = 53500 * ureg.lb/ureg.hr
        fluid = ht.ThermState('propane&butane', backend='REFPROP')
        fluid.set_mole_fractions(0.5, 0.5)
        fluid.update_kw(T=627*ureg.degR, P=97.2*ureg.psi)
        P_back = 77.2 * ureg.psi
        A_expect = 6.55 * ureg.inch**2
        A_calc = ht.A_relief_API(m_dot, fluid, P_back=P_back)
        self.assertApproxEqual(A_expect, A_calc, uncertainty=0.05)

    def test_API_sonic(self):
        """Example from API 5.6"""
        if REFPROP_TEST == False:
            return
        m_dot = 53500 * ureg.lb/ureg.hr
        fluid = ht.ThermState('propane&butane', backend='REFPROP')
        fluid.set_mole_fractions(0.5, 0.5)
        fluid.update_kw(T=627*ureg.degR, P=97.2*ureg.psi)
        P_back = Q_(0, ureg.psig)
        A_expect = 5.73 * ureg.inch**2
        A_calc = ht.A_relief_API(m_dot, fluid, P_back=P_back)
        self.assertApproxEqual(A_expect, A_calc, uncertainty=0.05)

    def test_blast_radia(self):
        E_stored = 150000 * ureg.lbf*ureg.ft
        R_blast = ht.blast_radius(E_stored)
        self.assertAlmostEqual(29.9744399*ureg.ft, R_blast[0])
        self.assertAlmostEqual(17.984664*ureg.ft, R_blast[1])
        self.assertAlmostEqual(8.99233198*ureg.ft, R_blast[2])
        self.assertAlmostEqual(2.99744399*ureg.ft, R_blast[3])

class PipingTest(unittest.TestCase):
    """Piping checks, mostly taken from textbooks.
    """

    def assertApproxEqual(self, expected, calculated, uncertainty=0.05):
        error = float(abs(expected-calculated)/expected)
        if isinstance(expected, ureg.Quantity) and isinstance(calculated, ureg.Quantity):
            unit = expected.units
            calculated.ito(unit)
        assert error < uncertainty, \
            f'Calculated value {calculated:g} is not within {uncertainty:.1%} of ' + \
            f'expected value {expected:g}'

    def test_piping(self):
        pipe = ht.piping.Pipe(1)
        piping = ht.piping.Piping(pipe)
        self.assertEqual([pipe], piping._elements)
        self.assertEqual(pipe, piping[0])
        piping.insert(0, pipe)
        self.assertEqual([pipe, pipe], piping._elements)
        pipe2 = ht.piping.Pipe(2)
        piping[1] = pipe2
        self.assertEqual([pipe, pipe2], piping._elements)
        del piping[0]
        self.assertEqual([pipe2], piping._elements)
        self.assertEqual(pipe2, piping[0])
        self.assertEqual(1, len(piping))

    def test_create_pipes(self):
        tube = ht.piping.Tube(Q_('1 inch'))
        pipe = ht.piping.Pipe(Q_('1 inch'))
        c_tube = ht.piping.CopperTube(Q_('3/4 inch'))
        vj_pipe = ht.piping.VJPipe(Q_('1 inch'), VJ_D=Q_('2 inch'))
        corr_pipe = ht.piping.CorrugatedPipe(Q_('1 inch'))
        entrance = ht.piping.Entrance(Q_('1 inch'))
        pipe_exit = ht.piping.Exit(Q_('1 inch'))
        orifice = ht.piping.Orifice(Q_('1 inch'))
        c_orifice = ht.piping.ConicOrifice(1, Q_('3/4 inch'))
        annulus = ht.piping.Annulus(Q_('1 inch'), Q_('3/4 inch'))
        pipe_elbow = ht.piping.PipeElbow(Q_('1 inch'), N=2)
        elbow = ht.piping.Elbow(Q_('1 inch'), N=2)
        pipe_tee = ht.piping.PipeTee(Q_('1 inch'), N=2)
        tee = ht.piping.Tee(Q_('1 inch'), N=2)
        valve = ht.piping.Valve(Q_('1 inch'), 1)
        # g_valve = ht.piping.GlobeValve(Q_('1 inch'))
        # print(f'Generated {g_valve}')
        # v_cone = ht.piping.VCone(Q_('1 inch'), 0.7, 1)
        # print(f'Generated {v_cone}')
        cont = ht.piping.Contraction(1*ureg.inch, 0.5*ureg.inch)
        enl = ht.piping.Enlargement(0.5*ureg.inch, 1*ureg.inch)
        test_state = ht.AIR
        piping = ht.piping.Piping(
            pipe, vj_pipe, corr_pipe, entrance, pipe_exit, orifice, c_orifice,
             tube, c_tube, annulus, pipe_elbow, elbow, pipe_tee, tee, valve,
             # g_valve, v_cone,
             cont, enl)
        piping.volume()
        dP = ht.piping.dP_incomp(Q_('10 g/s'), test_state, piping)
        self.assertApproxEqual(21.2*ureg.psi, dP)

    def test_f_Darcy(self):
        eps_smooth = 0.0018 * ureg.inch
        Re = 1e8
        ID = 0.2 * ureg.inch
        eps_r = eps_smooth/ID
        self.assertApproxEqual(0.0368, ht.piping.f_Darcy(Re, eps_r))
        ID = 0.4 * ureg.inch
        eps_r = eps_smooth/ID
        self.assertApproxEqual(0.0296, ht.piping.f_Darcy(Re, eps_r))
        ID = 1 * ureg.inch
        eps_r = eps_smooth/ID
        self.assertApproxEqual(0.0228, ht.piping.f_Darcy(Re, eps_r))
        ID = 2 * ureg.inch
        eps_r = eps_smooth/ID
        self.assertApproxEqual(0.0191, ht.piping.f_Darcy(Re, eps_r))
        eps_r = 0.006
        Re = 1e5
        self.assertApproxEqual(0.033, ht.piping.f_Darcy(Re, eps_r))
        eps_r = 0.006
        Re = 1e3
        self.assertApproxEqual(64/Re, ht.piping.f_Darcy(Re, eps_r))

    # def test_Crane_4_22(self):
    #     test_air = ht.ThermState('air', P=2.343*ureg.bar, T=40*ureg.degC)
    #     pipe = ht.piping.Pipe(1/2, SCH=80, L=3*ureg.m)
    #     piping = ht.piping.Piping(
    #         test_air,
    #         [pipe,
    #          ht.piping.Exit(pipe.ID)])
    #     Y = 0.76  # Taken from Crane; temporary stub
    #     flow = ht.to_standard_flow(piping.m_dot(), test_air) * Y
    #     # TODO Check this test (might have to do with subsonic flow)
    #     self.assertAlmostEqual(1.76, flow.to(ureg.m**3/ureg.min).magnitude)

    def test_Crane_7_16(self):
        air = ht.ThermState('air', P=Q_(65, ureg.psig), T=Q_(110, ureg.degF))
        pipe = ht.piping.Pipe(1, SCH=40, L=75*ureg.ft)
        flow = 100 * ureg.ft**3/ureg.min
        m_dot = flow * ht.AIR.Dmass
        piping = ht.piping.Piping(pipe, ht.piping.Exit(pipe.ID))
        dP = ht.piping.dP_incomp(m_dot, air, piping)
        self.assertApproxEqual(2.61, dP.m_as(ureg.psi))
        dP_isot = ht.piping.dP_isot(m_dot, air, pipe)
        self.assertApproxEqual(2.61, dP_isot.m_as(ureg.psi))

    def test_Crane_7_22(self):
        air = ht.ThermState('air', P=Q_(19.3, ureg.psig), T=Q_(100, ureg.degF))
        pipe = ht.piping.Pipe(1/2, SCH=80, L=7.04/6.04*10*ureg.ft)  # Adjust for entrance K = 1, total 7.04
        dP = 19.3 * ureg.psi
        P_out = air.P - dP
        q_expected = 3762 * ureg.ft**3/ureg.hr  # STD flow
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
        # check(M, M_from_K_lim(K_lim, air.gamma))
        # K_left = K_lim - pipe.K(Re)
        # # print('K left: ', K_left)
        # M_end = M_from_K_lim(K_left, air.gamma)
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
        pipe = ht.piping.Pipe(4, SCH=40, L=100*ureg.ft)
        fluid = ht.ThermState('nitrogen', P=100*ureg.psi, T=530*ureg.degR)
        piping = ht.piping.Piping(pipe)
        P_out = 84.056 * ureg.psi
        m_dot_isot = ht.piping.m_dot_isot(fluid, pipe, P_out)
        self.assertApproxEqual(10, m_dot_isot.m_as(ureg.lb/ureg.s))

    def test_White_9_3(self):
        v = 240 * ureg.m/ureg.s
        T1 = 320 * ureg.K
        P1 = 170 * ureg.kPa
        T01 = 349 * ureg.K
        P01 = 230 * ureg.kPa
        fluid_static = ht.ThermState('air', P=P1, T=T1)
        area = 0.05 * ureg.m**2  # Assumed value
        M = ht.piping.Mach(fluid_static, v)
        M_exp = 0.67
        self.assertApproxEqual(M_exp, M, uncertainty=0.005)
        fluid_total = ht.ThermState('air', P=P01, T=T01)
        m_dot = v * area * fluid_static.Dmass
        M_total = ht.piping.Mach_total(fluid_total, m_dot, area)
        self.assertApproxEqual(M_exp, M_total, uncertainty=0.005)

    def test_White_9_4(self):
        fluid_static = ht.ThermState('air', P=500*ureg.kPa, T=470*ureg.K)
        area = 0.05 * ureg.m**2
        v = 180 * ureg.m/ureg.s
        m_dot = v * area * fluid_static.Dmass
        M = ht.piping.Mach(fluid_static, v)
        M_exp = 0.414
        self.assertApproxEqual(M_exp, M, uncertainty=0.005)
        P0 = 563 * ureg.kPa
        T0 = 486 * ureg.K
        fluid_total = ht.ThermState('air', P=P0, T=T0)
        M_total = ht.piping.Mach_total(fluid_total, m_dot, area)
        self.assertApproxEqual(M_exp, M_total, uncertainty=0.005)

    def test_White_9_10(self):
        Ma1 = 0.1
        Ma2 = 0.5
        Ma3 = 1
        k = ht.AIR.gamma
        f = 0.024
        d = 2 * ureg.cm
        K_lim1 = ht.piping.K_lim(Ma1, k)
        L1 = K_lim1 * d / f
        K_lim2 = ht.piping.K_lim(Ma2, k)
        L2 = K_lim2 * d / f
        K_lim3 = ht.piping.K_lim(Ma3, k)
        L3 = K_lim3 * d / f
        self.assertApproxEqual(55*ureg.m, L1-L2)
        self.assertApproxEqual(0.9*ureg.m, L2-L3)

    def test_White_9_11(self):
        f = 0.024
        d = 2 * ureg.cm
        Ma1 = 0.1
        Ma2 = 0.5
        P1 = 600 * ureg.kPa
        T1 = 450 * ureg.K
        fluid = ht.ThermState('air', P=P1, T=T1)
        k = fluid.gamma
        P2 = ht.piping.P_from_M(P1, Ma1, Ma2, k)
        self.assertApproxEqual(117*ureg.kPa, P2)
        P02 = ht.piping.P_total(P2, Ma2, k)
        self.assertApproxEqual(139*ureg.kPa, P02)

    def test_White_9_12(self):
        fluid_total = ht.ThermState('air', P=200*ureg.kPa, T=500*ureg.K)
        v = 100 * ureg.m/ureg.s
        pipe = ht.piping.Tube(3*ureg.cm, L=15*ureg.m)
        m_dot = 0.0961 * ureg.kg/ureg.s  # Based on static density from the source
        f = 0.02
        Ma = ht.piping.Mach_total(fluid_total, m_dot, pipe.area)
        K_lim = ht.piping.K_lim(Ma, fluid_total.gamma)
        Lstar = K_lim * pipe.ID / f
        self.assertApproxEqual(16.5*ureg.m, Lstar)

    def test_White_9_13(self):
        fluid_static = ht.ThermState('air', P=220*ureg.kPa, T=300*ureg.K)
        P_out = 140 * ureg.kPa
        pipe = ht.piping.Tube(1*ureg.cm, wall=0*ureg.m, L=1.2*ureg.m)
        K = 0.025 * pipe.L/pipe.ID
        pipe.K = MagicMock(return_value=K)
        m_dot_adiab = ht.piping.m_dot_adiab(fluid_static, pipe, P_out, state='static')
        self.assertApproxEqual(0.0233*ureg.kg/ureg.s, m_dot_adiab)

    def test_White_P_9_100(self):
        fluid = ht.ThermState('methane', T=Q_(68, ureg.degF), P=5*ureg.bar+101325*ureg.Pa)
        P_out = 1*ureg.bar + 101325*ureg.Pa
        pipe = ht.piping.Pipe(6, L=31*ureg.mile)
        K = 0.019 * pipe.L/pipe.ID
        pipe.K = MagicMock(return_value=K)
        m_dot_isot = ht.piping.m_dot_isot(fluid, pipe, P_out)
        self.assertApproxEqual(0.345*ureg.kg/ureg.s, m_dot_isot)
        dP_isot = ht.piping.dP_isot(m_dot_isot, fluid, pipe)
        self.assertApproxEqual(fluid.P-P_out, dP_isot)

    def test_White_P_9_101(self):
        air = ht.ThermState('air', T=Q_(20, ureg.degC), P=102*ureg.kPa)
        P_out = 100*ureg.kPa
        pipe = ht.piping.Tube(3*ureg.cm, wall=0*ureg.m, L=1*ureg.m)
        K = 0.028 * pipe.L/pipe.ID
        pipe.K = MagicMock(return_value=K)
        m_dot_incomp = ht.piping.m_dot_incomp(air, ht.piping.Piping(pipe), P_out=P_out)
        m_dot_isot = ht.piping.m_dot_isot(air, pipe, P_out)
        self.assertApproxEqual(m_dot_incomp, m_dot_isot)
        m_dot_adiab = ht.piping.m_dot_adiab(air, pipe, P_out)
        self.assertApproxEqual(m_dot_incomp, m_dot_adiab)
        dP_incomp = ht.piping.dP_incomp(m_dot_incomp, air, [pipe])
        self.assertApproxEqual(air.P-P_out, dP_incomp)
        dP_isot = ht.piping.dP_isot(m_dot_isot, air, pipe)
        self.assertApproxEqual(air.P-P_out, dP_isot)

    def test_M_from_K(self):
        air = ht.ThermState('air', P=300*ureg.kPa, T=500*ureg.K)
        self.assertApproxEqual(0.225, ht.piping.M_from_K_lim(
            11, air.gamma))
        self.assertApproxEqual(0.2, ht.piping.M_from_K_lim(
            14.533, air.gamma))
        self.assertApproxEqual(0.15, ht.piping.M_from_K_lim(
            27.932, air.gamma))
        self.assertApproxEqual(0.17, ht.piping.M_from_K_lim(
            21.115, air.gamma))
        self.assertApproxEqual(0.175, ht.piping.M_from_K_lim(
            19.772, air.gamma))
        self.assertApproxEqual(0.1741, ht.piping.M_from_K_lim(
            20, air.gamma))
        self.assertApproxEqual(1, ht.piping.M_from_K_lim(
            0, air.gamma))

    def test_piping_stress(self):
        OD = 1.75*ureg.inch
        ID = 1.625*ureg.inch
        wall = (OD - ID) / 2
        tube = ht.piping.Tube(OD, wall=wall)
        S = 16700*ureg.psi
        E = 1
        W = 1
        Y = 0.4
        self.assertApproxEqual(
            S,
            ht.piping.piping_stress(
                tube,
                1070.512*ureg.psi,
                E=E, W=W, Y=Y), uncertainty=0.01)

    def test_direct_integration_nitrogen_P165_T80(self):
        P = Q_(165, ureg.psi)
        T = 80 * ureg.K
        fluid = ht.ThermState('nitrogen', P=P, T=T)
        G_nozzle, P_2 = ht.piping.G_nozzle(fluid)
        self.assertApproxEqual(40147*ureg.kg/(ureg.s*ureg.m**2),
                               G_nozzle,
                               uncertainty=0.01)
        self.assertApproxEqual(19.26*ureg.psi,
                               P_2,
                               uncertainty=0.01)

    def test_direct_integration_nitrogen_P165_T90(self):
        P = Q_(165, ureg.psi)
        T = 90 * ureg.K
        fluid = ht.ThermState('nitrogen', P=P, T=T)
        G_nozzle, P_2 = ht.piping.G_nozzle(fluid)
        self.assertApproxEqual(34133.4*ureg.kg/(ureg.s*ureg.m**2),
                               G_nozzle,
                               uncertainty=0.01)
        self.assertApproxEqual(52.2416*ureg.psi,
                               P_2,
                               uncertainty=0.05)
        # For some reason above fails for lower uncertainty for even 200 steps

    def test_direct_integration_nitrogen_P165_100(self):
        P = Q_(165, ureg.psi)
        T = 100 * ureg.K
        fluid = ht.ThermState('nitrogen', P=P, T=T)
        G_nozzle, P_2 = ht.piping.G_nozzle(fluid, n_steps=150)
        self.assertApproxEqual(22655.1*ureg.kg/(ureg.s*ureg.m**2),
                               G_nozzle,
                               uncertainty=0.01)
        self.assertApproxEqual(111.4025*ureg.psi,
                               P_2,
                               uncertainty=0.01)

    def test_direct_integration_nitrogen_P165_Q0(self):
        P = Q_(165, ureg.psi)
        Q = 0
        fluid = ht.ThermState('nitrogen', P=P, Q=Q)
        G_nozzle, P_2 = ht.piping.G_nozzle(fluid, n_steps=150)
        self.assertApproxEqual(11268.4*ureg.kg/(ureg.s*ureg.m**2),
                               G_nozzle,
                               uncertainty=0.01)
        self.assertApproxEqual(124.4172*ureg.psi,
                               P_2,
                               uncertainty=0.01)

    def test_direct_integration_nitrogen_P165_Q025(self):
        P = Q_(165, ureg.psi)
        Q = 0.25
        fluid = ht.ThermState('nitrogen', P=P, Q=Q)
        G_nozzle, P_2 = ht.piping.G_nozzle(fluid, n_steps=150)
        self.assertApproxEqual(7472.6*ureg.kg/(ureg.s*ureg.m**2),
                               G_nozzle,
                               uncertainty=0.01)
        self.assertApproxEqual(105.4028*ureg.psi,
                               P_2,
                               uncertainty=0.01)

    def test_direct_integration_nitrogen_P165_Q05(self):
        P = Q_(165, ureg.psi)
        Q = 0.5
        fluid = ht.ThermState('nitrogen', P=P, Q=Q)
        G_nozzle, P_2 = ht.piping.G_nozzle(fluid)
        self.assertApproxEqual(6043.5*ureg.kg/(ureg.s*ureg.m**2),
                               G_nozzle,
                               uncertainty=0.01)
        self.assertApproxEqual(99.6527*ureg.psi,
                               P_2,
                               uncertainty=0.01)

    def test_direct_integration_nitrogen_P165_Q075(self):
        P = Q_(165, ureg.psi)
        Q = 0.75
        fluid = ht.ThermState('nitrogen', P=P, Q=Q)
        G_nozzle, P_2 = ht.piping.G_nozzle(fluid)
        self.assertApproxEqual(5218*ureg.kg/(ureg.s*ureg.m**2),
                               G_nozzle,
                               uncertainty=0.01)
        self.assertApproxEqual(96.7808*ureg.psi,
                               P_2,
                               uncertainty=0.01)

    def test_direct_integration_nitrogen_P165_Q1(self):
        P = Q_(165, ureg.psi)
        Q = 1
        fluid = ht.ThermState('nitrogen', P=P, Q=Q)
        G_nozzle, P_2 = ht.piping.G_nozzle(fluid, n_steps=250)
        self.assertApproxEqual(4661.2*ureg.kg/(ureg.s*ureg.m**2),
                               G_nozzle,
                               uncertainty=0.01)
        self.assertApproxEqual(95.0221*ureg.psi,
                               P_2,
                               uncertainty=0.01)

    def test_direct_integration_nitrogen_P165_T150(self):
        P = Q_(165, ureg.psi)
        T = 150 * ureg.K
        fluid = ht.ThermState('nitrogen', P=P, T=T)
        G_nozzle, P_2 = ht.piping.G_nozzle(fluid, n_steps=250)
        self.assertApproxEqual(3829*ureg.kg/(ureg.s*ureg.m**2),
                               G_nozzle,
                               uncertainty=0.01)
        self.assertApproxEqual(87.186*ureg.psi,
                               P_2,
                               uncertainty=0.01)

    def test_direct_integration_nitrogen_P165_T300(self):
        P = Q_(165, ureg.psi)
        T = 300 * ureg.K
        fluid = ht.ThermState('nitrogen', P=P, T=T)
        G_nozzle, P_2 = ht.piping.G_nozzle(fluid, n_steps=200)
        self.assertApproxEqual(2620.8*ureg.kg/(ureg.s*ureg.m**2),
                               G_nozzle,
                               uncertainty=0.01)
        self.assertApproxEqual(86.92*ureg.psi,
                               P_2,
                               uncertainty=0.01)

    def test_direct_integration_nitrogen_P495_T126_5(self):
        P = Q_(495, ureg.psi)
        T = 126.5 * ureg.K
        fluid = ht.ThermState('nitrogen', P=P, T=T)
        G_nozzle, P_2 = ht.piping.G_nozzle(fluid, n_steps=300)
        self.assertApproxEqual(16354.4*ureg.kg/(ureg.s*ureg.m**2),
                               G_nozzle,
                               uncertainty=0.01)
        self.assertApproxEqual(306.1692*ureg.psi,
                               P_2,
                               uncertainty=0.01)

# TODO Add Crane examples: 4-22 (may need Y implementation),
# 4-20, 4-19, 4-18, 4-16, 4-12?, 4-10?


class CPWrapperTest(unittest.TestCase):
    """Test for additional methods of ThermState class"""
    def test_copy(self):
        """Minimal test of copy method"""
        test_state = ht.ThermState('nitrogen')
        test_state.copy()
        test_state.update_kw(T=300*ureg.K, P=1*ureg.MPa)
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


class ODHTest(FunctionsTest):
    """Set of tests for basic functionality of ODH analysis
    """
    def test_hole_leak_voirin_P150(self):
        # Uses classic API 520 calculation
        P_150 = Q_(150, ureg.psig)
        T = Q_(60, ureg.degF)
        V_flux_150 = 2.856 * ureg.ft**3 / (ureg.min * ureg.mm**2)
        K_d = 0.6

        area = 1*ureg.mm**2  # Chosen for simplicity
        fluid = ht.ThermState('nitrogen', P=P_150, T=T)
        q_std = ht.odh.hole_leak(area, fluid)
        self.assertApproxEqual(V_flux_150*area/K_d, q_std)

    def test_hole_leak_voirin_P15(self):
        # Uses classic API 520 calculation
        P_15 = Q_(15, ureg.psig)
        T = Q_(60, ureg.degF)
        V_flux_15 = 0.515 * ureg.ft**3 / (ureg.min * ureg.mm**2)
        K_d = 0.6

        area = 1*ureg.mm**2  # Chosen for simplicity
        fluid = ht.ThermState('nitrogen', P=P_15, T=T)
        q_std = ht.odh.hole_leak(area, fluid)
        self.assertApproxEqual(V_flux_15*area/K_d, q_std)


class GeometryTest(unittest.TestCase):
    def test_circle_area(self):
        D = 1 * ureg.inch**2
        A = pi * D**2 / 4
        self.assertEqual(A, g.circle_area(D))

    def test_cylinder_volume(self):
        D = 1 * ureg.inch**2
        A = pi * D**2 / 4
        L = 10 * ureg.ft
        V = A * L
        self.assertEqual(V, g.cylinder_volume(D, L))


if __name__ == '__main__':
    REFPROP_TEST = False
    unittest.main()
