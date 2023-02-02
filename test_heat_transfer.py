import heat_transfer as ht
from heat_transfer import odh
import pprint
import unittest
import random
from math import pi
from unittest.mock import MagicMock

pp = pprint.PrettyPrinter()

u = ht.ureg
Q_ = ht.Q_


# Temporarily disabled as no REFPROP available on the machine
# class RefPropTest(unittest.TestCase):
#     """Verify fluid properties are calculated accurately.

#     Tests borrowed from pyrefprop package.
#     """

#     criteria = 0.00001  # NIST acceptance criteria

#     def assertNIST(self, nist, calc, criteria=None):
#         criteria = criteria or self.criteria
#         assert abs(nist-calc) / nist < criteria, \
#             f'Calculated value {calc} is not within {criteria:%} of NIST value ' + \
#             f'{nist}'

#     def test_air(self):
#         fluid = ht.ThermState('air.mix', backend='REFPROP')
#         self.assertNIST(fluid.M, 28.958600656)
#         # TODO Open an issue with CoolProp

#     def test_argon(self):
#         fluid = ht.ThermState('argon', backend='REFPROP')
#         P = 2 * 1000 * u.kPa
#         D = 15 / fluid.M * u.mol/u.L
#         fluid.update_kw(P=P, Dmolar=D)
#         self.assertNIST(637.377588657857, fluid.T.to(u.K).magnitude)

#     def test_r134a(self):
#         fluid = ht.ThermState('r134a', backend='REFPROP')
#         T = 400 * u.K
#         D = 50 / fluid.M * u.mol/u.L
#         fluid.update_kw(T=T, Dmolar=D)
#         self.assertNIST(1.45691892789737, fluid.P.to(u.kPa).magnitude/1000)

#     def test_oxygen(self):
#         fluid = ht.ThermState('oxygen', backend='REFPROP')
#         T = 100 * u.K
#         P = 1000 * u.kPa
#         fluid.update_kw(T=T, P=P)
#         self.assertNIST(153.886680663753,
#                         fluid.viscosity.to(u.uPa*u.s).magnitude)

#     def test_nitrogen(self):
#         fluid = ht.ThermState('nitrogen', backend='REFPROP')
#         T = 100 * u.K
#         fluid.update_kw(T=T, Q=0)
#         self.assertNIST(100.111748964945,
#                         fluid.conductivity.to(u.W/(u.m*u.K)).magnitude*
#                         1000)

#     def test_air_density(self):
#         fluid = ht.ThermState('air.mix', backend='REFPROP')
#         T = (((70 - 32) * 5 / 9) + 273.15) * u.K
#         P = 14.7 / 14.50377377 * (10**5) / 1000 * u.kPa
#         fluid.update_kw(P=P, T=T)
#         self.assertNIST(0.0749156384666842,
#                         fluid.Dmolar.to(u.mol/u.L).magnitude*fluid.M *
#                         0.062427974)

#     # def test_freon_mixture(self):
#     #     fluid = ht.ThermState('R32&R125', backend='REFPROP')
#     #     fluid.set_mole_fractions(0.3, 0.7)
#     #     P = 10 * 1000 * u.kPa
#     #     Smolar = 110 * u.J / (u.mol * u.K)
#     #     fluid.update_kw(P=P, Smolar=Smolar)
#     #     self.assertNIST(23643.993624382,
#     #                     fluid.Hmolar.to(u.J/u.mol).magnitude)
#     # TODO Figure out the issue with mixtures (possibly refprop_setref/ixflag thing)

#     # def test_ethane_butane(self):
#     #     fluid = ht.ThermState('ethane&butane', backend='REFPROP')
#     #     fluid.set_mole_fractions(0.5, 0.5)
#     #     Dmolar = 30 * 0.45359237 / 0.028316846592 / fluid.M * u.mol / u.L
#     #     Hmolar = 283 * 1.05435026448889 / 0.45359237 * fluid.M * u.J / u.mol
#     #     fluid.update_kw(Dmolar=Dmolar, Hmolar=Hmolar)
#     #     self.assertNIST(298.431320311048,
#     #                     fluid.T.to(u.degF).magnitude)
#     # TODO Figure out the issue with mixtures (possibly refprop_setref/ixflag thing)

#     def test_ammonia_water(self):
#         fluid = ht.ThermState('ammonia&water', backend='REFPROP')
#         fluid.set_mole_fractions(0.4, 0.6)
#         T = (((300 - 32) * 5 / 9) + 273.15) * u.K
#         P = 10000 / 14.50377377 * (10**5) / 1000 * u.kPa
#         fluid.update_kw(T=T, P=P)
#         self.assertNIST(5536.79144924071,
#                         fluid.speed_sound.to(u.m/u.s).magnitude *
#                         1000 / 25.4 / 12)

#     def test_octane(self):
#         fluid = ht.ThermState('octane', backend='REFPROP')
#         T = (100+273.15) * u.K
#         fluid.update_kw(T=T, Q=0)
#         self.assertNIST(319.167499870568,
#                         fluid.latent_heat.to(u.J/u.g).magnitude)

#     def test_R410A_mole_fraction(self):
#         fluid = ht.ThermState('R410A.MIX', backend='REFPROP')
#         self.assertNIST(0.697614699375863,
#                         fluid.mole_fractions[0])

#     def test_R410A_mass_fraction(self):
#         fluid = ht.ThermState('R410A.MIX', backend='REFPROP')
#         self.assertNIST(0.5,
#                         fluid.mass_fractions[0])


class FunctionsTest(unittest.TestCase):
    def assertApproxEqual(self, data, calc, uncertainty=0.1):
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

    # Temporarily disabled as no REFPROP available
    # def test_API_subsonic(self):
    #     """Example from API 5.6"""
    #     m_dot = 53500 * u.lb/u.hr
    #     fluid = ht.ThermState('propane&butane', backend='REFPROP')
    #     fluid.set_mole_fractions(0.5, 0.5)
    #     fluid.update_kw(T=627*u.degR, P=97.2*u.psi)
    #     P_back = 77.2 * u.psi
    #     A_expect = 6.55 * u.inch**2
    #     A_calc = ht.A_relief_API(m_dot, fluid, P_back=P_back)
    #     self.assertApproxEqual(A_expect, A_calc, uncertainty=0.05)

    # def test_API_sonic(self):
    #     """Example from API 5.6"""
    #     m_dot = 53500 * u.lb/u.hr
    #     fluid = ht.ThermState('propane&butane', backend='REFPROP')
    #     fluid.set_mole_fractions(0.5, 0.5)
    #     fluid.update_kw(T=627*u.degR, P=97.2*u.psi)
    #     P_back = Q_(0, u.psig)
    #     A_expect = 5.73 * u.inch**2
    #     A_calc = ht.A_relief_API(m_dot, fluid, P_back=P_back)
    #     self.assertApproxEqual(A_expect, A_calc, uncertainty=0.05)


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
        piping = ht.piping.Piping(
            vj_pipe, entrance, pipe_exit, orifice, c_orifice,
             annulus, pipe_tee, tee, valve,
             # g_valve, v_cone,
             cont, enl)
        piping.volume()
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
        piping = ht.piping.Piping(pipe, ht.piping.Exit(pipe.ID))
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
        piping = ht.piping.Piping(pipe)
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
        self.assertApproxEqual(0.0233*u.kg/u.s, m_dot_adiab)

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
        P_out = 100*u.kPa
        pipe = ht.piping.Tube(3*u.cm, wall=0*u.m, L=1*u.m)
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

    def test_piping_stress(self):
        OD = 1.75*u.inch
        ID = 1.625*u.inch
        wall = (OD - ID) / 2
        tube = ht.piping.Tube(OD, wall=wall)
        S = 16700*u.psi
        E = 1
        W = 1
        Y = 0.4
        self.assertApproxEqual(
            S,
            ht.piping.piping_stress(
                tube,
                1070.512*u.psi,
                E=E, W=W, Y=Y), uncertainty=0.01)

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


class ODHTest(FunctionsTest):
    """Set of tests for basic functionality of ODH analysis
    """
    def test_hole_leak_Crane(self):
        LHe = ht.ThermState('helium')
        LHe.update('P', ht.P_NTP, 'Q', Q_('0'))
        Portable_500L = odh.Source('Portable 500L', LHe, Q_(500, u.L))
        Portable_500L.dewar_insulation_failure(1054*u('ft^3/min'))
        fluid = ht.ThermState('water', P=Q_(2, u.psig), T=Q_(60, u.degF))
        tube = ht.piping.Pipe(3, SCH=80)
        d_hole = 2 * u.inch
        area = 3.14159 * d_hole**2 / 4
        self.assertApproxEqual(131*u.gal/u.min,
                               odh.hole_leak(tube, area, fluid),
                               uncertainty=0.2)
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
        m_dot = ht.odh.hole_leak(tube, area, fluid)
        self.assertApproxEqual(V_flux_150*area, m_dot)
        fluid = ht.ThermState('nitrogen', P=P_15, T=T)
        m_dot = ht.odh.hole_leak(tube, area, fluid)
        self.assertApproxEqual(V_flux_15*area, m_dot)


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
        C_t_bar = ht.C_t_bar_cyl(Pr)
        for Ra, Nu in figure_4_17_HHT:
            self.assertApproxEqual(Nu, ht.Nu_hcyl(Pr, Ra, C_t_bar), 0.1)

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
        L = R_D * ID * angle
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


if __name__ == '__main__':
    print("NO REFPROP TESTS PERFORMED!")
    unittest.main()
