import heat_transfer as ht
import pprint
import unittest

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
        fluid = ht.ThermState('air.mix', backend='REFPROP')
        self.assertNIST(fluid.M, 28.958600656)
        # TODO Open an issue with CoolProp

    def test_argon(self):
        fluid = ht.ThermState('argon', backend='REFPROP')
        P = 2 * 1000 * ureg.kPa
        D = 15 / fluid.M * ureg.mol/ureg.L
        fluid.update_kw(P=P, Dmolar=D)
        self.assertNIST(637.377588657857, fluid.T.to(ureg.K).magnitude)

    def test_r134a(self):
        fluid = ht.ThermState('r134a', backend='REFPROP')
        T = 400 * ureg.K
        D = 50 / fluid.M * ureg.mol/ureg.L
        fluid.update_kw(T=T, Dmolar=D)
        self.assertNIST(1.45691892789737, fluid.P.to(ureg.kPa).magnitude/1000)

    def test_oxygen(self):
        fluid = ht.ThermState('oxygen', backend='REFPROP')
        T = 100 * ureg.K
        P = 1000 * ureg.kPa
        fluid.update_kw(T=T, P=P)
        self.assertNIST(153.886680663753,
                        fluid.viscosity.to(ureg.uPa*ureg.s).magnitude)

    def test_nitrogen(self):
        fluid = ht.ThermState('nitrogen', backend='REFPROP')
        T = 100 * ureg.K
        fluid.update_kw(T=T, Q=0)
        self.assertNIST(100.111748964945,
                        fluid.conductivity.to(ureg.W/(ureg.m*ureg.K)).magnitude*
                        1000)

    def test_air_density(self):
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
        fluid = ht.ThermState('ammonia&water', backend='REFPROP')
        fluid.set_mole_fractions(0.4, 0.6)
        T = (((300 - 32) * 5 / 9) + 273.15) * ureg.K
        P = 10000 / 14.50377377 * (10**5) / 1000 * ureg.kPa
        fluid.update_kw(T=T, P=P)
        self.assertNIST(5536.79144924071,
                        fluid.speed_sound.to(ureg.m/ureg.s).magnitude *
                        1000 / 25.4 / 12)

    def test_octane(self):
        fluid = ht.ThermState('octane', backend='REFPROP')
        T = (100+273.15) * ureg.K
        fluid.update_kw(T=T, Q=0)
        self.assertNIST(319.167499870568,
                        fluid.latent_heat.to(ureg.J/ureg.g).magnitude)

    def test_R410A_mole_fraction(self):
        fluid = ht.ThermState('R410A.MIX', backend='REFPROP')
        self.assertNIST(0.697614699375863,
                        fluid.mole_fractions[0])

    def test_R410A_mass_fraction(self):
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
        T1 = 27 * ureg.degC
        T2 = -183 * ureg.degC
        A = 0.05 * ureg.m**2
        heat_flow = abs(ht.rad_hl(T1, eps, T2, eps)) * A
        self.assertAlmostEqual(0.23054, heat_flow.to(ureg.W).magnitude, 5)

    def test_rad_hl_2(self):
        eps = 0.55
        T1 = 27 * ureg.degC
        T2 = -183 * ureg.degC
        A = 0.05 * ureg.m**2
        heat_flow = abs(ht.rad_hl(T1, eps, T2, eps)) * A
        self.assertAlmostEqual(8.65727, heat_flow.to(ureg.W).magnitude, 5)

    def test_rad_hl_3(self):
        eps = 0.8
        T1 = 327 * ureg.degC
        T2 = 127 * ureg.degC
        eps_b = 0.05
        heat_flow = abs(ht.rad_hl(T1, eps, T2, eps,
                                  baffles={'N': 1, 'eps': eps_b}))
        self.assertAlmostEqual(145.7373409,
                               heat_flow.to(ureg.W/ureg.m**2).magnitude, 5)

    def test_rad_hl_4(self):
        T1 = 327 * ureg.degC
        eps_1 = 0.8
        T2 = 127 * ureg.degC
        eps_2 = 0.4
        eps_b = 0.05
        heat_flow = abs(ht.rad_hl(T1, eps_1, T2, eps_2,
                                  baffles={'N': 1, 'eps': eps_b}))
        self.assertAlmostEqual(141.3739475,
                               heat_flow.to(ureg.W/ureg.m**2).magnitude, 5)

    def test_API_subsonic(self):
        """Example from API 5.6"""
        m_dot = 53500 * ureg.lb/ureg.hr
        fluid = ht.ThermState('propane&butane', backend='REFPROP')
        fluid.set_mole_fractions(0.5, 0.5)
        fluid.update_kw(T=627*ureg.degR, P=97.2*ureg.psi)
        P_back = 77.2 * ureg.psi
        # print(fluid.M, fluid.compressibility_factor, fluid.gamma)
        A_expect = 6.55 * ureg.inch**2
        A_calc = ht.A_relief_API(m_dot, fluid, P_back=P_back)
        self.assertApproxEqual(A_expect, A_calc, uncertainty=0.05)

    def test_API_sonic(self):
        """Example from API 5.6"""
        m_dot = 53500 * ureg.lb/ureg.hr
        fluid = ht.ThermState('propane&butane', backend='REFPROP')
        fluid.set_mole_fractions(0.5, 0.5)
        fluid.update_kw(T=627*ureg.degR, P=97.2*ureg.psi)
        P_back = 0 * ureg.psig
        # print(fluid.M, fluid.compressibility_factor, fluid.gamma)
        A_expect = 5.73 * ureg.inch**2
        A_calc = ht.A_relief_API(m_dot, fluid, P_back=P_back)
        self.assertApproxEqual(A_expect, A_calc, uncertainty=0.05)


class PipingTest(unittest.TestCase):
    """Simple piping check to see if basic functionality works.

    Doesn't check for correctness yet."""

    def assertApproxEqual(self, data, calc, uncertainty=0.1):
        assert abs(data-calc) / data < uncertainty, \
            f'Calculated value {calc} is not within {uncertainty:.1%} of data \
            value {data}'

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
            test_state,
            [pipe, vj_pipe, corr_pipe, entrance, pipe_exit, orifice, c_orifice,
             tube, c_tube, annulus, pipe_elbow, elbow, pipe_tee, tee, valve,
             # g_valve, v_cone,
             cont, enl])
        piping.volume()
        self.assertApproxEqual(21.2*ureg.psi, piping.dP(Q_('10 g/s')))

        # Test the solver inside Piping.m_dot
        # for i in range(30):
        #     m_dot = random.random()*10**(random.randrange(-1, 5)) * ureg.g/ureg.s
        #     dP_calc = piping.dP(m_dot)
        #     dP = dP_calc
        #     m_dot_calc = piping.m_dot(P_out=ht.AIR.P-dP)
        #     self.assertApproxEqual(m_dot, m_dot_calc)

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
        air = ht.ThermState('air', P=65*ureg.psig, T=110*ureg.degF)
        pipe = ht.piping.Pipe(1, SCH=40, L=75*ureg.ft)
        flow = 100 * ureg.ft**3/ureg.min
        m_dot = flow * ht.AIR.Dmass
        piping = ht.piping.Piping(
            air,
            [pipe, ht.piping.Exit(pipe.ID)])
        self.assertApproxEqual(2.61, piping.dP(m_dot).m_as(ureg.psi))

    def test_Rennels_4_5(self):
        pipe = ht.piping.Pipe(4, SCH=40, L=100*ureg.ft)
        fluid = ht.ThermState('nitrogen', P=100*ureg.psi, T=530*ureg.degR)
        piping = ht.piping.Piping(fluid, [pipe])
        P_out = 84.056 * ureg.psi
        self.assertApproxEqual(10, piping.m_dot(P_out).m_as(ureg.lb/ureg.s))

    def test_direct_integration(self):
        fluid = ht.ThermState('helium', P=200*ureg.psi, T=7*ureg.K)
        self.assertApproxEqual(13920*ureg.kg/(ureg.s*ureg.m**2),
                               ht.piping.G_nozzle(fluid, P_out=1*ureg.atm),
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


if __name__ == '__main__':
    unittest.main()
