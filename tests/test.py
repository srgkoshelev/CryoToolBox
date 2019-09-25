import heat_transfer as ht
from scipy.integrate import quad


Test_State = ht.ThermState('helium')
#Test_State.update('PT_INPUTS', ht.Q_('49.17 psi'), ht.Q_('11.029 degR'))
P_SHI = ht.Q_('100 psi')
T_SHI = 2.6579 * P_SHI.to(ht.ureg.psi).magnitude**0.3653 * ht.ureg.degR #Bruce S. formula
Test_State.update('P', P_SHI, 'T', T_SHI)
Test_State.update('T', ht.Q_('4.2 K'), 'Q', ht.Q_('1'))
print(ht.to_scfma(ht.Q_('1 g/s'), Test_State))
print(Test_State.Prandtl)
print(Test_State.Dmolar)
print(Test_State.Cpmass)
print(Test_State._AbstractState.first_partial_deriv(ht.CP.iHmass, ht.CP.iT, ht.CP.iP))
print(Test_State.first_partial_deriv('Hmass', 'T', 'P'))
print(Test_State.specific_heat_input.to(ht.ureg.BTU/ht.ureg.lb))
@ht.ureg.wraps(ht.ureg.BTU/ht.ureg.lb, ht.ureg.psi)
def theta_bruce(P):
    return 0.5724 * P**0.6813
print(theta_bruce(P_SHI))
print('\nCalculating evaporation heat')
Test_State.update('T', ht.Q_('4.2 K'), 'Q', ht.Q_('0'))
Hmass_liq = Test_State.Hmass
print(Hmass_liq)
print(Test_State.specific_heat_input)
Test_State.update('T', ht.Q_('4.2 K'), 'Q', ht.Q_('1'))
Hmass_vap = Test_State.Hmass
print(Hmass_vap)
Hmass_evap = Hmass_vap - Hmass_liq
print(Hmass_evap)
print(Test_State.specific_heat_input)
#print(ht.max_theta(Fluid))
#print(ht.spec_heat(ht.Air))
#Test_pipe = ht.piping.Pipe(1/8, L=ht.ureg('1 m'))
#print(Test_pipe)
#Test_piping = ht.piping.Piping(ht.Air, [Test_pipe])
#print(Test_piping.m_dot(P_out = ht.piping.ureg('1 psi')))
#print ("""100 SCFM of air is equivalent to {:.3g} of Nitrogen flow for P = 1 atm
#       and T = 38 C.""".format(ht.from_scfma(100*ht.ureg('ft^3/min'), {'fluid':'air', 'P':1*ht.ureg.atm, 'T':38*ht.ureg.degC})))
#print ("CGA S-1.3 Formula from 6.1.4 a) gives 0.0547 kg/s for the same air capacity.")
#Re = ht.Re(Fluid, ht.ureg('1g/s'), ht.ureg('5 mm'))
#print(Re)
#theta_temp = ht.theta_temp(ht.ureg('100 K'), ht.ureg('300 K'), ht.ureg('77 K'))
#Bi = ht.Bi(ht.ureg('1 W/(m*K)'), ht.ureg('1 cm'), ht.ureg('10 W/(m**2*K)'))
#Fo = ht.Fo_cyl(theta_temp, Bi)
#print(f'Biot number is: {Bi:.3}')
#print(f'Fourier number is: {Fo:.3}')
#G10_sc = [-2.4083, 7.6006, -8.2982, 7.3301, -4.2386, 1.4294, -0.24396, 0.015236, 0]
#G10_tc = [-4.1236, 13.788, -26.088, 26.272, -14.663, 4.4954, -0.6905, 0.0397, 0] #normal direction
#print(ht.nist_curve_fit(300, G10_tc))
#print(quad(lambda x: ht.nist_curve_fit(x, G10_tc ), 77, 300)[0]/(77-300))
#
#print('\nTesting invert dP calc')
#m_dot = ht.ureg('1000 g/s')
#P_out = ht.ureg('0 psig')
#Test_piping.P_in(m_dot, P_out)
##for p in range(1,100,10):
##    m_dot = ht.ureg('1 g/s')
##    P_test = ht.Q_(p, ht.ureg.psig)
##    Test_piping.init_cond['fluid'] = 'helium'
##    Test_piping.init_cond['P'] = P_test
##    T_test = ht.max_theta(Test_piping.init_cond)
##    Test_piping.init_cond['T'] = T_test
##    print(Test_piping.init_cond['P'].to(ht.ureg.psig), Test_piping.init_cond['T'].to(ht.ureg.K))
##    print (Test_piping.dP(m_dot))
##
##    
##
##if __name__ == "__main__":
##        print (Ra().to_base_units())
##        print (gamma())
##        print (rp_init({'fluid':'helium', 'T':Q_(20,ht.ureg.degC), 'P':Q_(101325, ht.ureg.Pa)}))
##        print (rp_init({'fluid':'helium', 'T':Q_(4.2,ht.ureg.K), 'P':Q_(101325, ht.ureg.Pa)}))
##        print (satp(Q_(101325, ht.ureg.Pa), [1])['t'])
##        print ('Decorator test:', satp(Q_(101325, ht.ureg.Pa), [1]))
##        print(tc_304(150*ht.ureg.K))
##        Leak = tc_304(150*ht.ureg.K)*3.14159*0.125*ht.ureg.inch*0.035*ht.ureg.inch/(1*ht.ureg.ft)*300*ht.ureg.K
##        print(Leak)
##        print(Leak.to(ht.ureg.W))
##        print((Leak/(7*ht.ureg('kJ/kg'))).to(ht.ureg.g/ht.ureg.s))
##        print(therm_exp(ht.ureg('4.5 K'))*ht.ureg('20 ft').to(ht.ureg.inch))
