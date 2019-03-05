import heat_transfer as ht
from heat_transfer import piping


Fluid = ht.Air
print (Fluid)
print(ht.max_theta(Fluid))
print(ht.spec_heat(ht.Air))
Test_pipe = piping.Pipe(1, L=ht.ureg('1 m'))
print(Test_pipe)
Test_piping = piping.Piping(ht.Air, Test_pipe)
print(Test_piping.m_dot(P_out = piping.ureg('1 psi')))
print ("""100 SCFM of air is equivalent to {:.3g} of Nitrogen flow for P = 1 atm
       and T = 38 C.""".format(ht.from_scfma(100*ht.ureg('ft^3/min'), {'fluid':'air', 'P':1*ht.ureg.atm, 'T':38*ht.ureg.degC})))
print ("CGA S-1.3 Formula from 6.1.4 a) gives 0.0547 kg/s for the same air capacity.")
#
#    
#
#if __name__ == "__main__":
#        print (Ra().to_base_units())
#        print (gamma())
#        print (rp_init({'fluid':'helium', 'T':Q_(20,ht.ureg.degC), 'P':Q_(101325, ht.ureg.Pa)}))
#        print (rp_init({'fluid':'helium', 'T':Q_(4.2,ht.ureg.K), 'P':Q_(101325, ht.ureg.Pa)}))
#        print (satp(Q_(101325, ht.ureg.Pa), [1])['t'])
#        print ('Decorator test:', satp(Q_(101325, ht.ureg.Pa), [1]))
#        print(tc_304(150*ht.ureg.K))
#        Leak = tc_304(150*ht.ureg.K)*3.14159*0.125*ht.ureg.inch*0.035*ht.ureg.inch/(1*ht.ureg.ft)*300*ht.ureg.K
#        print(Leak)
#        print(Leak.to(ht.ureg.W))
#        print((Leak/(7*ht.ureg('kJ/kg'))).to(ht.ureg.g/ht.ureg.s))
#        print(therm_exp(ht.ureg('4.5 K'))*ht.ureg('20 ft').to(ht.ureg.inch))
