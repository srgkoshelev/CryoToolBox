import heat_transfer as ht


print (ht.Air)
#
#    
#
#if __name__ == "__main__":
#        print (Ra().to_base_units())
#        print (gamma())
#        print (rp_init({'fluid':'helium', 'T':Q_(20,ureg.degC), 'P':Q_(101325, ureg.Pa)}))
#        print (rp_init({'fluid':'helium', 'T':Q_(4.2,ureg.K), 'P':Q_(101325, ureg.Pa)}))
#        print (satp(Q_(101325, ureg.Pa), [1])['t'])
#        print ('Decorator test:', satp(Q_(101325, ureg.Pa), [1]))
#        print(tc_304(150*ureg.K))
#        Leak = tc_304(150*ureg.K)*3.14159*0.125*ureg.inch*0.035*ureg.inch/(1*ureg.ft)*300*ureg.K
#        print(Leak)
#        print(Leak.to(ureg.W))
#        print((Leak/(7*ureg('kJ/kg'))).to(ureg.g/ureg.s))
#        print(therm_exp(ureg('4.5 K'))*ureg('20 ft').to(ureg.inch))
