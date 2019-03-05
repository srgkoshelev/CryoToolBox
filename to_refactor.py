#!python3
#From piping.Piping
    def K(self):
        if len(self) > 0:
            K0 = 0*ureg.dimensionless
            A0 = self[0].Area
            #ID_prev = self[0].ID #ID of the previous piping section; used for sudden contraction and enlargement calculations
            for section in self:
#                if ID_prev < section.ID: #Sudden enlargement
#                    K0 += (1-beta(ID_prev, section.ID)**2)**2/beta(ID_prev, section.ID)**4*(A0/section.Area)**2
#                    logger.debug('Enlargement: {:.3g} -> {:.3g}'.format(ID_prev, section.ID))
#                if ID_prev > section.ID: #Sudden contraction
#                    K0 += 0.5*(1-beta(ID_prev, section.ID))/2**0.25*(A0/section.Area)**2
#                    logger.debug('Contraction: {:.3g} -> {:.3g}'.format(ID_prev, section.ID))
                K0 += section.K*(A0/section.Area)**2
#                ID_prev = section.ID
            return (K0, A0)
        else:
            logger.error('Piping has no elements! Use Piping.add to add sections to piping.')

# I think this has been used for Nu, Gr and other routines that may need
# refactoring too

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

def tc_304(T):
    Coefs = [-1.4087, 1.3982, 0.2543, -0.6260, 0.2334, 0.4256, -0.4658, 0.1650, -0.0199]
    log10_k = 0
    for ind, coef in enumerate(Coefs):
        log10_k += log10(T.magnitude)**ind*coef
    return 10**log10_k*ureg('(W/(m*K))')

def polyval(p,x):
    """Simple evaluation of a polynomial at a specific value using Horner's method
    """
    res = 0
    for coef in p:
        res = res*x + coef
    return res

def therm_exp(T, material='304'):
    """Linear expansion, (L-L_293)/L_293, dimensionless
    """
    LinearExpansion = {'304':((1.7127E-8, -2.0261E-5, 9.2683E-3, -3.9811E-1, -2.9554E2,),(23*ureg.K, -300.04), (ureg('4 K'), ureg('300 K'))),
     } #Data from NIST website
    (Coeffs, (T_low, low_val), (T_min, T_max)) = LinearExpansion[material]
    if T<T_min or T>T_max:
        raise ValueError('NIST model is only defined for T=({T_min:.3~}..{T_max:.3~} for {material}: {T:.3~}'.format(T_min=T_min, T_max=T_max, material=material, T=T))
    if T<T_low:
        res = low_val
    else:
        res = polyval(Coeffs, T.to(ureg.K).magnitude)
    return res*1e-5

    

if __name__ == "__main__":
        print (Ra().to_base_units())
        print (gamma())
        print (rp_init({'fluid':'helium', 'T':Q_(20,ureg.degC), 'P':Q_(101325, ureg.Pa)}))
        print (rp_init({'fluid':'helium', 'T':Q_(4.2,ureg.K), 'P':Q_(101325, ureg.Pa)}))
        print (satp(Q_(101325, ureg.Pa), [1])['t'])
        print ('Decorator test:', satp(Q_(101325, ureg.Pa), [1]))
        print(tc_304(150*ureg.K))
        Leak = tc_304(150*ureg.K)*3.14159*0.125*ureg.inch*0.035*ureg.inch/(1*ureg.ft)*300*ureg.K
        print(Leak)
        print(Leak.to(ureg.W))
        print((Leak/(7*ureg('kJ/kg'))).to(ureg.g/ureg.s))
        print(therm_exp(ureg('4.5 K'))*ureg('20 ft').to(ureg.inch))





#def Q_vac (Pipe, External, P_fr):  #TODO Clean/rewrite if needed
#        """
#        Calculate required relief capacity due to Vacuum Loss.
#        Based on CGA S-1.3 2008 6.2.2. F = 1.
#
#        """
#        (x, M, D_fluid) = ht.rp_init(Pipe['fluid'])
#        P = Pipe['fluid']['P']
#        P_crit = rp.critp(x)['pcrit']*kPa #Critical pressure 
#        T = Pipe['fluid']['T']
#
#        C_coef = 356*(K**0.5*kg/m**3)
#        if P <= P_crit: 
#                D_vap = rp.satp(P_fr/kPa, x)['Dvap']*mol/L
#                Z = rp.therm2(T/K, D_vap/(mol/L), x)['Z'] #Compressibility factor
#                r = (rp.flsh("PQ",P/kPa,1,x)['h']*J/(mol) - rp.flsh("PQ",P/kPa,0,x)['h']*J/(mol))/M # latent heat
#                G_i = 241*(922*K-T)/(C_coef*r)*(Z*T/(M/u))**0.5
#
#        else:  #If supercritical flow need to use specific input calculation
#                D_vap = rp.flsh('TP', T/K, P/kPa, x)['Dvap']*mol/L
#                Z = rp.therm2(T/K, D_vap/(mol/L), x)['Z'] #Compressibility factor
#                theta = (rp.therm3(T/K, D_vap/(mol/L), x)['spht']*J/mol)/M
#                G_i = 241*(922*K-T)/(C_coef*theta)*(Z*T/(M/u))**0.5
#        G_i.display_unit = 'K*m3/kJ'
#
#        Surface = make_surface(Pipe)
#        if Pipe['Orientation'] == 'Horizontal': #make it a proper function in future
#                Case = {'convection':'free', 'body':'cyl_hor'}
#        elif Pipe['Orientation'] == 'Vertical':
#                Case = {'convection':'free', 'body':'cyl_vert'}
#        h = ht.Nu(External, Surface, Case)['h']
#        Nuss = ht.Nu(External, Surface, Case)['Nu']
#        A_surf = Surface['Dim']*Surface['Dim_sec']
#        T_ext = External['T']
#        q_conv = h*(T_ext-T) #"W/m^2 convective heat load"
#        q_conv.display_unit = 'W/m2'
#        q_conv_disp = q_conv*A_surf
#        q_conv_disp.display_unit = 'W'
#        q_rad = ht.rad_hl(T_hot = T_ext, T_cold = T)['q0']
#        q_rad_disp = q_rad*A_surf
#        q_rad_disp.display_unit = 'W'
#        U_coef = (q_conv+q_rad)/(T_ext-T) #convert(W/m^2-K,kJ/hr-m^2-K)
#        U_coef.display_unit = 'kJ/(hr*m2*K)'
#        Q_a = 0.383*(328*K - T)/(922*K-T)*G_i*U_coef*A_surf # "m^3/hr of air"
#        Q_a.display_unit = 'ft3/min'
#
#        Pipe.update({'q_rad':q_rad_disp, 'q_conv':q_conv_disp})
#
#        return Q_a


#def Q_air (Pipe, External, P_fr): #Required relief capacity due to air condensation  #TODO Clean/rewrite if needed
#        """
#        Calculation of required relief capacity due to Air Condensation.
#        Based on CGA S-1.3 2008 6.2.2. F = 1.
#        U is calculated from condensation heat load.
#        """
#        Diam = OD(Pipe)
#        L = Pipe['L']
#
#        (x, M, D_fluid) = ht.rp_init(Pipe['fluid'])
#        P = Pipe['fluid']['P']
#        P_crit = rp.critp(x)['pcrit']*kPa #Critical pressure; 
#        T = Pipe['fluid']['T']
#
#        C_coef = C_gas_flow(Pipe['fluid'])
#        if P <= P_crit:
#                D_vap = rp.satp(P_fr/kPa, x)['Dvap']*mol/L
#                Z = rp.therm2(T/K, D_vap/(mol/L), x)['Z'] #Compressibility factor
#                r = (rp.flsh("PQ",P/kPa,1,x)['h']*J/(mol) - rp.flsh("PQ",P/kPa,0,x)['h']*J/(mol))/M # latent heat
#                G_i = 241*(922*K-T)/(C_coef*r)*(Z*T/(M/u))**0.5
#
#        else:  #If supercritical flow need to use specific input calculation
#                D_vap = rp.flsh('TP', T/K, P/kPa, x)['Dvap']*mol/L
#                Z = rp.therm2(T/K, D_vap/(mol/L), x)['Z'] #Compressibility factor
#                theta = (rp.therm3(T/K, D_vap/(mol/L), x)['spht']*J/mol)/M
#                G_i = 241*(922*K-T)/(C_coef*theta)*(Z*T/(M/u))**0.5
#        G_i.display_unit = 'K*m3/kJ'
#
#
#        A_surf = pi*Diam*L
#        flux = 0.6*W/cm**2 #Heat flux for  MLI insulated LHe tank from Lehman, Zahn
#        q_cond = A_surf*flux
#        q_cond.display_unit = 'W'
#        T_cold = Pipe['fluid']['T'] #Temperature of the cold surface used for air condensation
#        T_ext = External['T']
#        U_coef = flux/(T_ext-T_cold) #convert(W/m^2-K,kJ/hr-m^2-K)
#        U_coef.display_unit = 'kJ/(hr*m2*K)'
#
#        Q_a = 0.383*(328*K - T)/(922*K-T)*G_i*U_coef*A_surf # "m^3/hr of air"
#        Q_a.display_unit = 'ft3/min'
#
#        Pipe.update({'q_cond':q_cond, })
#
#        return Q_a



