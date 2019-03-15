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

def Re (M_dot, Fluid_data, Dim):
        """
        Reynolds number.
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

def f_friction(M_dot, pipe, Fluid_data):
        """
        Calculate friction coefficient for pressure drop calculation. 
        Based on Handbook of Hydraulic Resistance by I.E. Idelchik.
        More accurate value using Re.
        """
        Re_num = Re(M_dot, Fluid_data, pipe.ID())
        mult = 1 #Default value for multiplier
        if Re_num < 2000:
            return 64/Re_num*mult
        elif Re_num > 4000:
            return 1/(1.8*log10(Re_num)-1.64)**2*mult
        else:
            #For transitional region the highest of 2 regimes is used.
            return max(64/Re_num*mult, 1/(1.8*log10(Re_num)-1.64)**2*mult) 



        
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

    
#' Functions for basic dimensionless quantities (Re number is a part of piping module).
def Pr (Fluid_data = {'fluid':'air', 'P':Q_(101325,ureg.Pa), 'T':Q_(15,ureg.degC)}):
        """
        Calculate Prandtl number for given fluid
        Properties of fluid are obtained from Refprop package
        """
        (fluid, T_fluid, P_fluid) = unpack_fluid(Fluid_data)
        (x, M, D_fluid) = rp_init(Fluid_data)
        fluid_trans_prop = trnprp(T_fluid, D_fluid, x)
        mu = fluid_trans_prop['eta']*ureg('uPa*s') #dynamic viscosity
        k = fluid_trans_prop['tcx']*ureg('W/(m*K)') #thermal conductivity
        nu = mu/(D_fluid*M) #kinematic viscosity
        Pr = mu*flsh ("TP", T_fluid, P_fluid, x)['cp']/(k*M)
        return Pr

def Gr (Fluid_data = {'fluid':'air', 'P':Q_(101325,ureg.Pa), 'T':Q_(15,ureg.degC)}, Surface_data = {'Dim':Q_(1.315,ureg.inch), 'T':Q_(100,ureg.K)}):
        """
        Calculate Grashof number for given fluid and surface
        Properties of fluid are obtained from Refprop package
        """
        (fluid, T_fluid, P_fluid) = unpack_fluid(Fluid_data)
        T_surf = Surface_data['T']
        Dim = Surface_data['Dim']
        (x, M, D_fluid) = rp_init(Fluid_data)
        fluid_trans_prop = trnprp(T_fluid, D_fluid, x)
        mu_fluid = fluid_trans_prop['eta']*ureg('uPa*s') #dynamic viscosity
        nu_fluid = mu_fluid/(D_fluid*M) #kinematic viscosity
        beta_exp = 1/(T_fluid.to('K')) #volumetric thermal expansion coefficient
        Gr = ureg.g_0*Dim**3*beta_exp*(T_fluid-T_surf)/nu_fluid**2 #Grashof number
        return Gr

def Ra (Fluid_data = {'fluid':'air', 'P':Q_(101325,ureg.Pa), 'T':Q_(15,ureg.degC)}, Surface_data = {'Dim':Q_(1.315,ureg.inch), 'T':Q_(100,ureg.K)}):
        """
        Calculate Rayleigh number for given fluid and surface
        Properties of fluid are obtained from Refprop package
        """
        return  Gr(Fluid_data, Surface_data)*Pr(Fluid_data) 


#' Nu routine calculates both Nu number and convective heat transfer coefficient.
def  heat_trans_coef (k_fluid = 0.02647*ureg('W/(m*K)'), Nu = 4, Dim = Q_(1.315,ureg.inch), Case = {'convection':'free', 'body':'cyl_hor'}):
        """
        Calculate heat transfer coefficient for Nu routine
        Cases like external flow and pipe inside a pipe have different equations for Nu number - to be implemented
        """
        h = k_fluid*Nu/Dim #convective heat transfer coefficient
        return h



def Nu (Fluid_data = {'fluid':'air', 'P':Q_(101325,ureg.Pa), 'T':Q_(15,ureg.degC)}, Surface_data = {'Dim':10*ureg('ft'), 'T':Q_(100,ureg.K), 'Dim_sec':Q_(1.315,ureg.inch)}, Case = {'convection':'free', 'body':'cyl_hor'}):
        """
        Calculate Nusselt number for given Fluid/Surface and specific case.
        Only natural convection currently supported.
        Properties of fluid are obtained from Refprop package
        """
        Prandtl = Pr(Fluid_data)
        Convection = not (Case['convection'] == 'free' or Case['convection'] == 'natural')
        if not Convection:
                Rayleigh = Ra(Fluid_data, Surface_data)
        elif Case['convection'] == 'forced':
                raise BaseException ('{} convection is not implemented. Possible uses: \'natural\', \'free\''.format(Case['convection']))

        C_l = 0.671/(1+(0.492/Prandtl)**(9/16))**(4/9) #Handbook of heat transfer, Rohsenow, Hartnet, Cho (4.13)
        if Case['body'] == 'cyl_hor' and not Convection:
                Nu_T = 0.772*C_l*Rayleigh**(1/4) #Handbook of heat transfer, Rohsenow, Hartnet, Cho (4.45)
                f = 1-0.13/Nu_T**0.16
                Nu_l = 2*f/log(1+2*f*Nu_T)
                C_t = 0.0002*log(Prandtl)**3 - 0.0027*log(Prandtl)**2 + 0.0061*log(Prandtl) + 0.1054
                Nu_t = 0.103*Rayleigh**(1/3)
                Nu = (Nu_l**10 + Nu_t**10)**(1/10) #Nu number, Handbook of heat transfer, Rohsenow, Hartnet, Cho

        elif Case['body'] == 'cyl_vert' and not Convection:
                Dim = Surface_data['Dim']
                Dim_sec = Surface_data['Dim_sec']
                C_t_vert = (0.13*Prandtl**0.22)/(1+0.61*Prandtl**0.81)**0.42 #Handbook of heat transfer, Rohsenow, Hartnet, Cho (4.24)
                Nu_T = C_l * Rayleigh**0.25
                Nu_l_plate = 2/log(1+2/Nu_T) #Handbook of heat transfer, Rohsenow, Hartnet, Cho (4.33)
                zeta = 1.8*Dim/(Dim_sec*Nu_T)  #Handbook of heat transfer, Rohsenow, Hartnet, Cho (4.44)
                Nu_l = zeta/(log(1+zeta)*Nu_l_plate)
                Nu_t = C_t_vert*Rayleigh**(1/3)/(1+1.4e9*Prandtl/Rayleigh)
                Nu = (Nu_l**6 + Nu_t**6)**(1/6)

        else:
                raise ValueError ("Mode \'{}\' is not supported!".format(Mode))

        (fluid, T_fluid, P_fluid) = unpack_fluid(Fluid_data)
        Dim = Surface_data['Dim']
        (x, M, D_fluid) = rp_init(Fluid_data)
        fluid_trans_prop = trnprp(T_fluid, D_fluid, x)
        k_fluid = fluid_trans_prop['tcx']*ureg('W/(m*K)') #thermal conductivity
        h = heat_trans_coef (k_fluid, Nu, Dim, Case)    
        return {'Nu':Nu, 'h':h}


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



