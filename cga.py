#!python3
from math import *
from heat_transfer import functions as ht
ureg = ht.ureg
Q_ = ureg.Quantity


def max_theta(Fluid_data, step = 0.01):
        """Calculate tepmerature at which sqrt(v)/SHI is max for safety calculations (CGA S-1.3 2008)
        SHI - specific heat input v*(dh/dv)|p
        step - temperature step
        """
        x, _, _ = ht.rp_init(Fluid_data)
        T_start = ht.satp(Q_(101325, ureg.Pa),x)['t']  #Starting temperature - liquid temperature for atmospheric pressure. Vacuum should be handled separately.
        T_end = 300*ureg.K
        theta = ht.ureg('0 mol**1.5/(J*l**0.5)')
        T = T_start
        P = Fluid_data['P']
        while T <= T_end:
            D_vap = ht.flsh('TP', T, P, x)['Dvap']
            theta_new = (D_vap**0.5)/ht.therm3(T, D_vap, x)['spht'] 
            if theta_new > theta:
                theta = theta_new
            else:
                break #Function has only one maximum
            T += step*ureg.K
        return T


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


def to_scfma (M_dot, Fluid_data):
        """
        Convert mass flow rate into SCFM of air.
        Assumption: Flow through a relief device with invariant Area/discharge coefficient (KA).
        """
        (x_fluid, M_fluid, D_fluid) = ht.rp_init(Fluid_data)
        T_fluid = Fluid_data['T']
        Z_fluid = ht.therm2(T_fluid, D_fluid, x_fluid)['Z'] #Compressibility factor
        (x_air, M_air, D_air) = ht.rp_init(ht.Air)
        T_air = ht.Air['T']
        Z_air = ht.therm2(T_air, D_air, x_air)['Z'] #Compressibility factor
        Q_air =  M_dot*C_gas_const(ht.Air)/(D_air*M_air*C_gas_const(Fluid_data))*(T_fluid*Z_fluid*M_air/(M_fluid*T_air*Z_air))**0.5
        Q_air.ito(ureg.ft**3/ureg.min)
        return Q_air

def from_scfma (Q_air, Fluid_data):
        """
        Convert SCFM of air into mass flow of specified fluid.
        """
        (x_fluid, M_fluid, D_fluid) = ht.rp_init(Fluid_data)
        T_fluid = Fluid_data['T']
        Z_fluid = ht.therm2(T_fluid, D_fluid, x_fluid)['Z'] #Compressibility factor
        (x_air, M_air, D_air) = ht.rp_init(ht.Air)
        T_air = ht.Air['T']
        Z_air = ht.therm2(T_air, D_air, x_air)['Z'] #Compressibility factor
        M_dot =  Q_air/(C_gas_const(ht.Air))*(D_air*M_air*C_gas_const(Fluid_data))/(T_fluid*Z_fluid*M_air/(M_fluid*T_air*Z_air))**0.5
        M_dot.ito(ureg.kg/ureg.s)
        return M_dot


def C_gas_const (Fluid_data):
        """
        Constant for gal or vapor which is the function of the ratio of specific heats k = Cp/Cv. ASME VIII.1-2015 pp. 423-424.
        """
        k = ht.k(Fluid_data)
        C = 520*(k*(2/(k+1))**((k+1)/(k-1)))**0.5*ureg('lb/(hr*lbf)*(degR)^0.5')
        return C

def is_super_crit(Fluid_data):
    x, M, D_fluid = ht.rp_init(Fluid_data)
    P_crit = ht.critp(x)['pcrit']
    P = Fluid_data['P']
    return P > P_crit

def spec_heat(Fluid_data):
    """
    Calculate specific heat for given Fluid State.
    For subcritical flow returns latent heat of evaporation,
    for supercritical: specific heat input
    """
    x, M, D_fluid = ht.rp_init(Fluid_data)
    _, T, P = ht.unpack_fluid(Fluid_data)
    Z = ht.therm2(T, D_fluid, x)['Z'] #Compressibility factor
    C = C_gas_const(Fluid_data)
    if is_super_crit(Fluid_data):
        theta = ht.therm3(T, D_fluid, x)['spht']/M
        spec_heat = theta #For simplification
    else:
        h_vapor = ht.flsh("PQ",P,1,x)['h']/M #Saturated vapor enthalpy
        h_liquid = ht.flsh("PQ",P,0,x)['h']/M #Saturated liquid enthalpy
        r = h_vapor - h_liquid
        spec_heat = r #For simplification
    return spec_heat


if __name__ == "__main__":
    print ("100 SCFM of air is equivalent to {:.3g} of Nitrogen flow for P = 1 atm and T = 38 C.".format(from_scfma(100*ureg('ft^3/min'), {'fluid':'air', 'P':1*ureg.atm, 'T':38*ureg.degC})))
    print ("CGA S-1.3 Formula from 6.1.4 a) gives 0.0547 kg/s for the same air capacity.")
    Fluid = ht.Air
    Fluid['P'] = ht.ureg('690 kPa')
    print(max_theta(Fluid))
    print(spec_heat(Fluid))
