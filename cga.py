#!python3
from math import *
from pyrefprop import refprop as rp
from heat_transfer import functions as ht
from heat_transfer.piping import *
ureg = ht.ureg
Q_ = ureg.Quantity


def max_theta(Fluid_data, step = 0.01):
        """Calculate tepmerature at which sqrt(v)/SHI is max for safety calculations (CGA S-1.3 2008)
        SHI - specific heat input v*(dh/dv)|p
        step - temperature step
        """
        x, __ = ht.rp_init(Fluid_data)
        T_start = rp.satp(Q_(101325, ureg.Pa),x)['t']*ureg.K #Starting temperature - liquid temperature for atmospheric pressure. Vacuum should be handled separately.
        T_end = 300*ureg.K
        theta = 0
        T = T_start
        while T <= T_end:
                D_vap = ht.flsh('TP', T, P, x)['Dvap']
                theta_new = (D_vap**0.5)/rp.therm3(T, D_vap, x)['spht'] #TODO need wrapper
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
#        (x, M, D_fluid) = rp_init(Pipe['fluid'])
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
#        (x, M, D_fluid) = rp_init(Pipe['fluid'])
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
        (x_fluid, M_fluid, D_fluid) = rp_init(Fluid_data)
        T_fluid = Fluid_data['T']
        Z_fluid = therm2(T_fluid, D_fluid, x_fluid)['Z'] #Compressibility factor
        (x_air, M_air, D_air) = rp_init(Air)
        T_air = Air['T']
        Z_air = therm2(T_air, D_air, x_air)['Z'] #Compressibility factor
        Q_air =  M_dot*C_gas_const(Air)/(D_air*M_air*C_gas_const(Fluid_data))*(T_fluid*Z_fluid*M_air/(M_fluid*T_air*Z_air))**0.5
        Q_air.ito(ureg.ft**3/ureg.min)
        return Q_air

def from_scfma (Q_air, Fluid_data):
        """
        Convert SCFM of air into mass flow of specified fluid.
        """
        (x_fluid, M_fluid, D_fluid) = rp_init(Fluid_data)
        T_fluid = Fluid_data['T']
        Z_fluid = therm2(T_fluid, D_fluid, x_fluid)['Z'] #Compressibility factor
        (x_air, M_air, D_air) = rp_init(Air)
        T_air = Air['T']
        Z_air = therm2(T_air, D_air, x_air)['Z'] #Compressibility factor
        M_dot =  Q_air/(C_gas_const(Air))*(D_air*M_air*C_gas_const(Fluid_data))/(T_fluid*Z_fluid*M_air/(M_fluid*T_air*Z_air))**0.5
        M_dot.ito(ureg.kg/ureg.s)
        return M_dot


def C_gas_const (Fluid_data):
        """
        Constant for gal or vapor which is the function of the ratio of specific heats k = Cp/Cv. ASME VIII.1-2015 pp. 423-424.
        """
        k = ht.k(Fluid_data)
        C = 520*(k*(2/(k+1))**((k+1)/(k-1)))**0.5*ureg('lb/(hr*lbf)*(degR)^0.5')
        return C

#def make_sections(Piping, Section_start): #TODO Where is this used? Do I need it still?
#        """
#        Separate Piping into independent sections for relief capacity calculation
#        """
#        Sections = []
#        for i, ind in enumerate(Section_start):
#                if i == 0:
#                        start = 0
#                end = ind
#                Sections.append(Piping[start:end])
#                start = ind
#        Sections.append(Piping[start:])
#        return Sections



if __name__ == "__main__":
    print ("100 SCFM of air is equivalent to {:.3g} of Nitrogen flow for P = 1 atm and T = 38 C.".format(from_scfma(100*ureg('ft^3/min'), {'fluid':'air', 'P':1*ureg.atm, 'T':38*ureg.degC})))
    print ("CGA S-1.3 Formula from 6.1.4 a) gives 0.0547 kg/s for the same air capacity.")
