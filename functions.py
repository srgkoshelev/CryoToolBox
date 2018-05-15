from math import log, log10
from pyrefprop import refprop as rp
from pint import UnitRegistry
import logging
from functools import wraps
import sys, os

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.WARNING)
logger = logging.getLogger(__name__)

ureg = UnitRegistry(autoconvert_offset_to_baseunit = True)
Q_ = ureg.Quantity
__location__ = os.path.dirname(os.path.abspath(__file__))
ureg.load_definitions(os.path.join(__location__, 'pint definitions.txt'))
sigma = Q_('stefan_boltzmann_constant')

#flsh = ureg.wraps (None, (None, ureg.K, ureg.kPa, None))(rp.flsh)
def flsh(routine, var1, var2, x, kph=1):
    """
    Multifunctional wrapper for refprop flsh function
    """
    if routine == 'TP':   #temperature; pressure
        var1_unitless = var1.to(ureg.K).magnitude
        var2_unitless = var2.to(ureg.kPa).magnitude
    if routine == 'TD':   #temperature; Molar Density
        var1_unitless = var1.to(ureg.K).magnitude
        var2_unitless = var2.to(ureg.mol/ureg.L).magnitude
    if routine == 'TH':   #temperature; enthalpy
        var1_unitless = var1.to(ureg.K).magnitude
        var2_unitless = var2.to(ureg.J/ureg.mol).magnitude
    if routine == 'TS':   #temperature; entropy
        var1_unitless = var1.to(ureg.K).magnitude
        var2_unitless = var2.to(ureg.J/(ureg.mol*ureg.K)).magnitude
    if routine == 'TE':   #temperature; internal energy
        var1_unitless = var1.to(ureg.K).magnitude
        var2_unitless = var2.to(ureg.J/ureg.mol).magnitude
    if routine == 'PD':   #pressure; molar density
        var1_unitless = var1.to(ureg.kPa).magnitude
        var2_unitless = var2.to(ureg.mol/ureg.L).magnitude
    if routine == 'PH':   #pressure; enthalpy
        var1_unitless = var1.to(ureg.kPa).magnitude
        var2_unitless = var2.to(ureg.J/ureg.mol).magnitude
    if routine == 'PS':   #pressure; entropy
        var1_unitless = var1.to(ureg.kPa).magnitude
        var2_unitless = var2.to(ureg.J/(ureg.mol*ureg.K)).magnitude
    if routine == 'PE':   #pressure; internal energy
        var1_unitless = var1.to(ureg.kPa).magnitude
        var2_unitless = var2.to(ureg.J/ureg.mol).magnitude
    if routine == 'HS':   #enthalpy; entropy
        var1_unitless = var1.to(ureg.J/ureg.mol).magnitude
        var2_unitless = var2.to(ureg.J/(ureg.mol*ureg.K)).magnitude
    if routine == 'ES':   #internal energy; entropy
        var1_unitless = var1.to(ureg.J/ureg.mol).magnitude
        var2_unitless = var2.to(ureg.J/(ureg.mol*ureg.K)).magnitude
    if routine == 'DH':   #molar density; enthalpy
        var1_unitless = var1.to(ureg.mol/ureg.L).magnitude
        var2_unitless = var2.to(ureg.J/ureg.mol).magnitude
    if routine == 'DS':   #molar density; entropy
        var1_unitless = var1.to(ureg.mol/ureg.L).magnitude
        var2_unitless = var2.to(ureg.J/(ureg.mol*ureg.K)).magnitude
    if routine == 'DE':   #molar density; internal energy
        var1_unitless = var1.to(ureg.mol/ureg.L).magnitude
        var2_unitless = var2.to(ureg.J/ureg.mol).magnitude
    if routine == 'TQ':   #temperature; vapour quality
        var1_unitless = var1.to(ureg.K).magnitude
        var2_unitless = var2
    if routine == 'PQ':   #pressure; vapour qaulity
        var1_unitless = var1.to(ureg.kPa).magnitude
        var2_unitless = var2

    refprop_output = rp.flsh(routine, var1_unitless, var2_unitless, x, kph=1)
    return rp_out_unit(refprop_output)




Outputs = {'t':ureg.K, 'p':ureg.kPa, 'D':ureg.mol/ureg.L, 'Dliq':ureg.mol/ureg.L, 'Dvap':ureg.mol/ureg.L,
           'x':None, 'xliq':None, 'xvap':None, 'q':None, 'e':ureg.J/ureg.mol, 'h':ureg.J/ureg.mol,
           's':ureg.J/ureg.mol*ureg.K, 'cv':ureg.J/(ureg.mol*ureg.K), 'cp':ureg.J/(ureg.mol*ureg.K),
           'w':ureg.m/ureg.s, 'hfmix':None, 'kph':None, 'hrf':None, 'hfld':None, 'nc':None,
           'xkappa':None, 'beta':None, 'xisenk':None, 'xkt':None, 'betas':None, 'bs':None, 'xkkt':None,
           'thrott':None, 'pint':None, 'spht':ureg.J/ureg.mol, 'wmm':ureg.g/ureg.mol, 'ttrp':ureg.K,
           'tnbpt':ureg.K, 'tcrit':ureg.K, 'pcrit':ureg.kPa, 'Dcrit':ureg.mol/ureg.L, 'zcrit':None,
           'acf':None, 'dip':ureg.debye, 'Rgas':ureg.J/(ureg.mol*ureg.K), 'icomp':None,
           }

def rp_value(value, name):
    """Prepare the quantity for passing to rp function by obtaining a dimensionless value.
    """
    if Outputs[name]:
        return value.to(Outputs[name]).magnitude
    else:
        return value

def rp_out_unit (rp_output):
    """ Add units to output of a refprop function.
    """
    Output_w_units = {}
    for quantity, value in rp_output.items():
        try:
            unit = Outputs[quantity]
            if unit:
                Output_w_units[quantity] = value*unit
            else:
                Output_w_units[quantity] = value
        except KeyError:
            logger.warning('Quantity is missing from unit list, please update: {}'.format(quantity))
    return Output_w_units

def rp_unitize(*names):
    """Decorator for refprop functions.
    Given the string names of the variables convert input to rp accepted dimensionless values.
    Results are parsed to assign units.
    Use:
    @rp_unitize('p', 'x', 'kph')
    def satp(p, x, kph=2):
        return rp.satp(p, x, kph)
    """
    def rp_decorate(rp_func):
        def rp_wrapper(*rp_args, **rp_kwargs):
            rp_input_values = list(rp_args) + list(rp_kwargs.values())
            rp_input = []
            for value,name in zip(rp_input_values,names):
                rp_input.append(rp_value(value,name))
            return rp_out_unit(rp_func(*rp_input))
        return rp_wrapper
    return rp_decorate

@rp_unitize('p', 'x', 'kph')
def satp(p, x, kph=2):
    return rp.satp(p, x, kph)

@rp_unitize('t', 'D', 'x')
def therm3(t, D, x):
    return rp.therm3(t, D, x)

@rp_unitize('icomp')
def info(icomp=1):
    return rp.info(icomp)

def latent_heat(Fluid_data):
    """Calculate latent heat/specific heat input for given conditions.
    """
    x,M,D = rp_init(Fluid_data)
    _,T,P = unpack_fluid(Fluid_data)
    crit_props = info()
    T_crit = crit_props['tcrit']
    D_crit = crit_props['Dcrit']
    P_crit = flsh('TD', T_crit, D_crit, x)['p']
    if T < T_crit and P < P_crit: #Two phase region only possible below critical point
        props = flsh('TP', T, P, x)
        quality = props['q']
        if quality < 1: #For 2 phase region (including subcooled liquid) using latent heat of evaporation
            props = satp(P, x)
            Dliq = props['Dliq']
            Dvap = props['Dvap']
            h_liq = flsh('TD', T, Dliq, x)['h']
            h_vap = flsh('TD', T, Dvap, x)['h']
            L = (h_vap - h_liq)/M
        else:
            L = therm3(T,D,x)['spht']/M #Specific heat input
    else:
        L = therm3(T,D,x)['spht']/M #Specific heat input
    return L



trnprp = ureg.wraps (None, (ureg.K, ureg('mol/L'), None))(rp.trnprp)


def gamma (Fluid_data = {'fluid':'air', 'P':Q_(101325,ureg.Pa), 'T':Q_(15,ureg.degC)}):
    """
    Calculate gamma (k) coefficient for the fluid
    gamma = k = Cp/Cv
    """
    (fluid, T_fluid, P_fluid) = unpack_fluid(Fluid_data)
    (x, M, D_fluid) = rp_init(Fluid_data)
    fluid_prop = flsh ("TP", T_fluid, P_fluid, x)
    Cp = fluid_prop['cp']
    Cv = fluid_prop['cv']
    return Cp/Cv

def pack_fluid (fluid, T_fluid = Q_(15, ureg.degC), P_fluid = Q_(101325, ureg.Pa)):
        return {'fluid':fluid, 'P':P_fluid, 'T':T_fluid}

def unpack_fluid (Fluid_data = {'fluid':'air', 'P':Q_(101325,ureg.Pa), 'T':Q_(15,ureg.degC)}):
        fluid = Fluid_data.get('fluid')
        T_fluid = Fluid_data.get('T')
        P_fluid = Fluid_data.get('P')
        return (fluid, T_fluid, P_fluid)

def rp_init (Fluid_data = {'fluid':'air', 'P':Q_(101325,ureg.Pa), 'T':Q_(15,ureg.degC)}):
        """
        Shortcut initialization of refprop
        Returns (x, M, D_fluid)
        """
        (fluid, T_fluid, P_fluid) = unpack_fluid(Fluid_data)
        prop = rp.setup('def', fluid)
        x = prop.get('x', [1]) 
        M = rp.wmol(x)['wmix']*ureg('g/mol')
        if T_fluid and P_fluid:
            fluid_prop = flsh ("TP", T_fluid, P_fluid, x)
            D_fluid = fluid_prop['D'] 
            return (x, M, D_fluid)
        else: 
            return (x, M)


def Pr (Fluid_data = {'fluid':'air', 'P':Q_(101325,ureg.Pa), 'T':Q_(15,ureg.degC)}):
        """
        Calculate Prandtl number for given fluid
        Properties of fluid are obtained from Refprop package
        Units are handled by natu package
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
        Units are handled by natu package
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
        Units are handled by natu package
        """
        return  Gr(Fluid_data, Surface_data)*Pr(Fluid_data) 

def heat_trans_coef (k_fluid = 0.02647*ureg('W/(m*K)'), Nu = 4, Dim = Q_(1.315,ureg.inch), Case = {'convection':'free', 'body':'cyl_hor'}):
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
        Units are handled by natu package
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







def rad_hl (eps_cold = 0.55, eps_hot = 0.55, T_hot = 300*ureg.K, T_cold = 77*ureg.K, F1_2 = 1, eps_baffle = 0.02, N_baffles = 5):
        """
        Calculate radiative heat load including reduction due to baffles
        Based on Kaganer "Thermal insulation in cryogenic engineering", p. 42.
        F1_2 = F_cold/F_hot
        """
        
        Eps_mut = 1/(1/eps_cold + F1_2*(1/eps_hot-1)) #Mutual emissivity
        q0 = Eps_mut*sigma*(T_hot**4 - T_cold**4)*F1_2
        Eps_baffle_mut = eps_baffle/(2-eps_baffle)
        eta = (1+N_baffles*Eps_mut/Eps_baffle_mut)**(-1)
        q_baffle = eta*q0
        return {'q0':q0.to(ureg.W/ureg.m**2), 'q_baffle':q_baffle.to(ureg.W/ureg.m**2), 'eta':eta}

def tc_304(T):
    Coefs = [-1.4087, 1.3982, 0.2543, -0.6260, 0.2334, 0.4256, -0.4658, 0.1650, -0.0199]
    log10_k = 0
    for ind, coef in enumerate(Coefs):
        log10_k += log10(T.magnitude)**ind*coef
    return 10**log10_k*ureg('(W/(m*K))')




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
