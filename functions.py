from math import log
from pyrefprop import refprop as rp
from pint import UnitRegistry

ureg = UnitRegistry(autoconvert_offset_to_baseunit = True)
Q_ = ureg.Quantity
ureg.load_definitions('D:/Personal/Python repo/pint definitions.txt')
#ureg.auto_reduce_dimensions = True
#from natu.units import *
#from natu.units import kPa, uPa, kJ

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
    Outputs = [('t', ureg.K), ('p', ureg.kPa), ('D', ureg.mol/ureg.L), ('Dliq', ureg.mol/ureg.L), ('Dvap', ureg.mol/ureg.L), ('x', None), ('xliq', None), ('xvap', None), ('q', None), ('e', ureg.J/ureg.mol), ('h', ureg.J/ureg.mol), ('s', ureg.J/(ureg.mol*ureg.K)), ('cv', ureg.J/(ureg.mol*ureg.K)), ('cp', ureg.J/(ureg.mol*ureg.K)), ('w', ureg.m/ureg.s)]
    Output_units = {}
    for output in Outputs:
        if output[0] in refprop_output:
            if output[1]:
                var = refprop_output[output[0]]*output[1]
            else:
                var = refprop_output[output[0]]
            Output_units[output[0]] = var
    return Output_units


















trnprp = ureg.wraps (None, (ureg.K, ureg('mol/L'), None))(rp.trnprp)
satp = ureg.wraps(None, (ureg.kPa, None))(rp.satp) 


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
            D_fluid = fluid_prop['Dvap'] #currently supporting only vapor phase
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
        return {'q0':q0, 'q_baffle':q_baffle, 'eta':eta}


if __name__ == "__main__":
        print (Ra().to_base_units())
        print (gamma())
        print (rp_init({'fluid':'helium', 'T':Q_(20,ureg.degC), 'P':Q_(101325, ureg.Pa)}))
        print (rp_init({'fluid':'helium', 'T':Q_(4.2,ureg.K), 'P':Q_(101325, ureg.Pa)}))
        print (satp(Q_(101325, ureg.Pa), [1])['Dliq']*ureg('mol/L'))

