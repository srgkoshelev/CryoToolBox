
from math import log
from pyrefprop import refprop as rp
from natu.units import *
from natu.units import kPa, uPa, kJ




def Pr (Fluid_data = {'fluid':'air', 'P':101325*Pa, 'T':38*degC}):
	#Function for calculating Prandtl number for given fluid
	#Properties of fluid are obtained from Refprop package
	#Units are handled by natu package
	fluid = Fluid_data['fluid']
	Temperature = Fluid_data['T']
	Pressure = Fluid_data['P']
	prop = rp.setup('def', fluid)
	M = rp.wmol(prop['x'])['wmix']*g/mol
	fluid_prop = rp.flsh ("TP", Temperature/K, Pressure/kPa, prop['x'])
	Density = fluid_prop['Dvap']*mol/L #currently supporting only vapor phase
	fluid_trans_prop = rp.trnprp(Temperature/K, Density/(mol/L), prop['x'])
	mu = fluid_trans_prop['eta']*uPa*s #dynamic viscosity
	k = fluid_trans_prop['tcx']*W/(m*K) #thermal conductivity
	nu = mu/(Density*M) #kinematic viscosity
	Pr = mu*fluid_prop['cp']*(J/(mol*K))/(k*M)
	return Pr

def Gr (Fluid_data = {'fluid':'air', 'P':101325*Pa, 'T':38*degC}, Surface_data = {'Dim':1.315*inch, 'T':100*K}):
	#Function for calculating Grashof number for given fluid and surface
	#Properties of fluid are obtained from Refprop package
	#Units are handled by natu package
	fluid = Fluid_data['fluid']
	T_fluid = Fluid_data['T']
	P_fluid = Fluid_data['P']
	T_surf = Surface_data['T']
	Dim = Surface_data['Dim']
	prop = rp.setup('def', fluid)
	M = rp.wmol(prop['x'])['wmix']*g/mol
	fluid_prop = rp.flsh ("TP", T_fluid/K, P_fluid/kPa, prop['x'])
	D_fluid = fluid_prop['Dvap']*mol/L #currently supporting only vapor phase
	fluid_trans_prop = rp.trnprp(T_fluid/K, D_fluid/(mol/L), prop['x'])
	mu_fluid = fluid_trans_prop['eta']*uPa*s #dynamic viscosity
	nu_fluid = mu_fluid/(D_fluid*M) #kinematic viscosity
	nu_fluid.display_unit = 'm2/s'
	beta_exp = 1/(T_fluid/K)*(1/K) #volumetric thermal expansion coefficient
	Gr = g_0*Dim**3*beta_exp*(T_fluid-T_surf)/nu_fluid**2 #Grashof number
	return Gr

def Ra (Fluid_data = {'fluid':'air', 'P':101325*Pa, 'T':38*degC}, Surface_data = {'Dim':1.315*inch, 'T':100*K}):
	#Function for calculating Rayleigh number for given fluid and surface
	#Properties of fluid are obtained from Refprop package
	#Units are handled by natu package
	return 	Gr(Fluid_data, Surface_data)*Pr(Fluid_data) 

def heat_trans_coef (k_fluid = 0.02647*W/(m*K), Nu = 4, Dim = 1.315*inch, Case = {'convection':'free', 'body':'cyl_hor'}):
	#Calculating heat transfer coefficient for Nu routine
	#Cases like external flow and pipe inside a pipe have different equations for Nu number - to be implemented
	h = k_fluid*Nu/Dim #convective heat transfer coefficient
	return h



def Nu (Fluid_data = {'fluid':'air', 'P':101325*Pa, 'T':38*degC}, Surface_data = {'Dim':10*ft, 'T':100*K, 'Dim_sec':1.315*inch}, Case = {'convection':'free', 'body':'cyl_hor'}):
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

	fluid = Fluid_data['fluid']
	T_fluid = Fluid_data['T']
	P_fluid = Fluid_data['P']
	Dim = Surface_data['Dim']
	prop = rp.setup('def', fluid)
	fluid_prop = rp.flsh ("TP", T_fluid/K, P_fluid/kPa, prop['x'])
	D_fluid = fluid_prop['Dvap']*mol/L #currently supporting only vapor phase
	fluid_trans_prop = rp.trnprp(T_fluid/K, D_fluid/(mol/L), prop['x'])
	k_fluid = fluid_trans_prop['tcx']*W/(m*K) #thermal conductivity
	h = heat_trans_coef (k_fluid, Nu, Dim, Case)	
	h.display_unit = 'W/(m2*K)'
	return {'Nu':Nu, 'h':h}







def rad_hl (eps_cold = 0.55, eps_hot = 0.55, T_hot = 300*K, T_cold = 77*K, F1_2 = 1, eps_baffle = 0.02, N_baffles = 5):
	#F1_2 shows what part of cold surface is being affected by radiation; see Kaganer, p. 42
	Eps_mut = eps_cold*eps_hot/(eps_cold + eps_hot - eps_cold*eps_hot) #Mutual (term?) epsilon, assuming surfaces are equal; otherwise use: Eps_mut = 1/(1/eps_cold + F_cold/F_hot*(1/eps_hot-1))
	q0 = Eps_mut*sigma*(T_hot**4 - T_cold**4)*F1_2
	q0.display_unit = 'W/m2'
	Eps_baffle_mut = eps_baffle/(2-eps_baffle)
	eta = (1+N_baffles*Eps_mut/Eps_baffle_mut)**(-1)
	q_baffle = eta*q0
	return {'q0':q0, 'q_baffle':q_baffle}