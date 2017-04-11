#!python3
from math import *
from pyrefprop import refprop as rp
from natu import config; config.simplification_level=0
from natu.units import *
from natu.units import kPa, uPa, kJ

import sys
sys.path.append(r'D:\Personal\Python repo')
from heat_transfer import functions as ht
from heat_transfer.piping import *
u=g/mol #numerical equivalent of atomic mass unit


def max_theta(P, x, step = 0.01):
	"""Calculate tepmerature at which sqrt(v)/SHI is max for safety calculations (CGA S-1.3 2008)
	SHI - specific heat input v*(dh/dv)|p
	P - flow rating pressure
	x - flow composition (refprop array)
	step - temperature step
	"""
	T_start = rp.satp(101.325,x)['t'] #Starting temperature - liquid temperature for atmospheric pressure. Vacuum should be handled separately.
	T_end = 300 #K
	theta = 0
	T = T_start
	T_max = 0
	while T <= T_end:
		D_vap = rp.flsh('TP', T, P/kPa, x)['Dvap']
		theta_new = (D_vap**0.5)/rp.therm3(T, D_vap, x)['spht']
		if theta_new > theta:
			theta = theta_new
			T_max = round(T,3)
		else:
			break #Function has only one maximum
		T += step
	return T_max*K


def Q_vac (Pipe, External, P_fr):
	"""
	Calculate required relief capacity due to Vacuum Loss.
	Based on CGA S-1.3 2008 6.2.2. F = 1.

	"""
	fluid = Pipe['fluid']['fluid']
	prop = rp.setup('def', fluid)
	x = prop.get('x', [1])
	M = rp.wmol(x)['wmix']*g/mol
	P = Pipe['fluid']['P']
	P_crit = rp.critp(x)['pcrit']*kPa #Critical pressure 
	T = Pipe['fluid']['T']

	C_coef = 356*(K**0.5*kg/m**3)
	if P <= P_crit: 
		D_vap = rp.satp(P_fr/kPa, x)['Dvap']*mol/L
		Z = rp.therm2(T/K, D_vap/(mol/L), x)['Z'] #Compressibility factor
		r = (rp.flsh("PQ",P/kPa,1,x)['h']*J/(mol) - rp.flsh("PQ",P/kPa,0,x)['h']*J/(mol))/M # latent heat
		G_i = 241*(922*K-T)/(C_coef*r)*(Z*T/(M/u))**0.5

	else:  #If supercritical flow need to use specific input calculation
		D_vap = rp.flsh('TP', T/K, P/kPa, x)['Dvap']*mol/L
		Z = rp.therm2(T/K, D_vap/(mol/L), x)['Z'] #Compressibility factor
		theta = (rp.therm3(T/K, D_vap/(mol/L), x)['spht']*J/mol)/M
		G_i = 241*(922*K-T)/(C_coef*theta)*(Z*T/(M/u))**0.5
	G_i.display_unit = 'K*m3/kJ'

	Surface = make_surface(Pipe)
	if Pipe['Orientation'] == 'Horizontal': #make it a proper function in future
		Case = {'convection':'free', 'body':'cyl_hor'}
	elif Pipe['Orientation'] == 'Vertical':
		Case = {'convection':'free', 'body':'cyl_vert'}
	h = ht.Nu(External, Surface, Case)['h']
	Nuss = ht.Nu(External, Surface, Case)['Nu']
	A_surf = Surface['Dim']*Surface['Dim_sec']
	T_ext = External['T']
	q_conv = h*(T_ext-T) #"W/m^2 convective heat load"
	q_conv.display_unit = 'W/m2'
	q_conv_disp = q_conv*A_surf
	q_conv_disp.display_unit = 'W'
	q_rad = ht.rad_hl(T_hot = T_ext, T_cold = T)['q0']
	q_rad_disp = q_rad*A_surf
	q_rad_disp.display_unit = 'W'
	U_coef = (q_conv+q_rad)/(T_ext-T) #convert(W/m^2-K,kJ/hr-m^2-K)
	U_coef.display_unit = 'kJ/(hr*m2*K)'
	Q_a = 0.383*(328*K - T)/(922*K-T)*G_i*U_coef*A_surf # "m^3/hr of air"
	Q_a.display_unit = 'ft3/min'

	Pipe.update({'q_rad':q_rad_disp, 'q_conv':q_conv_disp})

	return Q_a


def Q_air (Pipe, External, P_fr): #Required relief capacity due to air condensation
	"""
	Calculation of required relief capacity due to Air Condensation.
	Based on CGA S-1.3 2008 6.2.2. F = 1.
	U is calculated from condensation heat load.
	"""
	Diam = OD(Pipe)
	L = Pipe['L']

	fluid = Pipe['fluid']['fluid']
	prop = rp.setup('def', fluid)
	x = prop.get('x', [1])
	M = rp.wmol(x)['wmix']*g/mol
	P = Pipe['fluid']['P']
	P_crit = rp.critp(x)['pcrit']*kPa #Critical pressure; 
	T = Pipe['fluid']['T']

	C_coef = 356*(K**0.5*kg/m**3)
	if P <= P_crit:
		D_vap = rp.satp(P_fr/kPa, x)['Dvap']*mol/L
		Z = rp.therm2(T/K, D_vap/(mol/L), x)['Z'] #Compressibility factor
		r = (rp.flsh("PQ",P/kPa,1,x)['h']*J/(mol) - rp.flsh("PQ",P/kPa,0,x)['h']*J/(mol))/M # latent heat
		G_i = 241*(922*K-T)/(C_coef*r)*(Z*T/(M/u))**0.5

	else:  #If supercritical flow need to use specific input calculation
		D_vap = rp.flsh('TP', T/K, P/kPa, x)['Dvap']*mol/L
		Z = rp.therm2(T/K, D_vap/(mol/L), x)['Z'] #Compressibility factor
		theta = (rp.therm3(T/K, D_vap/(mol/L), x)['spht']*J/mol)/M
		G_i = 241*(922*K-T)/(C_coef*theta)*(Z*T/(M/u))**0.5
	G_i.display_unit = 'K*m3/kJ'


	A_surf = pi*Diam*L
	flux = 0.6*W/cm**2 #Heat flux for  MLI insulated LHe tank from Lehman, Zahn
	q_cond = A_surf*flux
	q_cond.display_unit = 'W'
	T_cold = Pipe['fluid']['T'] #Temperature of the cold surface used for air condensation
	T_ext = External['T']
	U_coef = flux/(T_ext-T_cold) #convert(W/m^2-K,kJ/hr-m^2-K)
	U_coef.display_unit = 'kJ/(hr*m2*K)'

	Q_a = 0.383*(328*K - T)/(922*K-T)*G_i*U_coef*A_surf # "m^3/hr of air"
	Q_a.display_unit = 'ft3/min'

	Pipe.update({'q_cond':q_cond, })

	return Q_a


def to_scfma (M_dot, Fluid_data):
	(x_fluid, M_fluid, D_fluid) = rp_init(Fluid_data)
	T_fluid = Fluid_data['T']
	Z_fluid = rp.therm2(T_fluid/K, D_fluid/(mol/L), x_fluid)['Z'] #Compressibility factor
	Air = {'fluid':'air', 'P':101325*Pa, 'T':38*degC}
	(x_air, M_air, D_air) = rp_init(Air)
	T_air = Air['T']
	Z_air = rp.therm2(T_air/K, D_air/(mol/L), x_air)['Z'] #Compressibility factor
	Q_air =  M_dot/(D_air*M_air)*(T_fluid*Z_fluid*M_air/(M_fluid*T_air*Z_air))**0.5
	Q_air.display_unit = 'ft3/min'
	return Q_air

def from_scfma (Q_air, Fluid_data):
	(x_fluid, M_fluid, D_fluid) = rp_init(Fluid_data)
	T_fluid = Fluid_data['T']
	Z_fluid = rp.therm2(T_fluid/K, D_fluid/(mol/L), x_fluid)['Z'] #Compressibility factor
	Air = {'fluid':'air', 'P':101325*Pa, 'T':38*degC}
	(x_air, M_air, D_air) = rp_init(Air)
	T_air = Air['T']
	Z_air = rp.therm2(T_air/K, D_air/(mol/L), x_air)['Z'] #Compressibility factor
	M_dot =  Q_air*(D_air*M_air)/(T_fluid*Z_fluid*M_air/(M_fluid*T_air*Z_air))**0.5
	M_dot.display_unit = 'kg/s'
	return M_dot