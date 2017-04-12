#python3
from natu.units import *
from natu.units import kPa, uPa, kJ
from math import pi, log10, sin, log
from pyrefprop import refprop as rp
if __name__ != "__main__":
	from .functions import rp_init

NPS_raw = [{'NPS': 0.125,
		'OD': 0.405,
		5:  0.035,
		10:  0.049,
		30:  0.057,
		40:  0.068,
		80:  0.095,
		},
		{'NPS': 0.25,
		'OD': 0.540,
		5:  0.049,
		10:  0.065,
		30:  0.073,
		40:  0.088,
		80:  0.119,
		},
		{'NPS': 0.375,
		'OD': 0.675,
		5:  0.049,
		10:  0.065,
		30:  0.073,
		40:  0.091,
		80:  0.126,
		},
		{'NPS': 0.5,
		'OD': 0.840,
		5:  0.065,
		10:  0.083,
		30:  0.095,
		40:  0.109,
		80:  0.147,
		},
		{'NPS': 0.75,
		'OD': 1.050,
		5:  0.065,
		10:  0.083,
		30:  0.095,
		40:  0.113,
		80:  0.154,
		},
		{'NPS': 1,
		'OD': 1.315,
		5:  0.065,
		10:  0.109,
		30:  0.114,
		40:  0.133,
		80:  0.179,
		},
		{'NPS': 1.125,
		'OD': 1.660,
		5:  0.065,
		10:  0.109,
		30:  0.117,
		40:  0.140,
		80:  0.191,
		},
		{'NPS': 1.5,
		'OD': 1.90,
		5:  0.065,
		10:  0.109,
		30:  0.125,
		40:  0.145,
		80:  0.200,
		},
		{'NPS': 2,
		'OD': 2.375,
		5:  0.065,
		10:  0.109,
		30:  0.125,
		40:  0.154,
		80:  0.218,
		},
		{'NPS': 2.5,
		'OD': 2.875,
		5:  0.083,
		10:  0.120,
		30:  0.188,
		40:  0.203,
		80:  0.276,
		},
		{'NPS': 3,
		'OD': 3.5,
		5:  0.083,
		10:  0.120,
		30:  0.188,
		40:  0.216,
		80:  0.300,
		},
		{'NPS': 3.5,
		'OD': 4.0,
		5:  0.083,
		10:  0.120,
		30:  0.188,
		40:  0.226,
		80:  0.318,
		},
	]

NPS_table = {}
for el in NPS_raw:
	NPS_table.update({el['NPS']:{'OD':el['OD']*inch, 5:el[5]*inch, 10:el[10]*inch, 20:el[30]*inch, 30:el[30]*inch, 40:el[40]*inch, 80:el[80]*inch, }})



def OD(Pipe):
	"""
	Return OD of the Pipe element based on NPS table
	"""
	if 'OD' in Pipe:
		return Pipe['OD']

	D_nom = Pipe['D_nom']
	piping_type = Pipe.get('type', 'pipe')
	if piping_type == 'pipe' or piping_type == 'NPS':
		OD_pipe = NPS_table[D_nom]['OD']
		Pipe.update({'OD':OD_pipe})
		return OD_pipe
	elif piping_type == 'tube':
		OD_pipe = D_nom*inch
		Pipe.update({'OD':OD_pipe})
		return OD_pipe
	else:
		raise BaseException ('Wrong piping type: {}.'.format(Pipe['type']))


def wall(Pipe):
	"""
	Return wall thickness of Pipe element based on NPS table
	"""
	if 'wall' in Pipe:
		return Pipe['wall']
	else:
		SCH = Pipe['SCH']
		D_nom = Pipe['D_nom']
		wall_thick = NPS_table[D_nom][SCH]
		Pipe.update({'wall':wall_thick})
	return wall_thick

def ID(Pipe):
	"""
	Return ID of the Pipe element based on NPS table
	"""
	if 'ID' in Pipe:
		return Pipe['ID']
	else:
		ID_pipe = OD(Pipe) - 2*wall(Pipe)
		Pipe.update({'ID':ID_pipe})
		return ID_pipe


def make_surface (Pipe, method = 'OD'):
	"""
	Make surface element for convection heat load calculation.
	Method determines which surface is considered. Orientation changes which dimension should be used for Nu  calculation. 
	"""
	T = Pipe['fluid']['T']
	if method == 'OD':
		Diam = OD(Pipe)
	elif method == 'VJ':
		Diam = VJOD(Pipe)
	elif method == 'average':
		Diam = (OD(Pipe) + VJOD(Pipe))/2

	if Pipe['Orientation'] == 'Horizontal':
		Dim = Diam
		Dim_sec = Pipe['L']
	elif Pipe['Orientation'] == 'Vertical':
		Dim = Pipe['L']
		Dim_sec = Diam
	return {'T':T, 'Dim':Dim, 'Dim_sec':Dim_sec}


def VJOD(Pipe):
	"""
	Return OD of Vacuum Jacket pipe (NPS)
	"""
	Jacket = {'D_nom':Pipe.get('VJ'), 'SCH':Pipe.get('VJSCH', 10)}
	return OD(Jacket)


#Hydraulic functions
def Re (M_dot = 0.01*kg/s, Fluid_data = {'fluid':'air', 'P':101325*Pa, 'T':38*degC}, Dim = 1.097*inch):
	"""
	Calculate Reynolds number for internal flow inside the pipe size of Dim.
	"""
	fluid = Fluid_data['fluid']
	T_fluid = Fluid_data['T']
	P_fluid = Fluid_data['P']


	(x, M, D_fluid) = rp_init(Fluid_data)
	fluid_trans_prop = rp.trnprp(T_fluid/K, D_fluid/(mol/L), x)
	mu_fluid = fluid_trans_prop['eta']*uPa*s #dynamic viscosity

	d = Dim
	A = pi*d**2/4
	rho_fluid = D_fluid*M
	w_flow = M_dot/(rho_fluid*A)
	return w_flow*d*rho_fluid/mu_fluid


#Pressure drops for different configurations of piping

def f_friction(M_dot = 0.01*kg/s, Pipe = {'D_nom':1, 'SCH':10, 'L':10*ft}, Fluid_data = {'fluid':'air', 'P':101325*Pa, 'T':38*degC}):
	"""
	Calculate friction coefficient for pressure drop calculation. Based on Handbook of Hydraulic Resistance by I.E. Idelchik.
	"""
	Reynolds = Re(M_dot, Fluid_data, ID(Pipe))
	if corrugated(Pipe):
		mult = 4 #Using 4x multiplicator compared to straight pipe
	else:
		mult = 1
	if Reynolds < 2000:
		return 64/Reynolds*mult
	elif Reynolds > 4000:
		return 1/(1.8*log10(Reynolds)-1.64)**2*mult
	else:
		# print ("Warning: Re = {:g}, transition flow. Maximum value between pure laminar or pure turbulent flow will be used".format(Reynolds))
		return max (64/Reynolds*mult, 1/(1.8*log10(Reynolds)-1.64)**2*mult)
	

def dp_pipe (M_dot = 0.01*kg/s, Pipe = {'D_nom':1, 'SCH':10, 'L':10*ft}, Fluid_data = {'fluid':'air', 'P':101325*Pa, 'T':38*degC}):
	"""
	Calculate pressure drop for a straight Pipe element. Works with hard pipe and corrugated hose (4x coefficient used). Calculation is based on Crane TP-410.
	"""
	fluid = Fluid_data['fluid']
	T_fluid = Fluid_data['T']
	P_fluid = Fluid_data['P']

	ID_pipe = ID(Pipe)
	L_pipe = Pipe['L']

	(x, M, D_fluid) = rp_init(Fluid_data)
	rho_fluid = D_fluid*M

	A = pi*ID_pipe**2/4
	w_flow = M_dot/(rho_fluid*A)
	w_flow.display_unit = 'm/s'

	f = f_friction(M_dot, Pipe, Fluid_data)

	delta_P = rho_fluid*f*L_pipe*w_flow**2/(2*ID_pipe)
	delta_P.display_unit = 'psi'

	if delta_P/P_fluid < 0.1:
		return delta_P
	elif delta_P/P_fluid < 0.4: #If pressure drop above 10% of input pressure, recalculating for average density
		D_fluid_in = D_fluid
		P_fluid_out = P_fluid - delta_P
		fluid_prop_out = rp.flsh ("TP", T_fluid/K, P_fluid_out/kPa, x)
		D_fluid_out = fluid_prop_out['Dvap']*mol/L #currently supporting only vapor phase
		rho_fluid = (D_fluid_in+D_fluid_out)/2*M
		w_flow = M_dot/(rho_fluid*A)
		f = f_friction(M_dot, Pipe, Fluid_data)
		delta_P = rho_fluid*f*L_pipe*w_flow**2/(2*ID_pipe)
		delta_P.display_unit = 'psi'
		return delta_P
	else:
		raise BaseException ('Pressure drop is {:.0%} which is greater than 40% recommended by Crane TP-410. Consider separating pipeline into sections.'.format(delta_P/P_fluid))


def corrugated (Pipe = {'Corrugated':True}):
	"""
	Determines whether piping is corrugated
	"""
	return Pipe.get('Corrugated', False)

def dp_elbow (M_dot = 0.01*kg/s, Elbow = {'R/D':1, 'D_nom':1, 'SCH':10, 'L':10*ft}, Fluid_data = {'fluid':'air', 'P':101325*Pa, 'T':38*degC}):
	"""
	Pressure drop in an elbow fitting. Based on Handbook of Hydraulic Resistance by I.E. Idelchik.
	"""
		
	delta = Elbow.get('Angle', 90*deg)
	if delta <= 70*deg:
		A1 = 0.9*sin(delta/rad)
	elif delta == 90*deg:
		A1 = 1
	elif delta >= 100*deg:
		A1 = 0.7+0.35*delta/(90*deg)
	else:
		raise BaseException ('Non standard angle is used: {:g}.'.format(delta))

	R_frac_D = Elbow.get('R/D', 1) #Most of the elbows we use have R/D = 1.5 but it is not covered by this calculation. Using conservative value instead
	if (R_frac_D >= 0.5) and (R_frac_D <= 1):
		B1 = 0.21*(R_frac_D)**(-2.5)
	elif R_frac_D >= 10:
		B1 = 0.21*(R_frac_D)**(-0.5)
	else:
		raise BaseException ('Non standard elbow R/D is used: {:g}.'.format(R_frac_D))

	a0_frac_b0 = Elbow.get('a0/b0', 1)
	if a0_frac_b0 == 1:
		C1 = 1
	elif a0_frac_b0 <=4:
		C1 = 0.85+0.125/a0_frac_b0
	else:
		C1 = 1.115-0.84/a0_frac_b0

	zeta = A1*B1*C1

	fluid = Fluid_data['fluid']
	T_fluid = Fluid_data['T']
	P_fluid = Fluid_data['P']

	ID_pipe = ID(Elbow)
	L_pipe = pi*R_frac_D*ID_pipe/2


	(x, M, D_fluid) = rp_init(Fluid_data)
	rho_fluid = D_fluid*M

	A = pi*ID_pipe**2/4
	w_flow = M_dot/(rho_fluid*A)

	delta_P_local = zeta*rho_fluid*w_flow**2/2

	delta_P_frict = dp_pipe (M_dot, Elbow, Fluid_data)

	return delta_P_frict+delta_P_local


#Simple colver for calculating flow
def derivative(f, x, h):
		return (f(x+h) - f(x-h)) / (2*h)  # might want to return a small non-zero if ==0


def solve(f, x0, h):
	lastX = x0
	nextX = lastX + 10* h  # different than lastX so loop starts OK
	i = 0
	i_max = 10
	while (abs(lastX - nextX) > h): 
		newY = f(nextX)                     
		# print ("f(", nextX, ") = ", newY)     # print out progress... again just debug
		lastX = nextX
		nextX = lastX - newY / derivative(f, lastX, h)  # update estimate using N-R
		i += 1
		if i > i_max:
			raise BaseException ("To many cycles! Possible instability.")
	return nextX

def flow_solve (f, x0, step = 1e-3*kg/s):
	"""
	Simpe iteration-based solver, backup if Newton solver fails.
	"""
	X = x0
	while f(X) < -0.5*psi:
		X += step
		# print (X, f(X))
	return X


def calculate_flow(Piping, P_in=1e-3*psig, P_out=0*psig, M_dot_0 = 1e-3*kg/s, M_step = 0.001*kg/s):
	def to_solve (M_dot):
		return dP_piping(M_dot, Piping)-(P_in-P_out)

	try:
		return solve(to_solve, M_dot_0, M_step)
	except:
		return flow_solve(to_solve, M_dot_0, M_step)


def dP_piping(M_dot, Piping):
	"""
	Calculate pressure drop for whole piping; fluid properties should be already initialized for this to work.
	Enthalpy is assumed constant. Each two straight pipes are connected by an elbow.
	"""
	(x, M, D_fluid) = rp_init(Piping[0]['fluid'])
	P_in = Piping[0]['fluid']['P']
	for sec_num, Pipe in enumerate (Piping):
		P = Pipe['fluid']['P']
		if 'h' not in Pipe['fluid']:
			#Calculating enthalpy for the first pipe
			T = Pipe['fluid']['T']
			h = rp.flsh('TP',T/K, P/kPa, x)['h']*(J/mol)
			Pipe['fluid'].update({'h':h})
		else:
			#Calculating new conditions for every pipe after the first
			h = Pipe['fluid']['h']
			T = rp.flsh('PH', P/kPa, h/(J/mol), x)['t']*K
			Pipe['fluid'].update({'T':T})
			
		delta_P = dp_pipe(M_dot, Pipe, Pipe['fluid'])
		if sec_num < len (Piping) - 1:
			if Pipe['Corrugated'] == False and Piping[sec_num+1]['Corrugated'] == False:
				delta_P += dp_elbow(M_dot, Pipe, Pipe['fluid'])
			if delta_P < P:
				P_next = P-delta_P
				h_next = h
				Piping[sec_num+1]['fluid'].update({'P':P_next, 'h':h_next})
			else:
				delta_P.display_unit = 'psi'
				raise BaseException ('Pressure drop {:.3} is greater than input pressure {:.3}! '.format(delta_P, P))
		else:
			P_out = P-delta_P
			Delta_P = P_in - P_out
			Delta_P.display_unit = 'psi'

	return Delta_P


if __name__ == "__main__":
	from functions import rp_init
	print("Testing piping module...")
	print ("D_nom    OD    SCH    wall")
	for D in [.125, .25, .375, .5, .75, 1, 1.125, 1.5, 2, 2.5, 3, 3.5]:
		for sch in [5, 10, 20, 40, 80]:
			Pipe = {'D_nom':D, 'SCH':sch}
			print (D, OD(Pipe), sch, wall(Pipe))

	Test_pipe = {'D_nom':1, 'SCH':10}

	print("Here is a dimensionless Reynolds number for 0.01 kg/s of standard air in 1 in SCH 10 pipe: {:g}. The pipe is stainless if you want to know.".format(Re()))
	print("Pressure drop for 10 ft pipe would be {:g}, while for corrugated hose would 4 times bigger: {:g}".format (dp_pipe(), dp_pipe(Pipe = {'Corrugated':True, 'D_nom':1, 'SCH':10, 'L':10*ft})))
	print ("And for a 90 deg elbow, pressure drop is {:g}.".format(dp_elbow()))

	Piping = [{'D_nom':1, 
			'SCH':10,
			'L':10*ft,
			'VJ':2.5,
			'Orientation':'Vertical', 
			'Corrugated':False, #Corrugated or straight
			'fluid': {'fluid':'air', 'P':101325*Pa, 'T':38*degC}
			},
			]
	print ("Flow for pressure drop for the straight pipe above: {:g}".format(calculate_flow(Piping, P_in=0.0482788*psig)))

