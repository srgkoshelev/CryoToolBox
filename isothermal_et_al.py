import heat_transfer as ht
from math import log
from scipy.optimize import root_scalar
from attrdict import AttrDict
ureg = ht.ureg
Q_ = ht.Q_
u = ureg

def check(expected, calculated, uncertainty=0.05, name=''):
    error = float(abs(expected-calculated)/expected)
    if isinstance(expected, u.Quantity) and isinstance(calculated, u.Quantity):
        unit = expected.units
        calculated.ito(unit)
    if error > uncertainty:
        if name:
            print(name)
        print(f"Expected: {expected}, Calculated: {calculated}, error: {error:.3g}")
        return False
    else:
        return True

print("Analyzing Rennel's 4.5 example by hand")
nitrogen = ht.ThermState('nitrogen', P=100*u.psi, T=520*u.degR)
print(nitrogen)
pipe = ht.piping.Pipe(4, L=100*u.ft)
print(pipe)
check(4.026*u.inch, pipe.ID)
check(0.088405*u.ft**2, pipe.area)
P_out = 84.056*u.psi
check(28.013, nitrogen.molar_mass.m_as(u.g/u.mole))  # Potential for errors
def simplified_isothermal(P1, P2, phi, f, Lm, Tavg, Z, Sg, d):
    Tb = 520
    Pb = 14.7
    q_h = 3.2308 * (Tb/Pb) * ((P1**2 - P2**2 - phi)/(f*Lm*Tavg*Z*Sg))**0.5 * d**2.5
    return q_h * u.ft**3/u.hr  # STD flow
f = 0.016466
Sg = nitrogen.molar_mass / ht.Air.molar_mass
q_h_std = simplified_isothermal(nitrogen.P.m_as(u.psi), P_out.m_as(u.psi), 0,
                                f, pipe.L.m_as(u.mile), nitrogen.T.m_as(u.degR),
                                nitrogen.Z, Sg, pipe.ID.m_as(u.inch))
m_dot_crane = ht.to_mass_flow(q_h_std, nitrogen)
print(f"Crane works: {check(10*u.lb/u.s, m_dot_crane)}")
def rennels_m_dot(A, P1, P2, R, T, f, L, D):
    A = pipe.area
    m_dot = A * ((P1**2 - P2**2)/(R*T*(2*log(P1/P2)+f*L/D)))**0.5
    return m_dot
def rennels_m_dot_true(fluid, pipe, P_out):
    A_ = pipe.area.m_as(ureg.m**2)
    P1_ = fluid.P.m_as(ureg.Pa)
    P2_ = P_out.m_as(ureg.Pa)
    R_ = fluid.specific_gas_constant.m_as(ureg.J/ureg.kg/ureg.K)
    T_ = fluid.T.m_as(ureg.K)
    L_ = pipe.L.m_as(ureg.m)
    D_ = pipe.ID.m_as(ureg.m)
    def to_solve(m_dot_):
        Re = ht.Re(fluid, m_dot_*ureg.kg/ureg.s, pipe.ID, pipe.area)
        m_dot = A_ * ((P1_**2 - P2_**2)/(R_*T_*(2*log(P1_/P2_)+pipe.K(Re))))**0.5
        return m_dot - m_dot_
    x0 = 0.01
    x1 = 0.5 * x0
    bracket = [1e-9, 1000]
    solution = root_scalar(to_solve, bracket=bracket)
    m_dot = solution.root * ureg.kg/ureg.s
    return m_dot
def rennels_m_dot_iter(fluid, pipe, P_out, max_iter=50, tol=1e-4, m_dot_guess=None):
    if max_iter <=0:
        raise ValueError(f'No convergence found.')
    A = pipe.area
    P1 = fluid.P
    P2 = P_out
    R = fluid.specific_gas_constant
    T = fluid.T
    L = pipe.L
    D = pipe.ID
    rho = fluid.Dmass
    if m_dot_guess is None:
        f_guess = pipe.f_T()
        K_guess = f_guess * L / D
        m_dot_guess = A * (2*(P1-P2)*rho/K_guess)**0.5
    Re = ht.Re(fluid, m_dot_guess, D, A)
    m_dot = A * ((P1**2 - P2**2)/(R*T*(2*log(P1/P2)+pipe.K(Re))))**0.5
    if abs(m_dot-m_dot_guess)/m_dot > tol:
        print(m_dot_guess.to_base_units(), m_dot.to_base_units())
        m_dot = rennels_m_dot_iter(fluid, pipe, P_out, max_iter-1, tol, m_dot)
    return m_dot
m_dot_rennels = rennels_m_dot(pipe.area, nitrogen.P, P_out,
                            nitrogen.specific_gas_constant, nitrogen.T,
                            f, pipe.L, pipe.ID)
print(f"Rennels m_dot works: {check(10*u.lb/u.s, m_dot_rennels)}")
m_dot_rennels_iter = rennels_m_dot_iter(nitrogen, pipe, P_out)
print(f"Iter Rennels m_dot works: {check(10*u.lb/u.s, m_dot_rennels_iter)}")
Re = ht.Re(nitrogen, 10*u.lb/u.s, pipe.ID, pipe.area)
check(f*pipe.L/pipe.ID, pipe.K(Re))

def Mach(fluid, m_dot, area):
    """Calculate Mach number for given conditions.

    Static pressure and temperature should be used.
    """
    # T = fluid.T
    # P = fluid.P
    # Z = fluid.Z
    # R_uni = ureg.R  # Universal gas constant
    # m = fluid.molar_mass
    # k = fluid.gamma
    # A = area
    # B = m_dot / A * (Z*R_uni/(k*m))**0.5
    # M = B * T**0.5 / P
    M = m_dot/(area*fluid.Dmass*fluid.speed_sound)
    return M.to_base_units()

def Mach_total(fluid, m_dot, area):
    """Calculate Mach number for total temperature and pressure.
    """
    k = float(fluid.gamma)
    M_ = float(Mach(fluid, m_dot, area))
    def M_sq_total(Msq, M_core, k):
        return M_core**2 * (1+Msq*(k-1)/2)**((k+1)/(k-1))

    def to_solve(Msq):
        return Msq - M_sq_total(Msq, M_, k)
    x0 = M_**2
    x1 = 0.9 * x0
    bracket = [0, 1]
    solution = root_scalar(to_solve, x0=x0, x1=x1,
                            )
    M_root = solution.root**0.5

    w = m_dot/(fluid.Dmass*area)
    T = fluid.T - w**2 / fluid.Cpmass
    a = (fluid.gamma*fluid.specific_gas_constant*T)**0.5
    M_simple = float(w / a)
    check(M_root, M_simple, name='Mach comparison, solver vs through T')
    if isinstance(M_root, complex):
        raise ValueError(f'No real solutions for Mach number found.')
    # print("Mach total", solution.iterations)
    return M_root
    # return M_simple

def K_limit(M, k):
    """Calculate max resistance coefficient of a pipe.
    """
    A = (1-M**2) / (k*M**2)
    B = (k+1)/(2*k)
    C = (k+1) * M**2
    D = 2 + (k-1) * M**2
    K = A + B * log(C/D)
    return K

def M_from_K_limit(K, k):
    K_ = float(K)
    def to_solve(M):
        return K_ - K_limit(M, k)
    x0 = 0.9
    x1 = 0.5
    bracket = [1e-15, 1]
    try:
        solution = root_scalar(to_solve, bracket=bracket)
    except ValueError:
        solution = root_scalar(to_solve, x0=x0, x1=x1)
    M = solution.root
    return M

def M_complex(M, k):
    return 1 + M**2 * (k-1) /2

def P_from_M(P1, M1, M2, k):
    PM1 = P1 * M1 * M_complex(M1, k)**0.5
    P2 = PM1 / (M2*M_complex(M2, k)**0.5)
    return P2

def P_total_from_static(P, M, k):
    return P * M_complex(M, k)**(k/(k-1))

def dP_isothermal(m_dot, fluid, pipe, P2):
    R = fluid.specific_gas_constant
    T = fluid.T
    P1 = fluid.P
    A = pipe.area
    Re = ht.Re(fluid, m_dot, pipe.ID, pipe.area)
    K = pipe.K(Re)

    def sq_diff(P2):
        return m_dot**2 * R * T / A**2 * (2*log(P1/P2)+K)
    P_out = (P1**2 - sq_diff(P2))**0.5

    def to_check(P2_):
        P1_ = P1.m_as(ureg.Pa)
        m_dot_ = m_dot.m_as(ureg.kg/ureg.s)
        R_ = R.m_as(ureg.J/ureg.kg/ureg.K)
        T_ = T.m_as(ureg.K)
        A_ = A.m_as(ureg.m**2)
        K_ = float(K)
        sq_diff = m_dot_**2 * R_ * T_ / A_**2 * (2*log(P1_/P2_)+K_)
        return (P1_**2 - sq_diff)**0.5
    check(P_out, to_check(P2.m_as(ureg.Pa))*ureg.Pa, uncertainty=0.0001)

    def to_solve(P2_):
        P1_ = P1.m_as(ureg.Pa)
        m_dot_ = m_dot.m_as(ureg.kg/ureg.s)
        R_ = R.m_as(ureg.J/ureg.kg/ureg.K)
        T_ = T.m_as(ureg.K)
        A_ = A.m_as(ureg.m**2)
        K_ = float(K)
        B = m_dot_**2 * R_ * T_ / A_**2
        sq_diff = B * (2*log(P1_/P2_)+K_)
        return P2_ - (P1_**2 - sq_diff)**0.5, 2*P2_-B/P2_, 2*B/P2_**2+2
    x0 = P1.m_as(ureg.Pa)
    x1 = 0.9 * x0
    bracket = [1e-9, x0]
    methods = ['secant', 'newton', 'halley']
    for method in methods:
        try:
            solution = root_scalar(to_solve, x0=x0, x1=x1, fprime=True,
                                    fprime2=True, bracket=bracket, method=method)
            P_2_root = solution.root * ureg.Pa
            # print(method, solution.iterations)
        except TypeError:
            print(f'{method} method failed for non-square solve')
            P_2_root = P_out
            check(P_out, P_2_root, uncertainty=0.0005)

    def to_solve_sq(P2sq_):
        P1_ = P1.m_as(ureg.Pa)
        m_dot_ = m_dot.m_as(ureg.kg/ureg.s)
        R_ = R.m_as(ureg.J/ureg.kg/ureg.K)
        T_ = T.m_as(ureg.K)
        A_ = A.m_as(ureg.m**2)
        K_ = float(K)
        B = m_dot_**2 * R_ * T_ / A_**2
        sq_diff = B * (2*log(P1_/(P2sq_)**0.5)+K_)
        return P2sq_ - P1_**2 + sq_diff, 1 - B/P2sq_, B/P2sq_**2
    x0 = P1.m_as(ureg.Pa)**2
    x1 = 0.9 * x0
    bracket = [1e-9, x0]
    for method in methods:
        try:
            solution = root_scalar(to_solve_sq, x0=x0, x1=x1, fprime=True,
                                    fprime2=True, bracket=bracket, method=method)
            P_2_root_sq = solution.root**0.5 * ureg.Pa
            print(method, solution.iterations)
        except TypeError:
            print(f'{method} method failed for square solve')
            P_2_root_sq = P_2_root
            check(P_out, P_2_root_sq, uncertainty=0.0005)

    P_final = P_2_root_sq

    return P1 - P_final
dP_rennels = dP_isothermal(10*u.lb/u.s, nitrogen, pipe, P_out)
print(f"Rennels dP works: {check(nitrogen.P-P_out, dP_rennels)}")
M = Mach_total(nitrogen, 10*u.lb/u.s, pipe.area)
# print('Mach at inlet: ', M)
K_lim = K_limit(M, nitrogen.gamma)
# print('K limit: ', K_lim)
# print('K pipe: ', pipe.K(Re).to_base_units())
check(M, M_from_K_limit(K_lim, nitrogen.gamma))
K_left = K_lim - pipe.K(Re)
# print('K left: ', K_left)
M_end = M_from_K_limit(K_left, nitrogen.gamma)
# print('M end: ', M_end)
P_static_end = P_from_M(nitrogen.P, M, M_end, nitrogen.gamma)
# print('P static end: ', P_static_end)
M_out = Mach_total(ht.ThermState('nitrogen', P=P_out, T=nitrogen.T), 10*u.lb/u.s, pipe.area)
P_static_out = P_from_M(nitrogen.P, M, M_out, nitrogen.gamma)
check(M_out, M_end)
check(P_static_end, P_static_out)
P_total_end = P_total_from_static(P_static_end, M, nitrogen.gamma)
# print('P total end: ', P_total_end)
print(f"Mach dP works: {check(P_out, P_total_end)}")

print()
print("Analyzing Crane 7-16")
air = ht.ThermState('air', P=65*u.psig, T=110*u.degF)
print(air)
pipe = ht.piping.Pipe(1, L=75*u.ft)
dP = 2.61 * u.psi
P_out = air.P - dP
q_expected = 100 * u.ft**3/u.minute  # STD flow
m_dot_expected = ht.to_mass_flow(q_expected, air)
Re = ht.Re(air, m_dot_expected, pipe.ID, pipe.area)
eps_r = 0.002
f = ht.piping.f_Darcy(Re, eps_r)
q_h_std = simplified_isothermal(air.P.m_as(u.psi), P_out.m_as(u.psi), 0,
                                f, pipe.L.m_as(u.mile), air.T.m_as(u.degR),
                                air.Z, 1, pipe.ID.m_as(u.inch))
m_dot_crane = ht.to_mass_flow(q_h_std, air)
print(f"Crane works: {check(m_dot_expected, m_dot_crane)}")
m_dot_rennels = rennels_m_dot(pipe.area, air.P, P_out,
                            air.specific_gas_constant, air.T,
                            f, pipe.L, pipe.ID)
print(f"Rennels m_dot works: {check(m_dot_expected, m_dot_rennels)}")
m_dot_rennels_iter = rennels_m_dot_iter(air, pipe, P_out)
print(f"Iter Rennels m_dot works: {check(m_dot_expected, m_dot_rennels_iter)}")
dP_rennels = dP_isothermal(m_dot_expected, air, pipe, P_out)
print(f"Rennels dP works: {check(dP, dP_rennels)}")
M = Mach_total(air, m_dot_expected, pipe.area)
# print('Mach at inlet: ', M)
K_lim = K_limit(M, air.gamma)
# print('K limit: ', K_lim)
# print('K pipe: ', pipe.K(Re).to_base_units())
check(M, M_from_K_limit(K_lim, air.gamma))
K_left = K_lim - pipe.K(Re)
# print('K left: ', K_left)
M_end = M_from_K_limit(K_left, air.gamma)
# print('M end: ', M_end)
M_out = Mach_total(ht.ThermState('air', P=P_out, T=air.T), m_dot_expected, pipe.area)
check(M_out, M_end)
P_static_end = P_from_M(air.P, M, M_end, air.gamma)
# print('P static end: ', P_static_end)
P_static_out = P_from_M(air.P, M, M_out, air.gamma)
check(P_static_end, P_static_out)
P_total_end = P_total_from_static(P_static_end, M, air.gamma)
# print('P total end: ', P_total_end)
print(f"Mach dP works: {check(P_out, P_total_end)}")

print("\nTerry's sonic test")
eps = 0.00015 * u.ft
pipe = ht.piping.Tube(10.42*u.inch, L=32*u.ft, eps=eps)
print(pipe)
m_dot_expected = 10 * u.kg/u.s
P_out = 14.7 * u.psi
dP_Erik = 12.23 * u.psi
P1_assumed = P_out + dP_Erik
helium = ht.ThermState('helium',P=P1_assumed, T=288.89*u.K)
print(helium)
Re = ht.Re(helium, m_dot_expected, pipe.ID, pipe.area)
f = ht.piping.f_Darcy(Re, eps/pipe.ID)
Sg = helium.molar_mass / ht.Air.molar_mass
q_h_std = simplified_isothermal(helium.P.m_as(u.psi), P_out.m_as(u.psi), 0,
                                f, pipe.L.m_as(u.mile), helium.T.m_as(u.degR),
                                helium.Z, Sg, pipe.ID.m_as(u.inch))
m_dot_crane = ht.to_mass_flow(q_h_std, air)
print(f"Crane works: {check(m_dot_expected, m_dot_crane)}")
m_dot_rennels = rennels_m_dot(pipe.area, helium.P, P_out,
                            helium.specific_gas_constant, helium.T,
                            f, pipe.L, pipe.ID)
print(f"Rennels m_dot works: {check(m_dot_expected, m_dot_rennels)}")
m_dot_rennels_iter = rennels_m_dot_iter(helium, pipe, P_out)
print(f"Iter Rennels m_dot works: {check(m_dot_expected, m_dot_rennels_iter)}")
dP_rennels = dP_isothermal(m_dot_expected, helium, pipe, P_out)
print(f"Rennels dP works: {check(dP_Erik, dP_rennels)}")
# print('Static Mach at inlet: ', Mach(helium, m_dot_expected, pipe.area))
# print('Static Mach at outlet: ', Mach(ht.ThermState('helium', P=P_out, T=helium.T), m_dot_expected, pipe.area))
# print('Max Mach isothermal applicability: ', (1/helium.gamma)**0.5)
M = Mach(helium, m_dot_expected, pipe.area)
# print('Mach at inlet: ', M)
K_lim = K_limit(M, helium.gamma)
# print('K limit: ', K_lim)
# print('K pipe: ', pipe.K(Re).to_base_units())
check(M, M_from_K_limit(K_lim, helium.gamma))
K_left = K_lim - pipe.K(Re)
# print('K left: ', K_left)
M_end = M_from_K_limit(K_left, helium.gamma)
# print('M end: ', M_end)
M_out = Mach(ht.ThermState('helium', P=P_out, T=helium.T), m_dot_expected, pipe.area)
check(M_out, M_end)
P_static_end = P_from_M(helium.P, M, M_end, helium.gamma)
# print('P static end: ', P_static_end)
P_static_out = P_from_M(helium.P, M, M_out, helium.gamma)
check(P_static_end, P_static_out)
P_total_end = P_total_from_static(P_static_end, M, helium.gamma)
# print('P total end: ', P_total_end)
print(f"Mach dP works: {check(P_out, P_total_end)}")

print("\n7-21 Crane test, coke")
coke = AttrDict()
Sg = 0.42
coke.P = 125 * u.psig
coke.T = 140 * u.degF
coke.Z = 1
coke.gamma = 1.4
coke.Dmass = ht.Air.Dmass * 0.42
coke.name = 'coke'
pipe = ht.piping.Pipe(3, L=20*u.ft)
coke.specific_gas_constant = ht.Air.specific_gas_constant * 0.42
coke.molar_mass = ht.Air.molar_mass * 0.42
P_out = 14.7 *u.psi
f = pipe.f_T()
m_dot_expected = ht.to_mass_flow(1028000*u.ft**3/u.hr, air)
q_h_std = simplified_isothermal(coke.P.m_as(u.psi), P_out.m_as(u.psi), 0,
                                f, pipe.L.m_as(u.mile), coke.T.m_as(u.degR),
                                coke.Z, Sg, pipe.ID.m_as(u.inch))
m_dot_crane = q_h_std * coke.Dmass
print(f"Pure Crane works: {check(m_dot_expected, m_dot_crane)}")
P_out_sonic = coke.P - 98.1*u.psi
q_h_std = simplified_isothermal(coke.P.m_as(u.psi), P_out_sonic.m_as(u.psi), 0,
                                f, pipe.L.m_as(u.mile), coke.T.m_as(u.degR),
                                coke.Z, Sg, pipe.ID.m_as(u.inch))
m_dot_crane = q_h_std * coke.Dmass
print(f"Crane with sonic pressure works: {check(m_dot_expected, m_dot_crane)}")
m_dot_rennels = rennels_m_dot(pipe.area, coke.P, P_out,
                            coke.specific_gas_constant, coke.T,
                            f, pipe.L, pipe.ID)
print(f"Rennels m_dot works: {check(m_dot_expected, m_dot_rennels)}")
# M = Mach(coke, m_dot_expected, pipe.area)
# print('Mach at inlet: ', M)
# K_lim = K_limit(M, coke.gamma)
# print('K limit: ', K_lim)
# print('K pipe: ', pipe.K(Re).to_base_units())
# check(M, M_from_K_limit(K_lim, coke.gamma))
# K_left = K_lim - pipe.K(Re)
# print('K left: ', K_left)
# M_end = M_from_K_limit(K_left, coke.gamma)
# print('M end: ', M_end)
# M_out = Mach(ht.ThermState('coke', P=P_out, T=coke.T), m_dot_expected, pipe.area)
# check(M_out, M_end)
# P_static_end = P_from_M(coke.P, M, M_end, coke.gamma)
# print('P static end: ', P_static_end)
# P_static_out = P_from_M(coke.P, M, M_out, coke.gamma)
# check(P_static_end, P_static_out)
# P_total_end = P_total_from_static(P_static_end, M, coke.gamma)
# print('P total end: ', P_total_end)
# print(f"Mach dP works: {check(P_out, P_total_end)}")

print()
print("Analyzing Crane 7-22")
air = ht.ThermState('air', P=19.3*u.psig, T=100*u.degF)
print(air)
pipe = ht.piping.Pipe(1/2, SCH=80, L=7.04/6.04*10*u.ft)  # Adjust for entrance K = 1, total 7.04
dP = 19.3 * u.psi
P_out = air.P - dP
q_expected = 3762 * u.ft**3/u.hr  # STD flow
m_dot_expected = ht.to_mass_flow(q_expected, air)
Re = ht.Re(air, m_dot_expected, pipe.ID, pipe.area)
f = 0.0275
q_h_std = simplified_isothermal(air.P.m_as(u.psi), P_out.m_as(u.psi), 0,
                                f, pipe.L.m_as(u.mile), air.T.m_as(u.degR),
                                air.Z, 1, pipe.ID.m_as(u.inch))
m_dot_crane = ht.to_mass_flow(q_h_std, air)
print(f"Crane works: {check(m_dot_expected, m_dot_crane)}")
m_dot_rennels = rennels_m_dot(pipe.area, air.P, P_out,
                            air.specific_gas_constant, air.T,
                            f, pipe.L, pipe.ID)
print(f"Rennels m_dot works: {check(m_dot_expected, m_dot_rennels)}")
m_dot_rennels_iter = rennels_m_dot_iter(air, pipe, P_out)
print(f"Iter Rennels m_dot works: {check(m_dot_expected, m_dot_rennels_iter)}")
dP_rennels = dP_isothermal(m_dot_expected, air, pipe, P_out)
print(f"Rennels dP works: {check(dP, dP_rennels)}")
M = Mach(air, m_dot_expected, pipe.area)
# print('Mach at inlet: ', M)
K_lim = K_limit(M, air.gamma)
# print('K limit: ', K_lim)
# print('K pipe: ', pipe.K(Re).to_base_units())
check(M, M_from_K_limit(K_lim, air.gamma))
K_left = K_lim - pipe.K(Re)
# print('K left: ', K_left)
M_end = M_from_K_limit(K_left, air.gamma)
# print('M end: ', M_end)
M_out = Mach(ht.ThermState('air', P=P_out, T=air.T), m_dot_expected, pipe.area)
check(M_out, M_end)
P_static_end = P_from_M(air.P, M, M_end, air.gamma)
# print('P static end: ', P_static_end)
P_static_out = P_from_M(air.P, M, M_out, air.gamma)
check(P_static_end, P_static_out)
P_total_end = P_total_from_static(P_static_end, M, air.gamma)
# print('P total end: ', P_total_end)
print(f"Mach dP works: {check(P_out, P_total_end)}")

print()
print("Analyzing White 9.12b")
air = ht.ThermState('air', P=300*u.kPa, T=500*u.K)
print(air)
v = 100 * u.m/u.s
pipe = ht.piping.Tube(3*u.cm, L=15*u.m)
P_out = None  # Calculate using mach number
m_dot_expected = v*pipe.area * air.Dmass
f = 0.02
# m_dot_rennels_iter = rennels_m_dot_iter(air, pipe, P_out)
# print(f"Iter Rennels m_dot works: {check(m_dot_expected, m_dot_rennels_iter)}")
# q_h_std = simplified_isothermal(air.P.m_as(u.psi), P_out.m_as(u.psi), 0,
#                                 f, pipe.L.m_as(u.mile), air.T.m_as(u.degR),
#                                 air.Z, 1, pipe.ID.m_as(u.inch))
# m_dot_crane = ht.to_mass_flow(q_h_std, air)
# print(f"Crane works: {check(m_dot_expected, m_dot_crane)}")
# m_dot_rennels = rennels_m_dot(pipe.area, air.P, P_out,
#                               air.specific_gas_constant, air.T,
#                               f, pipe.L, pipe.ID)
# print(f"Rennels m_dot works: {check(m_dot_expected, m_dot_rennels)}")

# dP_rennels = dP_isothermal(m_dot_expected, air, pipe, P_out)
# print(f"Rennels dP works: {check(dP, dP_rennels)}")
M = Mach(air, m_dot_expected, pipe.area)
# print('Mach at inlet: ', M)
K_lim = K_limit(M, air.gamma)
# print('K limit: ', K_lim)
# print('K pipe test: ', f*pipe.L/pipe.ID)
check(M, M_from_K_limit(K_lim, air.gamma))
K_left = K_lim - f*pipe.L/pipe.ID
# print('K left: ', K_left)
M_end = M_from_K_limit(K_left, air.gamma)
# print('M end: ', M_end)
M_out = M_from_K_limit(1, air.gamma)
print(f'Mach checks out: {check(M_out, M_end)}')

print()
print("Analyzing White 9.12c")
air = ht.ThermState('air', P=300*u.kPa, T=500*u.K)
print(air)
v = 100 * u.m/u.s
pipe = ht.piping.Tube(3*u.cm, L=30*u.m)
P_out = None
m_dot_expected = v*pipe.area * air.Dmass
f = 0.02
# q_h_std = simplified_isothermal(air.P.m_as(u.psi), P_out.m_as(u.psi), 0,
#                                 f, pipe.L.m_as(u.mile), air.T.m_as(u.degR),
#                                 air.Z, 1, pipe.ID.m_as(u.inch))
# m_dot_crane = ht.to_mass_flow(q_h_std, air)
# print(f"Crane works: {check(m_dot_expected, m_dot_crane)}")
# m_dot_rennels = rennels_m_dot(pipe.area, air.P, P_out,
#                               air.specific_gas_constant, air.T,
#                               f, pipe.L, pipe.ID)
# print(f"Rennels m_dot works: {check(m_dot_expected, m_dot_rennels)}")

# dP_rennels = dP_isothermal(m_dot_expected, air, pipe, P_out)
# print(f"Rennels dP works: {check(dP, dP_rennels)}")
M = Mach(air, m_dot_expected, pipe.area)
print('Mach at inlet: ', M)
check(0.225, M)
K_lim = K_limit(M, air.gamma)
# print('K limit: ', K_lim)
# print('K pipe: ', pipe.K(Re).to_base_units())
check(M, M_from_K_limit(K_lim, air.gamma))
K_left = K_lim - pipe.K(Re)
# print('K left: ', K_left)
M_end = M_from_K_limit(K_left, air.gamma)
# print('M end: ', M_end)
M_out = 0.174
print(f'Mach checks out: {check(M_out, M_end)}')
check(0.225, M_from_K_limit(11, air.gamma))
check(0.2, M_from_K_limit(14.533, air.gamma))
check(0.15, M_from_K_limit(27.932, air.gamma))
check(0.17, M_from_K_limit(21.115, air.gamma))
check(0.175, M_from_K_limit(19.772, air.gamma))
check(0.1741, M_from_K_limit(20, air.gamma))
