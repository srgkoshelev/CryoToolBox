"""Helpful thermodynamic functions.
"""

from math import log, log10, pi
from . import ureg, Q_
from . import Air
from .cp_wrapper import ThermState
from . import T_NTP, P_NTP
from . import cga
from . import logger
from scipy.interpolate import interp1d
from scipy.integrate import quad

E_TNT = Q_('4850 J/g')  # TNT equivalent Energy of Explosion (PNNL)
z_1 = Q_('200 ft')  # Scaled distance for debris and missile damage (PNNL)
z_2 = Q_('15 ft')  # Scaled distance for eardrum rupture (PNNL)
z_3 = Q_('6.7 ft')  # Scaled distance for lung damage (PNNL)

sigma = ureg.stefan_boltzmann_constant
# Basic thermodynamic functions


def to_scfma(M_dot_fluid, fluid):
    """
    Convert mass flow rate into equivalent flow of air.
    Flow through a relief device with invariant Area/discharge coefficient
    (KA).

    Parameters
    ----------
    M_dot_fluid : Quantity {mass: 1, time: -1}
        mass flow rate
    fluid : ThermState

    Returns
    -------
    ThermState
        volumetric air flow rate
    """
    C_fluid = fluid.C_gas_const
    C_air = Air.C_gas_const

    # Calculation
    M_dot_air = M_dot_fluid * C_air / C_fluid * Air.MZT / fluid.MZT
    Q_air = M_dot_air / Air.Dmass
    Q_air.ito(ureg.ft**3/ureg.min)
    return Q_air


def from_scfma(Q_air, fluid):
    """
    Convert volumetric air flow rate into equivalent mass flow of specified
    fluid. Flow through a relief device with invariant Area/discharge
    coefficient (KA).
    Invert function to to_scfma().

    Parameters
    ----------
    Q_air : Quantity {length: 3, time: -1}
        volumetric air flow rate
    fluid : ThermState

    Returns
    -------
    Quantity {mass: 1, time: -1}
        mass flow rate
    """
    C_fluid = fluid.C_gas_const
    C_air = Air.C_gas_const

    # Calculation
    M_dot_air = Q_air * Air.Dmass
    M_dot_fluid = M_dot_air * C_fluid / C_air * fluid.MZT / Air.MZT
    M_dot_fluid.ito(ureg.g/ureg.s)
    return M_dot_fluid


def to_standard_flow(flow_rate, fluid):
    '''
    Converting volumetric flow at certain conditions or mass flow to
    flow at NTP.
    '''
    fluid_NTP = ThermState(fluid.name)
    fluid_NTP.update('T', T_NTP, 'P', P_NTP)
    if flow_rate.dimensionality == ureg('kg/s').dimensionality:
        # mass flow, flow conditions are unnecessary
        q_std = flow_rate / fluid_NTP.Dmass
    elif flow_rate.dimensionality == ureg('m^3/s').dimensionality:
        # volumetric flow given, converting to standard pressure and
        # temperature
        if fluid.Dmass != -float('Inf')*ureg.kg/ureg.m**3:
            # By default ThermState is initialized with all fields == -inf
            q_std = flow_rate * fluid.Dmass / fluid_NTP.Dmass
        else:
            logger.warning('''Flow conditions for volumetric flow {:.3~}
                           are not set. Assuming standard flow at NTP.
                           '''.format(flow_rate))
            q_std = flow_rate
    else:
        logger.error('''Flow dimensionality is not supported: {:.3~}.
                       '''.format(flow_rate.dimensionality))
    q_std.ito(ureg.ft**3/ureg.min)
    return q_std


def to_mass_flow(Q_std, fluid):
    """
    Calculate mass flow for given volumetric flow at standard conditions.
    """
    fluid_NTP = ThermState(fluid.name)
    fluid_NTP.update('T', T_NTP, 'P', P_NTP)
    m_dot = Q_std * fluid_NTP.Dmass
    return m_dot.to(ureg.g/ureg.s)


def m_max(fluid, A):
    """Calculate max isentropic flow at sonic condition
    (9.46a, Fluid Mechanics, F. White, 2015)
    """
    C = fluid.C_gas_const
    P = fluid.P
    m_max_ = C * A * P * fluid.MZT
    return m_max_.to_base_units()


def PRV_flow(A, Kd, fluid):
    """Calculate mass flow through the relief valve based on
    BPVC VIII div. 1 UG-131 (e) (2).
    """
    W_T = m_max(fluid, A)  # Theoretical flow
    W_a = W_T * Kd  # Actual flow
    return W_a.to(ureg.g/ureg.s)


def A_relief_API(m_dot, fluid, *, P_back=P_NTP, K_d=0.975, K_b=1, K_c=1):
    """Calculate required relief area for given flow as per
    API 520 5.6.3/4

    Parameters
    ----------
    m_dot : Quantity {mass:1, time:-1}
        mass flow rate
    fluid : ThermState at relief conditions
    P_back : ureg.Quantity {length: -1, mass: 1, time: -2}
        backpressure
    K_d : discharge coefficient
        0.975 - when PRV installed with/without a rupture disk
        0.62 - for rupture disc only (see 5.11.1.1.2)
    K_c : capacity correction factor due to backpressure
        applies to balanced bellows valves only
    K_c : combination correction factor
        1 - for no rupture disc installed in combination
        0.9 - for rupture disc installed in combination
    """
    W = m_dot.m_as(ureg.lb/ureg.hr)
    P_1 = fluid.P.m_as(ureg.psi)
    P_2 = P_back.m_as(ureg.psi)
    k = fluid.gamma
    T = fluid.T.m_as(ureg.degR)
    Z = fluid.compressibility_factor
    M = fluid.M
    P_cf = P_1 * (2/(k+1))**(k/(k-1))
    if P_2 <= P_cf:
        critical = True  # For future verbose use
        C = fluid.C.m_as(ureg.lb/(ureg.hr*ureg.lbf)*(ureg.degR)**0.5)
        A = W / (C*K_d*P_1*K_b*K_c) * (T*Z/M)**0.5
    else:
        critical = False
        r = P_2 / P_1
        F_2_brack = (1-r**((k-1)/k)) / (1-r)
        F_2 = (k/(k-1) * r**(2/k) * F_2_brack)**0.5
        A = W / (735*F_2*K_d*K_c) * (T*Z/(M*P_1*(P_1-P_2)))**0.5
    return A * ureg.inch**2


def theta_heat(fluid, step=0.01):
    logger.warning('Deprecated. Use ht.cga.theta() instead.')
    return cga.theta(fluid, step)


def rad_hl(T_1, eps_1, T_2, eps_2, F1_2=1, baffles={'N': 0, 'eps': 0.02}):
    """
    Calculate radiative heat load including reduction due to baffles.
    Based on Kaganer "Thermal insulation in cryogenic engineering", p. 42.

    Parameters
    ----------
    eps_1 : float
        emissivity of the first surface
    eps_2 : float
        emissivity of the second surface
    T_1 : Quality {temperature: 1}
        temperature of the first surface
    T_2 : Quality {temperature: 1}
        temperature of the second surface
    F1_2 : float
        F1_2 = F_cold/F_hot
    baffles : dict
        N - number of baffles
        eps - emissivity of the baffle, assumed to be same on both sides

    Returns
    -------
    dict
            :q0: heat load without any baffles
            q_baffle : heat load with the baffles
            eta : effectiveness of the baffles
    """
    # TODO This function will be refactored
    N_baffles = baffles['N']
    eps_baffle = baffles['eps']

    eps_mut = 1/(1/eps_1 + F1_2*(1/eps_2-1))  # Mutual emissivity
    q0 = eps_mut*sigma*(T_2**4 - T_1**4)*F1_2
    eps_baffle_mut = eps_baffle/(2-eps_baffle)
    eta = (1+N_baffles*eps_mut/eps_baffle_mut)**(-1)
    q_baffle = eta*q0
    return q_baffle.to(ureg.W/ureg.m**2)


def Re(fluid, m_dot, D_H, A_cross):
    """
    Calculate Reynolds number.

    Parameters
    ----------
    fluid : ThermState object describing thermodynamic state (fluid, T, P)
    m_dot : mass flow
    D_H : hydraulic diameter of a pipe
    A_cross : cross-sectional area of a pipe

    Returns
    -------
    Reynolds number, dimensionless
    """
    rho = fluid.Dmass
    v = m_dot / (rho*A_cross)
    mu = fluid.viscosity
    Re_ = v * D_H * rho / mu
    return float(Re_)


def Pr(fluid):
    """
    Calculate Prandtl number.

    Parameters
    ----------
    fluid : ThermState object describing thermodynamic state (fluid, T, P)

    Returns
    -------
    Prandtl number, dimensionless
    """
    return fluid.Prandtl


def Gr(fluid, T_surf, L_surf):
    """
    Calculate Grashof number.

    Parameters
    ----------
    fluid : ThermState object describing thermodynamic state (fluid, T, P)
    T_surf : surface temperature
    L_surf : characteristic length

    Returns
    -------
    Grashof number, dimensionless
    """
    nu_fluid = fluid.viscosity/fluid.Dmass  # kinematic viscosity
    beta_exp = fluid.isobaric_expansion_coefficient
    Gr_ = ureg.g_0 * L_surf**3 * beta_exp * abs(T_surf-fluid.T) / nu_fluid**2
    return Gr_.to(ureg.dimensionless)


def Ra(fluid, T_surf, L_surf):
    """
    Calculate Rayleigh number.

    Parameters
    ----------
    fluid : ThermState object describing thermodynamic state (fluid, T, P)
    T_surf : surface temperature
    L_surf : characteristic length

    Returns
    -------
    Rayleigh number, dimensionless
    """
    return Gr(fluid, T_surf, L_surf)*fluid.Prandtl


def Nu_cyl_hor(fluid, T_cyl, D_cyl):
    """
    Calculate Nusselt number for vertical cylinder.
    Only natural convection currently supported.
    Based on Handbook of heat transfer by Rohsenow, Hartnet,
    Cho (HHT).

    Parameters
    ----------
    fluid : ThermState object describing thermodynamic state (fluid, T, P)
    T_cyl : surface temperature
    D_cyl : cylinder diameter

    Returns
    -------
    Nusselt number, dimensionless
    """
    Pr_ = fluid.Prandtl
    Ra_ = Ra(fluid, T_cyl, D_cyl)
    C_l = 0.671/(1+(0.492/Pr_)**(9/16))**(4/9)  # HHT (4.13)
    Nu_T = 0.772*C_l*Ra_**(1/4)  # HHT (4.45)
    f = 1-0.13/Nu_T**0.16
    Nu_l = 2*f/log(1+2*f*Nu_T)
    C_t = 0.0002*log(Pr_)**3 - 0.0027*log(Pr_)**2 + 0.0061*log(Pr_) + 0.1054
    Nu_t = 0.103*Ra_**(1/3)
    Nu_ = (Nu_l**10 + Nu_t**10)**(1/10)
    return Nu_.to(ureg.dimensionless)


def Nu_cyl_vert(fluid, T_cyl, D_cyl, L_cyl):
    """
    Calculate Nusselt number for vertical cylinder.
    Only natural convection currently supported.
    Based on Handbook of heat transfer by Rohsenow, Hartnet,
    Cho (HHT).

    Parameters
    ----------
    fluid : ThermState object describing thermodynamic state (fluid, T, P)
    T_cyl : surface temperature
    D_cyl : cylinder diameter
    L_cyl : cylinder length

    Returns
    -------
    Nusselt number, dimensionless
    """
    Pr_ = fluid.Prandtl
    Ra_ = Ra(fluid, T_cyl, D_cyl)
    C_l = 0.671/(1+(0.492/Pr_)**(9/16))**(4/9)  # HHT (4.13)
    C_t_vert = (0.13*Pr_**0.22)/(1+0.61*Pr_**0.81)**0.42  # HHT (4.24)
    Nu_T_plate = C_l * Ra_**0.25
    Nu_l_plate = 2 / log(1+2/Nu_T_plate)  # HHT (4.33)
    zeta = 1.8 * L_cyl / (D_cyl*Nu_T_plate)   # HHT (4.44)
    Nu_l = zeta / (log(1+zeta)*Nu_l_plate)
    Nu_t = C_t_vert*Ra_**(1/3)/(1+1.4e9*Pr_/Ra_)
    Nu_ = (Nu_l**6 + Nu_t**6)**(1/6)
    return Nu_.to(ureg.dimensionless)


def heat_trans_coef(fluid, Nu, L_surf):
    """
    Calculate heat transfer coefficient.

    Parameters
    ----------
    fluid : ThermState object describing thermodynamic state (fluid, T, P)
    Nu : Nusselt number
    L_surf : characteristic length:
        :Horizontal cylinder: L_surf = D_cyl
        :Vertical cylinder: L_surf = L_cyl

    Returns
    -------
    heat transfer coefficient
    """
    h = fluid.conductivity * Nu / L_surf
    return h.to(ureg.W/(ureg.m**2*ureg.K))


def Bi(k, L_c, h):
    """
    Calculate Biot number for a solid.

    Parameters
    ----------
    k : thermal conductivity of the solid
    L_c : characteristic length; L_c = V/A_s, where
        V : volume of the solid
        A_s : surface area of the solid
    h : heat transfer coefficient

    Returns
    -------
    Biot number, dimensionless
    """
    Bi_ = h * L_c / k
    return Bi_.to_base_units()


_zeta1_cyl_data = [0.1412, 0.1995, 0.2440, 0.2814, 0.3143, 0.3438, 0.3709,
                   0.3960, 0.4195, 0.4417, 0.5376, 0.6170, 0.6856, 0.7456,
                   0.8516, 0.9408, 1.0184, 1.0873, 1.1490, 1.2048, 1.2558,
                   1.5994, 1.7887, 1.9081, 1.9898, 2.0490, 2.0937, 2.1286,
                   2.1566, 2.1795, 2.2881, 2.3261, 2.3455, 2.3572, 2.3809]
# Table 5.1, Fundamentals of Heat and Mass Transfer, F. Incropera, 2006.


_C1_cyl_data = [1.0025, 1.0050, 1.0075, 1.0099, 1.0124, 1.0148, 1.0173, 1.0197,
                1.0222, 1.0246, 1.0365, 1.0483, 1.0598, 1.0712, 1.0932, 1.1143,
                1.1345, 1.1539, 1.1724, 1.1902, 1.2071, 1.3384, 1.4191, 1.4698,
                1.5029, 1.5253, 1.5411, 1.5526, 1.5611, 1.5677, 1.5919, 1.5973,
                1.5993, 1.6002, 1.6015]
# Table 5.1, Fundamentals of Heat and Mass Transfer, F. Incropera, 2006.


_Bi_data = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
            0.15, 0.20, 0.25, 0.30, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2, 3,
            4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 100]

_C1_cyl_fin = interp1d(_Bi_data, _C1_cyl_data)
# Linear interpolation for finite Biot numbers

_zeta1_cyl_fin = interp1d(_Bi_data, _zeta1_cyl_data)
# Linear interpolation for finite Biot numbers


def C1_cyl(Bi_):
    """
    Calculate first term C1 coefficient for infinite cylinder.

    Parameters
    ----------
    Bi_ : Biot number

    Returns
    -------
    C1 for infinite cylinder
    """
    if Bi_ > 100:
        C1 = 1.6018
        # Table 5.1, Fundamentals of Heat and Mass Transfer,
        # F. Incropera, 2006.
    else:
        C1 = _C1_cyl_fin(Bi_)
    return C1


def zeta1_cyl(Bi_):
    """
    Calculate first term zeta1 coefficient for infinite cylinder.

    Parameters
    ----------
    Bi_ : Biot number

    Returns
    -------
    zeta1 for infinite cylinder
    """
    if Bi_ > 100:
        zeta1 = 2.4050
        # Table 5.1, Fundamentals of Heat and Mass Transfer,
        # F. Incropera, 2006.
    else:
        zeta1 = _zeta1_cyl_fin(Bi_)
    return zeta1


def Fo_cyl(theta, Bi_):
    """
    Calculate Fourier number for infinite cylinder using approximate solution.
    Approximate solution is applicable when the solid has uniform temperature.

    Parameters
    ----------
    theta : dimensionless temperature difference
    Bi_ : Biot number

    Returns
    -------
    Fourier number, dimensionless
    """
    zeta1 = zeta1_cyl(Bi_)
    C1 = C1_cyl(Bi_)
    Fo_ = -1 / zeta1**2 * log(theta/C1)
    return Q_(Fo_, ureg.dimensionless)


def alpha(k, rho, C):
    """
    Calculate thermal diffusivity.

    Parameters
    ----------
    k : thermal conductivity of the solid
    rho : density of the solid
    C : specific heat capacity

    Returns
    -------
    thermal diffusivity
    """
    alpha_ = k / (rho*C)
    return alpha_.to(ureg.m**2/ureg.s)


def theta_temp(T, T_i, T_inf):
    """
    Calculate dimensionless temperature difference. Used for transient
    conduction and convection.

    Parameters
    ----------
    T : variable temperature of the solid
    T_i : initially uniform temperature of the solid
    T_inf : temperature of the medium

    Returns
    -------
    temperature difference, dimensionless
    """
    theta_temp_ = (T-T_inf) / (T_i-T_inf)
    return theta_temp_.to_base_units()


def _nist_log_fit(T, coefs):
    """
    Calculate NIST curve fit for given coefficients.
    https://trc.nist.gov/cryogenics/materials/materialproperties.htm

    Parameters
    ----------
    T : temperature, K
    coefs : coefficients from NIST cryo properties database

    Returns
    -------
    thermal property (e.g. thermal conductivity)
    """
    y = 0
    for ind, coef in enumerate(coefs):
        y += coef*log10(T)**ind
    return 10**y


def _nist_pow_fit(T, coefs):
    """
    Calculate NIST curve fit for given coefficients.
    https://trc.nist.gov/cryogenics/materials/materialproperties.htm

    Parameters
    ----------
    T : temperature, K
    coefs : coefficients from NIST cryo properties database

    Returns
    -------
    thermal property (e.g. linear expansion)
    """
    y = 0
    for ind, coef in enumerate(coefs):
        y += coef*T**ind
    return y


def _nist_Cu_log_fit(T, coefs):
    """
    Calculate NIST curve fit for copper thermal conductivity for
    given coefficients.
    https://trc.nist.gov/cryogenics/materials/materialproperties.htm

    Parameters
    ----------
    T : temperature, K
    coefs : coefficients from NIST cryo properties database

    Returns
    -------
    thermal property (e.g. thermal conductivity)
    """
    a, b, c, d, e, f, g, h, i = coefs
    y = (a + c*T**0.5 + e*T + g*T**1.5 + i*T**2) / \
        (1 + b*T**0.5 + d*T + f*T**1.5 + h*T**2)
    return 10**y


def _nist_quad(T1, T2, fun, coefs):
    """
    Calculate average value of the property for given temperature range.
    https://trc.nist.gov/cryogenics/materials/materialproperties.htm

    Parameters
    ----------
    T : Quantity, temperature
    coefs : coefficients from NIST cryo properties database

    Returns
    -------
    thermal property (e.g. thermal conductivity)
    """
    return quad(fun, T1, T2, args=coefs)[0] / (T2-T1)


def nist_property(material, prop, T1, T2=None, RRR_OFHC=None):
    """
    Calculate specific heat capacity using NIST properties database.
    https://trc.nist.gov/cryogenics/materials/materialproperties.htm

    Parameters
    ----------
    T1 : Quantity, temperature
    T2 : Quantity, temperature
        used for average value calculation

    Returns
    -------
    specific heat capacity
    """
    if material == 'OFHC':
        if RRR_OFHC is None:
            logger.warning('RRR for OFHC is not defined. Using RRR=100.')
            RRR_OHFC = 100
        RRR = str(RRR_OFHC)
        coefs = NIST_DATA[material][prop]['coefs'+RRR]
    else:
        coefs = NIST_DATA[material][prop]['coefs']
    fun = NIST_DATA[material][prop]['fun']
    unit = NIST_DATA[material][prop]['unit']
    eq_range = NIST_DATA[material][prop]['range']

    T1 = T1.to(ureg.K).magnitude
    if prop == 'LE':
        Tlow = NIST_DATA[material][prop]['Tlow']
        if T1 < Tlow:
            f = NIST_DATA[material][prop]['f']
            return f * unit
    if T1 < eq_range[0] or T1 > eq_range[1]:
        raise ValueError(f'Temperature is out of bounds: {T1} for'
                         f' {eq_range[0]}-{eq_range[1]}'
                         'limits.')
    if T2 is None:
        result = fun(T1, coefs)
    else:
        T2 = T2.to(ureg.K).magnitude
        result = _nist_quad(T1, T2, fun, coefs)
    return result * unit


# Temporary storage for NIST data
# Materials:
# 304SS - AISI 304 Stainless Steel
# 6061 - 6061-T6 Aluminum (UNS A96061)
# G10 - G10
# PTFE - PTFE/Teflon
# OFHC - Oxygen-free High thermal conductivity copper
#
# Properties
# TC - thermal conductivity, W/(m*K)
# SH - specific heat, W/(m*K)
# EC - expansion coefficient, 1/K
# LE - linear expansion, dimensionless
NIST_DATA = {
    # TODO Implement different TC, etc functions for the same material
    '304SS':
    {
        'TC': {'coefs': [-1.4087, 1.3982, 0.2543, -0.6260, 0.2334, 0.4256,
                         -0.4658, 0.1650, -0.0199],
               'range': (1, 300),
               'fun': _nist_log_fit,
               'unit': ureg.W/(ureg.m*ureg.K)},
        'LE': {'coefs': [-2.9554E2, -3.9811E-1, 9.2683E-3, -2.0261E-5,
                         1.7127E-8],
               'range': (4, 300),
               'fun': _nist_pow_fit,
               'Tlow': 23,
               'f': -300.04,
               'unit': 1e-5*ureg.m/ureg.m},
        'SH': {'coefs': [22.0061, -127.5528, 303.647, -381.0098, 274.0328,
                         -112.9212, 24.7593, -2.239153, 0],
               'range': (4, 300),
               'fun': _nist_log_fit,
               'unit': ureg.J/(ureg.kg*ureg.K)}
    },
    '6061':
    {
        'TC': {'coefs': [0.07918, 1.0957, -0.07277, 0.08084, 0.02803, -0.09464,
                         0.04179, -0.00571, 0],
               'range': (1, 300),
               'fun': _nist_log_fit,
               'unit': ureg.W/(ureg.m*ureg.K)},
        'SH': {'coefs': [46.6467, -314.292, 866.662, -1298.3, 1162.27, -637.795,
                         210.351, -38.3094, 2.96344],
               'range': (4, 300),
               'fun': _nist_log_fit,
               'unit': ureg.J/(ureg.kg*ureg.K)}
    },
    'G10':
    {
        'TC': {'coefs': [-4.1236, 13.788, -26.068, 26.272, -14.663, 4.4954,
                         -0.6905, 0.0397, 0],
               'range': (10, 300),
               'fun': _nist_log_fit,
               'unit': ureg.W/(ureg.m*ureg.K)},
        'SH': {'coefs': [-2.4083, 7.6006, -8.2982, 7.3301, -4.2386, 1.4294,
                         -0.24396, 0.015236, 0],
               'range': (4, 300),
               'fun': _nist_log_fit,
               'unit': ureg.J/(ureg.kg*ureg.K)}
    },
    'PTFE':
    {
        'TC': {'coefs': [2.7380, -30.677, 89.430, -136.99, 124.69, -69.556,
                         23.320, -4.3135, 0.33829],
               'range': (4, 300),
               'fun': _nist_log_fit,
               'unit': ureg.W/(ureg.m*ureg.K)},
        'SH': {'coefs': [31.88256, -166.51949, 352.01879, -393.44232, 259.98072,
                         -104.61429, 24.99276, -3.20792, 0.16503],
               'range': (4, 300),
               'fun': _nist_log_fit,
               'unit': ureg.J/(ureg.kg*ureg.K)}
    },
    'OFHC':
    {
        'EC': {'coefs': [-17.9081289, 67.131914, -118.809316, 109.9845997,
                         -53.8696089, 13.30247491, -1.30843441],
               'range': (4, 300),
               'fun': _nist_log_fit,
               'unit': 1e-6/ureg.K},
        'TC': {'coefs50': [1.8743, -0.41538, -0.6018, 0.13294, 0.26426, -0.0219,
                           -0.051276, 0.0014871, 0.003723],
               'coefs100': [2.2154, -0.47461, -0.88068, 0.13871, 0.29505,
                            -0.02043, -0.04831, 0.001281, 0.003207],
               'coefs150': [2.3797, -0.4918, -0.98615, 0.13942, 0.30475,
                            -0.019713, -0.046897, 0.0011969, 0.0029988],
               'coefs300': [1.357, 0.3981, 2.669, -0.1346, -0.6683, 0.01342,
                            0.05773, 0.0002147, 0],
               'coefs500': [2.8075, -0.54074, -1.2777, 0.15362, 0.36444,
                            -0.02105, -0.051727, 0.0012226, 0.0030964],
               'range': (4, 300),
               'fun': _nist_Cu_log_fit,
               'unit': ureg.W/(ureg.m*ureg.K)},
    },
}


def stored_energy(fluid, volume):
    """Calculate stored energy in volume using Baker equation."""
    P = fluid.P
    V = volume
    k = fluid.gamma
    E_stored = P * V / (k-1) * (1-(P_NTP/P)**((k-1)/k))
    return E_stored.to(ureg.lbf*ureg.ft)


def blast_radius(E_stored):
    """Calculate maximum distance for debris, eardrum rupture and
    lung damage based on PNNL paper."""
    W_TNT = E_stored / E_TNT  # Energy equivalent in TNT
    D_1 = z_1 * (W_TNT.to(ureg.kg).magnitude)**(1/3)
    D_2 = z_2 * (W_TNT.to(ureg.kg).magnitude)**(1/3)
    D_3 = z_3 * (W_TNT.to(ureg.kg).magnitude)**(1/3)
    return (D_1, D_2, D_3)
