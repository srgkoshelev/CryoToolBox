"""Helpful thermodynamic functions.
"""

from math import log, log10, pi
from . import ureg, Q_
from . import Air
from . import P_NTP
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


def to_scfma(M_dot_fluid, Fluid):
    """
    Convert mass flow rate into equivalent flow of air.
    Flow through a relief device with invariant Area/discharge coefficient
    (KA).

    Parameters
    ----------
    M_dot_fluid : Quantitiy {mass: 1, time: -1}
        mass flow rate
    Fluid : ThermState

    Returns
    -------
    ThermState
        volumetric air flow rate
    """
    C_fluid = Fluid.C_gas_constant
    C_air = Air.C_gas_constant

    # Calculation
    M_dot_air = M_dot_fluid * C_air / C_fluid * Air.MZT / Fluid.MZT
    Q_air = M_dot_air / Air.Dmass
    Q_air.ito(ureg.ft**3/ureg.min)
    return Q_air


def from_scfma(Q_air, Fluid):
    """
    Convert volumetric air flow rate into equivalent mass flow of specified
    fluid. Flow through a relief device with invariant Area/discharge
    coefficient (KA).
    Invert function to to_scfma().

    Parameters
    ----------
    Q_air : Quantity {length: 3, time: -1}
        volumetric air flow rate
    Fluid : ThermState

    Returns
    -------
    Quantitiy {mass: 1, time: -1}
        mass flow rate
    """
    C_fluid = Fluid.C_gas_constant
    C_air = Air.C_gas_constant

    # Calculation
    M_dot_air = Q_air * Air.Dmass
    M_dot_fluid = M_dot_air * C_fluid / C_air * Fluid.MZT / Air.MZT
    M_dot_fluid.ito(ureg.g/ureg.s)
    return M_dot_fluid


def theta_heat(Fluid, step=0.01):
    logger.warning('Deprecated. Use ht.cga.theta() instead.')
    return cga.theta(Fluid, step)


def rad_hl(eps_cold=0.55, eps_hot=0.55, T_hot=300*ureg.K, T_cold=77*ureg.K,
           F1_2=1, eps_baffle=0.02, N_baffles=5):
    """
    Calculate radiative heat load including reduction due to baffles.
    Based on Kaganer "Thermal insulation in cryogenic engineering", p. 42.

    Parameters
    ----------
    eps_cold : float
        emissivity of the cold surface
    eps_hot : float
        emissivity of the hot surface
    T_hot : Quality {temperature: 1}
        temperature of the hot surface
    T_cold : Quality {temperature: 1}
        temperature of the cold surface
    F1_2 : float
        F1_2 = F_cold/F_hot
    eps_baffle : float
        emissivity of the baffle, assumed to be same on both sides
    N_baffles : int
        number of baffles

    Returns
    -------
    dict
            :q0: heat load without any baffles
            q_baffle : heat load with the baffles
            eta : effectiveness of the baffles
    """
    # TODO This function will be refactored
    Eps_mut = 1/(1/eps_cold + F1_2*(1/eps_hot-1))  # Mutual emissivity
    q0 = Eps_mut*sigma*(T_hot**4 - T_cold**4)*F1_2
    Eps_baffle_mut = eps_baffle/(2-eps_baffle)
    eta = (1+N_baffles*Eps_mut/Eps_baffle_mut)**(-1)
    q_baffle = eta*q0
    return {'q0': q0.to(ureg.W/ureg.m**2),
            'q_baffle': q_baffle.to(ureg.W/ureg.m**2),
            'eta': eta}


def Re(Fluid, m_dot, D):
    """
    Calculate Reynolds number.

    Parameters
    ----------
    Fluid : ThermState object describing thermodynamic state (fluid, T, P)
    m_dot : mass flow
    D : characteristic length/hydraulic diameter

    Returns
    -------
    Reynolds number, dimensionless
    """
    A = pi * D**2 / 4
    w_flow = m_dot / (Fluid.Dmass*A)
    Re_ = w_flow * D * Fluid.Dmass / Fluid.viscosity
    return Re_.to_base_units()


def Pr(Fluid):
    """
    Calculate Prandtl number.

    Parameters
    ----------
    Fluid : ThermState object describing thermodynamic state (fluid, T, P)

    Returns
    -------
    Prandtl number, dimensionless
    """
    return Fluid.Prandtl


def Gr(Fluid, T_surf, L_surf):
    """
    Calculate Grashof number.

    Parameters
    ----------
    Fluid : ThermState object describing thermodynamic state (fluid, T, P)
    T_surf : surface temperature
    L_surf : characteristic length

    Returns
    -------
    Grashof number, dimensionless
    """
    nu_fluid = Fluid.viscosity/Fluid.Dmass  # kinematic viscosity
    beta_exp = Fluid.isobaric_expansion_coefficient
    Gr_ = ureg.g_0 * L_surf**3 * beta_exp * abs(T_surf-Fluid.T) / nu_fluid**2
    return Gr_.to(ureg.dimensionless)


def Ra(Fluid, T_surf, L_surf):
    """
    Calculate Rayleigh number.

    Parameters
    ----------
    Fluid : ThermState object describing thermodynamic state (fluid, T, P)
    T_surf : surface temperature
    L_surf : characteristic length

    Returns
    -------
    Rayleigh number, dimensionless
    """
    return Gr(Fluid, T_surf, L_surf)*Fluid.Prandtl


def Nu_cyl_hor(Fluid, T_cyl, D_cyl):
    """
    Calculate Nusselt number for vertical cylinder.
    Only natural convection currently supported.
    Based on Handbook of heat transfer by Rohsenow, Hartnet,
    Cho (HHT).

    Parameters
    ----------
    Fluid : ThermState object describing thermodynamic state (fluid, T, P)
    T_cyl : surface temperature
    D_cyl : cylinder diameter

    Returns
    -------
    Nusselt number, dimensionless
    """
    Pr_ = Fluid.Prandtl
    Ra_ = Ra(Fluid, T_cyl, D_cyl)
    C_l = 0.671/(1+(0.492/Pr_)**(9/16))**(4/9)  # HHT (4.13)
    Nu_T = 0.772*C_l*Ra_**(1/4)  # HHT (4.45)
    f = 1-0.13/Nu_T**0.16
    Nu_l = 2*f/log(1+2*f*Nu_T)
    C_t = 0.0002*log(Pr_)**3 - 0.0027*log(Pr_)**2 + 0.0061*log(Pr_) + 0.1054
    Nu_t = 0.103*Ra_**(1/3)
    Nu_ = (Nu_l**10 + Nu_t**10)**(1/10)
    return Nu_.to(ureg.dimensionless)


def Nu_cyl_vert(Fluid, T_cyl, D_cyl, L_cyl):
    """
    Calculate Nusselt number for vertical cylinder.
    Only natural convection currently supported.
    Based on Handbook of heat transfer by Rohsenow, Hartnet,
    Cho (HHT).

    Parameters
    ----------
    Fluid : ThermState object describing thermodynamic state (fluid, T, P)
    T_cyl : surface temperature
    D_cyl : cylinder diameter
    L_cyl : cylinder length

    Returns
    -------
    Nusselt number, dimensionless
    """
    Pr_ = Fluid.Prandtl
    Ra_ = Ra(Fluid, T_cyl, D_cyl)
    C_l = 0.671/(1+(0.492/Pr_)**(9/16))**(4/9)  # HHT (4.13)
    C_t_vert = (0.13*Pr_**0.22)/(1+0.61*Pr_**0.81)**0.42  # HHT (4.24)
    Nu_T_plate = C_l * Ra_**0.25
    Nu_l_plate = 2 / log(1+2/Nu_T_plate)  # HHT (4.33)
    zeta = 1.8 * L_cyl / (D_cyl*Nu_T_plate)   # HHT (4.44)
    Nu_l = zeta / (log(1+zeta)*Nu_l_plate)
    Nu_t = C_t_vert*Ra_**(1/3)/(1+1.4e9*Pr_/Ra_)
    Nu_ = (Nu_l**6 + Nu_t**6)**(1/6)
    return Nu_.to(ureg.dimensionless)


def heat_trans_coef(Fluid, Nu, L_surf):
    """
    Calculate heat transfer coefficient.

    Parameters
    ----------
    Fluid : ThermState object describing thermodynamic state (fluid, T, P)
    Nu : Nusselt number
    L_surf : characteristic length:
        :Horizontal cylinder: L_surf = D_cyl
        :Vertical cylinder: L_surf = L_cyl

    Returns
    -------
    heat transfer coefficient
    """
    h = Fluid.conductivity * Nu / L_surf
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


def nist_curve_fit(T, NIST_coefs):  # TODO make hidden
    """
    Calculate NIST curve fit for given coefficients.
    https://trc.nist.gov/cryogenics/materials/materialproperties.htm

    Parameters
    ----------
    T : temperature, K
    NIST_coefs : coefficients from NIST cryo properties database

    Returns
    -------
    thermal property (e.g. thermal conductivity)
    """
    y = 0
    for ind, coef in enumerate(NIST_coefs):
        # print('abcdefghi'[ind], coef) #TODO add DEBUG
        y += coef*log10(T)**ind
    return 10**y


def nist_quad(T1, T2, NIST_coefs):
    """
    Calculate average value of the property for given temperature range.
    https://trc.nist.gov/cryogenics/materials/materialproperties.htm

    Parameters
    ----------
    T : Quantity, temperature
    NIST_coefs : coefficients from NIST cryo properties database

    Returns
    -------
    thermal property (e.g. thermal conductivity)
    """
    return quad(nist_curve_fit, T1, T2, args=NIST_coefs)[0] / (T2-T1)


def nist_property(material, prop, T1, T2=None):
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
    if prop == 'TC':
        output_unit = ureg.W / (ureg.m*ureg.K)
    elif prop == 'HC':
        output_unit = ureg.J / (ureg.kg*ureg.K)
    else:
        raise NotImplementedError('Only thermal conductivity (TC) and'
                                  ' heat capacity (HC) currently implemented.')
    T1 = T1.to(ureg.K).magnitude
    property_data = NIST_DATA[material][prop]
    if T1 < property_data[1][0] or T1 > property_data[1][1]:
        raise ValueError(f'Temperature is out of bounds: {T1} for'
                         f' {property_data[1][0]}-{property_data[1][1]}'
                         'limits.')
    if T2 is None:
        result = nist_curve_fit(T1, property_data[0])
    else:
        T2 = T2.to(ureg.K).magnitude
        result = nist_quad(T1, T2, property_data[0])
    return result * output_unit


# Temporary storage for NIST data
NIST_DATA = {
    '304SS':
    {
        'TC': ([-1.4087, 1.3982, 0.2543, -0.6260, 0.2334, 0.4256, -0.4658,
                0.1650, -0.0199], (1, 300)),
    },
    'G10':
    {
        'TC': ([-4.1236, 13.788, -26.068, 26.272, -14.663, 4.4954, -0.6905,
                0.0397, 0], (4, 300)),
    }
}



def stored_energy(Piping):
    """Calculate stored energy in piping using Baker equation."""
    P = Piping.Fluid.P
    V = Piping.volume
    k = Piping.Fluid.gamma
    E_stored = P * V / (k-1) * (1-(P_NTP/P)**((k-1)/k))
    return E_stored.to(ureg.lbf*ureg.ft)


def blast_radius(E_stored):
    """Calculate maximum distance for debris, eardrum rupture and
    lung damage based on PNNL paper."""
    W_TNT = E_stored / E_TNT  # Energy equivalent in TNT
    D_1 = z_1 * (W_TNT.to(ureg.kg).magnitude)**0.5
    D_2 = z_2 * (W_TNT.to(ureg.kg).magnitude)**0.5
    D_3 = z_3 * (W_TNT.to(ureg.kg).magnitude)**0.5
    return (D_1, D_2, D_3)
