"""
Pressure relief calculations for CGA S-1.3.
2008 ed.
"""

from . import logger
from .std_conditions import ureg, Q_
from .cp_wrapper import CP_const_unit
from scipy import optimize
from scipy.integrate import quad
from math import e


def P_FR(P_set, factor=0.1):
    """Calculate flow rating pressure as per CGA S-1.3 2008 5.1.13.

    Parameters:
    -----------
    P_set : Set pressure
    factor : pressure increase factor as per BPVC VIII Div. 1, UG-125
        Common factors: 1.1, 1.21 (fire)

    Returns:
    --------
    P_FR : flow rating pressure
    """
    return ((1+factor) * P_set.to(ureg.psi)).to(ureg.psig)


def theta(fluid):
    """""Calculate temperature for flow capacity calculation per
    CGA S-1.3 2008 6.1.3.

    Parameters
    ----------
    fluid : ThermState
        Object describing thermodynamic state (fluid, T, P).

    Returns
    -------
    Quantity {temperature: 1}
    """
    temp_state = fluid.copy()
    if not fluid.is_super_critical:
        temp_state.update_kw(P=temp_state.P, Q=0)
        return temp_state.T
    T2_ = 1000      # K
    def func(T_):
        temp_state.update_kw(P=temp_state.P, T=T_*ureg.K)
        spec_vol = 1/temp_state.Dmass  # specific volume
        y = -(spec_vol**0.5)/temp_state.specific_heat_input
        return y.m_as(ureg('0 m/J * (m*kg)**0.5'))
    T_relief = optimize.fminbound(func, fluid.T_min.m_as(ureg.K), T2_)
    return T_relief * ureg.K


def calculate_fluid_FR(fluid, factor=0.1):
    """Calculate flow rating temperature and pressure.
    References: CGA S-1.3 2008 5.1.13 and 6.1.3.
    Parameters:
    -----------
    fluid : ThermState
        Fluid at design temperature and pressure.
    Returns:
    --------
    ThermState
        Fluid at flow rating pressure and temperature for flow capacity calculation.
    """
    P_flow = P_FR(fluid.P, factor=factor)
    fluid_FR = fluid.copy()
    fluid_FR.update_kw(P=P_flow, Hmass=fluid.Hmass)
    if not fluid_FR.is_super_critical:
        fluid_FR.update_kw(P=fluid_FR.P, Q=1)
    else:
        T_flow = theta(fluid_FR)
        fluid_FR.update_kw(P=P_flow, T=T_flow)
    return fluid_FR




def F(fluid_FR, fluid_PRD):
    """Calculate Correction factor F per 6.1.4.

    fluid_FR: fluid at flow rating temperature and pressure
    fluid_PRD: fluid at the inlet of PRD
    """
    if fluid_FR.name != fluid_PRD.name:
        raise TypeError('Both states have to be of the same fluid.')
    P_i = fluid_PRD.P
    v_i = 1 / fluid_PRD.Dmass  # Specific volume
    P = fluid_FR.P
    v = 1 / fluid_FR.v
    F_ = (P_i*v_i/(P*v))**0.5
    return F_.to(ureg.dimensionless)


def primary_insulated(fluid_FR, U, A, F=1, conservative=True):
    """Calculate minimum required flow capacity for primary PRD on insulated
    containers for liquefied compressed gases, refrigerated fluids, and
    refrigerated (cryogenic) fluids per 6.2.2.
    """
    T = fluid_FR.T.to(ureg.degR).magnitude
    G_i_ = G_i(fluid_FR, conservative=conservative)
    U = U.to(ureg.BTU/(ureg.hr*ureg.ft**2*ureg.degR)).magnitude
    A = A.to(ureg.ft**2).magnitude
    Q_a = (590-T) / (4*(1660-T)) * F * G_i_ * U * A
    return Q_a * ureg.ft**3 / ureg.min


def relief_fire_liquefied(fluid_FR, A, F=1):
    """Calculate required relief capacity for liquefied compressed gases,
    refrigerated fluids, and refrigerated (cryogenic) fluids in uninsulated and
    insulated containers.
    CGA S-1.3 2008 6.3.2."""
    return F * G_u(fluid_FR) * A.m_as(ureg.ft**2)**0.82 * ureg.ft**3/ureg.min


def G_i(fluid_FR, conservative=True):
    """Calculate gas factor for insulated containers per
    Notes to Table 1 and Table 2.

    fluid_FR: fluid at flow rating temperature and pressure
    """

    L_adj = _calculate_Gi_Gu_heat(fluid_FR)
    T = theta(fluid_FR)
    C = fluid_FR.C_gas_const
    TZM = 1 / fluid_FR.MZT  # sqrt(T*Z/M)
    # Conservative value for He is 52.5 for P <= 200 psig
    if conservative and fluid_FR.name.lower() == 'helium':
        return 52.5
    else:
        return _G_i_us(T, C, L_adj, TZM)

def G_u(fluid_FR):
    """Calculate gas factor for insulated containers per
    Notes to Table 1 and Table 2.

    fluid_FR: fluid at flow rating temperature and pressure
    """
    L_adj = _calculate_Gi_Gu_heat(fluid_FR)
    C = fluid_FR.C_gas_const
    TZM = 1 / fluid_FR.MZT  # sqrt(T*Z/M)
    return _G_u_us(C, L_adj, TZM)

def _calculate_Gi_Gu_heat(fluid_FR):
    """Calculate adjusted value for latent heat/specific heat input based on
    Notes for Tables 1 and 2 of CGA S-1.3 2008."""
    if fluid_FR.P >= fluid_FR.P_critical:
        L = fluid_FR.specific_heat_input  # L is replaced by theta
    elif fluid_FR.P >= 0.4*fluid_FR.P_critical:
        temp_state = fluid_FR.copy()
        temp_state.update_kw(P=temp_state.P, Q=0*ureg.dimensionless)
        v_l = 1/temp_state.Dmass
        temp_state.update_kw(P=temp_state.P, Q=1*ureg.dimensionless)
        v_g = 1/temp_state.Dmass
        L = fluid_FR.latent_heat * v_g/(v_g-v_l)
    else:
        L = fluid_FR.latent_heat
    return L

@ureg.wraps(None, (ureg.degR,
                   CP_const_unit['C_us'][1],
                   ureg.BTU/ureg.lb,
                   ureg.degR**0.5))
def _G_i_us(T, C, L, TZM):
    """Calculate G_i factor in US customary units.
    While the actual value has units, the quantity returned as dimensionless.
    CGA S-1.3 Notes for Table 1 and Table 2.
    """
    return 73.4 * (1660-T) / (C*L) * TZM


@ureg.wraps(None, (CP_const_unit['C_us'][1],
                   ureg.BTU/ureg.lb,
                   ureg.degR**0.5))
def _G_u_us(C, L, TZM):
    """Calculate G_i factor in US customary units.
    While the actual value has units, the quantity returned as dimensionless.
    CGA S-1.3 Notes for Table 1 and Table 2.
    """
    return 633_000/ (C*L) * TZM

def calculate_inlet_temp(fluid_FR, m_dot, pipe):
    """Calculate temperature at the relief inlet per 6.1.4 b)."""
    # Ts = fluid_FR.T.m_as(ureg.K)
    # D = pipe.OD.m_as(u.m)
    # L = pipe.L.m_as(u.m)
    # W = m_dot.m_as(u.kg/u.hr)
    print(spec_heat_ave.to(ureg.kJ/(ureg.kg*ureg.K)))

    # Ti_fire = 1192-(1192-Ts)/(e*1285*D*L/(W*Cp))

def _average_spec_heat(fluid, condition='fire'):
    """Calculate average specific heat for temperature range."""
    if condition == 'fire':
        Th_ = 922  # K, 6.1.4 b)
    else:
        Th_ = 328  # K, 6.1.4 b)
    def specific_heat(T_):
        """Calculate specific heat for given temperature."""
        fluid_temp = fluid.copy()
        T = T_ * ureg.K
        fluid_temp.update_kw(P=fluid.P, T=T)
        return fluid_temp.Cpmass.m_as(ureg.J/(ureg.kg*ureg.K))
    spec_heat_ave = quad(specific_heat, fluid.T.m_as(ureg.K), Th_)[0]
    spec_heat_ave *= ureg.J/(ureg.kg*(Th_*ureg.K-fluid.T))
    return spec_heat_ave
