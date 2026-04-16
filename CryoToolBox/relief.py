"""Relief and discharge calculations with backward-compatible aliases.

This module is the preferred home for safety and relief helpers. It gathers:

- API-style sizing and discharge calculations,
- CGA S-1.3 flow-rating and capacity helpers,
- nozzle/direct-integration helpers used in relief work.

Legacy imports from :mod:`functions`, :mod:`cga`, and :mod:`piping` remain
available through compatibility aliases in those modules.
"""

from math import e, pi
import warnings

from scipy import optimize
from scipy.integrate import quad

from . import logger
from .constants import AIR
from .cp_wrapper import CP_const_unit
from .std_conditions import P_NTP, Q_, ureg


def crit_discharge_m_dot(fluid, A):
    """Calculate critical discharge mass flow.

    Reference: White, *Fluid Mechanics*, eq. 9.46a.
    """
    C = fluid.C_gas_const
    P = fluid.P
    m_max_ = C * A * P * fluid.MZT
    return m_max_.to_base_units()


def PRV_flow(A, Kd, fluid):
    """Deprecated helper for actual discharge mass flow through a relief valve."""
    warnings.warn(
        '`PRV_flow()` is deprecated and may be removed in a future release. '
        'Use `crit_discharge_m_dot()` and apply the discharge coefficient in '
        'user code instead.',
        DeprecationWarning,
        stacklevel=2,
    )
    W_T = crit_discharge_m_dot(fluid, A)
    W_a = W_T * Kd
    return W_a.to(ureg.g / ureg.s)


def api_relief_area(m_dot, fluid, *, P_back=P_NTP, K_d=0.975, K_b=1, K_c=1):
    """Calculate required relief area as per API 520 5.6.3/4."""
    W = m_dot.m_as(ureg.lb / ureg.hr)
    P_1 = fluid.P.m_as(ureg.psi)
    P_2 = P_back.m_as(ureg.psi)
    k = fluid.gamma
    T = fluid.T.m_as(ureg.degR)
    Z = fluid.compressibility_factor
    M = fluid.M
    P_cf = P_1 * (2 / (k + 1)) ** (k / (k - 1))
    if P_2 <= P_cf:
        C = fluid.C.m_as(ureg.lb / (ureg.hr * ureg.lbf) * (ureg.degR) ** 0.5)
        A = W / (C * K_d * P_1 * K_b * K_c) * (T * Z / M) ** 0.5
    else:
        r = P_2 / P_1
        F_2_brack = (1 - r ** ((k - 1) / k)) / (1 - r)
        F_2 = (k / (k - 1) * r ** (2 / k) * F_2_brack) ** 0.5
        A = W / (735 * F_2 * K_d * K_c) * (T * Z / (M * P_1 * (P_1 - P_2))) ** 0.5
    return A * ureg.inch ** 2


def cga_flow_rating_pres(P_set, factor=0.1):
    """Calculate flow rating pressure as per CGA S-1.3 2008 5.1.13."""
    return ((1 + factor) * P_set.to(ureg.psi)).to(ureg.psig)


def cga_flow_calc_temp(fluid):
    """Calculate CGA temperature for flow capacity calculation."""
    temp_state = fluid.copy()
    if not fluid.is_super_critical:
        temp_state.update_kw(P=temp_state.P, Q=0)
        return temp_state.T
    T2_ = 1000

    def func(T_):
        temp_state.update_kw(P=temp_state.P, T=T_ * ureg.K)
        spec_vol = 1 / temp_state.Dmass
        y = -(spec_vol ** 0.5) / temp_state.specific_heat_input
        return y.m_as(ureg('0 m/J * (m*kg)**0.5'))

    T_relief = optimize.fminbound(func, fluid.T_min.m_as(ureg.K), T2_)
    return T_relief * ureg.K


def cga_flow_rating_state(fluid, factor=0.1):
    """Calculate flow-rating temperature and pressure for relief sizing."""
    P_flow = cga_flow_rating_pres(fluid.P, factor=factor)
    fluid_FR = fluid.copy()
    fluid_FR.update_kw(P=P_flow, Hmass=fluid.Hmass)
    if not fluid_FR.is_super_critical:
        fluid_FR.update_kw(P=fluid_FR.P, Q=1)
    else:
        T_flow = cga_flow_calc_temp(fluid_FR)
        fluid_FR.update_kw(P=P_flow, T=T_flow)
    return fluid_FR


def cga_factor_F(fluid_FR, fluid_i):
    """Calculate correction factor F per CGA S-1.3 6.1.4."""
    if fluid_FR.name != fluid_i.name:
        raise TypeError('Both states have to be of the same fluid.')
    Zi = fluid_i.Z
    Ti = fluid_i.T
    Z = fluid_FR.Z
    T = fluid_FR.T
    F_ = (Zi * Ti / (Z * T)) ** 0.5
    return float(F_)


def insul_prim_req_flow(fluid_FR, U, A, F=1, conservative=True):
    """Calculate minimum required primary PRD capacity for insulated vessels."""
    T = fluid_FR.T.to(ureg.degR).magnitude
    G_i_ = cga_factor_Gi(fluid_FR, conservative=conservative)
    U = U.to(ureg.BTU / (ureg.hr * ureg.ft ** 2 * ureg.degR)).magnitude
    A = A.to(ureg.ft ** 2).magnitude
    Q_a = (590 - T) / (4 * (1660 - T)) * F * G_i_ * U * A
    return Q_a * ureg.ft ** 3 / ureg.min


def fire_liq_req_flow(fluid_FR, A, F=1):
    """Calculate required relief capacity for fire exposure cases."""
    return F * cga_factor_Gu(fluid_FR) * A.m_as(ureg.ft ** 2) ** 0.82 * ureg.ft ** 3 / ureg.min


def cga_factor_Gi(fluid_FR, conservative=True):
    """Calculate the CGA ``G_i`` gas factor."""
    L_adj = _calculate_Gi_Gu_heat(fluid_FR)
    T = cga_flow_calc_temp(fluid_FR)
    C = fluid_FR.C_gas_const
    TZM = 1 / fluid_FR.MZT
    if conservative and fluid_FR.name.lower() == 'helium':
        return 52.5
    return _G_i_us(T, C, L_adj, TZM)


def cga_factor_Gu(fluid_FR):
    """Calculate the CGA ``G_u`` gas factor."""
    L_adj = _calculate_Gi_Gu_heat(fluid_FR)
    C = fluid_FR.C_gas_const
    TZM = 1 / fluid_FR.MZT
    return _G_u_us(C, L_adj, TZM)


def _calculate_Gi_Gu_heat(fluid_FR):
    if fluid_FR.P >= fluid_FR.P_critical:
        L = fluid_FR.specific_heat_input
    elif fluid_FR.P >= 0.4 * fluid_FR.P_critical:
        temp_state = fluid_FR.copy()
        temp_state.update_kw(P=temp_state.P, Q=0)
        v_l = 1 / temp_state.Dmass
        h_liq = temp_state.latent_heat
        temp_state.update_kw(P=temp_state.P, Q=1)
        v_g = 1 / temp_state.Dmass
        L = h_liq * v_g / (v_g - v_l)
    else:
        temp_state = fluid_FR.copy()
        temp_state.update_kw(P=temp_state.P, Q=0)
        L = temp_state.latent_heat
    return L


@ureg.wraps(None, (ureg.degR, CP_const_unit['C_us'][1], ureg.BTU / ureg.lb, ureg.degR ** 0.5))
def _G_i_us(T, C, L, TZM):
    return float(73.4 * (1660 - T) / (C * L) * TZM)


@ureg.wraps(None, (CP_const_unit['C_us'][1], ureg.BTU / ureg.lb, ureg.degR ** 0.5))
def _G_u_us(C, L, TZM):
    return float(633_000 / (C * L) * TZM)


def relief_inlet_temp(fluid_FR, m_dot, pipe, condition='fire'):
    """Calculate temperature at the relief inlet per CGA S-1.3 6.1.4 b)."""
    Ts_ = fluid_FR.T.m_as(ureg.K)
    if condition == 'fire':
        Th_ = 1192
        coef = 1285
    else:
        Th_ = 328
        coef = 134.8
    D = pipe.OD.m_as(ureg.m)
    L = pipe.L.m_as(ureg.m)
    W = m_dot.m_as(ureg.kg / ureg.hr)
    if 0 <= fluid_FR.Q < 1:
        fluid_FR = fluid_FR.copy()
        fluid_FR.update_kw(P=fluid_FR.P, Q=1)
    Cp = _average_spec_heat(fluid_FR).m_as(ureg.kJ / ureg.kg / ureg.K)
    Ti = Th_ - (Th_ - Ts_) / (e * coef * D * L / (W * Cp))
    Ti *= ureg.K
    return Ti


def _average_spec_heat(fluid, condition='fire'):
    if condition == 'fire':
        Th_ = 922
    else:
        Th_ = 328

    def specific_heat(T_):
        fluid_temp = fluid.copy()
        T = T_ * ureg.K
        fluid_temp.update_kw(P=fluid.P, T=T)
        return fluid_temp.Cpmass.m_as(ureg.J / (ureg.kg * ureg.K))

    spec_heat_ave = quad(specific_heat, fluid.T.m_as(ureg.K), Th_)[0]
    spec_heat_ave *= ureg.J / (ureg.kg * (Th_ * ureg.K - fluid.T))
    return spec_heat_ave


def equiv_orifice_ID(m_dot, dP, fluid=AIR):
    """Calculate ID for the equivalent square-edge orifice."""
    Cd = 0.61
    rho = fluid.Dmass
    ID = 2 * (m_dot / (pi * Cd * (2 * dP * rho) ** 0.5)) ** 0.5
    return ID.to(ureg.inch)


def crit_discharge_mass_flux(fluid, P_out=P_NTP, n_steps=100):
    """Calculate mass flux through a converging nozzle by direct integration."""
    fluid_temp = fluid.copy()
    P1_ = fluid.P.m_as(ureg.Pa)
    P2_ = P_out.m_as(ureg.Pa)
    S = fluid.Smass

    def v(P_):
        P = P_ * ureg.Pa
        fluid_temp.update_kw(P=P, Smass=S)
        rho = fluid_temp.Dmass.m_as(ureg.kg / ureg.m ** 3)
        return 1 / rho

    dP = (P1_ - P2_) / n_steps
    P_ = P1_
    G_max = 0
    while P_ > P2_:
        try:
            G = 1 / v(P_) * (-2 * quad(v, P1_, P_)[0]) ** 0.5
        except ValueError:
            P_ -= dP
            continue
        if G < G_max:
            break
        G_max = G
        P_ -= dP
    G_max *= ureg.kg / (ureg.s * ureg.m ** 2)
    return G_max.to_base_units()


# Preferred short descriptive names above; legacy compatibility aliases below.
m_max = crit_discharge_m_dot
A_relief_API = api_relief_area
P_FR = cga_flow_rating_pres
theta = cga_flow_calc_temp
cga_theta = cga_flow_calc_temp
calculate_fluid_FR = cga_flow_rating_state
F = cga_factor_F
G_i = cga_factor_Gi
G_u = cga_factor_Gu
primary_insulated = insul_prim_req_flow
relief_fire_liquefied = fire_liq_req_flow
calculate_inlet_temp = relief_inlet_temp
equivalent_orifice = equiv_orifice_ID
G_nozzle = crit_discharge_mass_flux


__all__ = [
    'A_relief_API',
    'F',
    'G_i',
    'G_nozzle',
    'G_u',
    'P_FR',
    'PRV_flow',
    'api_relief_area',
    'calculate_fluid_FR',
    'calculate_inlet_temp',
    'cga_factor_F',
    'cga_factor_Gi',
    'cga_factor_Gu',
    'cga_flow_calc_temp',
    'cga_flow_rating_pres',
    'cga_flow_rating_state',
    'cga_theta',
    'crit_discharge_m_dot',
    'crit_discharge_mass_flux',
    'equiv_orifice_ID',
    'equivalent_orifice',
    'fire_liq_req_flow',
    'insul_prim_req_flow',
    'm_max',
    'relief_inlet_temp',
    'primary_insulated',
    'relief_fire_liquefied',
    'theta',
]
