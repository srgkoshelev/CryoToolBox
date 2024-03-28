"""
Pressure relief calculations for CGA S-1.3.
2008 ed.
"""

from . import logger
from .std_conditions import ureg, Q_
from .cp_wrapper import CP_const_unit
from scipy import optimize


def P_fr(P_set, factor=1.1):
    """Calculate flow rating pressure as per CGA S-1.3 2008 5.1.13.

    Parameters:
    -----------
    P_set : Set pressure
    factor : pressure increase factor as per BPVC VIII Div. 1, UG-125
        Common factors: 1.1, 1.21 (fire)

    Returns:
    --------
    P_fr : flow rating pressure
    """
    return (factor * P_set.to(ureg.psi)).to(ureg.psig)


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


def G_i(fluid_FR, conservative=True):
    """Calculate gas factor for insulated containers per
    Notes to Table 1 and Table 2.

    fluid_FR: fluid at flow rating temperature and pressure
    """

    if fluid_FR.P >= fluid_FR.P_critical:
        L = fluid_FR.specific_heat_input  # L is replaced by theta
    elif fluid_FR.P >= 0.4*fluid_FR.P_critical:
        temp_state = fluid_FR.copy()
        temp_state.update_kw(P=temp_state.P, Q=0*ureg.dimensionless)
        v_l = 1/temp_state.Dmass
        temp_state.update_kw(P=temp_state.P, Q=1*ureg.dimensionless)
        v_g = 1/temp_state.Dmass
        L = fluid_FR.latent_heat * (v_g-v_l)/v_g
    else:
        L = fluid_FR.latent_heat
    L, T = theta(fluid_FR)
    C = fluid_FR.C_gas_const
    TZM = 1 / fluid_FR.MZT  # sqrt(T*Z/M)
    # Conservative value for He is 52.5 for P <= 200 psig
    if conservative and fluid_FR.name.lower() == 'helium':
        return 52.5
    else:
        return _G_i_us(T, C, L, TZM)


@ureg.wraps(None, (ureg.degR,
                   CP_const_unit['C_us'][1],
                   ureg.BTU/ureg.lb,
                   ureg.degR**0.5))
def _G_i_us(T, C, L, TZM):
    """Calculate G_i factor in US customary units.
    While the actual value has units, the quantity returned as dimensionless.
    """
    return 73.4 * (1660-T) / (C*L) * TZM
