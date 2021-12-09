"""
Pressure relief calculations for CGA S-1.3.
2008 ed.
"""

from . import ureg, Q_
from . import logger
from .cp_wrapper import CP_const_unit


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


def theta(fluid, step=0.01):
    """
    Calculate latent heat/specific heat input and temperature for flow
    capacity calculation per CGA S-1.3 2008 6.1.3.

    fluid : ThermState object describing thermodynamic state (fluid, T, P)
    step : temperature step
    returns : tuple (specific heat, temperature)
    """
    temp_state = fluid.copy()
    # Only working for pure fluids and pre-defined mixtures
    if fluid.is_super_critical:
        logger.warning(f'{fluid.name} is supercritical at {fluid.P:.3~g}'
                       f' and {fluid.T:.3g~}. Specific heat input will be used.')
        temp_state.update('P', fluid.P, 'T', fluid.T_min)
        T_start = temp_state.T
        T_end = 300*ureg.K
        cga_criterion = ureg('0 m/J * (m*kg)**0.5')
        # CGA criterion for calculating temperature of specific heat input
        T = T_start
        while T <= T_end:
            spec_vol = 1/temp_state.Dmass  # specific volume
            cga_criterion_new = (spec_vol**0.5)/temp_state.specific_heat_input
            if cga_criterion_new > cga_criterion:
                cga_criterion = cga_criterion_new
            else:
                break  # Function has only one maximum
            T += step*ureg.K
            temp_state.update('P', temp_state.P, 'T', T)
        return (temp_state.specific_heat_input, T)
    else:
        logger.warning(f'{fluid.name} is subcritical at {fluid.P:.3~g}. \
        Latent heat of evaporation will be used')
        latent_heat = fluid.latent_heat
        return (latent_heat, temp_state.T)  # saturation temperature


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
