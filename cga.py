"""
Pressure relief calculations for CGA S-1.3.
2008 ed.
"""

from . import ureg, Q_
from . import logger
from .cp_wrapper import CP_const_unit


def theta(Fluid, step=0.01):
    """
    Calculate latent heat/specific heat input and temperature for flow
    capacity calculation per CGA S-1.3 2008 6.1.3.

    :Fluid: ThermState object describing thermodynamic state (fluid, T, P)
    :step: temperature step
    :returns: tuple (specific heat, temperature)
    """
    TempState = Fluid.copy()
    # Only working for pure fluids and pre-defined mixtures
    if Fluid.is_super_critical:
        logger.warning(f'{Fluid.name} is supercritical at {Fluid.P:.3~g}'
                       f' and {Fluid.T:.3g~}. Specific heat input will be used.')
        TempState.update('P', Fluid.P, 'T', Fluid.T_min)
        T_start = TempState.T
        T_end = 300*ureg.K
        cga_criterion = ureg('0 m/J * (m*kg)**0.5')
        # CGA criterion for calculating temperature of specific heat input
        T = T_start
        while T <= T_end:
            spec_vol = 1/TempState.Dmass  # specific volume
            cga_criterion_new = (spec_vol**0.5)/TempState.specific_heat_input
            if cga_criterion_new > cga_criterion:
                cga_criterion = cga_criterion_new
            else:
                break  # Function has only one maximum
            T += step*ureg.K
            TempState.update('P', TempState.P, 'T', T)
        return (TempState.specific_heat_input, T)
    else:
        logger.warning(f'{Fluid.name} is subcritical at {Fluid.P:.3~g}. \
        Latent heat of evaporation will be used')
        latent_heat = Fluid.latent_heat
        return (latent_heat, TempState.T)  # saturation temperature


def F(Fluid_FR, Fluid_PRD):
    """Calculate Correction factor F per 6.1.4.

    Fluid_FR: Fluid at flow rating temperature and pressure
    Fluid_PRD: Fluid at the inlet of PRD
    """
    if Fluid_FR.name != Fluid_PRD.name:
        raise TypeError('Both states have to be of the same fluid.')
    P_i = Fluid_PRD.P
    v_i = 1 / Fluid_PRD.Dmass  # Specific volume
    P = Fluid_FR.P
    v = 1 / Fluid_FR.v
    F_ = (P_i*v_i/(P*v))**0.5
    return F_.to(ureg.dimensionless)


def primary_insulated(Fluid_FR, U, A, F=1, conservative=True):
    """Calculate minimum required flow capacity for primary PRD on insulated
    containers for liquefied compressed gases, refrigerated fluids, and
    refrigerated (cryogenic) fluids per 6.2.2.
    """
    T = Fluid_FR.T.to(ureg.degR).magnitude
    G_i_ = G_i(Fluid_FR, conservative=conservative)
    U = U.to(ureg.BTU/(ureg.hr*ureg.ft**2*ureg.degR)).magnitude
    A = A.to(ureg.ft**2).magnitude
    Q_a = (590-T) / (4*(1660-T)) * F * G_i_ * U * A
    return Q_a * ureg.ft**3 / ureg.min


def G_i(Fluid_FR, conservative=True):
    """Calculate gas factor for insulated containers per
    Notes to Table 1 and Table 2.

    Fluid_FR: Fluid at flow rating temperature and pressure
    """

    if Fluid_FR.P >= Fluid_FR.P_critical:
        L = Fluid_FR.specific_heat_input  # L is replaced by theta
    elif Fluid_FR.P >= 0.4*Fluid_FR.P_critical:
        TempState = Fluid_FR.copy()
        TempState.update_kw(P=TempState.P, Q=0*ureg.dimensionless)
        v_l = 1/TempState.Dmass
        TempState.update_kw(P=TempState.P, Q=1*ureg.dimensionless)
        v_g = 1/TempState.Dmass
        L = Fluid_FR.latent_heat * (v_g-v_l)/v_g
    else:
        L = Fluid_FR.latent_heat
    L, T = theta(Fluid_FR)
    C = Fluid_FR.C_gas_constant
    TZM = 1 / Fluid_FR.MZT  # sqrt(T*Z/M)
    # Conservative value for He is 52.5 for P <= 200 psig
    if conservative and Fluid_FR.name.lower() == 'helium':
        return 52.5
    else:
        return _G_i_us(T, C, L, TZM)

@ureg.wraps(None, (ureg.degR,
                   CP_const_unit['C_gas_constant'][1],
                   ureg.BTU/ureg.lb,
                   ureg.degR**0.5))
def _G_i_us(T, C, L, TZM):
    """Calculate G_i factor in US customary units.
    While the actual value has units, the quantity returned as dimensionless.
    """
    return 73.4 * (1660-T) / (C*L) * TZM
