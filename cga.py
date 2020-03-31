"""
Pressure relief calculations for CGA S-1.3.
2008 ed.
"""

from . import ureg, Q_
from . import logger
from .cp_wrapper import ThermState


def theta(Fluid, step=0.01):
    """
    Calculate latent heat/specific heat input and temperature for flow
    capacity calculation per CGA S-1.3 2008 6.1.3.

    :Fluid: ThermState object describing thermodynamic state (fluid, T, P)
    :step: temperature step
    :returns: tuple (specific heat, temperature)
    """
    TempState = ThermState(Fluid.name, backend=Fluid.backend)
    # Only working for pure fluids and pre-defined mixtures
    if Fluid.P > Fluid.P_critical:
        logger.warning(f'{Fluid.name} is supercritical at {Fluid.P:.3~g}. \
        Specific heat input will be used')
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
        TempState.update('P', Fluid.P, 'Q', Q_('0'))
        h_liquid = TempState.Hmass
        TempState.update('P', Fluid.P, 'Q', Q_('1'))
        h_vapor = TempState.Hmass
        latent_heat = h_vapor - h_liquid
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
    _F = (P_i*v_i/(P*v))**0.5
    return _F.to(ureg.dimensionless)


def primary_insulated(T, U, A, F=1):
    """Calculate minimum required flow capacity for primary PRD on insulated
    containers for liquefied compressed gases, refrigerated fluids, and
    refrigerated (cryogenic) fluids per 6.2.2.
    """
    Q_a = (Q_('590 degF')-T) / (4*(Q_('1660 degF')-T)) * F * G_i * U * A


def G_i(Fluid_FR):
    """Calculate gas factor for insulated containers per
    Notes to Table 1 and Table 2.

    Fluid_FR: Fluid at flow rating temperature and pressure
    """
    C = Fluid_FR.C_gas_constant
    TZM = 1 / Fluid_FR.MZT  # sqrt(T*Z/M)
    if Fluid_FR.P >= Fluid_FR.P_critical:
        L = Fluid_FR.specific_heat_input  # L is replaced by theta
    elif Fluid_FR.P >= 0.4*Fluid_FR.P_critical:
        # L =
        pass
    G_i = 73.4 * (Q_('1660 degF')-T) / (C*L) * TZM
