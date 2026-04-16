"""Flow-conversion helpers with backward-compatible aliases.

This module is the preferred home for unit-aware flow conversions. Legacy names
from :mod:`functions` remain available as aliases for compatibility.
"""

import logging

from .std_conditions import T_NTP, P_NTP, ureg
from .constants import AIR
from .cp_wrapper import ThermState

logger = logging.getLogger(__name__)


def to_equiv_air(m_dot_fluid, fluid):
    """Convert mass flow into equivalent standard flow of air.

    This is the forward conversion behind the legacy :func:`to_scfma` helper.
    """
    C_fluid = fluid.C_gas_const
    C_air = AIR.C_gas_const

    m_dot_air = m_dot_fluid * C_air / C_fluid * AIR.MZT / fluid.MZT
    q_air = m_dot_air / AIR.Dmass
    q_air.ito(ureg.ft**3 / ureg.min)
    return q_air


def from_equiv_air(q_air, fluid):
    """Convert equivalent standard air flow into fluid mass flow."""
    C_fluid = fluid.C_gas_const
    C_air = AIR.C_gas_const

    m_dot_air = q_air * AIR.Dmass
    m_dot_fluid = m_dot_air * C_fluid / C_air * fluid.MZT / AIR.MZT
    m_dot_fluid.ito(ureg.g / ureg.s)
    return m_dot_fluid


def to_std_flow(flow_rate, fluid):
    """Convert mass flow or actual volumetric flow to standard volumetric flow."""
    fluid_ntp = ThermState(fluid.name)
    fluid_ntp.update('T', T_NTP, 'P', P_NTP)
    if flow_rate.dimensionality == ureg('kg/s').dimensionality:
        q_std = flow_rate / fluid_ntp.Dmass
    elif flow_rate.dimensionality == ureg('m^3/s').dimensionality:
        if fluid.Dmass != -float('Inf') * ureg.kg / ureg.m**3:
            q_std = flow_rate * fluid.Dmass / fluid_ntp.Dmass
        else:
            logger.warning(
                'Flow conditions for volumetric flow %.3s are not set. '
                'Assuming standard flow at NTP.',
                flow_rate,
            )
            q_std = flow_rate
    else:
        raise ValueError(
            'Flow dimensionality is not supported: '
            f'{flow_rate.dimensionality}.'
        )
    q_std.ito(ureg.ft**3 / ureg.min)
    return q_std


def from_std_flow(q_std, fluid):
    """Convert standard volumetric flow to mass flow."""
    fluid_ntp = ThermState(fluid.name)
    fluid_ntp.update('T', T_NTP, 'P', P_NTP)
    m_dot = q_std * fluid_ntp.Dmass
    return m_dot.to(ureg.g / ureg.s)


# Compatibility aliases
to_scfma = to_equiv_air
from_scfma = from_equiv_air
to_standard_flow = to_std_flow
to_mass_flow = from_std_flow

__all__ = [
    'from_equiv_air',
    'from_scfma',
    'from_std_flow',
    'to_equiv_air',
    'to_mass_flow',
    'to_scfma',
    'to_standard_flow',
    'to_std_flow',
]
