"""
 `heat_transfer` is a Python module for simple heat transfer and hydraulics
calculation using [CoolProp](https://github.com/CoolProp/CoolProp) and
[Pint](https://github.com/hgrecco/pint) for unit handling.

Provides:
    1. Thermal properties via open source HEOS or REFPROP back-end
    2. Convinient calculation using a thermodynamic state `cp_wrapper.ThermState`
    3. Seemless use of units for all calculations
    4. Useful thermodynamic and hydraulic functions
    3. Relief valve sizing for CGA
"""

from pint import UnitRegistry
import logging
import logging.config
import os

# Setting up logging
__location__ = os.path.dirname(os.path.abspath(__file__))
logging.config.fileConfig(os.path.join(__location__, 'logging.ini'))
logger = logging.getLogger(__name__)

# Configuring units package:
ureg = UnitRegistry(autoconvert_offset_to_baseunit=True)
Q_ = ureg.Quantity
ureg.load_definitions(os.path.join(__location__, 'pint definitions.txt'))

# Setting units for "standard" flow
T_NTP = Q_(68, ureg.degF)  # Normal Temperature (NIST)
P_NTP = Q_(14.7, ureg.psi)  # Normal Pressure (NIST)

T_MSC = Q_(15, ureg.degC)  # Metric Standard Conditions (used by Crane TP-410)
P_MSC = Q_(101325, ureg.Pa)  # Metric Standard Conditions (used by Crane TP-410)

T_STD = Q_(60, ureg.degF)  # Standard conditions (BPVC VIII UG-129 (c))
P_STD = Q_(14.7, ureg.psi)  # Standard conditions (BPVC VIII UG-129 (c))

from .cp_wrapper import ThermState
# Default fluids
Air = ThermState('air')
Air.update('T', T_NTP, 'P', P_NTP)


from .functions import *
from . import piping
