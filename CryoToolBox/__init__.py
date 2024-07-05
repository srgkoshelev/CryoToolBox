"""
 `CryoToolBox` is a Python module for simple heat transfer and hydraulics
calculation using [CoolProp](https://github.com/CoolProp/CoolProp) and
[Pint](https://github.com/hgrecco/pint) for unit handling.

Provides:
    1. Thermal properties via open source HEOS or REFPROP back-end.
    2. Convenient calculation using a thermodynamic state `cp_wrapper.ThermState`.
    3. Seemless use of units for all calculations.
    4. Useful thermodynamic and hydraulic functions.
    5. Relief valve sizing functions.
    6. Oxygen deficiency analysis utilities.

.. include:: ./documentation.md
"""

import logging
import logging.config
import os

# Setting up logging
__location__ = os.path.dirname(os.path.abspath(__file__))
logging.config.fileConfig(os.path.join(__location__, 'logging.ini'))
logger = logging.getLogger(__name__)


from .std_conditions import ureg, T_NTP, P_NTP, P_MSC, T_MSC, P_STD, T_STD
from .functions import *
from . import piping, heated_pipe, line
