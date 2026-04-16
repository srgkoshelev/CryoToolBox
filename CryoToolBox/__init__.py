"""Cryogenic engineering utilities for thermodynamics, piping, and ODH.

`CryoToolBox` combines `CoolProp` for thermophysical properties and `Pint`
for unit handling to support cryogenic design calculations. The public API is
currently centered on:

1. :class:`cp_wrapper.ThermState` for fluid properties.
2. General thermodynamic and heat-transfer helpers from :mod:`functions`.
3. Piping geometry, loss coefficients, and flow calculations from
   :mod:`piping`.
4. Oxygen deficiency hazard analysis helpers from :mod:`odh`.

The package is still in its alpha phase, so compatibility is best-effort, but
the default `HEOS` backend is the supported path for general use.

.. include:: ./documentation.md
"""

from importlib.resources import files
import logging
import logging.config

# Setting up logging
__location__ = str(files(__package__))
logging.config.fileConfig(str(files(__package__).joinpath('logging.ini')))
logger = logging.getLogger(__name__)

# pdoc docformat
__docformat__ = "numpy"


from .std_conditions import Q_, ureg, T_NTP, P_NTP, P_MSC, T_MSC, P_STD, T_STD
from .cp_wrapper import ThermState, CryoToolBoxFutureWarning
from .functions import *
from . import piping
