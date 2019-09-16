from pint import UnitRegistry
import logging
import logging.config
import sys, os

#Setting up logging
__location__ = os.path.dirname(os.path.abspath(__file__))
logging.config.fileConfig(os.path.join(__location__, 'logging.ini'))
logger = logging.getLogger(__name__)

#' Configuring units package:
ureg = UnitRegistry(autoconvert_offset_to_baseunit = True)
Q_ = ureg.Quantity
ureg.load_definitions(os.path.join(__location__, 'pint definitions.txt'))


#Setting units for "standard" flow
T_NTP = Q_(68, ureg.degF) #Normal Temperature (NIST)
P_NTP = Q_(14.7, ureg.psi) #Normal Pressure (NIST)

T_MSC = Q_(15, ureg.degC) #Metric Standard Conditions (used by Crane TP-410)
P_MSC = Q_(101325, ureg.Pa) #Metric Standard Conditions (used by Crane TP-410)

#Default fluids
Air = {'fluid':'air', 'P':P_NTP, 'T':T_NTP}


from .cp_wrapper import *
from .functions import *
from . import piping
