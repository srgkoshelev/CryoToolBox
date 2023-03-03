from pint import UnitRegistry
import os

# Configuring units package:
ureg = UnitRegistry()
Q_ = ureg.Quantity
__location__ = os.path.dirname(os.path.abspath(__file__))
ureg.load_definitions(os.path.join(__location__, 'pint definitions.txt'))


# Setting units for "standard" flow
T_NTP = Q_(68, ureg.degF)  # Normal Temperature (NIST)
P_NTP = Q_(14.696, ureg.psi)  # Normal Pressure (NIST)

T_MSC = Q_(15, ureg.degC)  # Metric Standard Conditions (Crane TP-410)
P_MSC = Q_(101325, ureg.Pa)  # Metric Standard Conditions (Crane TP-410)

T_STD = Q_(60, ureg.degF)  # Standard conditions (BPVC VIII UG-129 (c))
P_STD = Q_(14.7, ureg.psi)  # Standard conditions (BPVC VIII UG-129 (c))
