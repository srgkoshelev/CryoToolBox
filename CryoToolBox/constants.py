from .std_conditions import P_NTP, T_NTP
from .cp_wrapper import ThermState
# Default fluids
AIR = ThermState('air', P=P_NTP, T=T_NTP)
