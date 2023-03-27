# Defining failure rates as per FESHM 4240 chapter
from .odh import Q_

TABLE_1 = {'Compressor':
           {'Leak': Q_('5 * 10^-6 / hr'),
            'Component rupture': Q_('3 * 10^-7 / hr')},
           'Dewar':
           {'Loss of vacuum': Q_('1 * 10^-6 / hr')},
           'Electrical Power Failure':
           {'Time rate': Q_('1 * 10^-4 / hr'),
            'Demand rate': Q_('3 * 10^-4')},
           'Fluid line':
           {'Leak': Q_('5 * 10^-7 /hr'),
            'Rupture': Q_('2 * 10^-8 / hr')},
           'Cryogenic magnet powered':
           {'Rupture': Q_('2 * 10 ^-7 / hr')},
           'Cryogenic magnet not powered':
           {'Rupture': Q_('2 * 10 ^-8 / hr')},
           'Header piping assembly':
           {'Rupture': Q_('1 * 10^-8 / hr')},
           'U-Tube change':
           {'Small event': Q_('3 * 10^-2'),
            'Large event': Q_('1 * 10^-3')},
}
TABLE_2 = {'Battery':
           {'No output': Q_('3 * 10^-6 /hr')},
           'Circuit Breaker':
           {'Failure to operate': Q_('1 * 10^-3'),
            'Premature transfer': Q_('1 * 10^-6 /hr')},
           'Diesel':
           {'Failure to start': Q_('3 * 10^-2'),
            'Failure to run (Emergency)': Q_('3 * 10^-3 / hr'),
            'Failure to run (Engine only)': Q_('3 * 10^-4 / hr')},
           'Fan':
           {'Failure to run': Q_('9 * 10^-6 / hr')},
           'Fuse':
           {'Premature open': Q_('1 * 10^-6 /hr'),
            'Failure to open': Q_('1 * 10^-5')},
           'Flange, reinforced gasket':
           {'Leak': {'Area': Q_('10 mm^2'), 'Failure rate': Q_('4 * 10^-7 / hr')},
            'Rupture': Q_('1 * 10 ^-9 / hr')},
           'Flange, soft gasket':
           {'Leak': {'Area': Q_('10 mm^2'), 'Failure rate': Q_('4 * 10^-7 / hr')},
            'Blowout': {'Area': Q_('1000 mm^2'), 'Failure rate': Q_('3 * 10 ^-8 / hr')},
            'Rupture': Q_('1 * 10 ^-9 / hr')},
           'Instrumentation':
           {'Failure to operate': Q_('1 * 10^-6 / hr'),
            'Shift': Q_('3 * 10^-5 / hr')},
           'Louver':
           {'Failure rate': Q_('3 * 10^-7 / hr')},
           'Piping':
           {'Small leak': {'Area': Q_('10 mm^2'), 'Failure rate': Q_('1 * 10^-9 / (m*hr)')},
            'Large leak': {'Area': Q_('1000 mm^2'), 'Failure rate': Q_('1 * 10^-10 / (m*hr)')}, #Only for pipes > 2"
            'Rupture': Q_('3 * 10^-11 /(m*hr)')},
           'Pipe weld':#Failure rates need to be multiplied my D/t, where D - diameter, t - wall thickness
           {'Small leak': {'Area': Q_('10 mm^2'), 'Failure rate': Q_('2 * 10^-11 / hr')},
            'Large leak': {'Area': Q_('1000 mm^2'), 'Failure rate': Q_('2 * 10^-12 / hr')}, #Only for pipes > 2"
            'Rupture': Q_('6 * 10^-13 / hr')},
           'Pump':
           {'Failure to start': Q_('1 * 10^-3'),
            'Failure to run, normal': Q_('3 * 10^-5 / hr'),
            'Failure to run, extreme': Q_('1 * 10^-3 / hr')},
           'Relay':
           {'Failure to energize': Q_('1 * 10^-4'),
            'Failure to close': Q_('3 * 10^-7 / hr'),
            'Short': Q_('1 * 10^-8 / hr'),
            'Open contact': Q_('1 * 10^-7 / hr')},
           'Solid State Device':
           {'HI PWR':
            {'Failure rate': Q_('3 * 10^-6 / hr'),
             'Short': Q_('1 * 10^-6 /hr')},
            'LOW PWR':
            {'Failure rate': Q_('1 * 10^-6 / hr'),
             'Short': Q_('1 * 10^-7 /hr')}},
           'Switch':
           {'Limit': Q_('3 * 10^-4'),
            'Torque': Q_('1 * 10^-4'),
            'Pressure': Q_('1 * 10^-4'),
            'Manual': Q_('1 * 10^-5'),
            'Short': Q_('1 * 10^-8 / hr')},
           'Transformer':
           {'Open': Q_('1 * 10^-6 / hr'),
            'Short': Q_('1 * 10^-6 / hr')},
           'Valve, motorized':
           {'Failure to operate': Q_('1 * 10^-3'),
            'Failure to remain open': Q_('1 * 10^-4'),
            'External leak': Q_('1 * 10^-8 /hr'),
            'Rupture': Q_('5 * 10^-10 / hr')},
           'Valve, solenoid':
           {'Failure to operate': Q_('1 * 10^-3')},
           'Valve, pneumatic':
           {'Failure to operate': Q_('3 * 10^-4'),
            'Failure to remain open': Q_('1 * 10^-4'),
            'External leak': Q_('1 * 10^-8 / hr'),
            'Rupture': Q_('5 * 10^-10 / hr')},
           'Valve, check':
           {'Failure to open': Q_('1 * 10^-4'),
            'Reverse leak': Q_('3 * 10^-7 / hr'),
            'External leak': Q_('1 * 10^-8 / hr'),
            'Rupture': Q_('5 * 10^-10 / hr')},
           'Orifice':
           {'Rupture': Q_('1 * 10^-8')},
           'Valve, manual':
           {'Failure to open': Q_('1 * 10^-4'),
            'External leak': Q_('1 * 10^-8 / hr'),
            'Rupture': Q_('5 * 10^-10 / hr')},
           'Valve, relief':
           {'Failure to open': Q_('1 * 10^-5'),
            'Premature open': Q_('1 * 10^-5 / hr')},
           'Vessel, pressure':
           {'Small leak': {'Area': Q_('10 mm^2'), 'Failure rate': Q_('8 * 10^-8 / hr')},
            'Failure': Q_('5 * 10^-9 / hr')},
           'Wire':
           {'Open': Q_('3 * 10^-6 / hr'),
            'Short to GND': Q_('3 * 10^-7 / hr'),
            'Short to PWR': Q_('1 * 10^-8 / hr')},
}
