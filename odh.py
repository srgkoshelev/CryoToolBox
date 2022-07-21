"""odh - Oxygen Deficiency Hazard analysis tool.

Based on Fermilab ES&H Manual chapter 4240 (see
https://eshq.fnal.gov/manuals/feshm/)

This package is maintained at https://github.com/srgkoshelev/ODH_analysis"""

import math
from scipy.special import binom
from . import ureg, Q_
from . import T_NTP, P_NTP
from .functions import ThermState
from .functions import to_standard_flow
from .piping import G_nozzle
from dataclasses import dataclass
import xlsxwriter
from enum import Enum, auto
# Loading FESHM 4240 Failure rates
from .FESHM4240_TABLES import TABLE_1, TABLE_2


# Probability of failure on demand and failure rate for ODH head+chassis
# Data from J. Anderson's analysis
PFD_ODH = Q_('2 * 10^-3')
lambda_ODH = 2.3e-6 / ureg.hr
TRANSFER_LINE_LEAK_AREA = Q_(10, ureg.mm**2)  # Taken for consistency
# Min required air intake from Table 6.1, ASHRAE 62-2001
ASHRAE_MIN_FLOW = 0.06 * ureg.ft**3/(ureg.min*ureg.ft**2)


class ODHError(Exception):
    pass


class ConstLeakTooBig(Exception):
    pass


class ConstLeakTooLong(Exception):
    pass

@ureg.check(None, None, '1/[time]', '[length]^3/[time]', '[time]', None)
@dataclass
class Leak:
    """Describe inert gas leak from a Source"""
    name: str
    fluid: ThermState
    failure_rate: ureg.Quantity
    q_std: ureg.Quantity
    tau: ureg.Quantity
    N_events: int
    _is_const = False


@ureg.check(None, None, '[length]^3/[time]', '[time]')
@dataclass
class ConstLeak:
    name: str
    fluid: ThermState
    q_std: ureg.Quantity
    tau: ureg.Quantity
    N_events = 1
    _is_const = True
    # TODO Remove once isinstance() issue resolved


class Source:
    """Source of inert gas

    Attributes
    ----------
    sol_PFD : float
        Probability of failure on demand (PFD) for isolation valve.
        If the source doesn't have isolating valve
        the probability is 1.
    """
    def __init__(self, name, fluid, volume, N=1, isol_valve=False):
        """Define the possible source of inert gas.

        Parameters
        ----------
        name : str
            Name of the source.
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        volume : ureg.Quantity {length: 3}
            Volume of the fluid stored.
        N : int
            Quantity of the sources if several similar sources exist,
            e.g. gas bottles.
        isol_valve : bool
            Denotes whether the source is protected using a normally closed
            valve.
        """
        self.name = name
        self.fluid = fluid
        self.leaks = []
        # Number of sources if multiple exist, e.g. gas cylinders
        # Increases probability of failure by N.
        self.N = N
        # Calculating volume at standard conditions
        temp_state = fluid.copy()
        temp_state.update('T', T_NTP, 'P', P_NTP)
        self.volume = volume*fluid.Dmass/temp_state.Dmass
        self.volume.ito(ureg.feet**3)
        # By default assume there is no isolation valve
        # that is used by ODH system
        self.isol_valve = isol_valve
        self.sol_PFD = (int(not isol_valve) or
                        TABLE_2['Valve, solenoid']['Failure to operate'])

    def pipe_failure(self, tube, N_welds, fluid=None, q_std_rupture=None):
        """Add pipe failure to the leaks dict.

        For a given tube calculate leak parameters as following:
        - failure rate
        - standard volumetric flow
        - duration of the release
        - number of possible similar events
        for all failure modes for piping and welds listed in Table 2 of
        FESHM 4240.

        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        tube : heat_transfer.piping.Tube
        N_welds : int
            Number of welds on the tube.
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid for the release.
        q_std_rupture : ureg.Quantity {length: 3, time: -1}
            Standard volumetric flow rate for flange rupture.
        """
        # If fluid not defined use fluid of the Source
        fluid = fluid or self.fluid
        # Failure rate coefficients; Piping failure rate is per unit of length,
        # weld is dependent on number of welds, pipe OD and wall thickness
        failure_rate_coeff = {'Piping': (tube.L, 1),
                              'Pipe weld': (tube.OD / tube.wall,
                                            N_welds)}
        # Piping and weld leaks as per Table 2
        for cause, fr_coefs in failure_rate_coeff.items():
            for mode in TABLE_2[cause].keys():
                if tube.D > 2 or mode != 'Large leak':  # Large leak only for D > 2"
                    name = f'{cause} {mode.lower()}: {tube}'
                    fr_coef, N_events = fr_coefs
                    if fr_coef.magnitude == 0:
                        raise ODHError('Failure rate should not be 0 '
                                       f'{cause} {mode}.')
                    if N_events == 0:
                        continue
                    if mode == 'Rupture':
                        failure_rate = fr_coef * TABLE_2[cause][mode]
                        if q_std_rupture is not None:
                            q_std = q_std_rupture
                        else:
                            area = tube.area
                            q_std = hole_leak(area, fluid)
                    else:
                        failure_rate = fr_coef * \
                            TABLE_2[cause][mode]['Failure rate']
                        area = TABLE_2[cause][mode]['Area']
                        if area > tube.area:
                            ODHError('Leak area cannot be larger'
                                     ' than pipe area.')
                        q_std = hole_leak(area, fluid)
                    self.add_leak(name, fluid, failure_rate, q_std, N_events)

    def dewar_insulation_failure(self, q_std):
        """Add dewar insulation failure to leaks dict.

        Store failure rate, flow rate and expected time duration of the
        failure event for the dewar insulation failure. Based on FESHM4240.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        q_std : ureg.Quantity {length: 3, time: -1}
            Standard volumetric flow rate of the relief for the case of
            dewar insulation failure.
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        """
        failure_rate = TABLE_1['Dewar']['Loss of vacuum']
        self.add_leak('Dewar insulation failure', 'fluid', failure_rate, q_std,
                      1)

    # def u_tube_failure(self, outer_tube, inner_tube, L, use_rate,
    #                    fluid=None, N=1):
    #     """Add U-Tube failure to leaks dict.

    #     Store failure rate, flow rate and expected time duration of the
    #     failure event for the dewar insulation failure. Based on FESHM4240.
    #     Failure modes are analyzed by `Volume.odh` method.

    #     Parameters
    #     ----------
    #     flow_rate : ureg.Quantity {mass: 1, time: -1} or {length: 3, time: -1}
    #         Relief flow rate for the case of dewar insulation failure.
    #     fluid : heat_transfer.ThermState
    #         Thermodynamic state of the fluid stored in the source.
    #     """
    #     # TODO Make areas adjustable, add info to docstring
    #     flow_path_cases = {'Small event': piping.Annulus(outer_tube.ID,
    #                                                         inner_tube.OD,
    #                                                         L=L),
    #                        'Large event': outer_tube}
    #     for mode in TABLE_1['U-Tube change']:
    #         flow_path = flow_path_cases[mode]
    #         name = f'U-Tube {mode.lower()}: {flow_path}'
    #         failure_rate = TABLE_1['U-Tube change'][mode] * \
    #             use_rate
    #         area = flow_path.area
    #         # TODO move this and gas leak check to separate method
    #         if area > outer_tube.area:
    #             logger.warning('Leak area cannot be larger'
    #                            ' than outer tube area.')
    #             continue
    #         # If fluid not defined use fluid of the Source
    #         fluid = fluid or self.fluid
    #         q_std = Source._leak_flow(flow_path, area, fluid)
#             self.failure_mode(name, failure_rate, q_std, N)

    def flange_failure(self, pipe, N=1, fluid=None, q_std_rupture=None):
        """Add reinforced or preformed gasket flange failure
        to leaks dict.

        Store failure rate, flow rate and expected time duration of
        the event for transfer line failure. Based on FESHM 4240.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        pipe : heat_transfer.Pipe
        N : int
            Number of reinforced seal connections on the pipe.
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        q_std_rupture : ureg.Quantity {length: 3, time: -1}
            Standard volumetric flow rate for flange rupture.
        """
        fluid = fluid or self.fluid
        # TODO Make leak and rupture areas adjustable, add info to docstring
        table = TABLE_2['Flange, reinforced gasket']
        # Leak case
        name = f'Flange leak: {pipe}'
        failure_rate = table['Leak']['Failure rate']
        area = table['Leak']['Area']
        q_std = hole_leak(area, fluid)
        self.add_leak(name, fluid, failure_rate, q_std, N)
        # Rupture
        name = f'Flange rupture: {pipe}'
        failure_rate = table['Rupture']
        if q_std_rupture is not None:
            q_std = q_std_rupture
        else:
            area = pipe.area
            q_std = hole_leak(area, fluid)
        self.add_leak(name, fluid, failure_rate, q_std, N)

    def valve_failure(self, pipe, N=1, fluid=None, q_std_rupture=None):
        """Add valve leak and rupture failure modes to leaks dict.

        Store failure rate, flow rate and expected time duration of
        the event for transfer line failure. Based on FESHM 4240.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        pipe : heat_transfer.Pipe
            Pipe/tube upstream the valve.
        N : int
            Number of valves
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        q_std_rupture : ureg.Quantity {length: 3, time: -1}
            Standard volumetric flow rate for flange rupture.
        """
        fluid = fluid or self.fluid
        table = TABLE_2['Valve, pneumatic']
        # Leak case
        name = f'Valve leak: {pipe}'
        failure_rate = table['External leak']
        # Using area from flange leak for consistency
        area = TABLE_2['Flange, soft gasket']['Leak']['Area']
        q_std = hole_leak(area, fluid)
        self.add_leak(name, fluid, failure_rate, q_std, N)
        # Rupture
        name = f'Valve rupture: {pipe}'
        failure_rate = table['Rupture']
        if q_std_rupture is not None:
            q_std = q_std_rupture
        else:
            area = pipe.area
            q_std = hole_leak(area, fluid)
        self.add_leak(name, fluid, failure_rate, q_std, N)

    def transfer_line_failure(self, pipe, fluid=None, q_std_rupture=None, N=1):
        """Add transfer line failure to leaks dict.

        For a given tube calculate leak parameters as following:
        - failure rate
        - standard volumetric flow
        - duration of the release
        - number of possible similar events
        for all failure modes for transfer line listed in Table 1 of
        FESHM 4240. Note that Rationale for Table 1: “Fermilab Equipment
        Failure Rate Estimates” clearly states that only bayonet failures
        are considered for these failure modes. Therefore only bayonet leak and
        blowout (rupture) are considered. Leak area is not defined in
        FESHM 4240 and is defined globally as `TRANSFER_LINE_LEAK_AREA`.

        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        pipe : heat_transfer.Pipe
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        q_std_rupture : ureg.Quantity {length: 3, time: -1}
            Standard volumetric flow rate for fluid line rupture.
        N : int
            Number of bayonets/soft seals on the transfer line.
        """
        # TODO This should be replaced by flange failure at some point
        # If fluid not defined use fluid of the Source
        fluid = fluid or self.fluid
        # Leak case
        name = f'Fluid line gasket leak: {pipe}'
        failure_rate = TABLE_1['Fluid line']['Leak']
        area = TRANSFER_LINE_LEAK_AREA
        q_std = hole_leak(area, fluid)
        self.add_leak(name, fluid, failure_rate, q_std, N)
        # Rupture case
        name = f'Fluid line gasket rupture: {pipe}'
        failure_rate = TABLE_1['Fluid line']['Rupture']
        if q_std_rupture is not None:
            q_std = q_std_rupture
        else:
            area = pipe.area
            q_std = hole_leak(area, fluid)
        self.add_leak(name, fluid, failure_rate, q_std, N)

    def line_failure(self, tube, N_welds, N_seals, N_valves, fluid=None,
                     q_std_rupture=None):
        """Adds pipe, pipe weld, seal/flange, and valve failures
        to the leaks dict.

        For a given tube calculate leak parameters as following:
        - failure rate
        - standard volumetric flow
        - duration of the release
        - number of possible similar events
        for all failure modes for piping and welds listed in Table 2 of
        FESHM 4240.

        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        tube : heat_transfer.piping.Tube
        N_welds : int
            Number of welds on the tube.
        N_seals : int
            Number of reinforced seal connections on the pipe.
        N_valves : int
            Number of valves
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid for the release.
        q_std_rupture : ureg.Quantity {length: 3, time: -1}
            Standard volumetric flow rate for flange rupture.
        """
        self.pipe_failure(tube, N_welds, fluid, q_std_rupture)
        if N_seals != 0:
            self.flange_failure(tube, N_seals, fluid, q_std_rupture)
        if N_valves != 0:
            self.valve_failure(tube, N_valves, fluid, q_std_rupture)

    def pressure_vessel_failure(self, relief_area=None, fluid=None):
        """Add pressure vessel failure to leaks dict.

        Store failure rate, flow rate and expected time duration of
        the event for transfer line failure. Based on FESHM 4240.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        relief_area : ureg.Quantity {length: 2}
            Vacuum jacket relief area if the vessel has one, None otherwise.
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        """
        # If fluid not defined use fluid of the Source
        fluid = fluid or self.fluid
        # Leak case
        name = 'Pressure vessel leak'
        area = TABLE_2['Vessel, pressure']['Small leak']['Area']
        failure_rate = TABLE_2['Vessel, pressure']['Small leak']['Failure rate']
        q_std = hole_leak(area, fluid)
        self.add_leak(name, fluid, failure_rate, q_std, 1)

        # Rupture case
        name = 'Pressure vessel rupture'
        if relief_area is None:
            # No vacuum jacket on the vessel
            area = 1000 * ureg.mm**2  # Equal to large leak of a pipe
        else:
            area = relief_area
        failure_rate = TABLE_2['Vessel, pressure']['Failure']
        q_std = hole_leak(area, fluid)
        self.add_leak(name, fluid, failure_rate, q_std, 1)

    def const_leak(self, name, fluid, q_std, N=1):
        """Add constant leak to leaks dict.

        Store flow rate and expected time duration of the
        failure event for general failure mode.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        name : str
            Name of the failure mode
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        q_std : ureg.Quantity {length: 3, time: -1}
            Standard volumetric flow rate.
        N : int
            Quantity of leaks.
        """
        N_events = N * self.N
        tau = self.volume/q_std
        self.leaks.append(ConstLeak(name, fluid, q_std*N_events, tau))

    def add_leak(self, name, fluid, failure_rate, q_std, N=1):
        """Add general leak to leaks dict.

        Store failure rate, flow rate and expected time duration of the
        failure event for general failure mode.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        name : str
            Name of the failure mode
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        failure rate : ureg.Quantity {time: -1}
            Failure rate of the failure mode,
            i.e. how often the failure occurs
        q_std : ureg.Quantity {length: 3, time: -1}
            Standard volumetric flow rate.
        N : int
            Quantity of similar failure modes.
        """
        N_events = N * self.N
        tau = self.volume/q_std
        total_failure_rate = N_events*failure_rate
        total_failure_rate.ito(1/ureg.hr)
        # self.leaks.append((name, total_failure_rate, q_std, tau.to(ureg.min),
        #                    N_events))
        self.leaks.append(Leak(name, fluid, total_failure_rate, q_std, tau.to(ureg.min),
                               N_events))

    @staticmethod
    def combine(name, sources):
        """Combine several ODH sources sharing volume.

        Can be used for failure modes affecting several sources in parallel.

        Parameters
        ----------
        name : str
            Name of the new combined source.
        sources : list of Source
            Sources connected together. If a source represents multiple
            identical sources, e.g., N > 1, the volumes are added together.

        Returns
        -------
        Source
            Combined source of inert gas.
        """
        fluid = ThermState(sources[0].fluid.name, T=T_NTP, P=P_NTP)
        if all([source.fluid.name == fluid.name for source in sources]):
            total_volume = sum([source.N*source.volume for source in sources])
            return Source(name, fluid, total_volume)
        else:
            ODHError('All volumes should contain the same fluid')
            return None

    def __repr__(self):
        return f'Source({self.name}, {self.fluid}, {self.volume}, N={self.N}, isol_valve={self.isol_valve})'

    def __str__(self):
        if self.N == 1:
            quant_str = ''
        else:
            quant_str = f'{self.N}x '
        return quant_str + f'{self.name}, ' + \
            f'{self.volume.to(ureg.ft**3):.3g~} ' + \
            f'of {self.fluid.name}'

    def print_leaks(self):
        """Print information on the leaks defined for the source."""
        print(f'Printing leaks for {self}')
        print()
        for leak in self.leaks:
            name, total_failure_rate, q_std, tau, N_events = leak
            print(name)
            print(f'Total failure rate: {total_failure_rate:.3g~}')
            print(f'Leak rate: {q_std:.3g~}')
            print(f'Event duration: {tau:.3g~}')
            print(f'Number of events: {N_events}')
            print()


@ureg.check(None, None, None, None, '1/[time]', None, '1/[time]', '1/[time]',
            None, None, '[length]^3/[time]', '[time]', '[length]^3/[time]',
            None, None)
@dataclass
class FailureMode:
    """Describes ODH failure modes, a combination of a leak and response to it."""
    name: str
    source: Source
    fluid: ThermState
    case: str
    phi: ureg.Quantity
    O2_conc: float
    leak_fr: ureg.Quantity
    P_i: ureg.Quantity
    F_i: float
    outage: bool
    q_leak: ureg.Quantity
    tau: ureg.Quantity
    Q_fan: ureg.Quantity
    N_fan: int
    N: int


@ureg.check(None, '[time]')
@dataclass
class BuildPower():
    pfd: float = TABLE_1['Electrical Power Failure']['Demand rate']
    max_outage: ureg.Quantity = float('inf') * ureg.hr


@ureg.check('[length]^3/[time]', None, '[time]', '1/[time]',
            '[length]^3/[time]', None)
@dataclass
class BuildVent():
    """
        Q_fan : ureg.Quantity {length: 3, time: -1}
            Volumetric flow of a single ODH fan installed in the volume.
        N_fans : int
            Number of fans installed.
        Test_period : ureg.Quantity {time: 1}
            Test period of the fans and louvers.
        lambda_fan : ureg.Quantity {time: -1}
            Failure rate of the fans in the building.
        min_vent : ureg.Quantity {length: 3, time: -1}
            Min volumetric flow required or present in the building.
            Default value is 0 CFM.
        indep_vent : bool
            `True` means min fresh air intake is independent from the fans,
            e.g., provided by HVAC system. `False` means min ventilation
            is achieved using the dedicated ODH fan, generally run at a lower
            capacity. For this case fan failure means no ventilation will be
            present in the building. Default value: `False`.
"""
    Q_fan: ureg.Quantity
    N_fans: int
    Test_period: ureg.Quantity
    lambda_fan: ureg.Quantity = TABLE_2['Fan']['Failure to run']
    min_vent: ureg.Quantity = 0*ureg.ft**3/ureg.min
    indep_vent: bool = False

    def const_vent(self, fail_case):
        if fail_case == FailCase.POWER:
            return 0 * ureg.ft**3/ureg.min
        elif fail_case == FailCase.ODH:
            return self.min_vent
        elif fail_case == FailCase.FAN:
            return float(self.indep_vent) * self.min_vent

    def fan_flowrates(self):
        """Calculate (Probability, flow) pairs for all combinations of fans
        working.

        All fans are expected to have same volume flow.
        """
        # TODO add fans with different volumetric rates (see report as well)
        fail_rate = self.lambda_fan
        fan_flowrates = []
        for m in range(self.N_fans+1):
            # Probability of exactly m units starting
            P_m_fan_work = prob_m_of_n(m, self.N_fans, self.Test_period,
                                       fail_rate)
            if m == 0:
                flowrate = self.const_vent(FailCase.FAN)
            else:
                flowrate = m * self.Q_fan
            fan_flowrates.append((P_m_fan_work, flowrate, m))
        return fan_flowrates


class FailCase(Enum):
    """Failure cases possible during a release event, i.e., power failure, ODH system failure, or fan failure."""
    POWER = auto()
    ODH = auto()
    FAN = auto()


class Volume:
    """Volume/building affected by inert gases."""
    def __init__(self, name, volume, build_power, build_vent):
        """Define a volume affected by inert gas release from  a `Source`.

        Parameters
        ----------
        name : str
            Name of the volume.
        volume : ureg.Quantity {length: 3}
            Volume of the building or part of the building.
        """
        self.name = name
        self.volume = volume
        self.vent = build_vent
        self.power = build_power
        self.PFD_ODH = PFD_ODH  # Default value for ODH system failure

    def odh(self, sources, power_outage=False):
        """Calculate ODH fatality rate for given `Source`s.

        For each leak of each source ODH conditions are analyzed and
        fatality rates are calculated. The results are collected in
        fail_modes list.

        Parameters
        ----------
        sources : list
            Sources affecting the volume.
        power_outage : bool
            Shows whether there is a power outage is in effect.
            Default is no outage.

        Returns
        -------
        None
        """
        self.fail_modes = []
        for source in sources:
            for leak in source.leaks:
                if not leak._is_const:
                    self._fatality_no_power(source, leak, power_outage)
                    if not power_outage:
                        self._fatality_no_odh(source, leak)
                        self._fatality_with_fans(source, leak)
                else:
                    self._fatality_const_leak(source, leak, power_outage)
        print(f'{self.name} failure modes analyzed: {len(self.fail_modes)}')
        print(f'{self.name} ODH class: {self.odh_class()}')
        print(f'{self.name} fatality rate: {self.phi:.2g}')

    def _fatality_no_power(self, source, leak, power_outage):
        """Calculate fatality rate in the case of power failure.

        Source isolation is the only protection for this case. No ventilation
        is available without electric power.

        Parameters
        ----------
        source : Source
        leak : Leak
            Leak failure rate, volumetric flow rate, event duration, and number
            of events.
        power_outage : bool
            Shows whether there is a power outage is in effect.

        Returns
        -------
        None
        """
        if power_outage:
            PFD_power = 1
            tau_event = min(leak.tau, float('inf') * ureg.hr)
        else:
            PFD_power = self.power.pfd
            tau_event = min(leak.tau, self.power.max_outage)
        PFD_isol_valve = source.sol_PFD
        P_case = PFD_power * PFD_isol_valve
        Q_fan = 0 * ureg.ft**3/ureg.min
        case = 'No power and isolation'
        self._add_failure_mode(case, P_case, tau_event, leak.fluid,
                               leak.failure_rate, source, leak, Q_fan, 0,
                               power_outage)

    def _fatality_no_odh(self, source, leak):
        """Calculate fatality rate in the case of ODH system failure.

        Power has to be  on for ODH system to fail. Source isolation is
        possible for this event. Only ventilation independent from ODH system
        is available, e.g., HVAC.

        Parameters
        ----------
        source : Source
        leak : Leak
            Leak failure rate, volumetric flow rate, event duration, and number
            of events.
        power_outage : bool
            Shows whether there is a power outage is in effect.

        Returns
        -------
        None
        """
        PFD_isol_valve = source.sol_PFD
        P_case = (1-self.power.pfd) * self.PFD_ODH * PFD_isol_valve
        Q_fan = self.vent.const_vent(FailCase.ODH)
        tau_event = min(leak.tau, self.vent.Test_period)
        case = 'ODH fail'
        self._add_failure_mode(case, P_case, tau_event, leak.fluid,
                               leak.failure_rate, source, leak, Q_fan, 0,
                               False)

    def _fatality_with_fans(self, source, leak):
        """Calculate fatality rate for any number of fans operational.

        Power has to be  on, and ODH system operational. Source isolation is
        possible for this event. If ventilation is independent from the ODH
        system, e.g., HVAC, it is used for the case of all fans failure at
        the time of the event.

        Parameters
        ----------
        source : Source
        leak : Leak
            Leak failure rate, volumetric flow rate, event duration, and number
            of events.

        Returns
        -------
        None
        """
        PFD_isol_valve = source.sol_PFD
        # Probability for this group of the cases
        P_group = (1-self.power.pfd) * (1-self.PFD_ODH) * PFD_isol_valve
        tau_event = min(leak.tau, self.vent.Test_period)
        case = 'Fan failure'
        for (P_fan, Q_fan, N_fans) in self.vent.fan_flowrates():
            P_case = P_group * P_fan  # Probability of the particular case
            self._add_failure_mode(case, P_case, tau_event, leak.fluid,
                                   leak.failure_rate, source, leak, Q_fan,
                                   N_fans, False)

    def _add_failure_mode(self, case, P_case, tau_event, fluid, failure_rate,
                          source, leak, Q_fan, N_fans, power_outage):
        """Add a failure mode.

        Parameters
        ----------
        case : str
            Description of the particular response case, e.g., no power.
        P_case : float
            Probability of a particular state of the system, assuming the leak
            has happened. Alternatively, probability of a constant leak or 1.
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        failure_rate : ureg.Quantity {time: -1}
            Failure rate of the leak. Alternatively, system failure rate for
            a constant leak.
        tau_event : ureg.Quantity {time: 1}
            Maximum duration of the event.
        source : Source
            Inert gas source.
        leak : Leak
            Leak failure rate, volumetric flow rate, event duration, and number
            of events.
        Q_fan : ureg.Quantity {length: 3, time: -1}
            Combined volumetric flow of ODH fans for the event.
        N_fans : int
            Number of fans operational during the event.
        power_outage : bool
            Shows whether there is a power outage is in effect.

        Returns
        -------
        None
        """
        P_i = failure_rate * P_case
        O2_conc = conc_vent(self.volume, leak.q_std, Q_fan, tau_event)
        F_i = self._fatality_prob(O2_conc)
        phi_i = P_i*F_i
        self.fail_modes.append(
            FailureMode(leak.name, source, fluid, case, phi_i, O2_conc,
                        failure_rate, P_i, F_i, power_outage, leak.q_std,
                        tau_event, Q_fan, N_fans, leak.N_events)
        )

    def _fatality_const_leak(self, source, leak, power_outage):
        """Calculate fatality rates for constant leak case.

        Calculate fatality rates for the case of ODH system responding and
        fans powered but some of the fans failing on demand.
        See wiki for further explanation.
        Adds calculation results to the fail_modes list.

        Parameters
        ----------
        source : Source
        leak : Leak
            Leak failure rate, volumetric flow rate, event duration, and number
            of events.
        power_outage : bool
            Shows whether there is a power outage is in effect.

        Returns
        -------
        None
        """

        PFD_isol_valve = source.sol_PFD
        if power_outage:
            t_max_outage = self.volume / leak.q_std * math.log(21/18)
            if self.power.max_outage > t_max_outage:
                raise ConstLeakTooLong(
                    f'Constant leak {leak.name} from {source.name} is not '
                    f'mitigated during power outage in {self.name}.'
                    ' Time to reach 18 % oxygen concentration: '
                    f'{t_max_outage:.3g~}. Power outage time: '
                    f'{self.power.max_outage:.3g~}.'
                )

        # TODO Replace hard coded value with changeable
        # Constant leak and power failure
        l_power = TABLE_1['Electrical Power Failure']['Time rate']
        failure_rate = l_power
        # Constant leak is assumed to be hazardous only when
        # the isolation valve fails
        P_case = PFD_isol_valve
        Q_fan = 0 * ureg.ft**3/ureg.min  # No ventilation without power
        tau_event = min(leak.tau, self.power.max_outage)  # And Test_period?
        case = 'No power and isolation'
        self._add_failure_mode(case, P_case, tau_event, leak.fluid,
                               failure_rate, source, leak, Q_fan, 0, False)

        # Constant leak and ODH system failure
        failure_rate = lambda_ODH
        P_case = PFD_isol_valve
        tau_event = min(leak.tau, self.vent.Test_period)
        # Ventilation is available on ODH system failure
        case = 'No ODH'
        self._add_failure_mode(case, P_case, tau_event, leak.fluid,
                               failure_rate, source, leak,
                               self.vent.const_vent(FailCase.ODH),
                               0, False)

        # Fan failure
        # If at least 1 fan is down, all fans are considered non-operational
        failure_rate = TABLE_2['Fan']['Failure to run']
        P_case = PFD_isol_valve
        tau_event = min(leak.tau, self.vent.Test_period)
        case = 'No fan'
        self._add_failure_mode(case, P_case, tau_event, leak.fluid,
                               failure_rate, source, leak,
                               self.vent.const_vent(FailCase.FAN),
                               0, False)

        # Fans operational
        # If fatality is possible with just one fan, raise error
        Q_1fan = self.vent.Q_fan
        # Event can't go undetected for longer than the test period
        tau = min(leak.tau, self.vent.Test_period)
        O2_conc_1fan = conc_vent(self.volume, leak.q_std, Q_1fan, tau)
        F_1fan = self._fatality_prob(O2_conc_1fan)
        if F_1fan > 0:
            raise ConstLeakTooBig('Constant leak creates ODH with at least 1 '
                                  f'fan operational: {leak}')

    def _fatality_prob(self, O2_conc):
        """Calculate fatality probability for given oxygen concentration.

        The equation is fitted from the FESHM 4240 plot.

        Parameters
        ----------
        O2_conc : float
            Oxygen concentration.

        Returns
        -------
        float
            Fatality rate.
        """
        if O2_conc >= 0.18:  # Lowest oxygen concentration above 18%
            Fi = 0
        elif O2_conc <= 0.088:  # 8.8% of oxygen is assumed to be 100% fatal
            Fi = 1
        else:
            # Fi formula, reverse engineered using 8.8% and 18% thresholds
            Fi = 10**(6.5-76*O2_conc)
        return Fi

    def odh_class(self):
        """Calculate ODH class as defined in FESHM 4240.

        Returns
        -------
        int
            ODH class.
        """
        if self.phi < 1e-7/ureg.hr:
            return 0
        elif self.phi < 1e-5/ureg.hr:
            return 1
        elif self.phi < 1e-3/ureg.hr:
            return 2
        else:
            # TODO add a custom exception for ODH > 2
            print('ODH fatality rate is too high. Please, check calculations')
            return None

    @property
    def phi(self):
        return sum((fm.phi for fm in self.fail_modes))

    def report(self, brief=True, sens=None):
        """Print a report for failure modes and effects.

        The report is sorted by fatality rate descending."""
        self.fail_modes.sort(key=lambda x: x.phi, reverse=True)
        sens = sens or 0/ureg.hr
        title = f'ODH report for {self}'
        padding = len(title) + 10
        print('#'*padding)
        print(title)
        print('-'*padding)
        if brief:
            print('Printing brief ODH report')
            print(f'Only leaks with Fatality rate > {sens} are shown')
        for f_mode in self.fail_modes:
            if f_mode.phi >= sens or not brief:
                print()
                print(f' Source:               {f_mode.source.name}')
                print(f' Failure:              {f_mode.name}')
                print(f' Fatality rate:        {f_mode.phi.to(1/ureg.hr):.2~}')
                print(f' Building is powered:  {not f_mode.outage}')
                print(f' Oxygen concentration: {f_mode.O2_conc:.0%}, '
                      f'{f_mode.O2_conc/0.21:.0%} percent of norm')
                print(f' Leak failure rate:    {f_mode.leak_fr:.3g~}')
                print(' ODH protection PFD:    '
                      f'{(f_mode.P_i/f_mode.leak_fr).to(ureg.dimensionless):.2~}')
                print(f' Total failure rate:   {f_mode.P_i.to(1/ureg.hr):.2~}')
                print(f' Leak rate:            {f_mode.q_leak:.2~}')
                print(f' Event duration:       {f_mode.tau:.2~}')
                print(f' Fans working:         {f_mode.N_fan}')
                print(f' Fan rate:             {f_mode.Q_fan:.2~}')
                print(f' Fatality prob:        {f_mode.F_i:.0%}')

    def report_short_table(self, sens=None):
        """Prepare a short table for failure modes and effects.

        The report is sorted by fatality rate descending."""
        self.fail_modes.sort(key=lambda x: x.phi, reverse=True)
        sens = sens or 0/ureg.hr
        table = [
            ["Failure mode", "Fans on", "O_2", "Duration, min", "\\phi_i"],
            None]
        for f_mode in self.fail_modes:
            if f_mode.phi <= sens:
                break
            row = []
            row.append(f'{f_mode.source.name} {f_mode.name}')
            row.append(f'{f_mode.N_fan}')
            row.append(f'{float(f_mode.O2_conc):.0%}')
            row.append(f'{f_mode.tau.m_as(ureg.min):,.1f}')
            row.append(f'{f_mode.phi.m_as(1/ureg.hr):.2}')
            table.append(row)
        return table

    def report_table(self, filename='ODH_report', sort='name'):
        """Make a table with the calculation results."""
        table = []
        header = ['Source', 'Failure', 'Fluid',
                  'Event failure rate, 1/hr', '# of',
                  'Total failure rate, 1/hr', 'Leak rate, SCFM',
                  '# fans working', 'Fan rate, SCFM', 'Event duration, min',
                  'Oxygen concentration', 'Fatality prob', 'Case',
                  'Case prob ', 'Fatality rate, 1/hr']
        # 'Total failure rate', 'ODH protection PFD', 'Building is powered'
        table.append(header)
        if sort == 'name':
            self.fail_modes.sort(key=lambda x: x.source.name)
        elif sort == 'phi':
            self.fail_modes.sort(key=lambda x: x.phi, reverse=True)
        else:
            ODHError(f'Sort option {sort} not supported.')
        for f_mode in self.fail_modes:
            tau = f_mode.tau.m_as(ureg.min)
            if tau == float('inf'):
                # Handling infinite release
                tau = str(tau)
            table.append([
                f_mode.source.name,
                f_mode.name,
                str(f_mode.fluid),
                (f_mode.leak_fr/f_mode.N).m_as(1/ureg.hr),
                f_mode.N,
                f_mode.leak_fr.m_as(1/ureg.hr),
                f_mode.q_leak.m_as(ureg.ft**3/ureg.min),
                f_mode.N_fan,
                f_mode.Q_fan.m_as(ureg.ft**3/ureg.min),
                tau,
                f_mode.O2_conc,
                f_mode.F_i,
                f_mode.case,
                f_mode.P_i/f_mode.leak_fr,
                f_mode.phi.m_as(1/ureg.hr)])
        wb_options = {}
        filename += '.xlsx'
        with xlsxwriter.Workbook(filename, wb_options) as workbook:
            worksheet = workbook.add_worksheet()
            col_width = [len(x) for x in table[0]]
            for row_n, row in enumerate(table):
                for col_n, data in enumerate(row):
                    worksheet.write(row_n, col_n, data)
                    if col_n in (0, 1, 2, 11, 12):
                        # For source names, failure names
                        # and 'Total failure rate'
                        col_width[col_n] = max(col_width[col_n], len(str(data)))
            header_format = workbook.add_format({'bold': True,
                                                 'font_size': 12,
                                                 'bottom': 3})
            sci_format = workbook.add_format({'num_format': '0.00E+00'},)
            percent_format = workbook.add_format({'num_format': '0.0%'},)
            number_format = workbook.add_format({'num_format': '0'},)
            worksheet.set_row(0, None, header_format)
            worksheet.set_column(3, 3, None, sci_format)  # Event fr
            worksheet.set_column(5, 5, None, sci_format)  # Total event fr
            worksheet.set_column(6, 6, None, sci_format)  # Leak rate
            worksheet.set_column(9, 9, None, sci_format)  # Duration
            worksheet.set_column(10, 10, None, percent_format)  # O2_conc
            worksheet.set_column(12, 14, None, sci_format)  # Fatality rate
            # Writing total/summary
            N_rows = len(table)
            N_cols = len(table[0])
            worksheet.write(N_rows+1, N_cols-2, 'Total fatality rate, 1/hr')
            worksheet.write(N_rows+1, N_cols-1,
                            self.phi.m_as(1/ureg.hr))
            worksheet.write(N_rows+2, N_cols-2, 'ODH class')
            worksheet.write(N_rows+2, N_cols-1, self.odh_class(),
                            number_format)
            # Autofit column width
            for col_n, width in enumerate(col_width):
                adj_width = width - 0.005 * width**2
                worksheet.set_column(col_n, col_n, adj_width)
            # Adding usability
            worksheet.conditional_format(
                1, N_cols-1, N_rows-1, N_cols-1,
                {'type': '3_color_scale', 'min_color': '#008000',
                 'max_color': '#FF0000'})
            worksheet.freeze_panes(1, 0)

    def info(self):
        print(f'ODH volume: {self.name}, {self.volume.to(ureg.ft**3):~}')
        print(f'Building power probability of failure: {self.power.pfd:.1e~}')
        print(f'Building power max outage period: {self.power.max_outage:.3g~}')
        print(f'Fan flow rate: {self.vent.Q_fan:.0f~}')
        print(f'Number of fans: {self.vent.N_fans:.0f}')
        print(f'ODH system test period: {self.vent.Test_period:.0f~}')
        print(f'Fan failure rate: {self.vent.lambda_fan:.1e~}')
        print(f'Building fresh air supply rate: {self.vent.min_vent:.3g~}')

    # def source_safe(self, source, escape = True):
    #    """
    #    Estimate the impact of the Source volume on oxygen concetration. Smaller sources might not be able to drop oxygen concentration to dangerous levels.
    #    """
    #    if escape == True:  # if mixed air is allowed to escape within considered volume
    #        O2_conc = 0.21*self.volume/(self.volume+source.volume)
    #    else:  # worst case; inert gas is trapped and expells the air outside the considered volume
    #        O2_conc = 0.21*(1-source.volume/self.volume)
    #    return self._fatality_prob(O2_conc) == 0


def prob_m_of_n(m, n, T, l):
    """Calculate the probability of m out of n units working.

    Calculation is done using binomial distribution.

    Parameters
    ----------
    m : int
        Number of units working.
    n : int
        Total number of units.
    T : ureg.Quantity {time: 1}
        Test period
    l : ureg.Quantity {time: -1}
        Failure rate (\\lambda) of a fan

    Returns
    -------
    float
        Probability of m out of n units working.
    """
    PFD_one_unit = l*T
    m_of_n = binom(n, m) * (PFD_one_unit)**(n-m) * (1-PFD_one_unit)**m
    return m_of_n.to_base_units()


def conc_vent(V, R, Q, t):
    """Calculate the oxygen concentration at the end of the event.

    As defined by FESHM 4240 6.1.A, Cases A, B, and C.

    Parameters
    ----------
    V : ureg.Quantity {length: 3}
        Volume of the confined space.
    R : ureg.Quantity {length: 3, time: -1}
        Volumetric spill rate into confined space.
    Q : ureg.Quantity {length: 3, time: -1}
        Volumetric ventilation rate of fan(s); positive value corresponds
        to blowing air into the confined space, negative - drawing contaminated
        air outside.
    t : ureg.Quantity {time: 1}
        time, beginning of release is at `t` = 0.

    Returns
    -------
    float
        Oxygen concentration.
    """
    if Q > 0:
        C = 0.21/(Q+R) * (Q+R*math.e**-((Q+R)/V*t))
    elif abs(Q) <= R:
        C = 0.21*math.e**-(R/V*t)
    elif abs(Q) > R:
        C = 0.21*(1-R/abs(Q)*(1-math.e**-(abs(Q)*t/V)))
    return C


def conc_inf(V, R, Q):
    """Calculate the final oxygen concentration for continuous flow.

    Equivalent to conc_vent(V, R, Q, float('inf')).

    Parameters
    ----------
    V : ureg.Quantity {length: 3}
        Volume of the confined space.
    R : ureg.Quantity {length: 3, time: -1}
        Volumetric spill rate into confined space.
    Q : ureg.Quantity {length: 3, time: -1}
        Volumetric ventilation rate of fan(s); positive value corresponds
        to blowing air into the confined space, negative - drawing contaminated
        air outside.

    Returns
    -------
    float
        Oxygen concentration.
    """
    if Q > 0:
        C = 0.21/(Q+R)*Q
    elif abs(Q) <= abs(R):
        C = 0
    elif abs(Q) > abs(R):
        C = 0.21*(1-R/abs(Q))
    return C


def conc_after(V, C_e, Q, t, t_e):
    """Calculate the oxygen concentration in the confined volume after
    the release has ended.

    As defined by FESHM 4240 6.1.A, Case D.

    Parameters
    ----------
    V : ureg.Quantity {length: 3}
        Volume of the confined space.
    C_e : float
        Oxygen concentration at the end of the release.
    Q : ureg.Quantity {length: 3, time: -1}
        Volumetric ventilation rate of fan(s); positive value corresponds
        to blowing air into the confined space, negative - drawing contaminated
        air outside.
    t : ureg.Quantity {time: 1}
        time, beginning of release is at `t`=0.
    t_e : ureg.Quantity {time: 1}
        time when release ended.

    Returns
    -------
    float
        Oxygen concentration.
    """
    C = 0.21-(0.21-C_e)*math.e**-(abs(Q)/V*(t-t_e))
    return C


def print_result(*Volumes):
    """Print the results of the ODH analysis for a volume.

    If several volumes given (in case of interlapping volumes) the worst case
    will be printed.
    """
    max_phi = -1/ureg.hr
    for volume in Volumes:
        if volume.phi > max_phi:
            max_volume = volume
    line_1 = '# Fatality rate for {} is {:.1e}  # '.format(max_volume, volume.phi)
    pad = len(line_1)
    line_2 = '# Recommended ODH class {}'.format(max_volume.odh_class()).ljust(pad-1)+'#'
    print('#'*pad)
    print(line_1)
    print(line_2)
    print('#'*pad)


def hole_leak(area, fluid, P_out=P_NTP):
    """Calculate leak flow rate through a hole in a pipe.

    Discharge flow through a hole is calculated using Darcy equation for
    compressible fluid using expansion factor Y Crane TP 410 2013 eq. 4-15.
    Resistance coefficient K of the orifice taken from Rennels, Pipe Flow,
    2012, eq. 13.5.
    For this calculation the gas is assumed to have no pressure loss on
    entry. This makes analysis simple and conservative.

    Parameters
    ----------
    tube : heat_transfer.Tube or None
        Tube upstream of the leak orifice. None if leaking from a vessel.
    area : ureg.Quantity {length: 2}
        Area of the leak.
    fluid : heat_transfer.ThermState
        Thermodynamic state of the fluid stored in the source.

    Returns
    -------
    ureg.Quantity {length: 3, time: -1}
        Standard volumetric flow at Normal Temperature and Pressure.

    """
    G_max = G_nozzle(fluid, P_out=P_out)[0]
    m_dot = G_max * area
    q_std = to_standard_flow(m_dot, fluid)
    return q_std
