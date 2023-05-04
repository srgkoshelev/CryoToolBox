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


class ConstLeakNoVent(Exception):
    pass

@ureg.check(None, '1/[time]', None, '[length]^3/[time]', '[time]', None)
@dataclass
class Leak:
    """Describe inert gas leak from a Source"""
    name: str
    failure_rate: ureg.Quantity
    fluid: ThermState
    q_std: ureg.Quantity
    tau: ureg.Quantity
    N_events: int
    _is_const = False


@ureg.check(None, None, '[length]^3/[time]', None)
@dataclass
class ConstLeak:
    name: str
    fluid: ThermState
    q_std: ureg.Quantity
    N_events: int
    _is_const = True
    # TODO Remove once isinstance() issue resolved


class Source:
    """Source of inert gas
    """
    def __init__(self, name, fluid, volume, N=1, PFD_isol_valve=1):
        """Define the possible source of inert gas.

        Parameters
        ----------
        name : str
            Name of the source.
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        volume : ureg.Quantity, [length]^3
            Volume of the fluid stored.
        N : int
            Number of the sources if several similar sources exist,
            e.g. gas bottles.
        PFD_isol_valve : float
            Probability of failure on demand of source isolation valve.
            If no isolation valve will trigger on ODH alarm, the probability
            is 1 (default).
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
        self.volume.ito_base_units()
        # By default assume there is no isolation valve
        # that is used by ODH system
        self.PFD_isol_valve = PFD_isol_valve

    def add_pipe_failure(self, tube, fluid, q_std_rupture=None, N_welds=1):
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
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid for the release.
        q_std_rupture : ureg.Quantity, [length]^3/[time]
            Standard volumetric flow rate for pipe rupture.
        N_welds : int
            Number of welds on the tube.
        """
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
                    self.add_failure_mode(name, failure_rate, fluid, q_std,
                                          N_events)

    def add_dewar_insulation_failure(self, q_std):
        """Add dewar insulation failure to leaks dict.

        Store failure rate, flow rate and expected time duration of the
        failure event for the dewar insulation failure. Based on FESHM4240.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        q_std : ureg.Quantity, [length]^3/[time]
            Standard volumetric flow rate of the relief for the case of
            dewar insulation failure.
        """
        failure_rate = TABLE_1['Dewar']['Loss of vacuum']
        self.add_failure_mode('Dewar insulation failure', failure_rate,
                              self.fluid, q_std, 1)

    def add_flange_failure(self, pipe, fluid, q_std_rupture=None, N=1,
                           gasket_type='reinforced', blowout_area=None):
        """Add reinforced or preformed gasket flange failure
        to leaks dict.

        Store failure rate, flow rate and expected time duration of
        the event for transfer line failure. Based on FESHM 4240.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        pipe : heat_transfer.Pipe
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        q_std_rupture : ureg.Quantity, [length]^3/[time]
            Standard volumetric flow rate for flange rupture.
        N : int
            Number of reinforced seal connections on the pipe.
        gasket_type : str
            Type of the gasket in the flange, 'soft' or 'reinforced'.
        blowout_area : Quantity [length^2]
            Leak area in case of the soft gasket blowout. Not applicable for reinforced gaskets.
        """
        # TODO Make leak and rupture areas adjustable, add info to docstring
        if gasket_type == 'reinforced':
            table = TABLE_2['Flange, reinforced gasket']
        if gasket_type == 'soft':
            table = TABLE_2['Flange, soft gasket']
        # Leak case
        name = f'Flange {gasket_type} gasket leak: {pipe}'
        failure_rate = table['Leak']['Failure rate']
        area = table['Leak']['Area']
        q_std = hole_leak(area, fluid)
        self.add_failure_mode(name, failure_rate, fluid, q_std, N)
        # Blowout for soft gasket
        if gasket_type == 'soft':
            name = f'Flange soft gasket blowout: {pipe}'
            failure_rate = table['Blowout']['Failure rate']
            area = blowout_area or table['Blowout']['Area']
            area = min(area, pipe.area)
            q_std = hole_leak(area, fluid)
            self.add_failure_mode(name, failure_rate, fluid, q_std, N)
        # Rupture
        name = f'Flange {gasket_type} gasket rupture: {pipe}'
        failure_rate = table['Rupture']
        if q_std_rupture is not None:
            q_std = q_std_rupture
        else:
            area = pipe.area
            q_std = hole_leak(area, fluid)
        self.add_failure_mode(name, failure_rate, fluid, q_std, N)

    def add_valve_failure(self, pipe, fluid, q_std_rupture=None, N=1):
        """Add valve leak and rupture failure modes to leaks dict.

        Store failure rate, flow rate and expected time duration of
        the event for transfer line failure. Based on FESHM 4240.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        pipe : heat_transfer.Pipe
            Pipe/tube upstream the valve.
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        q_std_rupture : ureg.Quantity, [length]^3/[time]
            Standard volumetric flow rate for flange rupture.
        N : int
            Number of valves
        """
        table = TABLE_2['Valve, pneumatic']
        # Leak case
        name = f'Valve leak: {pipe}'
        failure_rate = table['External leak']
        # Using area from flange leak for consistency
        area = TABLE_2['Flange, soft gasket']['Leak']['Area']
        q_std = hole_leak(area, fluid)
        self.add_failure_mode(name, failure_rate, fluid, q_std, N)
        # Rupture
        name = f'Valve rupture: {pipe}'
        failure_rate = table['Rupture']
        if q_std_rupture is not None:
            q_std = q_std_rupture
        else:
            area = pipe.area
            q_std = hole_leak(area, fluid)
        self.add_failure_mode(name, failure_rate, fluid, q_std, N)

    def add_line_failure(self, pipe, fluid, *, N_welds, N_reinforced, N_soft, N_valves,
                         q_std_rupture=None, blowout_area=None):
        """Add leaks for pipe, weld, flange, and valve failures.

        Store failure rate, flow rate and expected time duration of
        the event for pipe, weld, flange, and valve failures. Based on
        FESHM 4240. Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        pipe : heat_transfer.Pipe
            Pipe/tube upstream the valve.
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        N_welds : int
            Number of welds on the line.
        N_reinforced : int
            Number of flanges with reinforced, preformed or metal gaskets on the line.
        N_soft : int
            Number of flanges with soft gaskets on the line.
        N_valves : int
            Number of valves on the line.
        q_std_rupture : ureg.Quantity, [length]^3/[time]
            Standard volumetric flow rate for flange rupture.
        blowout_area : Quantity [length^2]
            Leak area in case of the soft gasket blowout. Not applicable for reinforced gaskets.
        """
        if N_welds > 0:
            self.add_pipe_failure(pipe, fluid, q_std_rupture, N_welds)
        if N_reinforced > 0:
            self.add_flange_failure(pipe, fluid, q_std_rupture, N=N_reinforced, gasket_type='reinforced')
        if N_soft > 0:
            self.add_flange_failure(pipe, fluid, q_std_rupture, N=N_soft, gasket_type='soft', blowout_area=blowout_area)
        if N_valves > 0:
            self.add_valve_failure(pipe, fluid, q_std_rupture, N=N_valves)

    # def transfer_line_failure(self, pipe, fluid, q_std_rupture=None, N=1):
    #     """Add transfer line failure to leaks dict.

    #     For a given tube calculate leak parameters as following:
    #     - failure rate
    #     - standard volumetric flow
    #     - duration of the release
    #     - number of possible similar events
    #     for all failure modes for transfer line listed in Table 1 of
    #     FESHM 4240. Note that Rationale for Table 1: “Fermilab Equipment
    #     Failure Rate Estimates” clearly states that only bayonet failures
    #     are considered for these failure modes. Therefore only bayonet leak and
    #     blowout (rupture) are considered. Leak area is not defined in
    #     FESHM 4240 and is defined globally as `TRANSFER_LINE_LEAK_AREA`.

    #     Failure modes are analyzed by `Volume.odh` method.

    #     Parameters
    #     ----------
    #     pipe : piping.PipingElement
    #     fluid : cp_wrapper.ThermState
    #         Thermodynamic state of the fluid stored in the source.
    #     q_std_rupture : ureg.Quantity {length: 3, time: -1}
    #         Standard volumetric flow rate for fluid line rupture.
    #     N : int
    #         Number of bayonets/soft seals on the transfer line.
    #     """
    #     # TODO This should be replaced by flange failure at some point
    #     # Leak case
    #     name = f'Fluid line gasket leak: {pipe}'
    #     failure_rate = TABLE_1['Fluid line']['Leak']
    #     area = TRANSFER_LINE_LEAK_AREA
    #     q_std = hole_leak(area, fluid)
    #     self.add_failure_mode(name, failure_rate, q_std, N)
    #     # Rupture case
    #     name = f'Fluid line gasket rupture: {pipe}'
    #     failure_rate = TABLE_1['Fluid line']['Rupture']
    #     if q_std_rupture is not None:
    #         q_std = q_std_rupture
    #     else:
    #         area = pipe.area
    #         q_std = hole_leak(area, fluid)
    #     self.add_failure_mode(name, failure_rate, q_std, N)

    def add_pressure_vessel_failure(self, relief_area=None):
        """Add pressure vessel failure to leaks dict.

        Store failure rate, flow rate and expected time duration of
        the event for transfer line failure. Based on FESHM 4240.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        relief_area : ureg.Quantity, [length]^2
            Vacuum jacket relief area if the vessel has one, None otherwise.
        """
        # Leak case
        name = 'Pressure vessel leak'
        area = TABLE_2['Vessel, pressure']['Small leak']['Area']
        failure_rate = TABLE_2['Vessel, pressure']['Small leak']['Failure rate']
        q_std = hole_leak(area, self.fluid)
        self.add_failure_mode(name, failure_rate, self.fluid, q_std, 1)

        # Rupture case
        name = 'Pressure vessel rupture'
        if relief_area is None:
            # No vacuum jacket on the vessel
            area = 1000 * ureg.mm**2  # Equal to large leak of a pipe
        else:
            area = relief_area
        failure_rate = TABLE_2['Vessel, pressure']['Failure']
        q_std = hole_leak(area, self.fluid)
        self.add_failure_mode(name, failure_rate, self.fluid, q_std, 1)

    def add_const_leak(self, name, fluid, q_std, N=1):
        """Add constant leak to leaks dict.

        Store flow rate and expected time duration of the
        failure event for general failure mode.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        name : str
            Name of the failure mode
        fluid : ThermState
            Thermodynamic state of the fluid at the leak location.
        q_std : ureg.Quantity, [length**3/time]
            Standard volumetric flow rate.
        N : int
            Number of leaks.
        """
        if not fluid:
            fluid = self.fluid
        # N_events = N * self.N
        self.leaks.append(ConstLeak(name, fluid, q_std*N_events, 1))

    def add_failure_mode(self, name, failure_rate, fluid, q_std, N=1):
        """Add general failure mode to leaks dict.

        Store failure rate, flow rate and expected time duration of the
        failure event for general failure mode.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        name : str
            Name of the failure mode
        failure rate : ureg.Quantity, 1/[time]
            Failure rate of the failure mode,
            i.e. how often the failure occurs
        fluid : ThermState
            Fluid state at the release location.
        q_std : ureg.Quantity, [length]^3/[time]
            Standard volumetric flow rate.
        N : int
            Number of similar failure modes.
        """
        N_events = N * self.N
        tau = self.volume/q_std
        total_failure_rate = N_events*failure_rate
        total_failure_rate.ito(1/ureg.hr)
        # self.leaks.append((name, total_failure_rate, q_std, tau.to(ureg.min),
        #                    N_events))
        self.leaks.append(Leak(name, total_failure_rate, fluid, q_std,
                               tau.to(ureg.min), N_events))

    @staticmethod
    def combine(name, fluid, sources, N=1, PFD_isol_valve=1):
        """Combine several ODH sources sharing volume.

        Can be used for failure modes affecting several sources in parallel.

        Parameters
        ----------
        name : str
            Name of the new combined source.
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        sources : list of Source
            Sources connected together.
        N : int
            Number of the sources if several similar sources exist,
            e.g. gas bottles.
        PFD_isol_valve : float
            Probability of failure on demand of source isolation valve.
            If no isolation valve will trigger on ODH alarm, the probability
            is 1 (default).

        Returns
        -------
        Source
            Combined source of inert gas.
        """
        if not all([source.fluid.name == fluid.name for source in sources]):
            ODHError('All volumes should contain the same fluid.')
            return None
        if sources[0].fluid.name != fluid.name:
            ODHError('New source should have the same fluid as combined.')
            return None
        volume_NTP = sum([source.N*source.volume for source in sources])
        fluid_NTP = ThermState(fluid.name, P=P_NTP, T=T_NTP)
        phys_volume = volume_NTP*fluid_NTP.Dmass / fluid.Dmass
        return Source(name, fluid, phys_volume, N=N, PFD_isol_valve=PFD_isol_valve)

    def __repr__(self):
        return f'Source({self.name!r}, {self.fluid!r}, {self.volume!r}, N={int(self.N)}, PFD_isol_valve={float(self.PFD_isol_valve)})'

    def __str__(self):
        return (f'{self.name}, '
            f'{self.volume:,.0f~} '
            f'of {self.fluid.name}')

    def print_leaks(self):
        """Print information on the leaks defined for the source."""
        print(f'Printing leaks for {self}')
        print()
        for leak in self.leaks:
            print(leak.name)
            # TODO This will break for constant leak
            print(f'Total failure rate: {leak.failure_rate:.3g~}')
            print(f'Leak rate: {leak.q_std:.3g~}')
            print(f'Event duration: {leak.tau:.3g~}')
            print(f'Number of events: {leak.N_events}')
            print()


@ureg.check(None, None, None, '1/[time]', None, '1/[time]', '1/[time]',
            None, None, '[length]^3/[time]', '[time]',
            '[length]^3/[time]', None, None, None, None)
@dataclass
class FailureMode:
    """Describes ODH failure modes, a combination of a leak and response to it."""
    name: str
    source: Source
    fluid: ThermState
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
    is_const: bool
    scenario: str


@ureg.check(None, None, '1/[time]', '[time]')
@dataclass
class BuildPower:
    """
    Dataclass representing the electrical power failure scenario of a building.

    Attributes
    ----------
    power_pfd : float
        The probability of failure demand rate.
    backup_pfd : float
        The probability of backup system failure.
    lambda_power : ureg.Quantity, 1/[time]
        The rate of power failure occurrences.
    max_outage : ureg.Quantity, [time]
        The maximum duration of an outage.

    Properties
    ----------
    pfd : float
        The total probability of a power failure, taking into account both
        the power probability of failure on demand and the backup system
        reliability.

    Methods
    -------
    None
    """
    power_pfd: float = TABLE_1['Electrical Power Failure']['Demand rate']
    backup_pfd: float = 1
    lambda_power: ureg.Quantity = TABLE_1['Electrical Power Failure']['Time rate']
    max_outage: ureg.Quantity = float('inf') * ureg.hr

    @property
    def pfd(self):
        return self.power_pfd * self.backup_pfd

    def info(self):
        return (f'Power provided with {self.pfd:.3g~} unavailability and '
                f'{self.max_outage:.3g~} max outage period.'
                f'\nPower failure rate is {self.lambda_power:.3e~}.')

@ureg.check('[length]^3/[time]', None, '[time]', '1/[time]',
            '[length]^3/[time]', None)
@dataclass
class BuildVent:
    """
    Dataclass representing the ventilation system of a building.

    Attributes
    ----------
    Q_fan : ureg.Quantity, [length]^3/[time]
        Volumetric flow of a single ODH fan installed in the volume.
    N_fans : int
        Number of fans installed.
    Test_period : ureg.Quantity, [time]
        Test period of the fans and dampers.
    lambda_fan : ureg.Quantity, 1/[time]
        Failure rate of the fans in the building.
    const_vent : ureg.Quantity, [length]^3/[time]
        Min volumetric flow required or present in the building.
    PFD_ODH : float
        The default value for the oxygen deficiency hazard (ODH) system failure.

    Properties
    ----------
    None

    Methods
    -------
    fan_flowrates() -> List[Tuple[float, ureg.Quantity, int]]:
        Calculate (probability, flow, number of fans) tuples for all
        combinations of fans working. All fans are expected to have the same
        volumetric flow.
    """

    Q_fan: ureg.Quantity
    N_fans: int
    Test_period: ureg.Quantity
    lambda_fan: ureg.Quantity = TABLE_2['Fan']['Failure to run']
    const_vent: ureg.Quantity = 0*ureg.ft**3/ureg.min
    PFD_ODH: float = PFD_ODH  # Default value for ODH system failure

    def fan_flowrates(self):
        # TODO add fans with different volumetric rates (see report as well)
        fail_rate = self.lambda_fan
        fan_flowrates = []
        for m in range(self.N_fans+1):
            # Probability of exactly m units starting
            P_m_fan_work = prob_m_of_n(m, self.N_fans, self.Test_period,
                                       fail_rate)
            flowrate = self.Q_fan*m
            if flowrate == Q_('0 m**3/min'):
                flowrate = self.const_vent
            fan_flowrates.append((P_m_fan_work, flowrate, m))
        return fan_flowrates

    def info(self):
        if self.Q_fan > 0*ureg.L/ureg.s:
            vent_action = 'supplying clean air to'
        else:
            vent_action = 'drawing contaminated air from'
        result = [f'Ventilation consists of {self.N_fans} fans {abs(self.Q_fan):,.0f~} each\n  {vent_action} the volume.']
        if self.const_vent > 0*ureg.L/ureg.s:
            vent_action = 'supplying clean air to'
        else:
            vent_action = 'drawing contaminated air from'
        result.append(f'Constant ventilation in the building is {abs(self.const_vent):,.1f~}\n  {vent_action} the volume.')
        result.append(f'ODH system is tested every {self.Test_period:.1f~} and fan failure rate is {self.lambda_fan:.3e~}.')
        result.append(f'ODH system unavailability is {self.PFD_ODH:.3e}.')
        return '\n'.join(result)


class Volume:
    """
    A class to represent a building or a part of a building that is affected by
    inert gas leaks.

    Parameters
    ----------
    name : str
        The name of the volume.
    volume : ureg.Quantity, [length]^3
        The volume of the building or part of the building.

    Attributes
    ----------
    name : str
        The name of the volume.
    volume : ureg.Quantity, [length]^3
        The volume of the building or part of the building.
    vent : BuildVent
        The ventilation system of the building or part of the building.
    power : BuildPower
        The power system of the building or part of the building.

    Methods
    -------
    info()
        Display information on the volume.
    odh(sources, power_outage=False)
        Calculate ODH fatality rate for given `Source`s.
    _fatality_no_power(source, leak, power_outage)
        Calculate fatality rate in the event of power failure.
    _fatality_no_odh(source, leak)
        Calculate fatality rate in the event of ODH system failure.
    _fatality_with_fans(source, leak)
        Calculate fatality rate for any number of fans operational.
    _add_failure_mode(P_event, tau_event, failure_rate, source, leak, Q_fan,
        N_fans, power_outage)
        Add a failure mode.
    _fatality_const_leak(source, leak, power_outage)
        Calculate fatality rates for constant leak case.
    _fatality_prob(O2_conc)
        Calculate fatality probability for given oxygen concentration.
    """
    def __init__(self, name, volume, *, build_vent, build_power):
        """Define a volume affected by inert gas release from  a `Source`.

        Parameters
        ----------
        name : str
            Name of the volume.
        volume : ureg.Quantity, [length]^3
            Volume of the building or part of the building.
        """
        self.name = name
        self.volume = volume.to_base_units()
        self.vent = build_vent
        self.power = build_power

    def info(self):
        """
        Display information on the volume.
        """
        print(f'{self.name} with affected volume of {self.volume:,.0f~}.')
        print(self.power.info())
        print(self.vent.info())

    def odh(self, sources, power_outage=False):
        """Calculate ODH fatality rate for given `Source`s.

        For each leak of each source ODH conditions are analyzed and
        fatality rates are calculated. The results are collected in
        fail_modes list.

        Parameters
        ----------
        sources : List[Source]
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

    def _fatality_no_power(self, source, leak, power_outage):
        """Calculate fatality rate in the event of power failure.

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
            tau_event = leak.tau
            scenario = 'Planned power outage'
        else:
            PFD_power = self.power.pfd
            tau_event = min(leak.tau, self.power.max_outage)
            scenario = 'Power failure'
        PFD_isol_valve = source.PFD_isol_valve
        P_event = PFD_power * PFD_isol_valve
        Q_fan = 0 * ureg.ft**3/ureg.min
        self._add_failure_mode(P_event, tau_event, leak.failure_rate, source,
                               leak, Q_fan, 0, power_outage, scenario)

    def _fatality_no_odh(self, source, leak):
        """Calculate fatality rate in the event of ODH system failure.

        Power has to be  on for ODH system to fail. Source isolation is
        possible for this event. Only ventilation independent from ODH system
        is available, e.g., HVAC.

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
        PFD_isol_valve = source.PFD_isol_valve
        P_event = (1-self.power.pfd) * self.vent.PFD_ODH * PFD_isol_valve
        Q_fan = self.vent.const_vent
        tau_event = min(leak.tau, self.vent.Test_period)
        scenario = 'ODH system failure'
        self._add_failure_mode(P_event, tau_event, leak.failure_rate, source,
                               leak, Q_fan, 0, False, scenario)

    def _fatality_with_fans(self, source, leak):
        """Calculate fatality rate for any number of fans operational.

        Power has to be  on, and ODH system operational. Source isolation is
        possible for this event. Ventilation independent from ODH system,
        e.g., HVAC, is considered for all fans unavailable at the time of
        the event.

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
        PFD_isol_valve = source.PFD_isol_valve
        # Probability for this group of the events
        P_group = (1-self.power.pfd) * (1-self.vent.PFD_ODH) * PFD_isol_valve
        tau_event = min(leak.tau, self.vent.Test_period)
        for (P_fan, Q_fan, N_fans) in self.vent.fan_flowrates():
            P_event = P_group * P_fan  # Probability of the particular event
            if Q_fan == 0 * ureg.ft**3/ureg.min:
                Q_fan = self.vent.const_vent
            if N_fans == self.vent.N_fans:
                scenario = 'Normal operation'
            else:
                scenario = 'Fan failure'
            self._add_failure_mode(P_event, tau_event,
                                   leak.failure_rate, source, leak, Q_fan,
                                   N_fans, False, scenario)

    def _add_failure_mode(self, P_event, tau_event, failure_rate, source, leak,
                          Q_fan, N_fans, power_outage, scenario):
        """Add a failure mode.

        Parameters
        ----------
        P_event : float
            Probability of a particular state of the system, assuming the leak
            has happened. Alternatively, probability of a constant leak or 1.
        failure_rate : ureg.Quantity, 1/[time]
            Failure rate of the leak. Alternatively, system failure rate for
            a constant leak.
        tau_event : ureg.Quantity, [time]
            Maximum duration of the event.
        source : Source
            Inert gas source.
        leak : Leak
            Leak failure rate, volumetric flow rate, event duration, and number
            of events.
        fluid : ThermState
            Fluid conditions at the release location.
        Q_fan : ureg.Quantity, [length]^3/[time]
            Combined volumetric flow of ODH fans for the event.
        N_fans : int
            Number of fans operational during the event.
        power_outage : bool
            Shows whether there is a power outage is in effect.
        scenario: str
            Name of the response scenario, e.g., power failure, ODH system failure,
            normal operation.

        Returns
        -------
        None
        """
        P_i = failure_rate * P_event
        O2_conc = conc_vent(self.volume, leak.q_std, Q_fan, tau_event)
        F_i = self._fatality_prob(O2_conc)
        phi_i = P_i*F_i
        self.fail_modes.append(
            FailureMode(leak.name, source, leak.fluid, phi_i, O2_conc,
                        failure_rate, P_i, F_i, power_outage, leak.q_std,
                        tau_event, Q_fan, N_fans, leak.N_events, leak._is_const,
                        scenario)
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

        PFD_isol_valve = source.PFD_isol_valve
        if power_outage:
            tau_outage = source.volume/leak.q_std  # No refills during outage
            Q_fan_outage = 0 * ureg.ft**3/ureg.min  # No ventilation without power
            O2_conc_outage = conc_vent(self.volume, leak.q_std, Q_fan_outage, tau_outage)
            if O2_conc_outage < 0.195:
                raise ConstLeakNoVent(f'Constant leak {leak.name} from '
                                    f'{source.name} is not mitigated during '
                                    f'power outage in {self.name}.')

        # TODO Replace hard coded value with changeable
        # Constant leak and power failure
        failure_rate = self.power.lambda_power * self.power.backup_pfd
        # Constant leak is assumed to be hazardous only when
        # the isolation valve fails
        P_event = PFD_isol_valve
        Q_fan = 0 * ureg.ft**3/ureg.min  # No ventilation without power
        tau_event = min(self.vent.Test_period, self.power.max_outage)
        scenario = 'Const leak and power failure'
        self._add_failure_mode(P_event, tau_event, failure_rate,
                               source, leak, Q_fan, 0, False,
                               scenario)

        # Constant leak and ODH system failure
        failure_rate = lambda_ODH
        P_event = PFD_isol_valve
        # Event can't go undetected for longer than the test period
        tau_event = self.vent.Test_period
        scenario = 'Const leak and ODH system failure'
        # Ventilation is available on ODH system failure
        self._add_failure_mode(P_event, tau_event, failure_rate,
                               source, leak, self.vent.const_vent, 0, False,
                               scenario)

        # Fan failure
        # If at least 1 fan is down, all fans are considered non-operational
        failure_rate = TABLE_2['Fan']['Failure to run']
        P_event = PFD_isol_valve
        # Event can't go undetected for longer than the test period
        tau_event = self.vent.Test_period
        scenario = 'Const leak and 1 fan failure'
        self._add_failure_mode(P_event, tau_event, failure_rate,
                               source, leak, self.vent.const_vent, 0, False,
                               scenario)

        # Fans operational
        # If fatality is possible with just one fan, raise error
        Q_1fan = self.vent.Q_fan
        tau = math.inf * ureg.s
        O2_conc_1fan = conc_vent(self.volume, leak.q_std, Q_1fan, tau)
        if O2_conc_1fan < 0.195:
            raise ConstLeakTooBig('Constant leak reduces O2 concentration to '
                                  f'{O2_conc_1fan:.1%} with at least 1 fan '
                                  f'operational: {leak}')

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

    @property
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
            raise ODHError('ODH fatality rate is too high. Please, check calculations.')

    @property
    def phi(self):
        return sum((fm.phi for fm in self.fail_modes))

    # def report(self, brief=True, sens=None):
    #     """Print a report for failure modes and effects.

    #     The report is sorted by fatality rate descending."""
    #     self.fail_modes.sort(key=lambda x: x.phi, reverse=True)
    #     sens = sens or 0/ureg.hr
    #     title = f'ODH report for {self}'
    #     padding = len(title) + 10
    #     print('#'*padding)
    #     print(title)
    #     print('-'*padding)
    #     if brief:
    #         print('Printing brief ODH report')
    #         print(f'Only leaks with Fatality rate > {sens} are shown')
    #     for f_mode in self.fail_modes:
    #         if f_mode.phi >= sens or not brief:
    #             print()
    #             print(f' Source:               {f_mode.source.name}')
    #             print(f' Failure:              {f_mode.name}')
    #             print(f' Fatality rate:        {f_mode.phi.to(1/ureg.hr):.2~}')
    #             print(f' Building is powered:  {not f_mode.outage}')
    #             print(f' Oxygen concentration: {f_mode.O2_conc:.0%}, '
    #                   f'{f_mode.O2_conc/0.21:.0%} percent of norm')
    #             print(f' Leak failure rate:    {f_mode.leak_fr:.3g~}')
    #             print(' ODH protection PFD:    '
    #                   f'{(f_mode.P_i/f_mode.leak_fr).to(ureg.dimensionless):.2~}')
    #             print(f' Total failure rate:   {f_mode.P_i.to(1/ureg.hr):.2~}')
    #             print(f' Leak rate:            {f_mode.q_leak:.2~}')
    #             print(f' Event duration:       {f_mode.tau:.2~}')
    #             print(f' Fans working:         {f_mode.N_fan}')
    #             print(f' Fan rate:             {f_mode.Q_fan:.2~}')
    #             print(f' Fatality prob:        {f_mode.F_i:.0%}')

    def make_doc_table(self, sens=None):
        """Make a short table for failure modes and effects to be included to a document.

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
            row.append(f'{f_mode.O2_conc:.0%}')
            row.append(f'{f_mode.tau.m_as(ureg.min):,.1f}')
            row.append(f'{f_mode.phi.m_as(1/ureg.hr):.2e}')
            table.append(row)
        return table

    def make_excel_table(self, filename='ODH_report'):
        """Make a table with the calculation results."""
        table = []
        header = ['Source', 'Failure', 'Fluid', 'Event failure rate, 1/hr',
                  '# of', 'Total failure rate, 1/hr', 'Leak rate, SCFM',
                  '# fans working', 'Fan rate, SCFM', 'Event duration, min',
                  'Oxygen concentration', 'Fatality prob', 'Scenario', 'Scenario prob',
                  'Fatality rate, 1/hr']
        # 'Total failure rate', 'ODH protection PFD', 'Building is powered'
        table.append(header)
        self.fail_modes.sort(key=lambda x: x.phi, reverse=True)
        for f_mode in self.fail_modes:
            tau = f_mode.tau.m_as(ureg.min)
            if tau == float('inf'):
                # Handling infinite release
                tau = str(tau)
            if f_mode.is_const:
                event_fr = f_mode.leak_fr
            else:
                event_fr = f_mode.leak_fr/f_mode.N
            table.append([
                f_mode.source.name,
                f_mode.name,
                str(f_mode.fluid),
                event_fr.m_as(1/ureg.hr),
                f_mode.N,
                f_mode.leak_fr.m_as(1/ureg.hr),
                f_mode.q_leak.m_as(ureg.ft**3/ureg.min),
                f_mode.N_fan,
                f_mode.Q_fan.m_as(ureg.ft**3/ureg.min),
                tau,
                f_mode.O2_conc,
                f_mode.F_i,
                f_mode.scenario,
                f_mode.P_i/f_mode.leak_fr,
                f_mode.phi.m_as(1/ureg.hr)])
        wb_options = {}
        filename += '.xlsx'
        with xlsxwriter.Workbook(filename, wb_options) as workbook:
            worksheet = workbook.add_worksheet()
            for row_n, row in enumerate(table):
                for col_n, data in enumerate(row):
                    worksheet.write(row_n, col_n, data)
            header_format = workbook.add_format({'bold': True,
                                                 'font_size': 12,
                                                 'bottom': 3})
            sci_format = workbook.add_format({'num_format': '0.00E+00'},)
            # flow_format = workbook.add_format({'num_format': '#'},)
            percent_format = workbook.add_format({'num_format': '0%'},)
            number_format = workbook.add_format({'num_format': '0'},)
            worksheet.set_row(0, None, header_format)
            worksheet.set_column(3, 3, None, sci_format)
            worksheet.set_column(5, 5, None, sci_format)
            # worksheet.set_column(5, 5, None, flow_format)
            worksheet.set_column(6, 6, None, sci_format)
            worksheet.set_column(9, 9, None, sci_format)
            worksheet.set_column(10, 10, None, percent_format)
            worksheet.set_column(11, 11, None, sci_format)
            worksheet.set_column(13, 14, None, sci_format)
            # Writing total/summary
            N_rows = len(table)
            N_cols = len(table[0])
            worksheet.write(N_rows+1, N_cols-2, 'Total fatality rate, 1/hr')
            worksheet.write(N_rows+1, N_cols-1,
                            self.phi.m_as(1/ureg.hr))
            worksheet.write(N_rows+2, N_cols-2, 'ODH class')
            worksheet.write(N_rows+2, N_cols-1, self.odh_class,
                            number_format)
            worksheet.autofit()
            worksheet.conditional_format(
                1, N_cols-1, N_rows-1, N_cols-1,
                {'type': '3_color_scale', 'min_color': '#008000',
                 'max_color': '#FF0000'})
            worksheet.freeze_panes(1, 0)

    def __str__(self):
        return (f'Volume: {self.name}, {self.volume.to(ureg.ft**3):~}')


def prob_m_of_n(m, n, T, l):
    """Calculate the probability of m out of n units working.

    Calculation is done using binomial distribution.

    Parameters
    ----------
    m : int
        Number of units working.
    n : int
        Total number of units.
    T : ureg.Quantity, [time]
        Test period
    l : ureg.Quantity, 1/[time]
        Failure rate (\\lambda) of a fan

    Returns
    -------
    float
        Probability of m out of n units working.
    """
    PFD_one_unit = l*T/2
    m_of_n = binom(n, m) * (PFD_one_unit)**(n-m) * (1-PFD_one_unit)**m
    return m_of_n


def conc_vent(V, R, Q, t):
    """Calculate the oxygen concentration at the end of the event.

    As defined by FESHM 4240 6.1.A, Cases A, B, and C.

    Parameters
    ----------
    V : ureg.Quantity, [length]^3
        Volume of the confined space.
    R : ureg.Quantity, [length]^3/[time]
        Volumetric spill rate into confined space.
    Q : ureg.Quantity, [length]^3/[time]
        Volumetric ventilation rate of fan(s); positive value corresponds
        to blowing air into the confined space, negative - drawing contaminated
        air outside.
    t : ureg.Quantity, [time]
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
    return float(C)


def conc_final(V, R, Q):
    """Calculate the final oxygen concentration for continuous flow.

    Equivalent to conc_vent(V, R, Q, float('inf')).

    Parameters
    ----------
    V : ureg.Quantity, [length]^3
        Volume of the confined space.
    R : ureg.Quantity, [length]^3/[time]
        Volumetric spill rate into confined space.
    Q : ureg.Quantity, [length]^3/[time]
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
    V : ureg.Quantity, [length]^3
        Volume of the confined space.
    C_e : float
        Oxygen concentration at the end of the release.
    Q : ureg.Quantity, [length]^3/[time]
        Volumetric ventilation rate of fan(s); positive value corresponds
        to blowing air into the confined space, negative - drawing contaminated
        air outside.
    t : ureg.Quantity, [time]
        time, beginning of release is at `t`=0.
    t_e : ureg.Quantity, [time]
        time when release ended.

    Returns
    -------
    float
        Oxygen concentration.
    """
    C = 0.21-(0.21-C_e)*math.e**-(abs(Q)/V*(t-t_e))
    return C


# def print_result(*Volumes):
#     """Print the results of the ODH analysis for a volume.

#     If several volumes given (in case of interlapping volumes) the worst case
#     will be printed.
#     """
#     max_phi = -1/ureg.hr
#     for volume in Volumes:
#         if volume.phi > max_phi:
#             max_volume = volume
#     line_1 = '# Fatality rate for {} is {:.1e}  # '.format(max_volume, volume.phi)
#     pad = len(line_1)
#     line_2 = '# Recommended ODH class {}'.format(max_volume.odh_class()).ljust(pad-1)+'#'
#     print('#'*pad)
#     print(line_1)
#     print(line_2)
#     print('#'*pad)


def hole_leak(area, fluid, P_out=P_NTP, Kd=0.62):
    """Calculate leak flow rate through a hole in a pipe.

    Discharge flow through a hole is calculated using Darcy equation for
    compressible fluid using expansion factor Y Crane TP 410 2013 eq. 4-15.
    Resistance coefficient K of the orifice taken from Rennels, Pipe Flow,
    2012, eq. 13.5.
    For this calculation the gas is assumed to have no pressure loss on
    entry. This makes analysis simple and conservative.

    Parameters
    ----------
    area : ureg.Quantity, [length]^2
        Area of the leak.
    fluid : heat_transfer.ThermState
        Thermodynamic state of the fluid stored in the source.
    P_out : ureg.Quantity, [mass]/([length]*[time]^2)
        Pressure at the exit of the hole.
    Kd : float, optional
        Discharge coefficient.
        For uneven holes 0.62 provides a good approximation.

    Returns
    -------
    ureg.Quantity, [length]^3/[time]
        Standard volumetric flow at Normal Temperature and Pressure.

    """
    m_dot = G_nozzle(fluid) * area * Kd
    q_std = to_standard_flow(m_dot, fluid)
    return q_std


def O2_sudden_release(release, volume, escape=True):
    """Calculate oxygen concentration after a sudden release.

    Parameters
    ----------
    release : Quantity [length*3]
        Standard volume of the inert gas released.
    volume : Quantity [length^3]
        Volume of the room, building, or area analyzed.
    escape : bool, optional
        If True, mixed air is allowed to escape from considered volume.
        If False, inert gas is trapped and expels the air outside the considered volume.
        Default is True.

    Returns
    -------
    float
        Resulting oxygen concentration.
    """
    if escape == True:
        O2_conc = 0.21*volume/(volume+release)
    else:
        O2_conc = 0.21*(1-release/volume)
    return O2_conc


def PFD_avg(l_du, l_dd, T_p, MTTR):
    """
    Calculate probability of failure on demand based on dangerous detected and dangerous undetected failure rates.

    The equation is taken from MSA Ultima X5000 safety manual.

    Parameters
    ----------
    l_du : Quantity [1/time]
        Dangerous undetected failure rate.
    l_dd : Quantity [1/time]
        Dangerous detected failure rate.
    T_p : Quantity [time]
        Proof test interval.
    MTTR : Quantity [time]
        Mean time to restoration.

    Returns
    -------
    float
        Probability of failure on demand.
    """
    PFD_du = l_du * (T_p/2 + MTTR)
    PFD_dd = l_dd * MTTR
    PFD = PFD_du + PFD_dd
    return float(PFD)
