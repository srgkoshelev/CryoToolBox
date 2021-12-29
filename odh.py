"""odh - Oxygen Deficiency Hazard analysis tool.

Based on Fermilab ES&H Manual chapter 4240 (see
https://eshq.fnal.gov/manuals/feshm/)

This package is maintained at https://github.com/srgkoshelev/ODH_analysis"""

import math
from scipy.special import binom
from . import ureg, Q_
from . import logger
from . import piping
from . import T_NTP, P_NTP
from .functions import ThermState
from .functions import to_standard_flow
from collections import namedtuple
import xlsxwriter
# Loading FESHM 4240 Failure rates
from .FESHM4240_TABLES import TABLE_1, TABLE_2


# Probability of failure on demand for main cases
PFD_ODH = Q_('2 * 10^-3')
# TODO Update to value from J. Anderson's document
TRANSFER_LINE_LEAK_AREA = Q_('10 mm^2')  # Taken for consistency
SHOW_SENS = 5e-8/ureg.hr
# Min required air intake from Table 6.1, ASHRAE 62-2001
ASHRAE_MIN_FLOW = 0.06 * ureg.ft**3/(ureg.min*ureg.ft**2)

failure_mode = namedtuple('Failure_mode', ['phi', 'source', 'name',
                                           'O2_conc', 'leak_fr', 'P_i',
                                           'F_i', 'outage', 'q_leak', 'tau',
                                           'Q_fan', 'N_fan', 'N'])


class ODHError(Exception):
    pass


class Source:
    """Source of inert gas

    Attributes
    ----------
    sol_PFD : float
        Probability of failure on demand (PFD) for solenoid valve.
        If the source doesn't have isolating solenoid valve
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
            solenoid valve.
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
        self.sol_PFD = (int(not isol_valve) or
                        TABLE_2['Valve, solenoid']['Failure to operate'])

    def pipe_failure(self, max_flow, tube, fluid=None, N_welds=1):
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
        max_flow : ureg.Quantity {mass: 1, time: -1} or {length: 3, time: -1}
            Max mass or volumetric flow through if limited,
            e.g. by compressor output.
        tube : heat_transfer.piping.Tube
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid for the release.
        N_welds : int
            Number of welds on the tube.
        """
        # If fluid not defined use fluid of the Source
        fluid = fluid or self.fluid
        # Failure rate coefficients; Piping failure rate is per unit of length,
        # weld is dependent on number of welds, pipe OD and wall thickness
        failure_rate_coeff = {'Piping': (tube.L, 1),
                              'Pipe weld': (tube.OD / tube.wall,
                                            N_welds)}
        # Piping and weld leaks as per Table 2
        for cause in ['Piping', 'Pipe weld']:
            for mode in TABLE_2[cause].keys():
                if tube.D > 2 or mode != 'Large leak':  # Large leak only for D > 2"
                    name = f'{cause} {mode.lower()}: {tube}, ' + \
                        f'{tube.L.to(ureg.ft):.3g~}'
                    fr_coef = failure_rate_coeff[cause][0]
                    N_events = failure_rate_coeff[cause][1]
                    if mode == 'Rupture':
                        failure_rate = fr_coef * TABLE_2[cause][mode]
                        # For rupture calculate flow through available
                        # pipe area
                        area = tube.area
                        q_std = to_standard_flow(max_flow, fluid)
                    else:
                        failure_rate = fr_coef * \
                            TABLE_2[cause][mode]['Failure rate']
                        area = TABLE_2[cause][mode]['Area']
                        if area > tube.area:
                            ODHError('Leak area cannot be larger'
                                     ' than pipe area.')
                        q_std = hole_leak(tube, area, fluid)
                    self.failure_mode(name, failure_rate, q_std, N_events)

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
        self.failure_mode('Dewar insulation failure', failure_rate, q_std, 1)

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

    def flange_failure(self, pipe, fluid=None, N=1):
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
        N : int
            Number of reinforced seal connections on the pipe.
        """
        # TODO Make leak and rupture areas adjustable, add info to docstring
        table = TABLE_2['Flange, reinforced gasket']
        area_cases = {
            'Leak': table['Leak']['Area'],
            'Rupture': pipe.area}
        for mode in table:
            name = f'Flange {mode.lower()}: {pipe}'
            if isinstance(table[mode], dict):
                failure_rate = table[mode]['Failure rate']
            else:
                failure_rate = table[mode]
            area = area_cases[mode]
            # TODO move this and gas leak check to separate method
            if area > pipe.area:
                logger.warning('Leak area cannot be larger'
                               ' than pipe area.')
                continue
            # If fluid not defined use fluid of the Source
            fluid = fluid or self.fluid
            q_std = hole_leak(pipe, area, fluid)
            self.failure_mode(name, failure_rate, q_std, N)

    def transfer_line_failure(self, pipe, fluid=None, N=1):
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
        N : int
            Number of bayonets/soft seals on the transfer line.
        """
        # TODO This should be replaced by flange failure at some point
        area_cases = {'Leak': TRANSFER_LINE_LEAK_AREA,
                      'Rupture': pipe.area}
        for mode in TABLE_1['Fluid line']:
            name = f'Fluid line {mode.lower()}: {pipe}'
            failure_rate = TABLE_1['Fluid line'][mode]
            area = area_cases[mode]
            if area > pipe.area:
                logger.warning('Leak area cannot be larger'
                               ' than pipe area.')
                continue
            # If fluid not defined use fluid of the Source
            fluid = fluid or self.fluid
            q_std = hole_leak(pipe, area, fluid)
            self.failure_mode(name, failure_rate, q_std, N)

    def pressure_vessel_failure(self, q_std_rupture, relief_area, fluid=None):
        """Add pressure vessel failure to leaks dict.

        Store failure rate, flow rate and expected time duration of
        the event for transfer line failure. Based on FESHM 4240.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        q_std_rupture : ureg.Quantity {length: 3, time: -1}
            Standard volumetric flow rate for pressure vessel rupture.
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
        q_std = hole_leak(None, area, fluid)
        self.failure_mode(name, failure_rate, q_std, 1)

        # Rupture case
        name = 'Pressure vessel rupture'
        if relief_area is None:
            # No vacuum jacket on the vessel
            area = 1000 * ureg.mm**2  # Equal to large leak of a pipe
        else:
            area = relief_area
        failure_rate = TABLE_2['Vessel, pressure']['Failure']
        q_std = hole_leak(None, area, fluid)
        self.failure_mode(name, failure_rate, q_std, 1)

    # def constant_leak(self, name, q_std, N=1):
    #     """Add constant leak to leaks dict.

    #     Store flow rate and expected time duration of the
    #     failure event for general failure mode.
    #     Failure modes are analyzed by `Volume.odh` method.

    #     Parameters
    #     ----------
    #     name : str
    #         Name of the failure mode
    #     q_std : ureg.Quantity {length: 3, time: -1}
    #         Standard volumetric flow rate.
    #     N : int
    #         Quantity of similar failure modes.
    #     """
    #     # Failure rate assumes the volume instantly refilled
    #     # after being completely emptied, and continues release
    #     # Failure rate for constant leak doesn't depend on N or self.N
    #     # Dividing by self.N*N to undo _make_leak multiplication
    #     failure_rate = q_std/(self.volume*self.N*N)
    #     self.failure_mode(name, failure_rate, N*q_std, N)

    def failure_mode(self, name, failure_rate, q_std, N=1):
        """Add general failure mode to leaks dict.

        Store failure rate, flow rate and expected time duration of the
        failure event for general failure mode.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        name : str
            Name of the failure mode
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
        self.leaks.append((name, total_failure_rate, q_std, tau.to(ureg.min), N_events))

    @staticmethod
    def combine(name, sources):
        """Combine several ODH sources sharing volume.

        Can be used for failure modes affecting several sources in parallel.

        Parameters
        ----------
        name : str
            Name of the new combined source.
        sources : list of Source
            Sources connected together.

        Returns
        -------
        Source
            Combined source of inert gas.
        """
        fluid = ThermState(sources[0].fluid.name, T=T_NTP, P=P_NTP)
        if all([source.fluid.name == fluid.name for source in sources]):
            total_volume = sum([source.volume for source in sources])
            return Source(name, fluid, total_volume)
        else:
            ODHError('All volumes should contain the same fluid')
            return None

    def __str__(self):
        return f'{self.name}, ' + \
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


class Volume:
    """Volume/building affected by inert gases."""
    def __init__(self, name, volume, *, Q_fan, N_fans, T_fan,
                 lambda_fan=TABLE_2['Fan']['Failure to run'],
                 vent_rate=0*ureg.ft**3/ureg.min):
        """Define a volume affected by inert gas release from  a `Source`.

        Parameters
        ----------
        name : str
            Name of the volume.
        volume : ureg.Quantity {length: 3}
            Volume of the building or part of the building.
        Q_fan : ureg.Quantity {length: 3, time: -1}
            Volumetric flow of a single ODH fan installed in the volume.
        N_fans : int
            Number of fans installed.
        T_fan : ureg.Quantity {time: 1}
            Test period of the fans.
        vent_rate : ureg.Quantity {length: 3, time: -1}
            Min volumetric flow required or present in the building.
        lambda_fan : ureg.Quantity {time: -1}
            Failure rate of the fans in the building.
        """
        self.name = name
        self.volume = volume
        self.vent_rate = vent_rate
        self.PFD_ODH = PFD_ODH  # Default value for ODH system failure
        # TODO Should be explicit/with default option
        self.lambda_fan = lambda_fan
        self.Q_fan = Q_fan
        self.N_fans = N_fans
        self.Test_period = T_fan
        # Calculate fan probability of failure
        self._fan_fail()
        # TODO should be external function; Don't need to keep fan info?

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
        """
        self.fail_modes = []
        # Probability of power failure in the building:
        # PFD_power if no outage, 1 if there is outage
        PFD_power_build = (power_outage or
                           TABLE_1['Electrical Power Failure']['Demand rate'])
        # Calculate fatality rates for each source
        for source in sources:
            for leak in source.leaks:
                leak_failure_rate = leak[0]
                if leak_failure_rate is not None:  # None for constant leak
                    self._fatality_no_response(source, leak, source.sol_PFD,
                                               PFD_power_build)
                    self._fatality_fan_powered(source, leak, source.sol_PFD,
                                               PFD_power_build)

    def _fatality_no_response(self, source, leak, sol_PFD,
                              PFD_power_build):
        """Calculate fatality rate in the volume for ODH protection failure.

        Calculate failure rate of leak occuring and no safety response
        occuring due to power failure and isolation solenoid failure,
        or power on and ODH system failure.
        O2 concentration is limited only by amount of inert gas the source has.
        Fans are not operational.
        Adds calculation results to the fail_modes list.

        Parameters
        ----------
        source : Source
        leak : tuple (str,
                      ureg.Quantity {time: -1},
                      ureg.Quantity {length: 3, time: -1},
                      ureg.Quantity {time: 1},
                      int)
            Leak failure rate, volumetric flow rate, event duration, and number
            of events.
        sol_PFD : float
            Probability of source solenoid failure.
        PFD_power_building : float
            Probability of power failure.
        """
        (failure_mode_name, leak_failure_rate, q_leak, tau, N) = leak
        P_no_response = float(PFD_power_build) * sol_PFD + \
            (1-PFD_power_build)*self.PFD_ODH
        P_i = leak_failure_rate * P_no_response
        Q_fan = self.vent_rate
        O2_conc = conc_vent(self.volume, q_leak, Q_fan, tau)
        F_i = self._fatality_prob(O2_conc)
        phi_i = P_i*F_i
        f_mode = failure_mode(phi_i, source, failure_mode_name, O2_conc,
                              leak_failure_rate, P_i, F_i,
                              PFD_power_build == 1, q_leak, tau, Q_fan, 0, N)
        self.fail_modes.append(f_mode)

    def _fatality_fan_powered(self, source, leak, sol_PFD, PFD_power_build):
        """Calculate fatality rates for fan failure on demand.

        Calculate fatality rates for the case of ODH system responding and
        fans powered but some of the fans failing on demand.
        See wiki for further explanation.
        Adds calculation results to the fail_modes list.

        Parameters
        ----------
        source : Source
        leak : tuple (str,
                      ureg.Quantity {time: -1},
                      ureg.Quantity {length: 3, time: -1},
                      ureg.Quantity {time: 1},
                      int)
            Leak failure rate, volumetric flow rate, event duration, and number
            of events.
        sol_PFD : float
            Probability of source solenoid failure.
        PFD_power_building : float
            Probability of power failure.
        """
        (failure_mode_name, leak_failure_rate, q_leak, tau, N) = leak
        for (P_fan, Q_fan, N_fan) in self.Fan_flowrates:
            # Probability of power on, ODH system working, and m number of fans
            # with flow rate Q_fan on.
            P_response = (1-PFD_power_build) * (1-self.PFD_ODH) * \
                sol_PFD * P_fan
            P_i = leak_failure_rate * P_response
            O2_conc = conc_vent(self.volume, q_leak, Q_fan, tau)
            F_i = self._fatality_prob(O2_conc)
            phi_i = P_i*F_i
            f_mode = failure_mode(phi_i, source, failure_mode_name, O2_conc,
                                  leak_failure_rate, P_i, F_i,
                                  PFD_power_build == 1, q_leak, tau, Q_fan,
                                  N_fan, N)
            self.fail_modes.append(f_mode)

    def _fan_fail(self):
        """Calculate (Probability, flow) pairs for all combinations of fans
        working.

        All fans are expected to have same volume flow.
        """
        # TODO add fans with different volumetric rates (see report as well)
        Fail_rate = self.lambda_fan
        Fan_flowrates = []
        for m in range(self.N_fans+1):
            # Probability of exactly m units starting
            P_m_fan_work = prob_m_of_n(m, self.N_fans, self.Test_period,
                                       Fail_rate)
            flowrate = self.Q_fan*m
            if flowrate == Q_('0 m**3/min'):
                flowrate = self.vent_rate
            Fan_flowrates.append((P_m_fan_work, flowrate, m))
        self.Fan_flowrates = Fan_flowrates

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
        sens = sens or SHOW_SENS
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
        sens = sens or SHOW_SENS
        table = [["Failure mode", "Fans on", "O_2", "Duration, min", "\\phi_i"]]
        table.append(None)
        for f_mode in self.fail_modes:
            if f_mode.phi >= sens:
                row = []
                row.append(f'{f_mode.source.name} {f_mode.name}')
                row.append(f'{f_mode.N_fan}')
                row.append(f'{f_mode.O2_conc:.0%}')
                row.append(f'{f_mode.tau.m_as(ureg.min):,.1f}')
                row.append(f'{f_mode.phi.m_as(1/ureg.hr):.2}')
                table.append(row)
        return table

    def report_table(self, filename='ODH_report'):
        """Make a table with the calculation results."""
        table = []
        header = ['Source', 'Failure', 'Event failure rate, 1/hr', '# of',
                  'Total failure rate, 1/hr', 'Leak rate, SCFM',
                  '# fans working', 'Fan rate, SCFM', 'Event duration, min',
                  'Oxygen concentration', 'Fatality prob', 'Case prob',
                  'Fatality rate, 1/hr']
        # 'Total failure rate', 'ODH protection PFD', 'Building is powered'
        table.append(header)
        self.fail_modes.sort(key=lambda x: x.source.name)
        for f_mode in self.fail_modes:
            table.append([
                f_mode.source.name,
                f_mode.name,
                (f_mode.leak_fr/f_mode.N).m_as(1/ureg.hr),
                f_mode.N,
                f_mode.leak_fr.m_as(1/ureg.hr),
                f_mode.q_leak.m_as(ureg.ft**3/ureg.min),
                f_mode.N_fan,
                f_mode.Q_fan.m_as(ureg.ft**3/ureg.min),
                f_mode.tau.m_as(ureg.min),
                f_mode.O2_conc,
                f_mode.F_i,
                f_mode.P_i/f_mode.leak_fr,
                f_mode.phi.m_as(1/ureg.hr)])
        filename += '.xlsx'
        with xlsxwriter.Workbook(filename) as workbook:
            header_format = workbook.add_format({'bold': True,
                                                 'font_size': 12,
                                                 'bottom': 3})
            worksheet = workbook.add_worksheet()
            col_width = [len(x) for x in table[0]]
            for row_n, row in enumerate(table):
                for col_n, data in enumerate(row):
                    worksheet.write(row_n, col_n, data)
                    if col_n in (0, 1, 10):
                        # For source names, failure names
                        # and 'Total failure rate'
                        col_width[col_n] = max(col_width[col_n], len(str(data)))
            sci_format = workbook.add_format({'num_format': '0.00E+00'},)
            flow_format = workbook.add_format({'num_format': '#'},)
            percent_format = workbook.add_format({'num_format': '0%'},)
            number_format = workbook.add_format({'num_format': '0'},)
            worksheet.set_row(0, None, header_format)
            worksheet.set_column(2, 2, None, sci_format)
            worksheet.set_column(4, 4, None, sci_format)
            worksheet.set_column(5, 5, None, flow_format)
            worksheet.set_column(8, 8, None, sci_format)
            worksheet.set_column(9, 9, None, percent_format)
            worksheet.set_column(10, 12, None, sci_format)
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

    def __str__(self):
        return (f'Volume: {self.name}, {self.volume.to(ureg.ft**3):.2~}')

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
    return m_of_n


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


def conc_final(V, R, Q):
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


def hole_leak(tube, area, fluid, P_out=P_NTP):
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
    P1 = fluid.P
    P2 = P_out
    if P2 >= P1:
        return 0*ureg.ft**3/ureg.min

    k = float(fluid.gamma)
    rc = (2/(k+1))**(k/(k-1))  # Critical pressure ratio
    if P2/P1 < rc:
        P2 = P1 * rc  # Choked flow

    dP = P1 - P2
    rho = fluid.Dmass
    A_hole = area
    beta = float((area/tube.area)**0.5)
    if fluid.phase in (0, 3, 6):
        # Fluid is a liquid
        Y = 1
    else:
        Y = 1 - (0.351 + 0.256*beta**4 + 0.93*beta**8) * (1-(P2/P1)**(1/k))
    if tube is None:
        K_orifice = 2.81  # Rennels, Pipe Flow, 13.2.3
    else:
        lam = 1 + 0.622*(1 - 0.215*beta**2 - 0.785*beta**5)
        K_orifice = lam**2 * (1 + 0.0696*(1-beta**5))

    q_local = Y * A_hole * ((2*dP)/(rho*K_orifice))**0.5
    q_std = to_standard_flow(q_local, fluid)
    return q_std

    #     """Calculate leak flow/release for a given piping element.

    #     For the full pipe rupture case it is usually assumed that the release
    #     area is equal to tube cross section area. For other leak cases, the
    #     hole in the piping is considered to be a square-edged orifice.

    #     Parameters
    #     ----------
    #     """
    #     d = (4*area/math.pi)**0.5  # diameter for the leak opening
    #     exit_ = piping.Exit(d)
    #     TempPiping = piping.Piping(fluid)
    #     TempPiping.add(
    #                    tube,
    #                    exit_,
    #     )
    #     if area != tube.area:
    #         Hole = piping.Orifice(d)
    #         TempPiping.insert(1, Hole)
    #     m_dot = TempPiping.m_dot(P_NTP)
    #     fluid_NTP = fluid.copy()
    #     fluid_NTP.update_kw(P=P_NTP, T=T_NTP)
    #     q_std = m_dot / fluid_NTP.Dmass
    #     return q_std
