"""Unit-aware wrapper layer around CoolProp and the optional HEPROP backend.

This module translates between Pint quantities and the lower-level property
backends so the rest of CryoToolBox can work with dimensional values while
keeping the backend-specific details in one place.
"""
import CoolProp.CoolProp as CP
from .heprop_state import HepropState
from math import inf
from .std_conditions import ureg, T_NTP, P_NTP, P_MSC, T_MSC, P_STD, T_STD
import warnings


class CryoToolBoxFutureWarning(FutureWarning):
    """FutureWarning subclass for user-facing API behavior changes."""


class BranchSelectionNotImplementedError(NotImplementedError):
    """Raised when a branch-selecting CoolProp update is requested."""


_HAS_QMASS = hasattr(CP, 'iQmass')


CP_const_unit = {
    'gas_constant': (CP.igas_constant, ureg.J / ureg.mol / ureg.K),
    'molar_mass': (CP.imolar_mass, ureg.kg / ureg.mol),
    # 'acentric_factor': (CP.iacentric_factor, ),
    # 'Dmolar_reducing': (CP.irhomolar_reducing, ),
    'Dmolar_critical': (CP.irhomolar_critical, ureg.mol / ureg.m**3),
    'T_reducing': (CP.iT_reducing, ureg.K),
    'T_critical': (CP.iT_critical, ureg.K),
    # 'Dmass_reducing': (CP.irhomass_reducing, ),
    'Dmass_critical': (CP.irhomass_critical, ureg.kg / ureg.m**3),
    'P_critical': (CP.iP_critical, ureg.Pa),
    'P_reducing': (CP.iP_reducing, ureg.Pa),
    'T_triple': (CP.iT_triple, ureg.K),
    'P_triple': (CP.iP_triple, ureg.Pa),
    'T_sat': (None, ureg.K),
    'P_sat': (None, ureg.Pa),
    'T_min': (CP.iT_min, ureg.K),
    'T_max': (CP.iT_max, ureg.K),
    'P_max': (CP.iP_max, ureg.Pa),
    'P_min': (CP.iP_min, ureg.Pa),
    # 'dipole_moment': (CP.idipole_moment, ),
    'T': (CP.iT, ureg.K),
    'P': (CP.iP, ureg.Pa),
    'Q': (CP.iQ, ureg.mol / ureg.mol),
    # 'Tau': (CP.iTau, ),
    # 'Delta': (CP.iDelta, ),
    'Dmolar': (CP.iDmolar, ureg.mol / ureg.m**3),
    'Hmolar': (CP.iHmolar, ureg.J / ureg.mol),
    'Smolar': (CP.iSmolar, ureg.J / ureg.mol / ureg.K),
    'Cpmolar': (CP.iCpmolar, ureg.J / ureg.mol / ureg.K),
    # 'Cp0molar': (CP.iCp0molar, ureg.J/ureg.mol/ureg.K),
    'Cvmolar': (CP.iCvmolar, ureg.J / ureg.mol / ureg.K),
    'Umolar': (CP.iUmolar, ureg.J / ureg.mol),
    # 'Gmolar': (CP.iGmolar, ),
    # 'Helmholtzmolar': (CP.iHelmholtzmolar, ),
    # 'Smolar_residual': (CP.iSmolar_residual, ),
    'Dmass': (CP.iDmass, ureg.kg / ureg.m**3),
    'Hmass': (CP.iHmass, ureg.J / ureg.kg),
    'Smass': (CP.iSmass, ureg.J / ureg.kg / ureg.K),
    'Cpmass': (CP.iCpmass, ureg.J / ureg.kg / ureg.K),
    # 'Cp0mass': (CP.iCp0mass, ureg.J/ureg.kg/ureg.K),
    'Cvmass': (CP.iCvmass, ureg.J / ureg.kg / ureg.K),
    'Umass': (CP.iUmass, ureg.J / ureg.kg),
    # 'Gmass': (CP.iGmass, ),
    # 'Helmholtzmass': (CP.iHelmholtzmass, ),
    'viscosity': (CP.iviscosity, ureg.Pa * ureg.s),
    'conductivity': (CP.iconductivity, ureg.W / ureg.m / ureg.K),
    'surface_tension': (CP.isurface_tension, ureg.N / ureg.m),
    'Prandtl': (CP.iPrandtl, None),
    'speed_sound': (CP.ispeed_sound, ureg.m / ureg.s),
    'isothermal_compressibility':
    (CP.iisothermal_compressibility, ureg.Pa**-1),
    'isobaric_expansion_coefficient':
    (CP.iisobaric_expansion_coefficient, ureg.K**-1),
    'isentropic_expansion_coefficient':
    (CP.iisentropic_expansion_coefficient, ureg.dimensionless),
    # 'fundamental_derivative_of_gas_dynamics':
    # (CP.ifundamental_derivative_of_gas_dynamics, ),
    # 'alphar': (CP.ialphar, ),
    # 'alpha0': (CP.ialpha0, ),
    # 'Bvirial': (CP.iBvirial, ),
    # 'Cvirial': (CP.iCvirial, ),
    # 'dBvirial_dT': (CP.idBvirial_dT, ),
    # 'dCvirial_dT': (CP.idCvirial_dT, ),
    'Z': (CP.iZ, ureg.dimensionless),
    # 'PIP': (CP.iPIP, ),
    # 'fraction_min': (CP.ifraction_min, ),
    # 'fraction_max': (CP.ifraction_max, ),
    # 'T_freeze': (CP.iT_freeze, ),
    # 'GWP20': (CP.iGWP20, ),
    # 'GWP100': (CP.iGWP100, ),
    # 'GWP500': (CP.iGWP500, ),
    # 'FH': (CP.iFH, ),
    # 'HH': (CP.iHH, ),
    # 'PH': (CP.iPH, ),
    # 'ODP': (CP.iODP, ),
    'Phase': (CP.iPhase, ureg.dimensionless),
    'C_us': (
        None,  # unit defined for external use
        ureg.lb / (ureg.hr * ureg.lbf) * (ureg.degR)**0.5),
    'C_m': (
        None,  # unit defined for external use
        ureg.s**2 / (ureg.hr * ureg.m) * (ureg.K)**0.5)
}

CP_const_unit['Qmass'] = (
    getattr(CP, 'iQmass', None),
    ureg.kg / ureg.kg,
)


def _input_parameter(name):
    """Return the CoolProp parameter ID for an update input."""
    try:
        parameter = CP_const_unit[name][0]
    except KeyError as exc:
        raise ValueError(f'Unsupported CoolProp input: {name!r}.') from exc
    if parameter is None:
        raise ValueError(
            f'{name!r} is not available as an input with this CoolProp version.'
        )
    return parameter


class ThermState:
    """Unit-aware wrapper around CoolProp's ``AbstractState``.

    Parameters
    ----------
    fluid : str
        CoolProp fluid string such as ``"nitrogen"`` or ``"helium"``.
    backend : str, default ``"HEOS"``
        Property backend. ``HEOS`` is the supported default. ``REFPROP`` and
        ``HEPROP`` require local third-party installations.
    **state_parameters
        Optional pair of state variables used to initialize the state, for
        example ``P=...`` and ``T=...``.

    Notes
    -----
    Public inputs and outputs are unit-aware quantities. Values are converted
    internally to the units expected by CoolProp or the optional HEPROP
    adapter.
    """

    def __init__(self, fluid, backend="HEOS", **state_parameters):
        """
        Available backends: HEOS, REFPROP, and HEPROP.

        HEOS is the default and supported alpha path.
        REFPROP and HEPROP are optional backends that require local
        third-party installations.
        """
        # TODO this should check for the fluid name too
        if backend == "HEPROP":
            if fluid != 'helium':
                raise ValueError(
                    f'Only helium fluid can be used with HEPROP backend: {fluid}'
                )
            try:
                self._AbstractState = HepropState()
            except OSError as exc:
                raise OSError(
                    'HEPROP backend is optional and requires a working '
                    '`HEPROPPATH` installation. Use the default HEOS backend '
                    'for the supported alpha path.'
                ) from exc
        else:
            self._AbstractState = CP.AbstractState(backend, fluid)
        if state_parameters:
            self.update_kw(**state_parameters)

    def update_kw(self, **state_parameters):
        """Update thermodynamic state using keyword arguments."""
        if len(state_parameters) != 2:
            raise TypeError(
                f'update_kw() takes 2 arguments '
                f'({len(state_parameters)} given)')
        (name1, value1), (name2, value2) = state_parameters.items()
        self.update(name1, value1, name2, value2)

    def update(self, name1, value1, name2, value2):
        """Update the state from two named thermodynamic inputs."""
        self.update_si(
            name1, self.prepare_input(name1, value1),
            name2, self.prepare_input(name2, value2))

    def update_si(self, name1, value1, name2, value2):
        """Update the state from raw CoolProp SI input values.

        CoolProp's low-level ``AbstractState`` API uses SI units for all
        Python backends, including the REFPROP wrapper. Use ``update`` or
        ``update_kw`` for Pint-aware inputs; use this method only when values
        are already in CoolProp's expected units such as Pa, K, kg/m^3,
        mol/m^3, J/kg, J/mol, J/(kg*K), or J/(mol*K).
        """
        input_pair, value1, value2 = CP.generate_update_pair(
            _input_parameter(name1), value1,
            _input_parameter(name2), value2)
        self._AbstractState.update(input_pair, value1, value2)

    def update_with_guesses(self, *args, **kwargs):
        """Placeholder for future branch-selecting saturation updates.

        CoolProp 8.0.0 can disambiguate multi-root saturation inputs through
        ``AbstractState.update_with_guesses`` and a temperature guess. CTB has
        not wrapped that API yet because it needs a unit-aware, backend-aware
        interface and clear tests around branch selection.
        """
        raise BranchSelectionNotImplementedError(
            'ThermState.update_with_guesses() is planned but not implemented. '
            'Use update()/update_kw() for unambiguous states; branch-selecting '
            'saturation flashes will be added in a future CTB release.'
        )

    def get_si(self, name):
        """Return a current-state property in raw CoolProp SI units.

        This is the output counterpart to :meth:`update_si`. Dimensionless
        values such as ``Q`` and ``Z`` are returned as plain floats.
        """
        if name == 'Qmass' and not _HAS_QMASS:
            return self.get_si('Q')
        try:
            parameter, unit = CP_const_unit[name]
        except KeyError as exc:
            raise ValueError(f'Unsupported CoolProp output: {name!r}.') from exc

        if parameter is not None:
            return self._AbstractState.keyed_output(parameter)

        value = getattr(self, name)
        if hasattr(value, 'to'):
            return value.to(unit).magnitude
        return float(value)

    @staticmethod
    def prepare_input(name, value):
        """Prepare value to input to CoolProp.

        Converts the value to units expected by CoolProp
        (see `CP_const_unit`).

        Parameters
        ----------
        name : str
            Parameter name.
        value : ureg.Quantity
            Parameter value.

        Returns
        -------
        float
        """
        if isinstance(value, (int, float)):
            return value
        unit = CP_const_unit[name][1]
        return value.to(unit).magnitude

    @property
    @ureg.wraps(CP_const_unit['T'][1], None)
    def T_critical(self):
        """Critical temperature."""
        return self._AbstractState.T_critical()

    @property
    @ureg.wraps(CP_const_unit['P'][1], None)
    def P_critical(self):
        """Critical pressure."""
        return self._AbstractState.p_critical()

    @property
    @ureg.wraps(CP_const_unit['Dmolar_critical'][1], None)
    def Dmolar_critical(self):
        """Critical molar density."""
        return self._AbstractState.rhomolar_critical()

    @property
    @ureg.wraps(CP_const_unit['Dmass_critical'][1], None)
    def Dmass_critical(self):
        """Critical mass density."""
        return self._AbstractState.rhomass_critical()

    @property
    @ureg.wraps(CP_const_unit['T'][1], None)
    def T(self):
        """Temperature of the current state."""
        return self._AbstractState.T()

    @property
    @ureg.wraps(CP_const_unit['Dmolar_critical'][1], None)
    def Dmolar(self):
        """Molar density of the current state."""
        return self._AbstractState.rhomolar()

    @property
    @ureg.wraps(CP_const_unit['Dmass_critical'][1], None)
    def Dmass(self):
        """Mass density of the current state."""
        return self._AbstractState.rhomass()

    @property
    @ureg.wraps(CP_const_unit['P'][1], None)
    def P(self):
        """Pressure of the current state."""
        return self._AbstractState.p()

    @property
    def Q(self):
        """Molar vapor quality for two-phase states."""
        return self._AbstractState.Q()

    @property
    def Qmass(self):
        """Mass-basis vapor quality for two-phase states."""
        if hasattr(self._AbstractState, 'Qmass'):
            return self._AbstractState.Qmass()
        return self.Q

    @property
    @ureg.wraps(CP_const_unit['molar_mass'][1], None)
    def molar_mass(self):
        """Fluid molar mass."""
        return self._AbstractState.molar_mass()

    @property
    @ureg.wraps(CP_const_unit['gas_constant'][1], None)
    def gas_constant(self):
        """Universal gas constant returned by the backend."""
        return self._AbstractState.gas_constant()

    @property
    def specific_gas_constant(self):
        """Specific gas constant derived from ``R/M``."""
        R_spec = self.gas_constant / self.molar_mass
        return R_spec

    @property
    # @ureg.wraps(CP_const_unit['Z'][1], None)
    def compressibility_factor(self):
        """Compressibility factor ``Z``."""
        return self._AbstractState.compressibility_factor()

    # Useful shorthand
    Z = compressibility_factor

    @property
    @ureg.wraps(CP_const_unit['Hmolar'][1], None)
    def Hmolar(self):
        """Molar enthalpy."""
        return self._AbstractState.hmolar()

    @property
    @ureg.wraps(CP_const_unit['Hmass'][1], None)
    def Hmass(self):
        """Specific enthalpy on a mass basis."""
        return self._AbstractState.hmass()

    @property
    @ureg.wraps(CP_const_unit['Smolar'][1], None)
    def Smolar(self):
        """Molar entropy."""
        return self._AbstractState.smolar()

    @property
    @ureg.wraps(CP_const_unit['Smass'][1], None)
    def Smass(self):
        """Specific entropy on a mass basis."""
        return self._AbstractState.smass()

    @property
    @ureg.wraps(CP_const_unit['Cpmolar'][1], None)
    def Cpmolar(self):
        """Constant-pressure heat capacity on a molar basis."""
        return self._AbstractState.cpmolar()

    @property
    @ureg.wraps(CP_const_unit['Cpmass'][1], None)
    def Cpmass(self):
        """Constant-pressure heat capacity on a mass basis."""
        return self._AbstractState.cpmass()

    @property
    @ureg.wraps(CP_const_unit['Cvmolar'][1], None)
    def Cvmolar(self):
        """Constant-volume heat capacity on a molar basis."""
        return self._AbstractState.cvmolar()

    @property
    @ureg.wraps(CP_const_unit['Cvmass'][1], None)
    def Cvmass(self):
        """Constant-volume heat capacity on a mass basis."""
        return self._AbstractState.cvmass()

    @property
    @ureg.wraps(CP_const_unit['viscosity'][1], None)
    def viscosity(self):
        """Dynamic viscosity."""
        return self._AbstractState.viscosity()

    @property
    @ureg.wraps(CP_const_unit['conductivity'][1], None)
    def conductivity(self):
        """Thermal conductivity."""
        return self._AbstractState.conductivity()

    @property
    @ureg.wraps(CP_const_unit['surface_tension'][1], None)
    def surface_tension(self):
        """Surface tension."""
        return self._AbstractState.surface_tension()

    @property
    @ureg.wraps(CP_const_unit['P'][1], None)
    def P_reducing(self):
        """Reducing pressure."""
        try:
            return self._AbstractState.keyed_output(CP.iP_reducing)
        except ValueError:
            # REFPROP does not expose reducing pressure for pure fluids;
            # for a pure fluid the reducing and critical states coincide.
            if len(self.mole_fractions) == 1:
                return self._AbstractState.p_critical()
            raise

    @property
    @ureg.wraps(CP_const_unit['P'][1], None)
    def P_triple(self):
        """Triple-point pressure."""
        return self._AbstractState.keyed_output(CP.iP_triple)

    @property
    @ureg.wraps(CP_const_unit['P'][1], None)
    def P_max(self):
        """Maximum pressure supported by the backend."""
        return self._AbstractState.pmax()

    @property
    @ureg.wraps(CP_const_unit['T'][1], None)
    def T_reducing(self):
        """Reducing temperature."""
        return self._AbstractState.T_reducing()

    @property
    @ureg.wraps(CP_const_unit['T'][1], None)
    def T_triple(self):
        """Triple-point temperature."""
        return self._AbstractState.Ttriple()

    @property
    @ureg.wraps(CP_const_unit['T'][1], None)
    def T_max(self):
        """Maximum temperature supported by the backend."""
        return self._AbstractState.Tmax()

    @property
    @ureg.wraps(CP_const_unit['T'][1], None)
    def T_min(self):
        """Minimum temperature supported by the backend."""
        return self._AbstractState.Tmin()

    @property
    @ureg.wraps(CP_const_unit['isothermal_compressibility'][1], None)
    def isothermal_compressibility(self):
        """Isothermal compressibility."""
        return self._AbstractState.isothermal_compressibility()

    @property
    @ureg.wraps(CP_const_unit['isobaric_expansion_coefficient'][1], None)
    def isobaric_expansion_coefficient(self):
        """Isobaric thermal expansion coefficient."""
        return self._AbstractState.isobaric_expansion_coefficient()

    @property
    @ureg.wraps(CP_const_unit['isentropic_expansion_coefficient'][1],
                None)
    def isentropic_expansion_coefficient(self):
        """Isentropic expansion coefficient."""
        return self._AbstractState.keyed_output(
            CP.iisentropic_expansion_coefficient)

    @property
    def Prandtl(self):
        """Prandtl number."""
        return self._AbstractState.Prandtl()

    @property
    @ureg.wraps(CP_const_unit['speed_sound'][1], None)
    def speed_sound(self):
        """Speed of sound at the current thermodynamic state."""
        return self._AbstractState.speed_sound()

    def first_partial_deriv(self, Of, Wrt, Constant):
        """Evaluate a first partial derivative exposed by CoolProp."""
        output_unit = CP_const_unit[Of][1] / CP_const_unit[Wrt][1]
        Of_CP_const = CP_const_unit[Of][0]
        Wrt_CP_const = CP_const_unit[Wrt][0]
        Constant_CP_const = CP_const_unit[Constant][0]
        result = self._AbstractState.first_partial_deriv(
            Of_CP_const, Wrt_CP_const, Constant_CP_const)
        return result * output_unit

    @property
    def phase(self):
        """Calculate the phase of the fluid.

        * 0: Subcritical liquid
        * 1: Supercritical (p > pc, T > Tc)
        * 2: Supercritical gas (p < pc, T > Tc)
        * 3: Supercritical liquid (p > pc, T < Tc)
        * 4: At the critical point.
        * 5: Subcritical gas.
        * 6: Twophase.
        * 7: Unknown phase
        """
        return self._AbstractState.phase()

    @property
    def phase_str(self):
        """Return string representation of the fluid phase."""
        phases = {
            0: 'Subcritical liquid',
            1: 'Supercritical (p > pc, T > Tc)',
            2: 'Supercritical gas (p < pc, T > Tc)',
            3: 'Supercritical liquid (p > pc, T < Tc)',
            4: 'At the critical point',
            5: 'Subcritical gas',
            6: 'Two-phase',
            7: 'Unknown phase',
        }
        return phases[self.phase]

    def set_mole_fractions(self, *fractions):
        """Set mole fractions for a mixture."""
        self._AbstractState.set_mole_fractions(fractions)

    def set_mass_fractions(self, *fractions):
        """Set mass fractions for a mixture."""
        self._AbstractState.set_mass_fractions(fractions)

    def set_volu_fractions(self, *fractions):
        """Set volume fractions for a mixture."""
        self._AbstractState.set_volu_fractions(fractions)

    @property
    def mole_fractions(self):
        """Mole fractions of a mixture."""
        return self._AbstractState.get_mole_fractions()

    @property
    def mass_fractions(self):
        """Mass fractions of a mixture."""
        return self._AbstractState.get_mass_fractions()

    @property
    def is_super_critical(self):
        """Return True if state is supercritical."""
        return self.phase in [1, 2, 3]

    @property
    def specific_heat_input(self):
        """
        Calculate Specific heat input, v * (dh/dv)|p.
        This function is not described in AbstractState class.
        """
        # Because specific volume is not available as function of
        # the AbstractState, density is used instead
        # The resulting function is: -Dmass*(dHmass/dDmass)|p
        specific_heat_input = getattr(self._AbstractState,
                                      'specific_heat_input', None)
        if specific_heat_input is not None:
            return specific_heat_input() * ureg.joule / ureg.kg
        return (-self.Dmass) * self.first_partial_deriv(
            'Hmass', 'Dmass', 'P')

    @property
    def gamma(self):
        """Calculate current-state gamma = k = Cp/Cv."""
        warnings.warn(
            'gamma behavior changed: it now returns current-state '
            'Cpmolar/Cvmolar. For an NTP value, calculate it explicitly from a '
            'standard-condition state.',
            CryoToolBoxFutureWarning,
            stacklevel=2,
        )
        return self.Cpmolar / self.Cvmolar

    @property
    def C_gas_const(self):
        """
        Constant for gas or vapor which is the function of the ratio of
        specific heats k = Cp/Cv. ASME VIII.1-2015 pp. 423-424.
        """
        standard_state = self.to_standard(conditions='NTP')
        k_ = standard_state.gamma
        root = (k_ * (2 / (k_ + 1))**((k_ + 1) / (k_ - 1)))**0.5
        R = ureg.R  # Universal gas constant
        C_ = root / R**0.5 * (ureg.g / ureg.mol)**0.5
        return C_

    @property
    def C(self):
        """A shortcut to C_gas_const"""
        return self.C_gas_const

    @property
    def C_us(self):
        """A shortcut to C_gas_const for US customary units"""
        us_unit = ureg.lb / (ureg.hr * ureg.lbf) * (ureg.degR)**0.5
        return self.C_gas_const.to(us_unit)

    @property
    def C_si(self):
        """A shortcut to C_gas_const for SI units"""
        si_unit = ureg.s**2 / (ureg.hr * ureg.m) * (ureg.K)**0.5
        return self.C_gas_const.to(si_unit)

    @property
    def name(self):
        """
        Return fluid name (backend dependent)
        """
        return self._AbstractState.name().lower()

    @property
    def backend(self):
        """
        Return backend name
        """
        return self._AbstractState.backend_name()

    def __str__(self):
        if self.Q == 0:
            return (f'saturated {self.name.lower()} liquid '
                    f'at P: {self.P.to(ureg.bar):.3g~}')
        elif self.Q == 1:
            return (f'saturated {self.name.lower()} vapor '
                    f'at P: {self.P.to(ureg.bar):.3g~}')
        elif self.Q > 0 and self.Q < 1:
            return (f'{self.name.lower()} '
                    f'at P: {self.P.to(ureg.bar):.3g~} and '
                    f'Q: {self.Q:.2f}')
        else:
            return (f'{self.name.lower()} at '
                    f'T: {self.T.to(ureg.K):.3g~} and '
                    f'P: {self.P.to(ureg.bar):.3g~}')

    @property
    def M(self):
        """Calculate relative molecular mass."""
        return self.molar_mass.m_as(ureg.g / ureg.mole)

    @property
    def MZT(self):
        """
        Calculate sqrt(M/(ZT)) a commonly used square root group for discharge
        flow calculation.
        """
        MZT_ = (self.M / (self.compressibility_factor * self.T))**0.5
        return MZT_

    @property
    @ureg.wraps(CP_const_unit['Hmass'][1], None)
    def latent_heat(self):
        """Calculate latent heat of evaporation for current quality."""
        warnings.warn(
            'latent_heat behavior changed: it now returns the remaining latent '
            'heat at the current quality (Q=0 -> full latent heat, Q=1 -> 0).',
            CryoToolBoxFutureWarning,
            stacklevel=2,
        )
        if not (0 <= self.Q <= 1):
            raise ValueError(
                f'Latent heat is undefined for Q={self.Q}, {self}.')
        try:
            h_vapor = self._AbstractState.saturated_vapor_keyed_output(CP.iHmass)
        except (AttributeError, NotImplementedError, ValueError):
            saturation_self = ThermState(self.name, backend=self.backend)
            saturation_self.update_si('P', self._AbstractState.p(), 'Q', 1)
            h_vapor = saturation_self.get_si('Hmass')
        return h_vapor - self._AbstractState.hmass()

    def copy(self):
        """Create a copy of current ThermState object."""
        temp_state = ThermState(self.name, backend=self.backend)
        density = self._AbstractState.rhomass()
        if density == -inf:
            return temp_state

        quality = self._AbstractState.Q()
        if 0 <= quality <= 1:
            temp_state.update_si('P', self._AbstractState.p(), 'Q', quality)
        else:
            temp_state.update_si('Dmass', density, 'T', self._AbstractState.T())
        return temp_state

    def to_standard(self, conditions='NTP'):
        """Create a copy of current ThermState object at standard conditions.
        """
        if conditions == 'NTP':
            P = P_NTP
            T = T_NTP
        elif conditions == 'MSC':
            P = P_MSC
            T = T_MSC
        elif conditions == 'STD':
            P = P_STD
            T = T_STD
        else:
            raise ValueError(f'Conditions {conditions!r} are not defined.')
        TempState = self.copy()
        TempState.update_kw(P=P, T=T)
        return TempState

    @property
    @ureg.wraps(CP_const_unit['P_sat'][1], None)
    def P_sat(self):
        """Saturation pressure at the current temperature."""
        try:
            P_sat = self._AbstractState.saturation_ancillary(
                CP.iP, 0, CP.iT, self._AbstractState.T())
            return P_sat
        except (AttributeError, NotImplementedError, ValueError):
            saturation_state = ThermState(self.name, backend=self.backend)
            saturation_state.update_si('Q', 0, 'T', self._AbstractState.T())
            return saturation_state._AbstractState.p()

    @property
    @ureg.wraps(CP_const_unit['T_sat'][1], None)
    def T_sat(self):
        """Saturation temperature at the current pressure."""
        try:
            T_sat = self._AbstractState.saturation_ancillary(
                CP.iT, 0, CP.iP, self._AbstractState.p())
            return T_sat
        except (AttributeError, NotImplementedError, ValueError):
            saturation_state = ThermState(self.name, backend=self.backend)
            saturation_state.update_si('P', self._AbstractState.p(), 'Q', 0)
            return saturation_state._AbstractState.T()
