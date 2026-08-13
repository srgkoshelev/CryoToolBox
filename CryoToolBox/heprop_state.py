"""Thin adapter around the optional external HEPROP helium property library."""

import numpy as np
import ctypes
import os
import CoolProp.CoolProp as CP

MOLAR_MASS = 0.004002602
T_CRITICAL = 5.1953
P_CRITICAL = 227462.3
DMOLAR_CRITICAL = 18130.0
DMASS_CRITICAL = 72.56717426
T_TRIPLE = 2.1768
P_TRIPLE = 5041.8
T_MAX = 1500.0
T_MIN = 0.8
P_MAX = 2.028e9
_QMASS_PARAMETER = getattr(CP, 'iQmass', None)

HEPROP_inputs = {
    'P': 1,
    'T': 2,
    'D': 3,
    'Dmass': 3,
    'V': 4,
    'S': 5,
    'Smass': 5,
    'H': 6,
    'Hmass': 6,
    'U': 7,
    'Umass': 7,
    'G': 8,
    'X': 9,
    'Q': 9,
    'dT': 10,
    'M': 11,
    '2': 12,
    'sL': 13,
    'sV': 14,
    'L': 15
}

HEPROP_input_pairs = {
    CP.QT_INPUTS: ('Q', 'T'),
    CP.PQ_INPUTS: ('P', 'Q'),
    CP.QSmolar_INPUTS: ('Q', 'Smolar'),
    CP.QSmass_INPUTS: ('Q', 'Smass'),
    CP.HmolarQ_INPUTS: ('Hmolar', 'Q'),
    CP.HmassQ_INPUTS: ('Hmass', 'Q'),
    CP.DmolarQ_INPUTS: ('Dmolar', 'Q'),
    CP.DmassQ_INPUTS: ('Dmass', 'Q'),
    CP.PT_INPUTS: ('P', 'T'),
    CP.DmassT_INPUTS: ('Dmass', 'T'),
    CP.DmolarT_INPUTS: ('Dmolar', 'T'),
    CP.HmolarT_INPUTS: ('Hmolar', 'T'),
    CP.HmassT_INPUTS: ('Hmass', 'T'),
    CP.SmolarT_INPUTS: ('Smolar', 'T'),
    CP.SmassT_INPUTS: ('Smass', 'T'),
    CP.TUmolar_INPUTS: ('T', 'Umolar'),
    CP.TUmass_INPUTS: ('T', 'Umass'),
    CP.DmassP_INPUTS: ('Dmass', 'P'),
    CP.DmolarP_INPUTS: ('Dmolar', 'P'),
    CP.HmassP_INPUTS: ('Hmass', 'P'),
    CP.HmolarP_INPUTS: ('Hmolar', 'P'),
    CP.PSmass_INPUTS: ('P', 'Smass'),
    CP.PSmolar_INPUTS: ('P', 'Smolar'),
    CP.PUmass_INPUTS: ('P', 'Umass'),
    CP.PUmolar_INPUTS: ('P', 'Umolar'),
    CP.HmassSmass_INPUTS: ('Hmass', 'Smass'),
    CP.HmolarSmolar_INPUTS: ('Hmolar', 'Smolar'),
    CP.SmassUmass_INPUTS: ('Smass', 'Umass'),
    CP.SmolarUmolar_INPUTS: ('Smolar', 'Umolar'),
    CP.DmassHmass_INPUTS: ('Dmass', 'Hmass'),
    CP.DmolarHmolar_INPUTS: ('Dmolar', 'Hmolar'),
    CP.DmassSmass_INPUTS: ('Dmass', 'Smass'),
    CP.DmolarSmolar_INPUTS: ('Dmolar', 'Smolar'),
    CP.DmassUmass_INPUTS: ('Dmass', 'Umass'),
    CP.DmolarUmolar_INPUTS: ('Dmolar', 'Umolar'),
}

if _QMASS_PARAMETER is not None:
    HEPROP_input_pairs.update({
        CP.QmassT_INPUTS: ('Qmass', 'T'),
        CP.PQmass_INPUTS: ('P', 'Qmass'),
        CP.QmassSmolar_INPUTS: ('Qmass', 'Smolar'),
        CP.QmassSmass_INPUTS: ('Qmass', 'Smass'),
        CP.HmolarQmass_INPUTS: ('Hmolar', 'Qmass'),
        CP.HmassQmass_INPUTS: ('Hmass', 'Qmass'),
        CP.DmolarQmass_INPUTS: ('Dmolar', 'Qmass'),
        CP.DmassQmass_INPUTS: ('Dmass', 'Qmass'),
    })


def _to_heprop_input(name, value):
    """Convert a CoolProp low-level input name/value to HEPROP form."""
    if name == 'Dmolar':
        return 'Dmass', value * MOLAR_MASS
    if name == 'Hmolar':
        return 'Hmass', value / MOLAR_MASS
    if name == 'Smolar':
        return 'Smass', value / MOLAR_MASS
    if name == 'Umolar':
        return 'Umass', value / MOLAR_MASS
    if name == 'Qmass':
        return 'Q', value
    return name, value


def _heprop_phase(heprop):
    """Return a CoolProp phase code for a raw HEPROP result table."""
    quality, pressure, temperature = heprop[0][:3]
    if quality <= 0 and pressure < P_CRITICAL:
        return CP.iphase_liquid
    if temperature > T_CRITICAL and pressure > P_CRITICAL:
        return CP.iphase_supercritical
    if temperature > T_CRITICAL and pressure < P_CRITICAL:
        return CP.iphase_supercritical_gas
    if temperature < T_CRITICAL and pressure > P_CRITICAL:
        return CP.iphase_supercritical_liquid
    if temperature == T_CRITICAL and pressure == P_CRITICAL:
        return CP.iphase_critical_point
    if quality >= 1:
        return CP.iphase_gas
    if 0 < quality < 1:
        return CP.iphase_twophase
    return CP.iphase_unknown


def _keyed_output_from_table(heprop, parameter):
    """Return a CoolProp-style keyed output from a raw HEPROP result table."""
    state = heprop[0]
    if _QMASS_PARAMETER is not None and parameter == _QMASS_PARAMETER:
        return state[0]
    if parameter == CP.isurface_tension and not 0 <= state[0] <= 1:
        raise ValueError(
            'Surface tension is only defined in the two-phase region.')
    outputs = {
        CP.iT: state[2],
        CP.iP: state[1],
        CP.iQ: state[0],
        CP.iDmolar: state[3] / MOLAR_MASS,
        CP.iDmass: state[3],
        CP.iHmolar: state[9] * MOLAR_MASS,
        CP.iHmass: state[9],
        CP.iSmolar: state[8] * MOLAR_MASS,
        CP.iSmass: state[8],
        CP.iCpmolar: state[14] * MOLAR_MASS,
        CP.iCpmass: state[14],
        CP.iCvmolar: state[15] * MOLAR_MASS,
        CP.iCvmass: state[15],
        CP.iviscosity: state[25],
        CP.iconductivity: state[26],
        CP.isurface_tension: heprop[1][29],
        CP.iPrandtl: state[27],
        CP.ispeed_sound: state[20],
        CP.iisothermal_compressibility: state[19],
        CP.iisobaric_expansion_coefficient: state[17] / state[2],
        CP.iZ: state[5],
        CP.iPhase: _heprop_phase(heprop),
        CP.iP_reducing: P_CRITICAL,
        CP.iP_triple: P_TRIPLE,
        CP.iP_max: P_MAX,
        CP.iT_reducing: T_CRITICAL,
        CP.iT_triple: T_TRIPLE,
        CP.iT_max: T_MAX,
        CP.iT_min: T_MIN,
        CP.iisentropic_expansion_coefficient:
            state[3] * state[20]**2 / state[1],
    }
    try:
        return outputs[parameter]
    except KeyError as exc:
        raise ValueError(
            f'HEPROP adapter does not implement keyed output {parameter!r}.'
        ) from exc


def _resolve_heprop_path():
    dll_path = os.getenv('HEPROPPATH')
    if not dll_path:
        raise OSError(
            'HEPROP backend is not configured. Set the HEPROPPATH '
            'environment variable to the HEPROP shared library.'
        )
    if not os.path.isfile(dll_path):
        raise OSError(f'HEPROP library not found at {dll_path}')
    return dll_path


def _init_he():
    """Load the configured HEPROP library and verify expected symbols."""
    dll_path = _resolve_heprop_path()
    try:
        # dll = ctypes.WinDLL(dll_path)
        dll = ctypes.CDLL(dll_path)
    except OSError as exc:
        raise OSError(f'Failed to load HEPROP library at {dll_path}') from exc

    for symbol in ('HEPROP', 'heprop_'):
        if hasattr(dll, symbol):
            return dll
    raise OSError(
        'Loaded HEPROP library does not expose the expected HEPROP or '
        'heprop_ symbols.'
    )


def _hecalc(j1, value1, j2, value2, unit, dll):
    """Internal wrapper around the HEPROP DLL entry points."""
    # No fallback mechanism, one of two expected at all times
    if hasattr(dll, 'HEPROP'):
        heprop_func = dll.HEPROP
    elif hasattr(dll, 'heprop_'):
        heprop_func = dll.heprop_
    else:
        raise OSError(
            'Loaded HEPROP library does not expose the expected HEPROP or '
            'heprop_ symbols.'
        )

    if isinstance(j1, str):
        j1_ptr = ctypes.c_int(HEPROP_inputs[j1])
    else:
        j1_ptr = ctypes.c_int(j1)
    value1_ptr = ctypes.c_double(value1)
    if isinstance(j2, str):
        j2_ptr = ctypes.c_int(HEPROP_inputs[j2])
    else:
        j2_ptr = ctypes.c_int(j2)
    value2_ptr = ctypes.c_double(value2)
    unit_ptr = ctypes.c_int(unit)
    PROP2 = np.zeros((3, 42), dtype=np.float64)
    ID = ctypes.c_int()

    try:
        heprop_func(ctypes.byref(j1_ptr), ctypes.byref(value1_ptr),
                    ctypes.byref(j2_ptr), ctypes.byref(value2_ptr),
                    ctypes.byref(unit_ptr),
                    PROP2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    ctypes.byref(ID))
    except Exception as exc:
        raise OSError('HEPROP backend call failed.') from exc
    return PROP2


class HepropState:
    """Expose the subset of the CoolProp ``AbstractState`` API used here.

    This adapter is only intended for the optional ``HEPROP`` backend in
    :class:`CryoToolBox.cp_wrapper.ThermState`.
    """

    def __init__(self):
        self.dll = _init_he()
        self._heprop = None

    def update(self, input_pair, value1, value2):
        """Update the cached state using the CoolProp AbstractState API."""
        try:
            name1, name2 = HEPROP_input_pairs[input_pair]
        except KeyError as exc:
            raise ValueError(
                f'HEPROP adapter does not support input pair {input_pair!r}.'
            ) from exc
        name1, value1 = _to_heprop_input(name1, value1)
        name2, value2 = _to_heprop_input(name2, value2)
        self._heprop = _hecalc(name1, value1, name2, value2, 1, self.dll)

    def keyed_output(self, parameter):
        """Return a CoolProp-style keyed output for the current state."""
        if self._heprop is None:
            return -np.inf
        return _keyed_output_from_table(self._heprop, parameter)

    def T_critical(self):
        """Critical temperature for helium."""
        return T_CRITICAL

    def p_critical(self):
        """Critical pressure for helium."""
        return P_CRITICAL

    def rhomolar_critical(self):
        """Critical molar density for helium."""
        return DMOLAR_CRITICAL

    def rhomass_critical(self):
        """Critical mass density for helium."""
        return DMASS_CRITICAL

    def T(self):
        """Temperature of the current state."""
        return self.keyed_output(CP.iT)

    def rhomolar(self):
        """Molar density of the current state."""
        return self.keyed_output(CP.iDmolar)

    def rhomass(self):
        """Mass density of the current state."""
        return self.keyed_output(CP.iDmass)

    def p(self):
        """Pressure of the current state."""
        return self.keyed_output(CP.iP)

    def Q(self):
        """Vapor quality."""
        return self.keyed_output(CP.iQ)

    def Qmass(self):
        """Mass-basis vapor quality."""
        return self.Q()

    def molar_mass(self):
        """Helium molar mass."""
        return MOLAR_MASS

    def gas_constant(self):
        """Universal gas constant used by HEPROP."""
        return 8.3144598

    def compressibility_factor(self):
        """Compressibility factor ``Z``."""
        return self._heprop[0][5]

    def hmolar(self):
        """Molar enthalpy."""
        return self._heprop[0][9] * self.molar_mass()

    def hmass(self):
        """Specific enthalpy on a mass basis."""
        return self._heprop[0][9]

    def smolar(self):
        """Molar entropy."""
        return self._heprop[0][8] * self.molar_mass()

    def smass(self):
        """Specific entropy on a mass basis."""
        return self._heprop[0][8]

    def cpmolar(self):
        """Constant-pressure heat capacity on a molar basis."""
        return self._heprop[0][14] * self.molar_mass()

    def cpmass(self):
        """Constant-pressure heat capacity on a mass basis."""
        return self._heprop[0][14]

    def cvmolar(self):
        """Constant-volume heat capacity on a molar basis."""
        return self._heprop[0][15] * self.molar_mass()

    def cvmass(self):
        """Constant-volume heat capacity on a mass basis."""
        return self._heprop[0][15]

    def viscosity(self):
        """Dynamic viscosity."""
        return self._heprop[0][25]

    def conductivity(self):
        """Thermal conductivity."""
        return self._heprop[0][26]

    def surface_tension(self):
        """Surface tension."""
        if not 0 <= self.Q() <= 1:
            raise ValueError(
                'Surface tension is only defined in the two-phase region.')
        return self._heprop[1][29]

    def p_reducing(self):
        """Reducing pressure for pure helium."""
        return P_CRITICAL

    def ptriple(self):
        """Triple-point pressure."""
        return P_TRIPLE

    def pmax(self):
        """Maximum pressure supported by HEPROP."""
        return P_MAX

    def T_reducing(self):
        """Reducing temperature."""
        return T_CRITICAL

    def Ttriple(self):
        """Triple-point temperature."""
        return T_TRIPLE

    def Tmax(self):
        """Maximum temperature supported by HEPROP."""
        return T_MAX

    def Tmin(self):
        """Minimum temperature supported by HEPROP."""
        return T_MIN

    def isothermal_compressibility(self):
        """Isothermal compressibility."""
        return self._heprop[0][19]

    def isobaric_expansion_coefficient(self):
        """Isobaric thermal expansion coefficient."""
        return self._heprop[0][17] / self._heprop[0][2]

    def isentropic_expansion_coefficient(self):
        """Return ``rho * speed_sound**2 / pressure``."""
        return self.rhomass() * self.speed_sound()**2 / self.p()

    def Prandtl(self):
        """Prandtl number."""
        return self._heprop[0][27]

    def speed_sound(self):
        """Speed of sound."""
        return self._heprop[0][20]

    def name(self):
        """Fluid name."""
        return 'helium'

    def backend_name(self):
        """Backend name used by the compatibility wrapper."""
        return 'HEPROP'

    def phase(self):
        """
        * 0: Subcritical liquid
        * 1: Supercritical (p > pc, T > Tc)
        * 2: Supercritical gas (p < pc, T > Tc)
        * 3: Supercritical liquid (p > pc, T < Tc)
        * 4: At the critical point.
        * 5: Subcritical gas.
        * 6: Two phase.
        * 7: Unknown phase
        """
        return _heprop_phase(self._heprop)

    def specific_heat_input(self):
        """Specific heat input term used by the CGA calculations."""
        return self._heprop[0][24]

    def saturated_liquid_keyed_output(self, parameter):
        """Return a keyed output at saturated liquid on current pressure."""
        heprop = _hecalc('P', self.p(), 'Q', 0, 1, self.dll)
        return _keyed_output_from_table(heprop, parameter)

    def saturated_vapor_keyed_output(self, parameter):
        """Return a keyed output at saturated vapor on current pressure."""
        heprop = _hecalc('P', self.p(), 'Q', 1, 1, self.dll)
        return _keyed_output_from_table(heprop, parameter)

    def saturation_ancillary(self, output, Q, given, value):
        """Approximate CoolProp's saturation_ancillary API for helium."""
        if given == CP.iT:
            heprop = _hecalc('T', value, 'Q', Q, 1, self.dll)
        elif given == CP.iP:
            heprop = _hecalc('P', value, 'Q', Q, 1, self.dll)
        else:
            raise ValueError(
                f'HEPROP adapter does not support saturation input {given!r}.'
            )
        return _keyed_output_from_table(heprop, output)

    def get_mole_fractions(self):
        """Return mole fractions for pure helium."""
        return [1.0]

    def get_mass_fractions(self):
        """Return mass fractions for pure helium."""
        return [1.0]
