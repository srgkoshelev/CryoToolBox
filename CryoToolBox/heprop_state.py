"""Thin adapter around the optional external HEPROP helium property library."""

import numpy as np
import ctypes
import os

HEPROP_inputs = {
    'P': 1,
    'T': 2,
    'D': 3,
    'V': 4,
    'S': 5,
    'Smass': 5,
    'H': 6,
    'Hmass': 6,
    'U': 7,
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
    idi = ID.value
    return PROP2


class HepropState:
    """Expose the subset of the CoolProp ``AbstractState`` API used here.

    This adapter is only intended for the optional ``HEPROP`` backend in
    :class:`CryoToolBox.cp_wrapper.ThermState`.
    """

    def __init__(self, **state_parameters):
        self.dll = _init_he()
        self._heprop = None

    def update(self, name1, value1, name2, value2):
        """Update the cached HEPROP state from two named inputs."""
        heprop = _hecalc(name1, value1, name2, value2, 1, self.dll)
        self._heprop = heprop

    def T_critical(self):
        """Critical temperature for helium."""
        return 5.1953

    def p_critical(self):
        """Critical pressure for helium."""
        return 227462.3

    def rhomolar_critical(self):
        """Critical molar density for helium."""
        return 18130.0

    def rhomass_critical(self):
        """Critical mass density for helium."""
        return 72.56717426

    def T(self):
        """Temperature of the current state."""
        return self._heprop[0][2]

    def rhomolar(self):
        """Molar density of the current state."""
        return self._heprop[0][3] / 0.004002602

    def rhomass(self):
        """Mass density of the current state."""
        return self._heprop[0][3]

    def p(self):
        """Pressure of the current state."""
        return self._heprop[0][1]

    def Q(self):
        """Vapor quality."""
        return self._heprop[0][0]

    def molar_mass(self):
        """Helium molar mass."""
        return 0.004002602

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
        return self._heprop[0][29]

    def p_reducing(self):
        """Reducing pressure placeholder for API compatibility."""
        return 0  #value to add, not existing in Coolprop

    def ptriple(self):
        """Triple-point pressure placeholder for API compatibility."""
        return 0  #value to add, not existing in Coolprop

    def pmax(self):
        """Maximum pressure supported by HEPROP."""
        return 2028000000

    def T_reducing(self):
        """Reducing temperature."""
        return 5.1953

    def Ttriple(self):
        """Triple-point temperature."""
        return 2.1768

    def Tmax(self):
        """Maximum temperature supported by HEPROP."""
        return 1500

    def Tmin(self):
        """Minimum temperature supported by HEPROP."""
        return 0.8

    def isothermal_compressibility(self):
        """Isothermal compressibility."""
        return self._heprop[0][19]

    def isobaric_expansion_coefficient(self):
        """Isobaric thermal expansion coefficient."""
        return self._heprop[0][17] / self._heprop[0][2]

    def isentropic_expansion_coefficient(self):
        """Isentropic expansion coefficient placeholder."""
        return 0  #value to add, I need to verify this part with Heprop

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
        if self._heprop[0][0] <= 0 and self._heprop[0][1] < self.p_critical():
            return 0
        elif self._heprop[0][2] > self.T_critical(
        ) and self._heprop[0][1] > self.p_critical():
            return 1
        elif self._heprop[0][2] > self.T_critical(
        ) and self._heprop[0][1] < self.p_critical():
            return 2
        elif self._heprop[0][2] < self.T_critical(
        ) and self._heprop[0][1] > self.p_critical():
            return 3
        elif self._heprop[0][2] == self.T_critical(
        ) and self._heprop[0][1] == self.p_critical():
            return 4
        elif self._heprop[0][0] >= 1:
            return 5
        elif self._heprop[0][0] > 0 and self._heprop[0][0] < 1:
            return 6
        else:
            return 7

    def specific_heat_input(self):
        """Specific heat input term used by the CGA calculations."""
        return self._heprop[0][24]
