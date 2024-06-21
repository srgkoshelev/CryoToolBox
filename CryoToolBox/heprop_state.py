import numpy as np
import ctypes


HEPAK_inputs = {
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

# Load the HEPAK DLL
def init_he():
    try:
        #dll = ctypes.WinDLL('./hepak.dll')
        dll = ctypes.CDLL('./hepak.dll')
        err='OK'
    except OSError as el:
        print("Helium dll not present in the folder") ###should raise an error
        dll = None
        err='No'
    return dll, err

def hecalc(j1, value1, j2, value2, unit, dll, err):
    if isinstance(j1, str):
            j1_ptr = ctypes.c_int(HEPAK_inputs[j1])
    else:
        j1_ptr = ctypes.c_int(j1)
    value1_ptr = ctypes.c_double(value1)
    if isinstance(j2, str):
        j2_ptr = ctypes.c_int(HEPAK_inputs[j2])
    else:
        j2_ptr = ctypes.c_int(j2)
    value2_ptr = ctypes.c_double(value2)
    unit_ptr = ctypes.c_int(unit)
    PROP2 = np.zeros((3,42), dtype=np.float64)
    ID = ctypes.c_int()
    if err == 'OK':
        try:
            dll.HEPROP(ctypes.byref(j1_ptr), ctypes.byref(value1_ptr), ctypes.byref(j2_ptr), ctypes.byref(value2_ptr), ctypes.byref(unit_ptr), PROP2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ctypes.byref(ID))
        except:
            dll.heprop_(ctypes.byref(j1_ptr), ctypes.byref(value1_ptr), ctypes.byref(j2_ptr), ctypes.byref(value2_ptr), ctypes.byref(unit_ptr), PROP2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ctypes.byref(ID))
        idi = ID.value
    return PROP2 if err == 'OK' else None


class HepropState:
    def __init__(self, **state_parameters):
        self.dll, self.err = init_he()
        self._heprop = None

    def update(self, name1, value1, name2, value2):
        hepak = hecalc(name1, value1, name2, value2, 1, self.dll, self.err)
        self._heprop = hepak

    def T_critical(self):
        return 5.1953

    def p_critical(self):
        return 227462.3

    def rhomolar_critical(self):
        return 18130.0

    def rhomass_critical(self):
        return 72.56717426

    def T(self):
        return self._heprop[0][2]

    def rhomolar(self):
        return self._heprop[0][3]/0.004002602

    def rhomass(self):
        return self._heprop[0][3]

    def p(self):
        return self._heprop[0][1]

    def Q(self):
        return self._heprop[0][0]

    def molar_mass(self):
        return 0.004002602

    def gas_constant(self):
        return 8.3144598

    def compressibility_factor(self):
        return self._heprop[0][5]

    def hmolar(self):
        return self._heprop[0][9]*self.molar_mass()

    def hmass(self):
        return self._heprop[0][9]

    def smolar(self):
        return self._heprop[0][8]*self.molar_mass()

    def smass(self):
        return self._heprop[0][8]

    def cpmolar(self):
        return self._heprop[0][14]*self.molar_mass()

    def cpmass(self):
        return self._heprop[0][14]

    def cvmolar(self):
        return self._heprop[0][15]*self.molar_mass()

    def cvmass(self):
        return self._heprop[0][15]

    def viscosity(self):
        return self._heprop[0][25]

    def conductivity(self):
        return self._heprop[0][26]

    def surface_tension(self):
        return self._heprop[0][29]

    def p_reducing(self):
        return 0 #value to add, not existing in Coolprop

    def ptriple(self):
        return 0 #value to add, not existing in Coolprop

    def pmax(self):
        return 2028000000

    def T_reducing(self):
        return 5.1953

    def Ttriple(self):
        return 2.1768

    def Tmax(self):
        return 1500

    def Tmin(self):
        return 0.8

    def isothermal_compressibility(self):
        return self._heprop[0][19]

    def isobaric_expansion_coefficient(self):
        return self._heprop[0][17]/self._heprop[0][2]

    def isentropic_expansion_coefficient(self):
        return 0 #value to add, I need to verify this part with Heprop

    def Prandtl(self):
        return self._heprop[0][27]

    def speed_sound(self):
        return self._heprop[0][20]

    def name(self):
        return 'helium'

    def backend_name(self):
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
        elif self._heprop[0][2] > self.T_critical() and self._heprop[0][1] > self.p_critical():
            return 1
        elif self._heprop[0][2] > self.T_critical() and self._heprop[0][1] < self.p_critical():
            return 2
        elif self._heprop[0][2] < self.T_critical() and self._heprop[0][1] > self.p_critical():
            return 3
        elif self._heprop[0][2] == self.T_critical() and self._heprop[0][1] == self.p_critical():
            return 4
        elif self._heprop[0][0] >= 1:
            return 5
        elif self._heprop[0][0] > 0 and self._heprop[0][0] < 1:
            return 6
        else:
            return 7

    def specific_heat_input(self):
        return self._heprop[0][24]
