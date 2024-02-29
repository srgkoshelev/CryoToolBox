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
        dll = ctypes.WinDLL('./hepak.dll')
        err='OK'
    except OSError as el:
        print("Helium dll not present in the folder")
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
        dll.HEPROP(ctypes.byref(j1_ptr), ctypes.byref(value1_ptr), ctypes.byref(j2_ptr), ctypes.byref(value2_ptr), ctypes.byref(unit_ptr), PROP2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ctypes.byref(ID))
        idi = ID.value
        # if idi < 0:
        #     print(f"error IDID = {idi}")
    return PROP2 if err == 'OK' else None


class HepropState:
    def __init__(self, **state_parameters):
        self.name = 'helium'
        # self.dll, self.err = init_he()
        self._heprop = None
        self._heprop_mol = None

    def T_critical(self):
        return 5.1953

    def p_critical(self):
        return 227462.3

    def Q(self):
        return self._heprop[0][0]
