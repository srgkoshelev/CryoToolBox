"""Wrapper module for CoolProp's AbstractState.

1. Provides unit support
2. Adds some small useful functions

"""
import CoolProp.CoolProp as CP
from .heprop_state import HepropState
from math import inf
from .std_conditions import ureg, T_NTP, P_NTP, P_MSC, T_MSC, P_STD, T_STD

CP_const_unit = {
    'gas_constant': (CP.igas_constant, ureg.J/ureg.mol/ureg.K),
    'molar_mass': (CP.imolar_mass, ureg.kg/ureg.mol),
    # 'acentric_factor': (CP.iacentric_factor, ),
    # 'Dmolar_reducing': (CP.irhomolar_reducing, ),
    'Dmolar_critical': (CP.irhomolar_critical, ureg.mol/ureg.m**3),
    'T_reducing': (CP.iT_reducing, ureg.K),
    'T_critical': (CP.iT_critical, ureg.K),
    # 'Dmass_reducing': (CP.irhomass_reducing, ),
    'Dmass_critical': (CP.irhomass_critical, ureg.kg/ureg.m**3),
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
    'Q': (CP.iQ, ureg.mol/ureg.mol),
    # 'Tau': (CP.iTau, ),
    # 'Delta': (CP.iDelta, ),
    'Dmolar': (CP.iDmolar, ureg.mol/ureg.m**3),
    'Hmolar': (CP.iHmolar, ureg.J/ureg.mol),
    'Smolar': (CP.iSmolar, ureg.J/ureg.mol/ureg.K),
    'Cpmolar': (CP.iCpmolar, ureg.J/ureg.mol/ureg.K),
    # 'Cp0molar': (CP.iCp0molar, ureg.J/ureg.mol/ureg.K),
    'Cvmolar': (CP.iCvmolar, ureg.J/ureg.mol/ureg.K),
    'Umolar': (CP.iUmolar, ureg.J/ureg.mol),
    # 'Gmolar': (CP.iGmolar, ),
    # 'Helmholtzmolar': (CP.iHelmholtzmolar, ),
    # 'Smolar_residual': (CP.iSmolar_residual, ),
    'Dmass': (CP.iDmass, ureg.kg/ureg.m**3),
    'Hmass': (CP.iHmass, ureg.J/ureg.kg),
    'Smass': (CP.iSmass, ureg.J/ureg.kg/ureg.K),
    'Cpmass': (CP.iCpmass, ureg.J/ureg.kg/ureg.K),
    # 'Cp0mass': (CP.iCp0mass, ureg.J/ureg.kg/ureg.K),
    'Cvmass': (CP.iCvmass, ureg.J/ureg.kg/ureg.K),
    'Umass': (CP.iUmass, ureg.J/ureg.kg),
    # 'Gmass': (CP.iGmass, ),
    # 'Helmholtzmass': (CP.iHelmholtzmass, ),
    'viscosity': (CP.iviscosity, ureg.Pa*ureg.s),
    'conductivity': (CP.iconductivity, ureg.W/ureg.m/ureg.K),
    'surface_tension': (CP.isurface_tension, ureg.N/ureg.m),
    'Prandtl': (CP.iPrandtl, None),
    'speed_sound': (CP.ispeed_sound, ureg.m/ureg.s),
    'isothermal_compressibility': (CP.iisothermal_compressibility,
                                   ureg.Pa**-1),
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
    'C_us': (None,  # unit defined for external use
             ureg.lb/(ureg.hr*ureg.lbf)*(ureg.degR)**0.5),
    'C_m': (None,  # unit defined for external use
            ureg.s**2/(ureg.hr*ureg.m)*(ureg.K)**0.5)
}

CP_inputs = {
    'QT': CP.QT_INPUTS,
    'PQ': CP.PQ_INPUTS,
    'QSmolar': CP.QSmolar_INPUTS,
    'QSmass': CP.QSmass_INPUTS,
    'HmolarQ': CP.HmolarQ_INPUTS,
    'HmassQ': CP.HmassQ_INPUTS,
    'DmolarQ': CP.DmolarQ_INPUTS,
    'DmassQ': CP.DmassQ_INPUTS,
    'PT': CP.PT_INPUTS,
    'DmassT': CP.DmassT_INPUTS,
    'DmolarT': CP.DmolarT_INPUTS,
    'HmolarT': CP.HmolarT_INPUTS,
    'HmassT': CP.HmassT_INPUTS,
    'SmolarT': CP.SmolarT_INPUTS,
    'SmassT': CP.SmassT_INPUTS,
    'TUmolar': CP.TUmolar_INPUTS,
    'TUmass': CP.TUmass_INPUTS,
    'DmassP': CP.DmassP_INPUTS,
    'DmolarP': CP.DmolarP_INPUTS,
    'HmassP': CP.HmassP_INPUTS,
    'HmolarP': CP.HmolarP_INPUTS,
    'PSmass': CP.PSmass_INPUTS,
    'PSmolar': CP.PSmolar_INPUTS,
    'PUmass': CP.PUmass_INPUTS,
    'PUmolar': CP.PUmolar_INPUTS,
    'HmassSmass': CP.HmassSmass_INPUTS,
    'HmolarSmolar': CP.HmolarSmolar_INPUTS,
    'SmassUmass': CP.SmassUmass_INPUTS,
    'SmolarUmolar': CP.SmolarUmolar_INPUTS,
    'DmassHmass': CP.DmassHmass_INPUTS,
    'DmolarHmolar': CP.DmolarHmolar_INPUTS,
    'DmassSmass': CP.DmassSmass_INPUTS,
    'DmolarSmolar': CP.DmolarSmolar_INPUTS,
    'DmassUmass': CP.DmassUmass_INPUTS,
    'DmolarUmolar': CP.DmolarUmolar_INPUTS,
}


class ThermState:
    def __init__(self, fluid, backend="HEOS", **state_parameters):
        """
        Available backends: HEOS (opensource), REFPROP.
        See http://www.coolprop.org/coolprop/REFPROP.html for details.
        """
        # TODO this should check for the fluid name too
        if backend == "HEPROP":
            self._AbstractState = HepropState
        else:
            self._AbstractState = CP.AbstractState(backend, fluid)
        if state_parameters:
            self.update_kw(**state_parameters)

    def update_kw(self, **state_parameters):
        """Update thermodynamic state using keyword arguments."""
        if len(state_parameters) != 2:
            raise TypeError(f'update_kw() takes 2 arguments \
            ({len(state_parameters)} given)')
        # TODO Code below can be simplified
        names = []
        values = []
        for name, value in state_parameters.items():
            names.append(name)
            values.append(value)
        self.update(names[0], values[0],
                    names[1], values[1])

    def update(self, name1, value1, name2, value2):
        if name1 > name2:  # Sorting inputs alphabetically
            name1, name2 = name2, name1
            value1, value2 = value2, value1
        CP_input_str = name1 + name2
        CP_value1 = self.prepare_input(name1, value1)
        CP_value2 = self.prepare_input(name2, value2)
        # TODO all this can be moved to HepropState
        # if self.name == 'Helium':
        #     hepak = hecalc(name1, CP_value1, name2, CP_value2, 1, self.dll, self.err)
        #     self._heprop = hepak
        #     try:
        #         hepak_satT = hecalc(1, hepak[0][1], 9, 0, 1, self.dll, self.err)
        #         hepak_satP = hecalc(2, hepak[0][2], 9, 0, 1, self.dll, self.err)
        #         self._heprop[1][1] = hepak_satP[0][1]
        #         self._heprop[1][2] = hepak_satT[0][2]
        #     except:
        #         self._heprop = hepak
        #     hepak_mol = hecalc(name1, CP_value1, name2, CP_value2, 3, self.dll, self.err)
        #     self._heprop_mol = hepak_mol
        #     try:
        #         self._AbstractState.update(
        #             CP_inputs[CP_input_str], CP_value1, CP_value2) ###update all properties in coolprop
        #     except:
        #         if hepak[0][1] < 5100:
        #             self._AbstractState.update(CP.PT_INPUTS, 5100, hepak[0][2]) ###update all properties in coolprop
        #         else:
        #             self._AbstractState.update(CP.PT_INPUTS, hepak[0][1], hepak[0][2]) ###update all properties in coolprop
        # TODO this will call the code above from the HepropState class
        self._AbstractState.update(
            CP_inputs[CP_input_str], CP_value1, CP_value2) ###update all properties in coolprop

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
            return value  # Dimensionless
        else:
            unit = CP_const_unit[name][1]
            return value.to(unit).magnitude

    @property
    @ureg.wraps(CP_const_unit['T'][1], None)
    def T_critical(self):
        return self._AbstractState.T_critical()

    @property
    @ureg.wraps(CP_const_unit['P'][1], None)
    def P_critical(self):
        return self._AbstractState.p_critical()


    @property
    @ureg.wraps(CP_const_unit['Dmolar_critical'][1], None)
    def Dmolar_critical(self):
        return self._AbstractState.rhomolar_critical()

    @property
    @ureg.wraps(CP_const_unit['Dmass_critical'][1], None)
    def Dmass_critical(self):
        return self._AbstractState.rhomass_critical()

    @property
    @ureg.wraps(CP_const_unit['T'][1], None)
    def T(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][2]
        # else:
        return self._AbstractState.T()

    @property
    @ureg.wraps(CP_const_unit['P'][1], None)
    def P_sat(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[1][1]
        # else:
        return CP.PropsSI('P','T',self._AbstractState.T(),'Q',0,self.name)

    @property
    @ureg.wraps(CP_const_unit['Dmolar_critical'][1], None)
    def Dmolar(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self.Dmass/self.molar_mass
        # else:
        return self._AbstractState.rhomolar()

    @property
    @ureg.wraps(CP_const_unit['Dmass_critical'][1], None)
    def Dmass(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][3]
        # else:
        return self._AbstractState.rhomass()

    @property
    @ureg.wraps(CP_const_unit['P'][1], None)
    def P(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][1]
        # else:
        return self._AbstractState.p()

    @property
    @ureg.wraps(CP_const_unit['T'][1], None)
    def T_sat(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[1][2]
        # else:
        return CP.PropsSI('T','P',self._AbstractState.p(),'Q',0,self.name)

    @property
    @ureg.wraps(CP_const_unit['Q'][1], None)
    def Q(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][0]
        # else:
        return self._AbstractState.Q()

    @property
    @ureg.wraps(CP_const_unit['molar_mass'][1], None)
    def molar_mass(self):
        return self._AbstractState.molar_mass()

    @property
    @ureg.wraps(CP_const_unit['gas_constant'][1], None)
    def gas_constant(self):
        return self._AbstractState.gas_constant()

    @property
    def specific_gas_constant(self):
        R_spec = self.gas_constant / self.molar_mass
        return R_spec

    @property
    # @ureg.wraps(CP_const_unit['Z'][1], None)
    def compressibility_factor(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     Z_ = self._heprop[0][5]
        # else:
        Z_ = (self.P * self.molar_mass / (self.Dmass*self.gas_constant*self.T)).m_as(ureg.dimensionless)
        return Z_
        # Temporarily unavailable function
        # return self._AbstractState.compressibility_factor()
    # Useful shorthand
    Z = compressibility_factor

    @property
    @ureg.wraps(CP_const_unit['Hmolar'][1], None)
    def Hmolar(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop_mol[0][9]
        # else:
        return self._AbstractState.hmolar()

    @property
    @ureg.wraps(CP_const_unit['Hmass'][1], None)
    def Hmass(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][9]
        # else:
        return self._AbstractState.hmass()

    @property
    @ureg.wraps(CP_const_unit['Smolar'][1], None)
    def Smolar(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop_mol[0][8]
        # else:
        return self._AbstractState.smolar()

    @property
    @ureg.wraps(CP_const_unit['Smass'][1], None)
    def Smass(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][8]
        # else:
        return self._AbstractState.smass()

    @property
    @ureg.wraps(CP_const_unit['Cpmolar'][1], None)
    def Cpmolar(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][14]
        # else:
        return self._AbstractState.cpmolar()

    @property
    @ureg.wraps(CP_const_unit['Cpmass'][1], None)
    def Cpmass(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][14]
        # else:
        return self._AbstractState.cpmass()

    @property
    @ureg.wraps(CP_const_unit['Cvmolar'][1], None)
    def Cvmolar(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][15]
        # else:
        return self._AbstractState.cvmolar()

    @property
    @ureg.wraps(CP_const_unit['Cvmass'][1], None)
    def Cvmass(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][15]
        # else:
        return self._AbstractState.cvmass()

    @property
    @ureg.wraps(CP_const_unit['viscosity'][1], None)
    def viscosity(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][25]
        # else:
        return self._AbstractState.viscosity()

    @property
    @ureg.wraps(CP_const_unit['conductivity'][1], None)
    def conductivity(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][26]
        # else:
        return self._AbstractState.conductivity()

    @property
    @ureg.wraps(CP_const_unit['surface_tension'][1], None)
    def surface_tension(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][29]
        # else:
        return self._AbstractState.surface_tension()

    @property
    @ureg.wraps(CP_const_unit['P'][1], None)
    def P_reducing(self):
        return self._AbstractState.p_reducing()

    @property
    @ureg.wraps(CP_const_unit['P'][1], None)
    def P_triple(self):
        return self._AbstractState.ptriple()

    @property
    @ureg.wraps(CP_const_unit['P'][1], None)
    def P_max(self):
        return self._AbstractState.pmax()

    @property
    @ureg.wraps(CP_const_unit['T'][1], None)
    def T_reducing(self):
        return self._AbstractState.T_reducing()

    @property
    @ureg.wraps(CP_const_unit['T'][1], None)
    def T_triple(self):
        return self._AbstractState.Ttriple()

    @property
    @ureg.wraps(CP_const_unit['T'][1], None)
    def T_max(self):
        return self._AbstractState.Tmax()

    @property
    @ureg.wraps(CP_const_unit['T'][1], None)
    def T_min(self):
        return self._AbstractState.Tmin()

    @property
    @ureg.wraps(CP_const_unit['isothermal_compressibility'][1], None)
    def isothermal_compressibility(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][19]
        # else:
        return self._AbstractState.isothermal_compressibility()

    @property
    @ureg.wraps(CP_const_unit['isobaric_expansion_coefficient'][1], None)
    def isobaric_expansion_coefficient(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][17]/self._heprop[0][2]
        # else:
        return self._AbstractState.isobaric_expansion_coefficient()

    @property
    @ureg.wraps(CP_const_unit['isentropic_expansion_coefficient'][1], None) ### not existing with coolprop issues
    def isentropic_expansion_coefficient(self):
        return self._AbstractState.isentropic_expansion_coefficient()

    @property
    def Prandtl(self):
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][27]
        # else:
        return self._AbstractState.Prandtl()

    @property
    @ureg.wraps(CP_const_unit['speed_sound'][1], None)
    def speed_sound(self):
        """Mass fractions of a mixture."""
        # if self.name == 'Helium' and self._heprop is not None:
        #     return self._heprop[0][20]
        # else:
        return self._AbstractState.speed_sound()

    def first_partial_deriv(self, Of, Wrt, Constant):
        output_unit = CP_const_unit[Of][1] / CP_const_unit[Wrt][1]
        Of_CP_const = CP_const_unit[Of][0]
        Wrt_CP_const = CP_const_unit[Wrt][0]
        Constant_CP_const = CP_const_unit[Constant][0]
        result = self._AbstractState.first_partial_deriv(
            Of_CP_const, Wrt_CP_const, Constant_CP_const)
        return result * output_unit

    @property
    @ureg.wraps(CP_const_unit['Phase'][1], None)
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
        if self.phase in [1, 2, 3, 4]:
            return True
        elif self.phase in [0, 4, 5, 6]:
            return False
        else:
            raise ValueError('Phase is unknown')

    @property
    def specific_heat_input(self):
        """
        Calculate Specific heat input, v * (dh/dv)|p.
        This function is not described in AbstractState class.
        """
        # Because specific volume is not available as function of
        # the AbstractState, density is used instead
        # The resulting function is: -Dmass*(dHmass/dDmass)|p
        return (-self.Dmass) * self.first_partial_deriv('Hmass', 'Dmass', 'P')

    @property
    def gamma(self):
        """Calculate gamma = k = Cp/Cv coefficient.
        """
        # To avoid real gas effects influencing Cp and Cv,
        # calculating at gamma at NTP
        T_current = self.T
        S_current = self.Smass
        self.update('P', P_NTP, 'T', T_NTP)
        _gamma = self.Cpmass / self.Cvmass
        self.update('T', T_current, 'Smass', S_current)
        return _gamma

    @property
    def C_gas_const(self):
        """
        Constant for gas or vapor which is the function of the ratio of
        specific heats k = Cp/Cv. ASME VIII.1-2015 pp. 423-424.
        """
        k_ = self.gamma.magnitude
        root = (k_ * (2/(k_+1))**((k_+1)/(k_-1)))**0.5
        R = ureg.R  # Universal gas constant
        C_ = root / R**0.5 * (ureg.g/ureg.mol)**0.5
        return C_

    @property
    def C(self):
        """A shortcut to C_gas_const"""
        return self.C_gas_const

    @property
    def C_us(self):
        """A shortcut to C_gas_const for US customary units"""
        us_unit = ureg.lb/(ureg.hr*ureg.lbf)*(ureg.degR)**0.5
        return self.C_gas_const.to(us_unit)

    @property
    def C_si(self):
        """A shortcut to C_gas_const for SI units"""
        si_unit = ureg.s**2/(ureg.hr*ureg.m)*(ureg.K)**0.5
        return self.C_gas_const.to(si_unit)

    @property
    def name(self):
        """
        Return fluid name (backend dependent)
        """
        return self._AbstractState.name()

    @property
    def backend(self):
        """
        Return backend name
        """
        return self._AbstractState.backend_name()

    def __str__(self):
        if self.Q == 0:
            return (f'saturated {self.name.lower()} liquid '
                    f'at P: {self.P.to(ureg.bar):.3g~}.')
        elif self.Q == 1:
            return (f'saturated {self.name.lower()} vapor '
                    f'at P: {self.P.to(ureg.bar):.3g~}.')
        elif self.Q > 0 and self.Q < 1:
            return (f'two-phase {self.name.lower()} '
                    f'at P: {self.P.to(ureg.bar):.3g~}.')
        else:
            return (f'{self.name.lower()} at '
                    f'T: {self.T.to(ureg.K):.3g~} and '
                    f'P: {self.P.to(ureg.bar):.3g~}.')

    @property
    def M(self):
        """Calculate relative molecular mass."""
        return self.molar_mass.m_as(ureg.g/ureg.mole)

    @property
    def MZT(self):
        """
        Calculate sqrt(M/(ZT)) a commonly used square root group for discharge
        flow calculation.
        """
        MZT_ = (self.M / (self.compressibility_factor*self.T))**0.5
        return MZT_

    @property
    def latent_heat(self):
        """Calculate latent heat of evaporation for current quality."""
        assert self.is_super_critical is False, (
            'Latent heat is only defined '
            'for subcritical phase')
        TempState = self.copy()
        TempState.update_kw(P=self.P, Q=0)
        h_liq = TempState.Hmass
        TempState = self.copy()
        TempState.update_kw(P=self.P, Q=1)
        h_gas = TempState.Hmass
        return h_gas - h_liq

    def copy(self):
        """Create a copy of current ThermState object."""
        TempState = ThermState(self.name, backend=self.backend)
        # If conditions are defined
        if self.Dmass != -inf*ureg.kg/ureg.m**3:
            TempState.update_kw(T=self.T, Smass=self.Smass)
        return TempState

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
