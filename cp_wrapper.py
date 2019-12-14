import CoolProp.CoolProp as CP
from . import ureg, T_NTP, P_NTP

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
    # 'surface_tension': (CP.isurface_tension, ),
    'Prandtl': (CP.iPrandtl, ureg.dimensionless),
    # 'speed_sound': (CP.ispeed_sound, ),
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
    # 'Phase': (CP.iPhase, ),
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
        self._AbstractState = CP.AbstractState(backend, fluid)
        if state_parameters:
            self.update_kw(**state_parameters)

    def update_kw(self, **state_parameters):
        """Update thermodynamic state using keyword arguments."""
        if len(state_parameters) != 2:
            raise TypeError(f'update_kw() takes 2 arguments \
            ({len(state_parameters)} given)')
        input_names = []
        input_values = []
        for name, value in state_parameters.items():
            input_names.append(name)
            input_values.append(value)
        self.update(input_names[0], input_values[0],
                    input_names[1], input_values[1])

    def update(self, input_name1, input_value1, input_name2, input_value2):
        if input_name1 < input_name2:  # Sorting inputs alphabetically
            CP_input_str = input_name1 + input_name2
            CP_input_value1 = input_value1.to(
                CP_const_unit[input_name1][1]).magnitude
            CP_input_value2 = input_value2.to(
                CP_const_unit[input_name2][1]).magnitude
        else:
            CP_input_str = input_name2 + input_name1
            CP_input_value1 = input_value2.to(
                CP_const_unit[input_name2][1]).magnitude
            CP_input_value2 = input_value1.to(
                CP_const_unit[input_name1][1]).magnitude
        self._AbstractState.update(
            CP_inputs[CP_input_str], CP_input_value1, CP_input_value2)

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
        return self._AbstractState.T()

    @property
    @ureg.wraps(CP_const_unit['Dmolar_critical'][1], None)
    def Dmolar(self):
        return self._AbstractState.rhomolar()

    @property
    @ureg.wraps(CP_const_unit['Dmass_critical'][1], None)
    def Dmass(self):
        return self._AbstractState.rhomass()

    @property
    @ureg.wraps(CP_const_unit['P'][1], None)
    def P(self):
        return self._AbstractState.p()

    @property
    @ureg.wraps(CP_const_unit['Q'][1], None)
    def Q(self):
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
    # @ureg.wraps(CP_const_unit['Z'][1], None)
    def compressibility_factor(self):
        Z_ = self.P * self.molar_mass / (self.Dmass*self.gas_constant*self.T)
        return Z_.to(ureg.dimensionless)
        # Temporarily unavailable function
        # return self._AbstractState.compressibility_factor()

    @property
    @ureg.wraps(CP_const_unit['Hmolar'][1], None)
    def Hmolar(self):
        return self._AbstractState.hmolar()

    @property
    @ureg.wraps(CP_const_unit['Hmass'][1], None)
    def Hmass(self):
        return self._AbstractState.hmass()

    @property
    @ureg.wraps(CP_const_unit['Smolar'][1], None)
    def Smolar(self):
        return self._AbstractState.smolar()

    @property
    @ureg.wraps(CP_const_unit['Smass'][1], None)
    def Smass(self):
        return self._AbstractState.smass()

    @property
    @ureg.wraps(CP_const_unit['Cpmolar'][1], None)
    def Cpmolar(self):
        return self._AbstractState.cpmolar()

    @property
    @ureg.wraps(CP_const_unit['Cpmass'][1], None)
    def Cpmass(self):
        return self._AbstractState.cpmass()

    @property
    @ureg.wraps(CP_const_unit['Cvmolar'][1], None)
    def Cvmolar(self):
        return self._AbstractState.cvmolar()

    @property
    @ureg.wraps(CP_const_unit['Cvmass'][1], None)
    def Cvmass(self):
        return self._AbstractState.cvmass()

    @property
    @ureg.wraps(CP_const_unit['viscosity'][1], None)
    def viscosity(self):
        return self._AbstractState.viscosity()

    @property
    @ureg.wraps(CP_const_unit['conductivity'][1], None)
    def conductivity(self):
        return self._AbstractState.conductivity()

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
        return self._AbstractState.isothermal_compressibility()

    @property
    @ureg.wraps(CP_const_unit['isobaric_expansion_coefficient'][1], None)
    def isobaric_expansion_coefficient(self):
        return self._AbstractState.isobaric_expansion_coefficient()

    @property
    @ureg.wraps(CP_const_unit['isentropic_expansion_coefficient'][1], None)
    def isentropic_expansion_coefficient(self):
        return self._AbstractState.isentropic_expansion_coefficient()

    @property
    @ureg.wraps(CP_const_unit['Prandtl'][1], None)
    def Prandtl(self):
        return self._AbstractState.Prandtl()

    def first_partial_deriv(self, Of, Wrt, Constant):
        output_unit = CP_const_unit[Of][1] / CP_const_unit[Wrt][1]
        Of_CP_const = CP_const_unit[Of][0]
        Wrt_CP_const = CP_const_unit[Wrt][0]
        Constant_CP_const = CP_const_unit[Constant][0]
        result = self._AbstractState.first_partial_deriv(
            Of_CP_const, Wrt_CP_const, Constant_CP_const)
        return result * output_unit

    @property
    def specific_heat_input(self):
        """
        Calculate Specific heat input, v * (dh/dv)|p.
        This function is not described in AbstractState class.
        """
        # Because specific volume is not available as function of
        # the AbstractState, density is used instead
        # The resulting function is: -Dmass*(dHmass/Dmass)|p
        return (-self.Dmass) * self.first_partial_deriv('Hmass', 'Dmass', 'P')

    @property
    def gamma(self):
        """
        Calculate gamma (k) coefficient.

        :returns: gamma = k = Cp/Cv
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
    def C_gas_constant(self):
        """
        Constant for gas or vapor which is the function of the ratio of
        specific heats k = Cp/Cv. ASME VIII.1-2015 pp. 423-424.
        """
        k_ = self.gamma
        C = 520*(k_*(2/(k_+1))**((k_+1)/(k_-1)))**0.5*ureg(
            'lb/(hr*lbf)*(mol*degR/g)^0.5')
        return C

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
        return f'{self.name.capitalize()} at T = {self.T.to(ureg.K):.3~g~} and \
        P = {self.P.to(ureg.psi):.3g~}.'

    @property
    def MZT(self):
        """
        Calculate sqrt(M/(ZT)) a commonly used square root group for discharge
        flow calculation.
        """
        MZT_ = (self.molar_mass / (self.compressibility_factor*self.T))**0.5
        return MZT_
