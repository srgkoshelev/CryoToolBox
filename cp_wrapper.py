# Pint wrapper for CoolProp AbstractState
# Abstract state is used as a data point of a process flow for which properties are calculated
import CoolProp.CoolProp as CP
from . import ureg

CP_const_unit = {
    'igas_constant': (CP.igas_constant, ureg.J/ureg.mol/ureg.K),
    'imolar_mass': (CP.imolar_mass, ureg.kg/ureg.mol),
    #'iacentric_factor': (CP.iacentric_factor, ),
    #'irhomolar_reducing': (CP._reducing, ),
    'irhomolar_critical': (CP.irhomolar_critical, ureg.mol/ureg.m**3),
    #'iT_reducing': (CP.iT_reducing, ), 
    'iT_critical': (CP.iT_critical, ureg.K),
    #'irhomass_reducing': (CP.irhomass_reducing, ),
    'irhomass_critical': (CP.irhomass_critical, ureg.kg/ureg.m**3),
    'iP_critical': (CP.iP_critical, ureg.Pa),
    #'iP_reducing': (CP.iP_reducing, ),
    #'iT_triple': (CP.iT_triple, ),
    #'iP_triple': (CP.iP_triple, ),
    #'iT_min': (CP.iT_min, ),
    #'iT_max': (CP.iT_max, ),
    #'iP_max': (CP.iP_max, ),
    #'iP_min': (CP.iP_min, ),
    #'idipole_moment': (CP.idipole_moment, ),
    'iT': (CP.iT, ureg.K),
    'iP': (CP.iP, ureg.Pa),
    'iQ': (CP.iQ, ureg.mol/ureg.mol),
    #'iTau': (CP.iTau, ),
    #'iDelta': (CP.iDelta, ),
    'iDmolar': (CP.iDmolar, ureg.mol/ureg.m**3),
    'iHmolar': (CP.iHmolar, ureg.J/ureg.mol),
    'iSmolar': (CP.iSmolar, ureg.J/ureg.mol/ureg.K),
    'iCpmolar': (CP.iCpmolar, ureg.J/ureg.mol/ureg.K),
    #'iCp0molar': (CP.iCp0molar, ureg.J/ureg.mol/ureg.K),
    'iCvmolar': (CP.iCvmolar, ureg.J/ureg.mol/ureg.K),
    'iUmolar': (CP.iUmolar, ureg.J/ureg.mol),
    #'iGmolar': (CP.iGmolar, ),
    #'iHelmholtzmolar': (CP.iHelmholtzmolar, ),
    #'iSmolar_residual': (CP.iSmolar_residual, ),
    'iDmass': (CP.iDmass, ureg.kg/ureg.m**3),
    'iHmass': (CP.iHmass, ureg.J/ureg.kg),
    'iSmass': (CP.iSmass, ureg.J/ureg.kg/ureg.K),
    'iCpmass': (CP.iCpmass, ureg.J/ureg.kg/ureg.K),
    #'iCp0mass': (CP.iCp0mass, ureg.J/ureg.kg/ureg.K),
    'iCvmass': (CP.iCvmass, ureg.J/ureg.kg/ureg.K),
    'iUmass': (CP.iUmass, ureg.J/ureg.kg),
    #'iGmass': (CP.iGmass, ),
    #'iHelmholtzmass': (CP.iHelmholtzmass, ),
    'iviscosity': (CP.iviscosity, ureg.Pa*ureg.s),
    'iconductivity': (CP.iconductivity, ureg.W/ureg.m/ureg.K),
    #'isurface_tension': (CP.isurface_tension, ),
    'iPrandtl': (CP.iPrandtl, ureg.dimensionless),
    #'ispeed_sound': (CP.ispeed_sound, ),
    #'iisothermal_compressibility': (CP.iisothermal_compressibility, ),
    #'iisobaric_expansion_coefficient': (CP.iisobaric_expansion_coefficient, ),
    #'iisentropic_expansion_coefficient': (CP.iisentropic_expansion_coefficient, ),
    #'ifundamental_derivative_of_gas_dynamics': (CP.ifundamental_derivative_of_gas_dynamics, ),
    #'ialphar': (CP.ialphar, ),
    #'ialpha0': (CP.ialpha0, ),
    #'iBvirial': (CP.iBvirial, ),
    #'iCvirial': (CP.iCvirial, ),
    #'idBvirial_dT': (CP.idBvirial_dT, ),
    #'idCvirial_dT': (CP.idCvirial_dT, ),
    'iZ': (CP.iZ, ureg.dimensionless),
    #'iPIP': (CP.iPIP, ),
    #'ifraction_min': (CP.ifraction_min, ),
    #'ifraction_max': (CP.ifraction_max, ),
    #'iT_freeze': (CP.iT_freeze, ),
    #'iGWP20': (CP.iGWP20, ),
    #'iGWP100': (CP.iGWP100, ),
    #'iGWP500': (CP.iGWP500, ),
    #'iFH': (CP.iFH, ),
    #'iHH': (CP.iHH, ),
    #'iPH': (CP.iPH, ),
    #'iODP': (CP.iODP, ),
    #'iPhase': (CP.iPhase, ),
}

input_pairs = {
    'QT_INPUTS': (CP.QT_INPUTS, ureg.dimensionless, ureg.K),
    'PQ_INPUTS': (CP.PQ_INPUTS, ureg.Pa, ureg.dimensionless),
    'QSmolar_INPUTS': (CP.QSmolar_INPUTS, ureg.dimensionless,  ureg.J/ureg.mol/ureg.K),
    'QSmass_INPUTS': (CP.QSmass_INPUTS, ureg.dimensionless,  ureg.J/ureg.kg/ureg.K),
    'HmolarQ_INPUTS': (CP.HmolarQ_INPUTS, ureg.J/ureg.mol, ureg.dimensionless),
    'HmassQ_INPUTS': (CP.HmassQ_INPUTS, ureg.J/ureg.kg, ureg.dimensionless),
    'DmolarQ_INPUTS': (CP.DmolarQ_INPUTS, ureg.mol/ureg.m**3, ureg.dimensionless),
    'DmassQ_INPUTS': (CP.DmassQ_INPUTS, ureg.kg/ureg.m**3, ureg.dimensionless),
    'PT_INPUTS': (CP.PT_INPUTS, ureg.Pa, ureg.K),
    'DmassT_INPUTS': (CP.DmassT_INPUTS, ureg.kg/ureg.m**3, ureg.K),
    'DmolarT_INPUTS': (CP.DmolarT_INPUTS, ureg.mol/ureg.m**3, ureg.K),
    'HmolarT_INPUTS': (CP.HmolarT_INPUTS, ureg.J/ureg.mol, ureg.K),
    'HmassT_INPUTS': (CP.HmassT_INPUTS, ureg.J/ureg.kg, ureg.K),
    'SmolarT_INPUTS': (CP.SmolarT_INPUTS, ureg.J/ureg.mol/ureg.K, ureg.K),
    'SmassT_INPUTS': (CP.SmassT_INPUTS, ureg.J/ureg.kg/ureg.K, ureg.K),
    'TUmolar_INPUTS': (CP.TUmolar_INPUTS, ureg.K,  ureg.J/ureg.mol),
    'TUmass_INPUTS': (CP.TUmass_INPUTS, ureg.K,  ureg.J/ureg.kg),
    'DmassP_INPUTS': (CP.DmassP_INPUTS, ureg.kg/ureg.m**3, ureg.Pa),
    'DmolarP_INPUTS': (CP.DmolarP_INPUTS, ureg.mol/ureg.m**3, ureg.Pa),
    'HmassP_INPUTS': (CP.HmassP_INPUTS, ureg.J/ureg.kg, ureg.Pa),
    'HmolarP_INPUTS': (CP.HmolarP_INPUTS, ureg.J/ureg.mol, ureg.Pa),
    'PSmass_INPUTS': (CP.PSmass_INPUTS, ureg.Pa,  ureg.J/ureg.kg/ureg.K),
    'PSmolar_INPUTS': (CP.PSmolar_INPUTS, ureg.Pa,  ureg.J/ureg.mol/ureg.K),
    'PUmass_INPUTS': (CP.PUmass_INPUTS, ureg.Pa,  ureg.J/ureg.kg),
    'PUmolar_INPUTS': (CP.PUmolar_INPUTS, ureg.Pa,  ureg.J/ureg.mol),
    'HmassSmass_INPUTS': (CP.HmassSmass_INPUTS, ureg.J/ureg.kg, ureg.J/ureg.kg/ureg.K),
    'HmolarSmolar_INPUTS': (CP.HmolarSmolar_INPUTS, ureg.J/ureg.mol, ureg.J/ureg.mol/ureg.K),
    'SmassUmass_INPUTS': (CP.SmassUmass_INPUTS, ureg.J/ureg.kg/ureg.K, ureg.J/ureg.kg),
    'SmolarUmolar_INPUTS': (CP.SmolarUmolar_INPUTS, ureg.J/ureg.mol/ureg.K, ureg.J/ureg.mol),
    'DmassHmass_INPUTS': (CP.DmassHmass_INPUTS, ureg.kg/ureg.m**3,  ureg.J/ureg.kg),
    'DmolarHmolar_INPUTS': (CP.DmolarHmolar_INPUTS, ureg.mol/ureg.m**3,  ureg.J/ureg.mol),
    'DmassSmass_INPUTS': (CP.DmassSmass_INPUTS, ureg.kg/ureg.m**3,  ureg.J/ureg.kg/ureg.K),
    'DmolarSmolar_INPUTS': (CP.DmolarSmolar_INPUTS, ureg.mol/ureg.m**3,  ureg.J/ureg.mol/ureg.K),
    'DmassUmass_INPUTS': (CP.DmassUmass_INPUTS, ureg.kg/ureg.m**3,  ureg.J/ureg.kg),
    'DmolarUmolar_INPUTS': (CP.DmolarUmolar_INPUTS, ureg.mol/ureg.m**3,  ureg.J/ureg.mol),
}

class State:
    def __init__(self, backend, fluid):
        self._AbstractState = CP.AbstractState(backend, fluid)

    def update(self, inputsstr, value1, value2):
        inputs = input_pairs[inputsstr]
        self._AbstractState.update(inputs[0], value1.to(inputs[1]).magnitude, value2.to(inputs[2]).magnitude)

    @property
    @ureg.wraps(CP_const_unit['iT_critical'][1], None)
    def T_critical(self):
        return self._AbstractState.T_critical()

    @property
    @ureg.wraps(CP_const_unit['iP_critical'][1], None)
    def p_critical(self):
        return self._AbstractState.p_critical()

    @property
    @ureg.wraps(CP_const_unit['irhomolar_critical'][1], None)
    def rhomolar_critical(self):
        return self._AbstractState.rhomolar_critical()

    @property
    @ureg.wraps(CP_const_unit['irhomass_critical'][1], None)
    def rhomass_critical(self):
        return self._AbstractState.rhomass_critical()

    @property
    @ureg.wraps(CP_const_unit['iT'][1], None)
    def T(self):
        return self._AbstractState.T()

    @property
    @ureg.wraps(CP_const_unit['irhomolar_critical'][1], None)
    def rhomolar(self):
        return self._AbstractState.rhomolar()

    @property
    @ureg.wraps(CP_const_unit['irhomass_critical'][1], None)
    def rhomass(self):
        return self._AbstractState.rhomass()

    @property
    @ureg.wraps(CP_const_unit['iP'][1], None)
    def p(self):
        return self._AbstractState.p()

    @property
    @ureg.wraps(CP_const_unit['iQ'][1], None)
    def Q(self):
        return self._AbstractState.Q()

    @property
    @ureg.wraps(CP_const_unit['imolar_mass'][1], None)
    def molar_mass(self):
        return self._AbstractState.molar_mass()

    @property
    @ureg.wraps(CP_const_unit['igas_constant'][1], None)
    def gas_constant(self):
        return self._AbstractState.gas_constant()

    @property
    @ureg.wraps(CP_const_unit['iZ'][1], None)
    def compressibility_factor(self):
        return self._AbstractState.compressibility_factor()

    @property
    @ureg.wraps(CP_const_unit['iHmolar'][1], None)
    def hmolar(self):
        return self._AbstractState.hmolar()

    @property
    @ureg.wraps(CP_const_unit['iHmass'][1], None)
    def hmass(self):
        return self._AbstractState.hmass()

    @property
    @ureg.wraps(CP_const_unit['iCpmolar'][1], None)
    def cpmolar(self):
        return self._AbstractState.cpmolar()

    @property
    @ureg.wraps(CP_const_unit['iCpmass'][1], None)
    def cpmass(self):
        return self._AbstractState.cpmass()

    @property
    @ureg.wraps(CP_const_unit['iCvmolar'][1], None)
    def cvmolar(self):
        return self._AbstractState.cvmolar()

    @property
    @ureg.wraps(CP_const_unit['iCvmass'][1], None)
    def cvmass(self):
        return self._AbstractState.cvmass()

    @property
    @ureg.wraps(CP_const_unit['iviscosity'][1], None)
    def viscosity(self):
        return self._AbstractState.viscosity()

    @property
    @ureg.wraps(CP_const_unit['iconductivity'][1], None)
    def conductivity(self):
        return self._AbstractState.conductivity()

    @property
    @ureg.wraps(CP_const_unit['iPrandtl'][1], None)
    def Prandtl(self):
        return self._AbstractState.Prandtl()

    def first_partial_deriv(self, Of, Wrt, Constant):
        output_unit = CP_const_unit[Of][1] / CP_const_unit[Wrt][1]
        Of_CP_const = CP_const_unit[Of][0]
        Wrt_CP_const = CP_const_unit[Wrt][0]
        Constant_CP_const = CP_const_unit[Constant][0]
        result = self._AbstractState.first_partial_deriv(Of_CP_const, Wrt_CP_const, Constant_CP_const)
        return result * output_unit
