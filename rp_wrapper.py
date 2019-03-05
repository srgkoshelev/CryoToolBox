#' % Utility functions for Refprop python package
#' % by Sergey Koshelev

#' Importing math functions, unit handling package, and logging. Refprop is used for material properties.
from math import log, log10
from pyrefprop import refprop as rp
from pint import UnitRegistry
import logging
import logging.config
from functools import wraps
#' Setting up logging configuration file:
import sys, os
if __name__ != '__main__':
    __location__ = os.path.dirname(os.path.abspath(__file__))
    # load the logging configuration:
    logging.config.fileConfig(os.path.join(__location__, 'logging.ini'))
else:
    # load the logging configuration:
    logging.config.fileConfig(os.path.join('./', 'logging.ini'))
logger = logging.getLogger(__name__)

#' Configuring units package:
ureg = UnitRegistry(autoconvert_offset_to_baseunit = True)
Q_ = ureg.Quantity
if __name__ != '__main__':
    ureg.load_definitions(os.path.join(__location__, 'pint definitions.txt'))
else:
    ureg.load_definitions(os.path.join('./', 'pint definitions.txt'))
sigma = Q_('stefan_boltzmann_constant')

#Setting units for "standard" flow
T_NTP = Q_(68, ureg.degF) #Normal Temperature (NIST)
P_NTP = Q_(14.7, ureg.psi) #Normal Pressure (NIST)

T_MSC = Q_(15, ureg.degC) #Metric Standard Conditions (used by Crane TP-410)
P_MSC = Q_(101325, ureg.Pa) #Metric Standard Conditions (used by Crane TP-410)

Air = {'fluid':'air', 'P':P_NTP, 'T':T_NTP}

#' Refprop package uses dimensionless values for its functions. A series of wrapper functions are used to allow for units seamless integration.
def flsh(routine, var1, var2, x, kph=1):
    """
    Multifunctional wrapper for refprop flsh function
    """
    if routine == 'TP':   #temperature; pressure
        var1_unitless = var1.to(ureg.K).magnitude
        var2_unitless = var2.to(ureg.kPa).magnitude
    if routine == 'TD':   #temperature; Molar Density
        var1_unitless = var1.to(ureg.K).magnitude
        var2_unitless = var2.to(ureg.mol/ureg.L).magnitude
    if routine == 'TH':   #temperature; enthalpy
        var1_unitless = var1.to(ureg.K).magnitude
        var2_unitless = var2.to(ureg.J/ureg.mol).magnitude
    if routine == 'TS':   #temperature; entropy
        var1_unitless = var1.to(ureg.K).magnitude
        var2_unitless = var2.to(ureg.J/(ureg.mol*ureg.K)).magnitude
    if routine == 'TE':   #temperature; internal energy
        var1_unitless = var1.to(ureg.K).magnitude
        var2_unitless = var2.to(ureg.J/ureg.mol).magnitude
    if routine == 'PD':   #pressure; molar density
        var1_unitless = var1.to(ureg.kPa).magnitude
        var2_unitless = var2.to(ureg.mol/ureg.L).magnitude
    if routine == 'PH':   #pressure; enthalpy
        var1_unitless = var1.to(ureg.kPa).magnitude
        var2_unitless = var2.to(ureg.J/ureg.mol).magnitude
    if routine == 'PS':   #pressure; entropy
        var1_unitless = var1.to(ureg.kPa).magnitude
        var2_unitless = var2.to(ureg.J/(ureg.mol*ureg.K)).magnitude
    if routine == 'PE':   #pressure; internal energy
        var1_unitless = var1.to(ureg.kPa).magnitude
        var2_unitless = var2.to(ureg.J/ureg.mol).magnitude
    if routine == 'HS':   #enthalpy; entropy
        var1_unitless = var1.to(ureg.J/ureg.mol).magnitude
        var2_unitless = var2.to(ureg.J/(ureg.mol*ureg.K)).magnitude
    if routine == 'ES':   #internal energy; entropy
        var1_unitless = var1.to(ureg.J/ureg.mol).magnitude
        var2_unitless = var2.to(ureg.J/(ureg.mol*ureg.K)).magnitude
    if routine == 'DH':   #molar density; enthalpy
        var1_unitless = var1.to(ureg.mol/ureg.L).magnitude
        var2_unitless = var2.to(ureg.J/ureg.mol).magnitude
    if routine == 'DS':   #molar density; entropy
        var1_unitless = var1.to(ureg.mol/ureg.L).magnitude
        var2_unitless = var2.to(ureg.J/(ureg.mol*ureg.K)).magnitude
    if routine == 'DE':   #molar density; internal energy
        var1_unitless = var1.to(ureg.mol/ureg.L).magnitude
        var2_unitless = var2.to(ureg.J/ureg.mol).magnitude
    if routine == 'TQ':   #temperature; vapor quality
        var1_unitless = var1.to(ureg.K).magnitude
        var2_unitless = var2
    if routine == 'PQ':   #pressure; vapor quality
        var1_unitless = var1.to(ureg.kPa).magnitude
        var2_unitless = var2

    refprop_output = rp.flsh(routine, var1_unitless, var2_unitless, x, kph=1)
    return rp_out_unit(refprop_output)



#Dictionary defining units used by refprop package
Outputs = {'t':ureg.K, 'p':ureg.kPa, 'D':ureg.mol/ureg.L,
           'Dliq':ureg.mol/ureg.L, 'Dvap':ureg.mol/ureg.L, 'x':None,
           'xliq':None, 'xvap':None, 'q':None, 'e':ureg.J/ureg.mol,
           'h':ureg.J/ureg.mol, 's':ureg.J/ureg.mol*ureg.K,
           'cv':ureg.J/(ureg.mol*ureg.K), 'cp':ureg.J/(ureg.mol*ureg.K),
           'w':ureg.m/ureg.s, 'hfmix':None, 'kph':None, 'hrf':None,
           'hfld':None, 'nc':None, 'xkappa':None, 'beta':None, 'xisenk':None,
           'xkt':None, 'betas':None, 'bs':None, 'xkkt':None, 'thrott':None,
           'pint':None, 'spht':ureg.J/ureg.mol, 'wmm':ureg.g/ureg.mol,
           'ttrp':ureg.K, 'tnbpt':ureg.K, 'tcrit':ureg.K, 'pcrit':ureg.kPa,
           'Dcrit':ureg.mol/ureg.L, 'zcrit':None, 'acf':None, 'dip':ureg.debye,
           'Rgas':ureg.J/(ureg.mol*ureg.K), 'icomp':None,
           'dpdD':ureg.kPa*ureg.L/ureg.mol,
           'd2pdD2':ureg.kPa*ureg.L**2/ureg.mol**2, 'dpdt':ureg.kPa/ureg.K,
           'dDdt':ureg.mol/(ureg.L*ureg.K), 'dDdp':ureg.mol/(ureg.L*ureg.kPa),
           'Z':None, 'hjt':ureg.K/ureg.kPa, 'A':ureg.J/ureg.mol,
           'G':ureg.J/ureg.mol, 'hmxnme':None,
           }

def rp_value(value, name):
    """
    Prepare the quantity for passing to rp function by obtaining a
    dimensionless value.
    """
    if Outputs[name]:
        return value.to(Outputs[name]).magnitude
    else:
        return value

def rp_out_unit (rp_output):
    """
    Add units to output of a refprop function.
    """
    Output_w_units = {}
    for quantity, value in rp_output.items():
        try:
            unit = Outputs[quantity]
            if unit:
                Output_w_units[quantity] = value*unit
            else:
                Output_w_units[quantity] = value
        except KeyError:
            logger.warning('''Quantity is missing from unit list, please update:
                           {}'''.format(quantity))
    return Output_w_units

def rp_unitize(*names):
    """
    Decorator for refprop functions.
    Given the string names of the variables convert input to rp accepted
    dimensionless values.
    Results are parsed to assign units.
    Use:
    @rp_unitize('p', 'x', 'kph')
    def satp(p, x, kph=2):
        return rp.satp(p, x, kph)
    """
    def rp_decorate(rp_func):
        def rp_wrapper(*rp_args, **rp_kwargs):
            rp_input_values = list(rp_args) + list(rp_kwargs.values())
            rp_input = []
            for value,name in zip(rp_input_values,names):
                rp_input.append(rp_value(value,name))
            return rp_out_unit(rp_func(*rp_input))
        return rp_wrapper
    return rp_decorate

@rp_unitize('p', 'x', 'kph')
def satp(p, x, kph=2):
    return rp.satp(p, x, kph)

@rp_unitize('t', 'D', 'x')
def therm2(t, D, x):
    return rp.therm2(t, D, x)

@rp_unitize('t', 'D', 'x')
def therm3(t, D, x):
    return rp.therm3(t, D, x)

@rp_unitize('x')
def critp(x):
    return rp.critp(x)

@rp_unitize('icomp')
def info(icomp=1):
    return rp.info(icomp)

trnprp = ureg.wraps (None, (ureg.K, ureg('mol/L'), None))(rp.trnprp)

#' The package uses unified way of keeping track of fluid state, specifying fluid, pressure, temperature. Helper functions to pack and unpack the fluid state.
def pack_fluid (fluid, T_fluid = Q_(15, ureg.degC), P_fluid = Q_(101325, ureg.Pa)):
        return {'fluid':fluid, 'P':P_fluid, 'T':T_fluid}

def unpack_fluid (Fluid_data = {'fluid':'air', 'P':Q_(101325,ureg.Pa),
                                'T':Q_(15,ureg.degC)}):
        fluid = Fluid_data.get('fluid')
        T_fluid = Fluid_data.get('T')
        P_fluid = Fluid_data.get('P')
        return (fluid, T_fluid, P_fluid)

#' To use refprop one needs first to set up the fluid. Following function simplifies initialization and also returns most commonly required properties when possible.
def rp_init (Fluid_data):
        """
        Shortcut initialization of refprop
        Returns (x, M, D_fluid)
        """
        (fluid, T_fluid, P_fluid) = unpack_fluid(Fluid_data)
        prop = rp.setup('def', fluid)
        x = prop.get('x', [1]) 
        M = rp.wmol(x)['wmix']*ureg('g/mol')
        if T_fluid and P_fluid:
            fluid_prop = flsh ("TP", T_fluid, P_fluid, x)
            D_fluid = fluid_prop['D'] 
            return (x, M, D_fluid)
        else: 
            return (x, M)
