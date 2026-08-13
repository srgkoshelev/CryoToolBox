"""Benchmark CryoToolBox's CoolProp integration paths.

Run from the repository root, for example:

    .venv313/bin/python benchmarks/coolprop_integration.py

The "legacy" rows use preserved CTB helper formulas that mirror the previous
implementation, allowing before/after comparisons from the current checkout.
"""

from __future__ import annotations

import platform
import statistics
import timeit
import warnings
from math import inf

import CoolProp
import CoolProp.CoolProp as CP

import CryoToolBox as ctb
from CryoToolBox.cp_wrapper import CryoToolBoxFutureWarning


warnings.filterwarnings('ignore', category=CryoToolBoxFutureWarning)


def _legacy_compressibility_factor(state):
    if state.backend == 'HEPROP':
        return state._AbstractState.compressibility_factor()
    return (state.P * state.molar_mass /
            (state.Dmass * state.gas_constant * state.T)).m_as(
                ctb.ureg.dimensionless)


def _legacy_copy(state):
    temp_state = ctb.ThermState(state.name, backend=state.backend)
    if state.Dmass != -inf * ctb.ureg.kg / ctb.ureg.m**3:
        temp_state.update_kw(T=state.T, Smass=state.Smass)
    return temp_state


def _legacy_latent_heat(state):
    if not (0 <= state.Q <= 1):
        raise ValueError(f'Latent heat is undefined for Q={state.Q}, {state}.')
    temp_state = _legacy_copy(state)
    temp_state.update_kw(P=state.P, Q=1)
    return temp_state.Hmass - state.Hmass


def _legacy_P_sat(state):
    temp_state = _legacy_copy(state)
    temp_state.update_kw(T=state.T, Q=0)
    return temp_state.P


def _legacy_T_sat(state):
    temp_state = _legacy_copy(state)
    temp_state.update_kw(P=state.P, Q=0)
    return temp_state.T


def _legacy_gamma(state):
    temp_state = _legacy_copy(state)
    temp_state.update_kw(P=ctb.P_NTP, T=ctb.T_NTP)
    return temp_state.Cpmass / temp_state.Cvmass


def best_time_us(stmt, number, repeat=5):
    """Return the best per-call time in microseconds."""
    timings = timeit.repeat(stmt, number=number, repeat=repeat)
    return min(timings) / number * 1e6


def row(label, before_us=None, after_us=None, note=''):
    """Print one benchmark row."""
    before = '' if before_us is None else f'{before_us:12.3f}'
    after = '' if after_us is None else f'{after_us:12.3f}'
    if before_us and after_us:
        speedup = f'{before_us / after_us:9.2f}x'
    else:
        speedup = ''
    print(f'{label:42s} {before:>12s} {after:>12s} {speedup:>10s} {note}')


def main():
    u = ctb.ureg
    state = ctb.ThermState('nitrogen', P=500 * u.kPa, T=300 * u.K)
    sat = ctb.ThermState('nitrogen', P=100 * u.kPa, Q=0.25)
    saturation_probe = ctb.ThermState('nitrogen', P=200 * u.kPa, T=90 * u.K)
    raw = CP.AbstractState('HEOS', 'Nitrogen')
    raw.update(CP.PT_INPUTS, 500_000.0, 300.0)

    print('CryoToolBox CoolProp integration benchmark')
    print(f'Python:   {platform.python_version()} ({platform.machine()})')
    print(f'CoolProp: {CoolProp.__version__}')
    print('')
    print(f'{"operation":42s} {"legacy us":>12s} {"current us":>12s} '
          f'{"speedup":>10s} note')
    print('-' * 92)

    row(
        'AbstractState.update(PT), raw CP',
        after_us=best_time_us(
            lambda: raw.update(CP.PT_INPUTS, 500_000.0, 300.0), 50_000),
    )
    row(
        'ThermState.update_kw(P,T)',
        after_us=best_time_us(
            lambda: state.update_kw(P=500 * u.kPa, T=300 * u.K), 20_000),
        note='Pint-aware inputs',
    )
    row(
        'ThermState.update_si(P,T)',
        after_us=best_time_us(
            lambda: state.update_si('P', 500_000.0, 'T', 300.0), 50_000),
        note='raw CoolProp units',
    )
    row(
        'ThermState.get_si(Dmass)',
        after_us=best_time_us(lambda: state.get_si('Dmass'), 100_000),
    )
    row(
        'P,T,Dmass,Hmass properties',
        after_us=best_time_us(
            lambda: (state.P, state.T, state.Dmass, state.Hmass), 20_000),
    )
    row(
        'P,T,Dmass,Hmass get_si',
        after_us=best_time_us(
            lambda: (state.get_si('P'), state.get_si('T'),
                     state.get_si('Dmass'), state.get_si('Hmass')),
            50_000,
        ),
    )
    row(
        'compressibility_factor',
        before_us=best_time_us(
            lambda: _legacy_compressibility_factor(state), 20_000),
        after_us=best_time_us(lambda: state.compressibility_factor, 100_000),
    )
    row(
        'latent_heat',
        before_us=best_time_us(lambda: _legacy_latent_heat(sat), 5_000),
        after_us=best_time_us(lambda: sat.latent_heat, 20_000),
    )
    row(
        'P_sat',
        before_us=best_time_us(lambda: _legacy_P_sat(saturation_probe), 5_000),
        after_us=best_time_us(lambda: saturation_probe.P_sat, 20_000),
        note='native ancillary',
    )
    row(
        'T_sat',
        before_us=best_time_us(lambda: _legacy_T_sat(saturation_probe), 5_000),
        after_us=best_time_us(lambda: saturation_probe.T_sat, 20_000),
        note='native ancillary',
    )
    row(
        'copy(), single-phase',
        before_us=best_time_us(lambda: _legacy_copy(state), 2_000),
        after_us=best_time_us(lambda: state.copy(), 10_000),
    )
    row(
        'gamma current-state',
        before_us=best_time_us(lambda: _legacy_gamma(state), 2_000),
        after_us=best_time_us(lambda: state.gamma, 20_000),
        note='behavior changed',
    )

    values = [
        _legacy_compressibility_factor(state),
        state.compressibility_factor,
        (_legacy_latent_heat(sat) - sat.latent_heat).m_as(u.J / u.kg),
        ((_legacy_P_sat(saturation_probe) - saturation_probe.P_sat) /
         saturation_probe.P_sat).m_as(u.dimensionless),
        ((_legacy_T_sat(saturation_probe) - saturation_probe.T_sat) /
         saturation_probe.T_sat).m_as(u.dimensionless),
    ]
    print('')
    print('Sanity checks:')
    print(f'  legacy/current Z: {values[0]:.12g} / {values[1]:.12g}')
    print(f'  latent_heat delta: {values[2]:.6g} J/kg')
    print(f'  P_sat relative delta: {values[3]:.6g}')
    print(f'  T_sat relative delta: {values[4]:.6g}')
    copy_timings = timeit.repeat(lambda: state.copy(), number=1000, repeat=5)
    copy_median_us = statistics.median(copy_timings) / 1000 * 1e6
    print(f'  timing spread note: median current copy timing is '
          f'{copy_median_us:.3f} us')


if __name__ == '__main__':
    main()
