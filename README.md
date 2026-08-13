# CryoToolBox

`CryoToolBox` is a Python package for cryogenic fluid, heat-transfer, and piping calculations built on top of [CoolProp](https://github.com/CoolProp/CoolProp) and [Pint](https://github.com/hgrecco/pint).

## Status

Version `0.2.0` is an early alpha. The API is still evolving and optional backends such as `REFPROP` and `HEPROP` may not be available in every environment.

## Installation

Fetch the source, create a virtual environment, install the package, and run
a quick smoke test:

```sh
git clone https://github.com/srgkoshelev/CryoToolBox.git
cd CryoToolBox
python3 -m venv .venv
source .venv/bin/activate
pip install .
python -c "import CryoToolBox as ctb; nitrogen = ctb.ThermState('nitrogen', P=500 * ctb.ureg.kPa, T=300 * ctb.ureg.K); print(nitrogen.Dmass.to(ctb.ureg.kg / ctb.ureg.meter**3))"
```

Optional extras:

```sh
pip install ".[docs]"
pip install ".[odh]"
```

To install directly from GitHub:

```sh
pip install git+https://github.com/srgkoshelev/CryoToolBox.git
```

## Quickstart

```python
import CryoToolBox as ctb

nitrogen = ctb.ThermState("nitrogen", P=500 * ctb.ureg.kPa, T=300 * ctb.ureg.K)
print(nitrogen.Dmass)
print(nitrogen.Dmass.to(ctb.ureg.lb / ctb.ureg.ft**3))
```

## Development Checks

Run the current test suite during development with:

```sh
.venv/bin/python -m unittest -q
```

The full suite exercises optional and dependency-sensitive behavior, so use the
installation smoke test above for a quick check that the package installed
correctly.

Build release artifacts with:

```sh
.venv/bin/python -m build --no-isolation
```

Run the CoolProp integration benchmark with:

```sh
.venv313/bin/python benchmarks/coolprop_integration.py
```

On a local arm64 Python 3.13 / CoolProp 8.0.0 environment, the benchmark
showed the new raw-SI paths and native CoolProp wrappers reducing representative
wrapper overhead:

- `ThermState.update_si(P,T)`: about `2.9 us`, close to raw
  `AbstractState.update(PT)` at about `2.8 us`.
- `P,T,Dmass,Hmass` via `get_si`: about `0.3 us`, compared with about
  `7.37 us` through unit-wrapped properties.
- `compressibility_factor`: about `0.04 us`, compared with about `36.7 us`
  for the preserved legacy formula.
- `latent_heat`, `P_sat`, and `T_sat`: about `38x` to `56x` faster than the
  preserved legacy copy/flash calculations.
- `copy()`: about `2.3x` faster for a representative single-phase state.

## Alpha Scope

Supported in alpha:

- `ThermState` with the default `HEOS` backend
- Flow-conversion helpers via `CryoToolBox.flow`, including `to_equiv_air`,
  `from_equiv_air`, `to_std_flow`, and `from_std_flow`
- Unit-aware thermodynamic helpers
- Relief and discharge helpers via `CryoToolBox.relief`
- Geometry helpers
- Piping and hydraulics calculations

Optional or experimental in alpha:

- `REFPROP`
- `HEPROP`
- ODH Excel export, which requires `CryoToolBox[odh]`

## Wishlist

- Add vectorized property-evaluation helpers for NumPy arrays. These should
  avoid per-point Pint object creation, accept raw SI arrays where appropriate,
  and return either raw arrays or unit-wrapped arrays at the boundary.
- Evaluate CoolProp 8's `SVDSBTL`/`fast_evaluate` table backends for large
  batch workloads. This should be opt-in because table creation has startup
  cost, supported fluids/backends differ from normal `HEOS`, and accuracy needs
  cryogenic edge-case validation.

## Documentation

Project documentation is published at <https://srgkoshelev.github.io/CryoToolBox/>.

## Development Notes

The package is intentionally notebook-friendly today, but the internal
structure is still evolving. A maintainer-facing review of the current layout
and the recommended `piping` refactor is available in
[STRUCTURE_AND_ROADMAP.md](/Users/skoshelev/repo/CryoToolBox/STRUCTURE_AND_ROADMAP.md).
