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

## Documentation

Project documentation is published at <https://srgkoshelev.github.io/CryoToolBox/>.

## Development Notes

The package is intentionally notebook-friendly today, but the internal
structure is still evolving. A maintainer-facing review of the current layout
and the recommended `piping` refactor is available in
[STRUCTURE_AND_ROADMAP.md](/Users/skoshelev/repo/CryoToolBox/STRUCTURE_AND_ROADMAP.md).
