# CryoToolBox

`CryoToolBox` is a Python package for cryogenic fluid, heat-transfer, and piping calculations built on top of [CoolProp](https://github.com/CoolProp/CoolProp) and [Pint](https://github.com/hgrecco/pint).

## Status

Version `0.2.0` is an early alpha. The API is still evolving and optional backends such as `REFPROP` and `HEPROP` may not be available in every environment.

## Installation

Create a virtual environment and install the package:

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install .
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

## Verification

Run the current test suite with:

```sh
.venv/bin/python -m unittest -q
```

Build release artifacts with:

```sh
.venv/bin/python -m build --no-isolation
```

## Alpha Scope

Supported in alpha:

- `ThermState` with the default `HEOS` backend
- Unit-aware thermodynamic helpers
- Geometry helpers
- Piping and hydraulics calculations

Optional or experimental in alpha:

- `REFPROP`
- `HEPROP`
- ODH Excel export, which requires `CryoToolBox[odh]`

## Documentation

Project documentation is published at <https://srgkoshelev.github.io/CryoToolBox/>.
