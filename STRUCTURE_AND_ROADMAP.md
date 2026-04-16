# CryoToolBox Structure And Roadmap

This note records the current package-structure review performed for the
`0.2.0` codebase after comparing the active branch with the remote GitHub
branches currently available through `origin` and after validating the package
in a working runtime with the full unit-test, doctest, and pdoc workflow.

## Current strengths

- The package already has clear domain groupings: thermodynamic state,
  general engineering functions, piping, geometry, and ODH analysis.
- The public API is convenient for exploratory work: `import CryoToolBox as ctb`
  exposes most day-to-day helpers immediately.
- Unit handling is a strong usability feature and is consistently present in
  the main user-facing calculations.

## Main structural issues

### 1. `piping.py` is doing too many jobs

`CryoToolBox/piping.py` currently combines:

- table loading and regex parsing,
- exceptions,
- geometric element classes,
- valve and fitting models,
- incompressible and compressible flow solvers,
- pressure-design helpers,
- stored-energy utilities, and
- a line-description parser.

This makes the file difficult to navigate, difficult to test in small pieces,
and costly to document because unrelated concerns sit in one namespace.

### 2. `functions.py` has become a catch-all module

`CryoToolBox/functions.py` currently mixes:

- flow-conversion helpers,
- relief-flow conversions,
- radiation and conduction helpers,
- transport-property dimensionless groups,
- natural-convection correlations,
- material-property curve fits, and
- stored-energy and blast-radius helpers.

The file name no longer communicates the contents well, and the mixed concerns
make both discoverability and documentation harder.

### 3. The top-level API is broader than it looks

`CryoToolBox/__init__.py` re-exports a large amount of functionality via
`from .functions import *`. That is convenient for notebooks, but it makes the
supported public API implicit and harder for collaborators to extend safely.

### 4. Relief logic is split between `functions.py`, `cga.py`, and `piping.py`

The relief-related functionality is spread across at least three modules:

- `functions.py` contains API-style relief sizing helpers such as
  `A_relief_API`, `PRV_flow`, and air-equivalent conversions.
- `cga.py` contains CGA S-1.3-specific relief logic.
- `piping.py` contains nozzle/choking helpers that are also relief-adjacent.

This split is understandable historically, but it is no longer obvious to a
user where relief functionality should live.

A minimal backward-compatible improvement is to expose a dedicated
`CryoToolBox.relief` facade now, while leaving implementations in place until a
larger reorganization is ready.

### 5. Flow-conversion helpers should live together

The following helpers are conceptually one group and should ultimately share a
module:

- `to_scfma` and `from_scfma`,
- `to_standard_flow` and `to_mass_flow`.

They are all unit-aware flow-conversion utilities rather than general
thermodynamic helpers. A dedicated `CryoToolBox.flow` facade is the least
disruptive near-term step.

Recommended future names:

- `to_scfma` -> `to_equiv_air`
- `from_scfma` -> `from_equiv_air`
- `to_standard_flow` -> `to_std_flow`
- `to_mass_flow` -> `from_std_flow`

The short legacy names should remain available through compatibility aliases.

### 6. `geometry.py` is small but still useful

The geometry module is not redundant by itself. It becomes valuable when it
remains narrowly scoped to shared, unit-aware primitives that are reused across
modules. The current issue is not that `geometry.py` exists; it is that its
role is not yet clearly bounded.

### 7. Table/data handling works, but can be made more robust

The current packaged YAML and Python-table approach is serviceable and tested,
but there are some long-term maintenance risks:

- table schemas are implicit rather than validated,
- reference metadata and provenance are not stored next to the values,
- package-resource access was filesystem-oriented before this review, and
- the ODH failure-rate tables live in executable Python instead of data files.

## Recommended package layout

The next structural refactor should preserve the current import paths while
moving implementation into a subpackage behind compatibility re-exports.

Suggested target layout:

```text
CryoToolBox/
  __init__.py
  thermo/
    __init__.py
    state.py
    reference_conditions.py
  heat_transfer/
    __init__.py
    radiation.py
    conduction.py
    convection.py
    transient.py
    materials.py
  relief/
    __init__.py
    api520.py
    cga.py
    nozzle.py
    equivalents.py
  flow/
    __init__.py
    conversions.py
  piping/
    __init__.py
    tables.py
    parsing.py
    elements.py
    flow_incompressible.py
    flow_compressible.py
    design.py
  safety/
    __init__.py
    odh.py
  geometry.py
  constants.py
```

Recommended responsibilities:

- `thermo/state.py`: `ThermState` and backend adapters.
- `thermo/reference_conditions.py`: shared standard conditions and unit setup.
- `heat_transfer/radiation.py`: radiation and baffle helpers.
- `heat_transfer/conduction.py`: steady conduction utilities.
- `heat_transfer/convection.py`: `Re`, `Pr`, `Gr`, `Ra`, `Nu_*`, and
  heat-transfer coefficients.
- `heat_transfer/transient.py`: Biot/Fourier helpers and lumped/transient
  approximations.
- `heat_transfer/materials.py`: NIST property fits and cryogenic materials data.
- `relief/api520.py`: API-style relief sizing helpers.
- `relief/cga.py`: CGA S-1.3 calculations now in `cga.py`.
- `relief/nozzle.py`: `G_nozzle`, choked-flow helpers, and discharge utilities.
- `relief/equivalents.py`: air-equivalent and standard-flow conversions.
- `flow/conversions.py`: standard-flow and air-equivalent flow conversions.
- `piping/tables.py`: YAML loading, schema validation, and dimensional lookup tables.
- `piping/parsing.py`: `LineContext`, regex patterns, `create_element`.
- `piping/elements.py`: `Tube`, `Pipe`, `CopperTube`, fittings, valves, beds.
- `piping/flow_incompressible.py`: `K_piping`, `dP_Darcy`, `dP_incomp`,
  `m_dot_incomp`.
- `piping/flow_compressible.py`: `dP_isot`, `m_dot_isot`, `dP_adiab`,
  `m_dot_adiab`, Mach helpers.
- `piping/design.py`: B31.3 wall-thickness and reinforcement calculations.
- `safety/odh.py`: ODH analysis models and reporting helpers.

Migration strategy:

1. Move implementations into the new submodules.
2. Move `functions.py` code into the new domain modules while keeping a
   compatibility facade.
3. Keep `CryoToolBox/piping.py` as a compatibility facade for one release.
4. Add tests for each extracted module before removing the facades.
5. Only then consider a public-API cleanup.

## Guidance for `geometry.py`

Keep `geometry.py`, but keep it disciplined.

Good candidates:

- circle, annulus, and cylinder primitives,
- shape factors reused by both piping and heat-transfer code,
- small dimensionally clear helpers with no domain-specific side effects.

Avoid putting these in `geometry.py`:

- pressure-drop correlations,
- piping standards logic,
- convenience wrappers that are only called once,
- parsing or table-lookup helpers.

If a helper is only used by one class and does not improve readability outside
that class, it should stay local to that module instead of being moved into
`geometry.py` just for DRY.

## Documentation findings from remote branches

Remote branches available under `origin` contain a long history of docstring
cleanup work, especially commits around the earlier NumPy-style doc migration.
Most of that historical work is already reflected in the current `develop`
branch. After comparing those branches with the validated `0.2.0` codebase,
the remaining gaps were mainly:

- undocumented public exceptions and compatibility wrappers,
- malformed or uneven docstrings in a few modules,
- missing maintainer-facing guidance on package boundaries.

Those are better addressed incrementally on the current branch than by trying
to merge an older documentation-heavy branch wholesale.

## Highest-priority usability improvements

1. Define the supported public API explicitly.
   Publish which names are stable at the package root and which are module-only.

2. Split `functions.py` into domain modules behind compatibility re-exports.
   This is the single biggest discoverability improvement after the docs pass.

3. Keep steering flow-conversion users toward `CryoToolBox.flow`.
   This creates a stable namespace before the deeper `functions.py` split.

4. Keep steering relief users toward `CryoToolBox.relief`.
   This creates a stable namespace before the deeper module split.

5. Split `piping.py` behind a compatibility layer.
   This offers the biggest maintainability improvement without forcing user
   code to change immediately.

6. Add focused examples for the three main workflows.
   These should cover thermodynamic state setup, piping pressure-drop analysis,
   and ODH/source modeling.

7. Expand tests around parsing and branch conditions.
   The line-description parsing and compressible-flow edge cases are the parts
   most likely to surprise users.

8. Document backend expectations clearly.
   `HEOS` should be the default documented path, while `REFPROP` and `HEPROP`
   should be marked as optional integrations with setup caveats.

9. Introduce lightweight validation around packaged tables.
   Required keys, units, and monotonic expectations should be tested so table
   updates cannot silently break calculations.

## Table recommendations

Near-term best practices:

1. Keep NPS and copper dimensions in data files, not code.
2. Add table schema checks in tests:
   required keys like `OD`, positive wall thicknesses, and valid schedules/types.
3. Store provenance in comments or adjacent docs:
   standard edition, source table, and any intentional omissions.
4. Prefer package-resource loading over direct filesystem assumptions.
   This review already updated the package-data access in code.
5. Keep all lookup-table loading and validation in one place.
   A future `piping.tables` module would simplify imports, testing, and schema
   evolution without changing the user-facing API.

Longer-term options:

1. Convert ODH failure-rate tables into YAML or JSON data files with explicit
   metadata, then load and unitize them centrally.
2. Introduce typed accessors or dataclasses for frequently used table rows.
3. Consider normalizing nominal-size keys to strings in raw data to avoid YAML
   float-key surprises across tools, then convert explicitly in the loader.

## Suggested near-term plan

For the next development cycle, the safest order is:

1. Freeze and document the current public API.
2. Add stable facades such as `flow` and `relief` before moving implementations.
3. Extract `functions.py` into domain modules behind compatibility re-exports.
4. Extract `piping` internals into submodules with no user-facing breakage.
5. Split the test suite by domain once the facades settle.
6. Add examples and tests around the extracted boundaries.
7. Revisit deeper API cleanup only after the modular split is stable.

## Test-plan notes

Current touched areas and coverage:

- `relief` facade: lightly covered by alias/export tests.
- flow conversions: only partially covered today.
- table loading: covered by import-time validation and integrity tests.

Easy tests that do not require new domain-specific reference data:

1. facade alias tests for `flow` and `relief`,
2. round-trip tests for inverse conversion helpers,
3. dimensionality/error-path tests for public conversion helpers,
4. packaging/import tests for data-backed modules.

Reasonable future split of the test suite:

- `tests/test_cp_wrapper.py`
- `tests/test_flow.py`
- `tests/test_relief.py`
- `tests/test_piping.py`
- `tests/test_odh.py`
- `tests/test_geometry.py`
