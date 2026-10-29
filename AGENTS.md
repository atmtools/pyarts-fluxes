# AGENTS.md

Guidance for AI coding agents working in this repository.

## Project overview

`pyarts-fluxes` (package name `FluxSimulator`) is a thin Python wrapper around
ARTS (via `pyarts`) for computing radiative fluxes, heating rates, and
radiances for single atmospheres and batches of atmospheres. Almost all logic
lives in two files:

- `src/FluxSimulator/_flux_simulator_module.py` — the `FluxSimulator` /
  `FluxSimulationConfig` classes: public API methods (single-profile and
  batch flux/radiance simulators, LUT generation, atmosphere prep helpers).
- `src/FluxSimulator/_flux_simulator_agendas.py` — ARTS `@arts_agenda`
  functions (gas scattering, surface, PSD/pnd, and `dobatch_calc_agenda_*`
  batch calculation agendas) used by the module via `self.ws.<agenda> = ...`.

`examples/` and `sims/` contain runnable scripts demonstrating the API and are
the closest thing to integration tests/documentation of expected usage.
`atmdata/` and `scattering_data/` hold sample ARTS XML input data.

## Environment

- Requires `pyarts >= 2.6.20` and `numpy >= 2.0.0` (see `pyproject.toml`).
- The system Python (`python3`) usually does **not** have `pyarts` installed.
  Use the `pyarts` conda environment for anything that imports the package:
  `conda run -n pyarts python ...`.
- Install editable: `python -m pip install --user -e .` (from within the
  `pyarts` env).

## Verification workflow (no formal test suite)

There is no `tests/` directory or CI test job. When changing code:

1. `python3 -m py_compile <changed files>` for a fast syntax check.
2. `conda run -n pyarts python -c "..."` to import the module and construct
   any new `@arts_agenda` functions/workspace variables — this catches ARTS
   API mistakes (unknown workspace variables/methods) that plain syntax
   checks miss. Agenda construction needs its referenced workspace variables
   to already exist (created via `ws.<Type>Create(name)`), so create minimal
   stand-ins when smoke-testing an agenda in isolation.
2. Prefer editing an existing example under `examples/` or `sims/` (or adding
   a small new one) over inventing throwaway scripts, since these double as
   living usage docs.
3. A full end-to-end run requires real ARTS catalog/LUT data and is slow;
   don't assume it's feasible in every session — say so explicitly if a
   runtime comparison can't be completed.

## Code conventions

- Batch methods (`flux_simulator_batch`, `radiance_simulator_batch`) mirror
  single-profile methods (`flux_simulator_single_profile`,
  `radiance_simulator_single_profile`) but drive ARTS via `DOBatchCalc` with
  per-profile data packed into `array_of_*` / `matrix_of_*` workspace
  variables and extracted per-index inside `dobatch_calc_agenda_*` functions
  using `ws.Extract(...)`.
- Results are returned as plain `dict`s; batch results use `array_of_<name>`
  keys holding one list entry per profile, single-profile results use bare
  `<name>` keys.
- Use `deepcopy`/`.copy()` when pulling values out of ARTS workspace
  variables into returned results to avoid aliasing ARTS-owned memory.
- Follow existing NumPy-style docstrings on public methods.
- Keep new agendas consistent with existing `dobatch_calc_agenda_*` ordering:
  extract atmosphere → surface → geography → sun → checks → RT call → (for
  flux agendas) reset large unneeded tensors to save memory.

## Repo memory

Session/repo notes from prior work (e.g. the batch radiance implementation
plan) are kept under the agent's `/memories/repo/` and `/memories/session/`
scopes, not in this file — check those before re-deriving ARTS API details
already verified in a previous session.
