---
mathjax: true
title: "rhoPimpleFluid"
description: "Transient, compressible, turbulent single-phase flow with energy equation on static or moving meshes using the PIMPLE algorithm."
tags: ["region-type", "fluid", "compressible", "turbulence", "moving-mesh", "heat-transfer"]
weight: 10
---

## Overview

`rhoPimpleFluid` solves the transient, compressible Navier-Stokes equations for
a Newtonian single-phase fluid with the **PIMPLE** (merged PISO–SIMPLE)
pressure-velocity coupling algorithm. It is a port of OpenFOAM's
`rhoPimpleFoam`: density follows from a run-time selectable thermophysical
model (`fluidThermo`), and an energy equation for sensible enthalpy or internal
energy is solved alongside momentum. It supports compressible turbulence
modelling (laminar, RANS, LES), a transonic pressure formulation, local time
stepping (LTS), moving meshes, MRF and `fvOptions`. It is suitable for
subsonic to supersonic gas flows, standalone or as the fluid side of a
multi-region problem.

## Governing Equations

### Continuity

$$
\frac{\partial \rho}{\partial t} + \nabla \cdot \left(\rho\mathbf{U}\right) = S_\rho
$$

### Momentum

$$
\frac{\partial \left(\rho\mathbf{U}\right)}{\partial t} + \nabla \cdot \left(\rho\mathbf{U}\mathbf{U}\right) - \nabla \cdot \mathbf{R}
= -\nabla p + \mathbf{S}_U
$$

The effective stress divergence is provided by the turbulence model via
`compressible::turbulenceModel::divDevRhoReff(U)`.

### Energy

The energy variable $e$ is either sensible enthalpy $h$ or sensible internal
energy $e$, selected by the `energy` keyword in `thermophysicalProperties`.
With $K = \tfrac{1}{2}|\mathbf{U}|^2$:

$$
\frac{\partial \left(\rho h\right)}{\partial t} + \nabla \cdot \left(\rho\mathbf{U}h\right) + \frac{\partial \left(\rho K\right)}{\partial t} + \nabla \cdot \left(\rho\mathbf{U}K\right) - \frac{\partial p}{\partial t} - \nabla \cdot \left(\alpha_\text{eff}\nabla h\right)
= S_h
$$

For internal energy the pressure-work term $-\partial p/\partial t$ is
replaced by $\nabla \cdot \left(p\,\mathbf{U}\right)$ (`div(phiv,p)`). The
$\partial p/\partial t$ term is only included if `dpdt` is enabled in the
thermophysical model (default `yes`). With active MRF zones the term
$\nabla \cdot (\phi_\text{MRF}\,p)$ is added.

### Equation of state

Density is closed by the thermophysical model, $\rho = \rho(p, T)$, with
compressibility $\psi = \partial\rho/\partial p$ used in the pressure
equation.

| Symbol | Description | Unit |
|--------|-------------|------|
| $\rho$ | Density | kg/m³ |
| $\mathbf{U}$ | Velocity | m/s |
| $p$ | (Static) pressure | Pa |
| $T$ | Temperature | K |
| $h$, $e$ | Sensible enthalpy / internal energy | J/kg |
| $K$ | Specific kinetic energy | m²/s² |
| $\psi$ | Compressibility | s²/m² |
| $\mathbf{R}$ | Effective stress tensor (viscous + turbulent) | Pa |
| $\alpha_\text{eff}$ | Effective thermal diffusivity (laminar + turbulent) | kg/(m·s) |
| $S_\rho$, $\mathbf{S}_U$, $S_h$ | Sources (from `fvOptions`) | — |

### Stress tensor (for FSI coupling)

After each pressure-correction stage the Cauchy stress tensor is updated:

$$
\boldsymbol{\sigma} = -p\,\mathbf{I} - \mathbf{R}^\text{dev}
$$

where $\mathbf{R}^\text{dev}$ is `turbulence::devRhoReff()`. Since $p$ is the
static pressure in Pa, no density scaling is applied (unlike `pimpleFluid`).

## Numerical Algorithm

The region implements the `regionType` PIMPLE interface
(`usesPIMPLE() = true`), so the outer PIMPLE loop is driven by
`multiRegionSystem`. Per time step:

1. **`correct()`** — reset the outer-iteration counter; if `correctPhi` is
   active, snapshot $\nabla\cdot(\rho\mathbf{U})$ for flux correction; if LTS is
   active, recompute the local time-step field `rDeltaT` (see below).
2. **Pre-predictor** — on the first outer iteration of the step (or every
   iteration with `moveMeshOuterCorrectors yes`): `mesh.controlledUpdate()`;
   if the mesh changed, update MRF, rebuild the flux from `rhoUf`, apply
   `CorrectPhi`, make the flux relative and correct `rhoUf`; optionally print
   the mesh Courant number.
3. **Momentum predictor** —
   - on the first outer iteration and unless `SIMPLErho` is set, solve the
     density transport equation;
   - assemble the momentum matrix (cached for the pressure corrector), relax,
     and solve it with $-\nabla p$ if `momentumPredictor yes`; update $K$;
   - assemble and solve the energy equation, then `thermo.correct()`.
4. **Pressure corrector (PISO loop, `nCorrectors`)** — for each corrector:
   update $\rho$ from thermo (unless `SIMPLErho`), build
   $\mathbf{H}/A$ and $\phi_{HbyA}$ (with `ddtCorr`, MRF and mesh-motion
   corrections, optional `consistent` SIMPLEC treatment), constrain pressure
   BCs, and solve the compressible pressure equation inside the
   non-orthogonal loop:

$$
\frac{\partial \rho}{\partial t} + \psi\,\frac{\partial p'}{\partial t} + \nabla\cdot\phi_{HbyA} - \nabla\cdot\left(\rho\, r_{AU}\,\nabla p\right) = S_p
$$

   With `transonic yes` the flux is split into a convective part
   $\phi_d = (\psi_f/\rho_f)\,\phi_{HbyA}$ treated implicitly
   ($\nabla\cdot(\phi_d\,p)$) and the pressure equation is relaxed for
   diagonal dominance. Afterwards: explicit pressure relaxation, velocity
   correction $\mathbf{U} = \mathbf{H}/A - r_{AU}\nabla p$, pressure limiting
   (`pressureControl`), `thermo.correctRho` with `rhoMin`/`rhoMax`, a second
   density solve, continuity-error report, $\rho$ reset to thermo, `rhoUf`
   correction and $\partial p/\partial t$ update.
5. **Turbulence correction** — `turbulence->correct()` if `pimple.turbCorr()`
   (see `turbOnFinalIterOnly`), then update of `sigma`.

Steps 2–5 are repeated for `nOuterCorrectors` outer iterations; step 4 for
`nCorrectors` inner iterations.

**Mesh motion corrector** (`meshMotionCorrector()`, called by
`multiRegionSystem::preSolve` and in the FSI coupling loop) performs
`controlledUpdate`, MRF update, flux reconstruction from `rhoUf` +
`CorrectPhi`, and `rhoUf` correction.

### Time-step control

- **Global time stepping** — `getMinDeltaT()` returns
  $\Delta t_\text{max} = \text{maxCo}/\text{Co}\;\Delta t$, where the Courant
  number is evaluated from the mass flux, $\text{Co} = \tfrac{1}{2}\max\left(\sum_f|\phi|/(\rho V)\right)\Delta t$.
  `maxCo` is read from `controlDict` (default `1.0`).
- **Local time stepping (LTS)** — enabled automatically when `fvSchemes`
  selects `ddtSchemes { default localEuler; }`. `getMinDeltaT()` then returns
  `GREAT` (the global step is left untouched) and `rDeltaT` is computed per
  cell from the Courant number (plus an acoustic estimate with `transonic`),
  then smoothed and damped, as in stock `rhoPimpleFoam`.

## Sub-Models

| Sub-model | Selection mechanism | Notes |
|-----------|---------------------|-------|
| Thermophysical model | `constant/thermophysicalProperties` | `fluidThermo::New`; `hePsiThermo` or `heRhoThermo`; energy must be `h` or `e` variant |
| Turbulence (laminar / RANS / LES) | `constant/turbulenceProperties` | `compressible::turbulenceModel` |
| Moving Reference Frame (MRF) | `constant/MRFProperties` | `IOMRFZoneList` |
| `fvOptions` | `system/fvOptions` (or `constant/fvOptions`) | Sources/constraints on $\rho$, $\mathbf{U}$, $h$/$e$, $p$ (e.g. `limitTemperature`) |
| Dynamic mesh | `constant/dynamicMeshDict` | Any `dynamicFvMesh`; enables `rhoUf` and `correctPhi` |

## Required Fields

| Field | Type | Description | Unit |
|-------|------|-------------|------|
| `p` | `volScalarField` | Pressure (read by thermo) | Pa |
| `T` | `volScalarField` | Temperature (read by thermo) | K |
| `U` | `volVectorField` | Velocity | m/s |
| Turbulence fields | — | As required by selected model (e.g. `k`, `omega`, `nut`, `alphat`) | — |

Additionally created at runtime:

| Field | Type | Description | Condition |
|-------|------|-------------|-----------|
| `rho` | `volScalarField` | Density (read if present, else from thermo) | always |
| `phi` | `surfaceScalarField` | Mass flux $\rho\mathbf{U}\cdot\mathbf{S}_f$ (read if present) | always |
| `K` | `volScalarField` | Kinetic energy $\tfrac12|\mathbf{U}|^2$ | always |
| `dpdt` | `volScalarField` | $\partial p/\partial t$ | always (zero if `dpdt no`; written only on dynamic meshes) |
| `kappa` | `volScalarField` | Effective thermal conductivity `kappaEff` | always (heat-transfer coupling) |
| `sigma` | `volSymmTensorField` | Cauchy stress tensor | always (FSI) |
| `rhoUf` | `surfaceVectorField` | Face momentum | dynamic mesh only |
| `rDeltaT` | `volScalarField` | Reciprocal local time step | LTS only |
| `uniform/cumulativeContErr` | `uniformDimensionedScalarField` | Cumulative continuity error | always |

## Case Setup

### Required dictionaries

```
constant/thermophysicalProperties   – thermo type, mixture properties
constant/turbulenceProperties       – turbulence model selection
system/fvSolution                   – PIMPLE sub-dict and linear solver settings
system/fvSchemes                    – discretisation schemes (localEuler → LTS)
```

### `multiRegionProperties` entry

```
regions
(
    ( gas (rhoPimpleFluid) )
);
```

### `controlDict` (singleRegionMode)

```
application         multiPhysicsFoam;
singleRegionMode    yes;
regionType          rhoPimpleFluid;
```

### `thermophysicalProperties` (from `aerofoilNACA0012`)

```
thermoType
{
    type            hePsiThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleInternalEnergy;   // or sensibleEnthalpy
}

mixture
{
    specie         { molWeight 28.9; }
    thermodynamics { Cp 1005; Hf 0; }
    transport      { mu 1.82e-05; Pr 0.71; }
}
```

### Key `fvSolution` settings

```
PIMPLE
{
    momentumPredictor        yes;
    nOuterCorrectors         1;
    nCorrectors              2;
    nNonOrthogonalCorrectors 1;

    transonic                no;     // implicit convective pressure term
    consistent               no;     // SIMPLEC-like correction
    SIMPLErho                no;     // skip rho transport, take rho from thermo
    turbOnFinalIterOnly      true;

    // Pressure limits (pressureControl): pMin/pMax or pMinFactor/pMaxFactor
    pMinFactor               0.1;
    pMaxFactor               2;

    // Density bounds for thermo.correctRho
    rhoMin                   0;      // default 0
    rhoMax                   1e5;    // default GREAT

    // Moving mesh
    correctPhi               yes;    // defaults to yes if mesh is dynamic
    checkMeshCourantNo       no;
    moveMeshOuterCorrectors  no;

    // LTS only (ddtSchemes default localEuler)
    maxCo                    0.5;    // default 0.8
    rDeltaTSmoothingCoeff    0.1;    // default 0.02
    rDeltaTDampingCoeff      0.9;    // default 1.0 (no damping)
    maxDeltaT                1;      // default GREAT
}
```

### Adjustable time stepping (`controlDict`)

```
adjustTimeStep  yes;
maxCo           0.5;        // default 1.0
maxDeltaT       1e-4;
```

## Coupling

| Interface type | Exchanged fields | Coupling modes |
|----------------|-----------------|----------------|
| `heatTransferInterface` | Temperature `T`, heat flux via `kappa` | partitioned (DNA) |
| `fsiInterface` | Displacement `D`, traction from `sigma` | partitioned (DNA) |
| `defaultInterface` | Any scalar (generic `regionCoupledScalarJump`/`Flux` BCs) | partitioned (DNA) |

For partitioned heat transfer, `setCoupledEqns` assembles a separate
temperature equation (independent of the `fieldName` argument):

$$
\frac{\partial \left(\rho c_{pv} T\right)}{\partial t} + \nabla\cdot\left(c_{pv,f}\,\phi\,T\right)
= \nabla\cdot\left(\kappa_\text{eff}\nabla T\right)
$$

with $c_{pv} = $ `thermo.Cpv()` and `kappa` refreshed from
`turbulence->kappaEff()` on each call. The matrix is registered as
`T<region>MeshrhoPimpleFluidTypeEqn`. The `regionCoupledTemperatureJump` /
`regionCoupledHeatFlux` BCs read the diffusivity field `kappa` by default.
`selfCoupledFieldNames()` is not overridden (empty list).

## Known Limitations

- **Coupled temperature equation is decoupled from the energy equation** —
  the DNA `T` equation in `setCoupledEqns` contains no pressure-work,
  kinetic-energy, MRF or `fvOptions` terms, and `postSolve()` is a no-op, so
  `thermo.correct()` is not called after the coupled `T` solve and `he` is
  not updated from the new `T`.
- **Heat-transfer and FSI coupling are untested** — no multi-region
  tutorial uses `rhoPimpleFluid`; only single-region cases exist.
- **Monolithic coupling not implemented** — `setCoupledEqns` only provides
  the partitioned `T` matrix; monolithic coupling is disabled in
  `multiRegionSystem`.
- **`solveRegion()` and `postSolve()` are no-ops** — all work happens inside
  the PIMPLE hooks.
- **`maxDeltaT` not applied in `getMinDeltaT()`** — `maxDeltaT` from
  `controlDict` is read but the returned value is only Courant-based.
- **Mesh motion corrector does not make the flux relative** — unlike
  `prePredictor`, `meshMotionCorrector()` does not call
  `fvc::makeRelative` after `CorrectPhi`.
- **Single phase, single species** — reacting or multi-species mixtures and
  buoyancy (gravity, `p_rgh`) are not handled; use `buoyantPimpleFluid` for
  buoyancy-driven flows.

## Tutorials

| Tutorial path | Description |
|---------------|-------------|
| `tutorials/singleRegion/rhoPimpleFluid/aerofoilNACA0012` | Compressible RANS (kOmegaSST) flow around a NACA0012 aerofoil, internal energy, `limitTemperature` fvOption |
| `tutorials/singleRegion/rhoPimpleFluid/obliqueShock` | Laminar supersonic oblique shock with `transonic` formulation and LTS (`localEuler`) |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
