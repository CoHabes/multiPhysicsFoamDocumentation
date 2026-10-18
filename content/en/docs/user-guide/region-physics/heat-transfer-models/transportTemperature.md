---
mathjax: true
title: "transportTemperature"
description: "Laminar advection-diffusion of temperature, stacked on top of an incompressible fluid region type that provides U and phi."
tags: ["region-type", "heat-transfer", "advection-diffusion", "add-on"]
weight: 10
---

## Overview

`transportTemperature` adds a temperature equation to a fluid region. It does
not solve for the flow itself: it is stacked on top of a flow region type in
the same region, e.g. `(fluid (pimpleFluid transportTemperature))`, and
advects `T` with the velocity `U` and face flux `phi` that the companion
region type owns. Material properties (`rho`, `cp`, `kappa`) are constant
(`kappa` may be a non-uniform field). There is no turbulent heat flux; the
turbulent variant `turbulentTransportTemperature` exists in the source tree
but is not compiled (see Known Limitations). The typical use is the fluid side
of conjugate heat transfer, coupled to `conductTemperature` through
`heatTransferInterface`.

## Governing Equations

### Energy (temperature)

$$
\rho\, c_p \left( \frac{\partial T}{\partial t} + \nabla \cdot \left( \phi\, T \right) \right)
= \nabla \cdot \left( \kappa\, \nabla T \right)
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $T$ | Temperature | K |
| $\phi$ | Volumetric face flux (`phi`) from the companion fluid region type | m³/s |
| $\rho$ | Density (constant) | kg/m³ |
| $c_p$ | Specific heat capacity (constant) | J/(kg K) |
| $\kappa$ | Thermal conductivity (uniform or field) | W/(m K) |

The convective term is written in conservative form $\nabla \cdot (\phi T)$
without a $-T\,\nabla\cdot\phi$ correction, so it assumes a divergence-free
(incompressible, volumetric) flux. There is no source term and no viscous
heating.

## Numerical Algorithm

`transportTemperature` has no solution loop of its own. `usesPIMPLE()` returns
`false`, and `solveRegion()`, `prePredictor()`, `momentumPredictor()`,
`pressureCorrector()`, `meshMotionCorrector()` and `postSolve()` are all
no-ops. Per time step:

1. `preSolve()`: `correct()` calls `kappa.correctBoundaryConditions()`.
2. The companion fluid region type solves its flow (e.g. the PIMPLE loop of
   `pimpleFluid`), which updates the shared `U` and `phi`.
3. Inside the DNA loop for `T`, `multiRegionSystem` calls
   `setCoupledEqns()`, which assembles the equation above as an
   `fvScalarMatrix` and registers it under
   `T<region>MeshtransportTemperatureTypeEqn`. The matrix is then relaxed
   (`relaxationFactors`) and solved, with optional Picard re-assembly
   controlled by the region's `PICARD` dictionary in `fvSolution`.
4. `getMinDeltaT()` returns `GREAT`, so this region never limits the time
   step. Time-step control comes from the companion fluid type (Courant
   number) or a fixed `deltaT`.

The temperature is therefore advanced once per time step, after the flow,
using the flux from the end of the flow solve (one-way, segregated coupling
from flow to temperature; `T` does not feed back into the flow).

`setCoupledEqns()` ignores its `fieldName` argument and always assembles the
`T` equation. `selfCoupledFieldNames()` is not overridden (returns an empty
list), so `T` is **only solved if it is listed as a partitioned-coupled field**
(a `DNA { T {...} }` entry in `multiRegionProperties` and `T` in an interface's
`coupledFields`).

## Sub-Models

None. `transportTemperature` does not construct `fvOptions`, a turbulence
model, or a thermophysical model.

## Required Fields

How each field is obtained:

| Field | Type | Description | Unit | Source |
|-------|------|-------------|------|--------|
| `U` | `volVectorField` | Velocity | m/s | Looked up from the mesh registry (provided by the companion fluid type); read from `0/<region>/U` (`MUST_READ`) if not registered yet |
| `phi` | `surfaceScalarField` | Volumetric face flux | m³/s | Looked up from the registry; otherwise read if a `phi` file exists, else created as `interpolate(U) & Sf` |
| `T` | `volScalarField` | Temperature | K | Looked up from the registry, otherwise read from `0/<region>/T` (`MUST_READ`, written) |
| `kappa` | `volScalarField` | Thermal conductivity | W/(m K) | Looked up from the registry, otherwise read from `0/<region>/kappa` if present, else created uniform from `transportProperties` (written) |

Because fields are looked up from the mesh registry first, the companion
region type must be listed **before** `transportTemperature` so that it
creates `U` and `phi`; every tutorial does this.

The coupling boundary conditions `regionCoupledTemperatureJump` and
`regionCoupledHeatFlux` look up the diffusivity field `kappa` by default (or
the name given by their `diffusivity` keyword), which matches the field this
region type registers.

## Case Setup

### Required dictionaries

```
constant/<region>/transportProperties   – rho, cp, kappa (all required)
system/<region>/fvSolution              – T solver, relaxation factor for T
system/<region>/fvSchemes               – ddt(T), div(phi,T), laplacian(kappa,T)
```

### `multiRegionProperties` entry

```
regions
(
    (fluid (pimpleFluid transportTemperature) )
    (solid (conductTemperature) )
);

DNA
{
    T
    {
        maxCoupleIter 20;

        residualControl
        {
            maxJumpRes 1e-07;
            outputJumpResField no;
            maxFluxRes 1e-07;
            outputFluxResField no;
        }
    }
}
```

### `transportProperties`

All three entries are read with `dimensionedScalar(name, dict)` and have no
defaults. They can share the dictionary with the companion fluid type's
entries (`transportModel`, `nu`):

```
rho             rho   [1 -3 0 0 0 0 0]  1;
cp              cp    [0 2 -2 -1 0 0 0] 500.0;
kappa           kappa [1 1 -3 -1 0 0 0] 5.0;
```

### Key `fvSolution` settings

```
solvers
{
    T
    {
        solver              GAMG;
        tolerance           1e-7;
        relTol              0;
        smoother            GaussSeidel;
    }
}

relaxationFactors
{
    T   1.0;
}
```

### Coupled `T` boundary condition (fluid side, partitioned)

```
interface
{
    type                    regionCoupledTemperatureJump;
    interfaceType           (heatTransferInterface);
    accType                 aitken;
    relax                   0.75;
    secondOrder             false;
    nonOrthCorr             false;
    value                   $internalField;
}
```

The solid side uses `regionCoupledHeatFlux` on its interface patch.

## Coupling

| Interface type | Exchanged fields | Coupling mode |
|----------------|-----------------|---------------|
| `heatTransferInterface` | `T` (jump on fluid side), heat flux $\kappa \nabla T \cdot \mathbf{n}$ (solid side) | partitioned (DNA) |

No interface type accesses `transportTemperature` directly. Coupling happens
through the `T` boundary conditions above, which read `kappa` from the region
registry. The partner region is `conductTemperature` in all tutorials,
possibly stacked with `solidStVK` for thermo-FSI.

Companion flow region types: all tutorials use `pimpleFluid`. Any flow type
that registers a volumetric `U` and `phi` on the same mesh should work in
principle (e.g. `icoFluid`), but this is untested. `rhoPimpleFluid`
registers a mass flux `phi` (`rho*U & Sf`), so it is dimensionally
inconsistent with the equation above.

## Known Limitations

- **Only solved as a DNA-coupled field** — `selfCoupledFieldNames()` is not
  overridden, so without a `DNA { T ... }` entry the temperature is never
  solved (in particular not in single-region mode). There is no
  `tutorials/singleRegion/transportTemperature` case.
- **Monolithic coupling not available** — monolithic coupling in
  `multiRegionSystem::solve()` is currently commented out.
- **Laminar heat conduction only** — no turbulent diffusivity
  ($\kappa_t$ / `Prt`). `turbulentTransportTemperature`
  (`src/regions/turbulentTransportTemperature`) adds this but is commented out
  in `src/multiRegionSystem/Make/files` and therefore not built.
- **Constant properties** — `rho` and `cp` are uniform constants; `kappa` is
  read once and only its boundary values are refreshed in `correct()`.
  No temperature-dependent properties.
- **One-way coupling** — `T` has no effect on the flow (no buoyancy). Use
  `buoyantPimpleFluid` for natural convection.
- **Incompressible flux assumed** — the convection term has no
  $\nabla\cdot\phi$ correction and expects a volumetric flux.
- **No `fvOptions`** — heat sources cannot be added.
- **No time-step constraint** — `getMinDeltaT()` returns `GREAT`; no
  diffusion-number based time-step limit.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/conjugateHeatTransfer/flowOverHeatedPlate` | Flow over a heated plate (pimpleFluid + transportTemperature / conductTemperature) |
| `tutorials/conjugateHeatTransfer/ShellAndTubeHeatExchanger` | Shell-and-tube heat exchanger with two fluid regions (inner, outer) and a solid wall |
| `tutorials/fluidStructureInteraction/heatedBeamInCrossFlow` | Thermo-FSI: heated beam in cross flow (fluid: pimpleFluid + transportTemperature; solid: solidStVK + conductTemperature; DNA on `D` and `T`) |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
