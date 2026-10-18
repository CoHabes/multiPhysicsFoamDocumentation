---
mathjax: true
title: "laplacianFluid"
description: "Transient pure-diffusion (heat conduction / passive scalar) solver on a static mesh, ported from laplacianFoam."
tags: ["region-type", "diffusion", "heat-transfer", "static-mesh"]
weight: 10
---

## Overview

`laplacianFluid` solves the transient diffusion equation for a scalar field
`T` on a **static** mesh. It is a direct port of the stock OpenFOAM
`laplacianFoam` solver into the `regionType` framework: no velocity field, no
pressure-velocity coupling, and no mesh motion are involved — only a single
diffusion equation with a constant or field diffusivity `DT`. Despite the
name `T`, the field is a generic passive scalar and can represent temperature
or any other diffusing quantity. It is suitable as a lightweight standalone
conduction region, or as one side of a multi-region coupling where only
diffusive transport of a scalar is required.

## Governing Equations

### Diffusion

$$
\frac{\partial T}{\partial t} - \nabla \cdot \left( D_T \, \nabla T \right) = S_T
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $T$ | Transported scalar (e.g. temperature) | K (or scalar-dependent) |
| $D_T$ | Diffusivity | m²/s |
| $S_T$ | Source term (from `fvOptions`) | [T]/s |

There is no continuity or momentum equation — `laplacianFluid` has no velocity
field and does not support advection.

## Numerical Algorithm

`laplacianFluid` does not run its own outer solve loop. `usesPIMPLE()` returns
`false`, and `solveRegion()`, `prePredictor()`, `momentumPredictor()`,
`pressureCorrector()`, and `meshMotionCorrector()` are all no-ops. Instead:

1. `setCoupledEqns("T")` assembles the diffusion equation as an
   `fvScalarMatrix` and registers it with `multiRegionSystem`.
2. The equation is solved through the standard coupled-matrix path
   (`multiRegionSystem::assembleAndSolveEqns<scalar>("T")`), the same
   mechanism used for DNA-coupled partitioned fields.
3. `getMinDeltaT()` returns `GREAT`, i.e. this region never constrains the
   global time step — `deltaT` must be set explicitly (no diffusion-number
   based adjustable time-stepping is implemented).
4. On write steps only, `postSolve()` recomputes the derived output fields
   `gradTx`, `gradTy`, `gradTz`, and `DTgradT = DT * grad(T)` if
   `writeGradT` is enabled (mirrors `laplacianFoam`'s `write.H`).

**Non-orthogonal correction caveat:** stock `laplacianFoam` wraps its solve in
a `while (simple.correctNonOrthogonal())` loop. The coupled-matrix solve path
used here has no equivalent non-orthogonal corrector loop, so results will
differ from `laplacianFoam` on non-orthogonal meshes. The closest substitute
is the region's Picard corrector (`maxCorr` / `convergenceTolerance` in
`fvSolution`), which re-invokes `setCoupledEqns()` per iteration.

## Sub-Models

| Sub-model | Selection mechanism | Notes |
|-----------|---------------------|-------|
| `fvOptions` | `system/fvOptions` | Explicit/implicit source terms for `T` |

No turbulence, transport-model, or MRF sub-models are used.

## Required Fields

| Field | Type | Description | Unit |
|-------|------|-------------|------|
| `T` | `volScalarField` | Transported scalar | — |
| `DT` | `volScalarField` | Diffusivity — read if present, otherwise created uniform from `transportProperties` | m²/s |

Additionally created at runtime if `writeGradT` is enabled (default `true`):

| Field | Type | Description |
|-------|------|-------------|
| `gradTx`, `gradTy`, `gradTz` | `volScalarField` | Components of $\nabla T$ |
| `DTgradT` | `volVectorField` | Diffusive flux $D_T \nabla T$ |

## Case Setup

### Required dictionaries

```
constant/transportProperties    – diffusivity (DT), if no DT field is supplied
system/fvSolution                – T solver settings, SIMPLE.writeGradT flag
system/fvSchemes                  – discretisation schemes
```

### `multiRegionProperties` entry

```
regions
(
    ( scalarRegion (laplacianFluid) )
);
```

### `controlDict` (singleRegionMode)

```
singleRegionMode    yes;
regionType          laplacianFluid;
```

### Key `fvSolution` settings

```
solvers
{
    T
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-06;
        relTol          0;
    }
}

SIMPLE
{
    writeGradT      true;   // write grad(T) components and DT*grad(T) at write time
}
```

### `transportProperties`

```
DT              4e-05;      // uniform diffusivity, used if no DT field is present
```

## Coupling

`laplacianFluid` exposes `T` via `selfCoupledFieldNames()`, so it is solved
standalone and can also participate as a pure-diffusion side in multi-region
interfaces (unlike `conductTemperature`, whose `T` is only solved through
partitioned coupling):

| Interface type | Exchanged fields | Coupling modes |
|----------------|-----------------|----------------|
| `heatTransferInterface` | `T`, heat flux `q` | partitioned (DNA) |
| `defaultInterface` | `T` | partitioned (DNA) |

## Known Limitations

- **Static mesh only** — `meshMotionCorrector()` is a no-op; moving/deforming
  meshes are not supported.
- **No non-orthogonal corrector loop** — unlike `laplacianFoam`, the coupled
  solve path used here does not iterate over non-orthogonal correctors within
  a single time step; see the numerical algorithm note above for the
  workaround.

## Tutorials

| Tutorial path | Description |
|---------------|-------------|
| `tutorials/singleRegion/laplacianFluid/flange` | Classic OpenFOAM "flange" thermal conduction case (single-region, fixed time step) |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
