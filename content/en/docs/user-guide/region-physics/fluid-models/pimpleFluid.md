---
title: "pimpleFluid"
mathjax: true
description: "Transient, incompressible, turbulent single-phase flow on static or moving meshes using the PIMPLE algorithm."
tags: ["fluid", "incompressible", "turbulence", "moving-mesh"]
weight: 10
---

## Overview

`pimpleFluid` solves the transient, incompressible Navier-Stokes equations for
a Newtonian single-phase fluid using the **PIMPLE** (merged PISO–SIMPLE)
pressure-velocity coupling algorithm. It supports turbulence modelling (laminar,
RANS, LES), arbitrary moving meshes (ALE formulation), multiple reference frames
(MRF), and run-time selectable source terms via `fvOptions`. It is the most
general-purpose fluid region type in multiPhysicsFoam and is suitable as the
fluid side in conjugate heat transfer, fluid-structure interaction, and
standalone flow problems.

## Governing Equations

### Continuity

$$
\nabla \cdot \mathbf{U} = 0
$$

### Momentum

$$
\frac{\partial \mathbf{U}}{\partial t} + \nabla\cdot\left(\mathbf{U}\mathbf{U}\right) - \nabla\cdot\mathbf{R}
 = -\nabla p_\text{kin} + \mathbf{S}_U
$$


| Symbol | Description | Unit |
|--------|-------------|------|
| $\mathbf{U}$ | Velocity | m/s |
| $p_\text{kin} = p/\rho$ | Kinematic pressure | m²/s² |
| $\mathbf{R}$ | Effective stress tensor (viscous + turbulent) | m²/s² |
| $\mathbf{S}_U$ | Momentum source (from `fvOptions`) | m/s² |

The effective stress tensor $\mathbf{R}$ is computed by the selected turbulence
model via `turbulence::divDevReff(U)`. For laminar flow this reduces to the
viscous stress $2\nu\,\text{dev}(\mathbf{D})$ where $\mathbf{D}$ is the
rate-of-strain tensor.

### Stress tensor (for FSI coupling)

After each pressure-correction step the Cauchy stress tensor is updated:

$$
\boldsymbol{\sigma} = \rho \left( -p_\text{kin}\,\mathbf{I} - \mathbf{R}^\text{dev} \right)
$$

This field (`sigma`) is used when `pimpleFluid` is coupled to a solid region
via `fsiInterface`.

## Numerical Algorithm

The PIMPLE loop consists of an outer iteration (SIMPLE-like relaxation) wrapping
an inner PISO loop:

1. **Pre-predictor** — mesh motion update (ALE); flux correction after mesh
   deformation (`correctPhi`); MRF zone update.
2. **Momentum predictor** — assemble and (optionally) solve the momentum
   equation to obtain a predicted velocity $\mathbf{U}^*$.
3. **Pressure correction (PISO inner loop)** — solve the pressure Poisson
   equation:

$$
\nabla \cdot \left( r_{AU}\,\nabla p_\text{kin} \right) = \nabla \cdot \hat{\mathbf{U}}
$$

   where $r_{AU} = (A_U)^{-1}$ and $\hat{\mathbf{U}} = \mathbf{H}(U)/A_U$.
   Velocity is then explicitly corrected:

$$
\mathbf{U} = \hat{\mathbf{U}} - r_{AU}\,\nabla p_\text{kin}
$$

4. **Non-orthogonal corrector sub-loop** — improves accuracy on non-orthogonal
   meshes (controlled by `nNonOrthogonalCorrectors`).
5. **Turbulence correction** — transport equations for the selected turbulence
   model are solved once per outer PIMPLE iteration (controlled by
   `turbOnFinalIterOnly` / `turbCorr`).

Steps 2–5 are repeated for `nOuterCorrectors` outer iterations; steps 3–4 are
repeated for `nCorrectors` inner iterations.

## Sub-Models

| Sub-model | Selection mechanism | Notes |
|-----------|---------------------|-------|
| Turbulence (laminar / RANS / LES) | `constant/turbulenceProperties` | Uses `incompressible::turbulenceModel` |
| Transport properties | `constant/transportProperties` | Kinematic viscosity $\nu$ via `singlePhaseTransportModel` |
| Moving Reference Frame (MRF) | `constant/MRFProperties` | `IOMRFZoneList`; multiple zones supported |
| `fvOptions` | `system/fvOptions` | Explicit/implicit source terms, porosity, momentum sources |

## Required Fields

| Field | Type | Description | Unit |
|-------|------|-------------|------|
| `U` | `volVectorField` | Velocity | m/s |
| `pKin` | `volScalarField` | Kinematic pressure $p/\rho$ | m²/s² |
| Turbulence fields | — | As required by selected model (e.g. `k`, `epsilon`, `omega`, `nut`) | — |

Additionally created at runtime:

| Field | Type | Description | Condition |
|-------|------|-------------|-----------|
| `phi` | `surfaceScalarField` | Face-normal flux | always |
| `Uf` | `surfaceVectorField` | Face velocity | moving mesh only |
| `sigma` | `volSymmTensorField` | Cauchy stress tensor | always (FSI) |

## Case Setup

### Required dictionaries

```
constant/transportProperties    – kinematic viscosity (nu)
constant/turbulenceProperties   – turbulence model selection
system/fvSolution               – PIMPLE sub-dict and linear solver settings
system/fvSchemes                – discretisation schemes
```

### `multiRegionProperties` entry

```
regions
(
    ( fluid (pimpleFluid) )
);
```

For combined fluid + passive scalar transport add further region types:

```
regions
(
    ( fluid (pimpleFluid transportTemperature) )
);
```

### `controlDict` (singleRegionMode)

When running as a standalone single-region case via `singleRegionMode`:

```
singleRegionMode    yes;
regionType          pimpleFluid;
```

### Key `fvSolution` settings

```
PIMPLE
{
    momentumPredictor       yes;
    nOuterCorrectors        1;      // outer SIMPLE-like iterations
    nCorrectors             2;      // inner PISO iterations
    nNonOrthogonalCorrectors 1;

    // Moving mesh
    correctPhi              yes;    // defaults to yes if mesh is dynamic
    checkMeshCourantNo      no;
    moveMeshOuterCorrectors no;

    // Pressure reference (closed domains)
    pKinRefCell             0;
    pKinRefValue            0.0;
}
```

### Adjustable time stepping (`controlDict`)

```
adjustTimeStep  yes;
maxCo           0.5;        // maximum Courant number
maxDeltaT       1e-3;       // upper bound on time-step size
```

## Coupling

`pimpleFluid` participates in the following multi-region interfaces:

| Interface type | Exchanged fields | Coupling modes |
|----------------|-----------------|----------------|
| `heatTransferInterface` | Temperature `T`, heat flux `q` | partitioned (DNA), monolithic |
| `fsiInterface` | Displacement `D`, traction `sigma·n` | partitioned (DNA) |
| `defaultInterface` | Any scalar/vector | partitioned (DNA) |

The stress tensor `sigma` is automatically maintained and made available to
`fsiInterface` for computing the traction boundary condition on the solid side.

## Known Limitations

- **Incompressible only** — density is constant (read from
  `transportProperties`). Compressibility, variable density, or
  buoyancy-driven flows are not supported; use `buoyantPimpleFluid` for
  the latter.
- **Newtonian fluid** — viscosity is treated as constant or temperature-
  independent. Non-Newtonian rheology is not supported via this region type.
- **Single phase** — two-phase or interface-tracking flows are not supported;
  use `interTrackFluid` instead.
- **No energy equation** — the fluid energy/temperature equation is not solved
  here. Pair with `transportTemperature` or `turbulentTransportTemperature` for
  heat transfer in the fluid.
- **Monolithic p-U coupling not implemented** — `setCoupledEqns` is a no-op.
  For monolithic pressure-velocity coupling use `pUCoupledIcoFluid`.

## Tutorials

| Tutorial path | Description |
|---------------|-------------|
| `tutorials/singleRegion/pimpleFluid/pitzDaily` | Backward-facing step flow (single-region, adjustable time-step) |
| `tutorials/singleRegion/pimpleFluid/mixerVesselAMI2D` | Rotating mixer with AMI interface (MRF / sliding mesh) |
| `tutorials/conjugateHeatTransfer/flowOverHeatedPlate` | Flow over a heated plate (pimpleFluid + conductTemperature) |
| `tutorials/conjugateHeatTransfer/ShellAndTubeHeatExchanger` | Shell-and-tube heat exchanger (pimpleFluid + conductTemperature) |
| `tutorials/fluidStructureInteraction/beamInCrossFlow` | Beam in cross-flow FSI (pimpleFluid + solidStVK) |
| `tutorials/fluidStructureInteraction/HronTurekFsi3` | Hron-Turek FSI benchmark (pimpleFluid + solidStVK) |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
