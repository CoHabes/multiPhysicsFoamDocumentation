---
mathjax: true
title: "icoFluid"
description: "Transient, incompressible, laminar Newtonian flow with density-weighted momentum and static pressure, on static or moving meshes, using a PISO corrector loop."
tags: ["region-type", "fluid", "incompressible", "laminar", "moving-mesh"]
weight: 10
---

## Overview

`icoFluid` solves the transient, incompressible Navier-Stokes equations for a
**laminar Newtonian** fluid. Unlike `pimpleFluid`, it works with the dynamic
viscosity $\mu$ and the density $\rho$ directly: the momentum equation is
density-weighted and the pressure `p` is the static pressure in Pa, not the
kinematic pressure. It runs a PISO pressure-corrector loop inside the global
PIMPLE outer loop of `multiRegionSystem`. It supports moving meshes (ALE) and
translational moving reference frames from the `movingReferenceFrame` library.
Several of its options (the `closedVolume` flux correction, the `space` patch
and the `PISO` pressure reference) match the conventions of the partitioned
two-phase (interface-tracking) coupling BCs. It does not provide turbulence
modelling or `fvOptions`.

## Governing Equations

### Continuity

$$
\nabla \cdot \mathbf{U} = 0
$$

### Momentum

$$
\frac{\partial (\rho \mathbf{U})}{\partial t} + \nabla \cdot \left( \rho_f \, \phi \, \mathbf{U} \right) - \nabla \cdot \left( \mu \nabla \mathbf{U} \right) - \nabla \cdot \left( \mu \, \text{dev}\!\left[ (\nabla \mathbf{U})^T \right] \right) + \mathbf{F}_\text{MRF}
= -\nabla p
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $\mathbf{U}$ | Velocity | m/s |
| $p$ | (Static) pressure | Pa |
| $\rho$ | Density, $\rho_f$ its face interpolate | kg/m³ |
| $\mu$ | Dynamic viscosity | Pa·s |
| $\phi$ | Volumetric face flux | m³/s |
| $\mathbf{F}_\text{MRF}$ | Frame-acceleration term added by translational MRF zones (`addFrameAcceleration`) | N/m³ |

The convection term uses the scheme name `div(phi,U)` explicitly, even though
the actual flux is $\rho_f\phi$.

### Stress tensor (for FSI coupling)

After each PISO corrector the Cauchy stress is updated:

$$
\boldsymbol{\sigma} = -p\,\mathbf{I} + \mu \left( \nabla \mathbf{U} + (\nabla \mathbf{U})^T \right)
$$

Because $p$ is already in Pa and $\mu$ is dynamic, no extra $\rho$ scaling is
applied (unlike `pimpleFluid`, which builds $\boldsymbol\sigma$ from $p_\text{kin}$ and
the turbulence model's $\mathbf{R}^\text{dev}$).

## Numerical Algorithm

`usesPIMPLE()` returns `true`. `solveRegion()` and `setCoupledEqns()` are
no-ops, so all work happens in the PIMPLE hooks called by
`regionTypeList::solvePIMPLE()`. That function runs `nOuterCorrectors` outer
iterations, read from the **top-level** `system/fvSolution` `PIMPLE` dictionary.
Each iteration calls the following hooks on every region:

1. **`prePredictor()`**
   - Once per time step (guarded by `myTimeIndex_`): compute the Courant
     number from `phi` (`CourantNo.H`, plus the stock `meshCourantNo.H` when
     the mesh is moving). If `adjustTimeStep` is on, change `deltaT` directly
     (`setDeltaT.H`):
     $\Delta t_\text{new} = \max\big(\min(f\,\Delta t,\ \texttt{maxDeltaT}),\ \texttt{minDeltaT}\big)$
     with $f = \min\big(\min(\texttt{maxCo}/Co,\ 1 + 0.1\,\texttt{maxCo}/Co),\ 1.2\big)$.
   - If the mesh is changing, run the flux correction `correctPhi.H`. It
     re-evaluates fixed-value `U` patches, then solves a `pcorr` Poisson
     equation with $(1/A)_f = \Delta t$ inside a non-orthogonal corrector loop,
     and subtracts its flux from `phi`.
2. **`momentumPredictor()`**
   - If the mesh is changing, make `phi` relative: `phi == phi - fvc::meshPhi(rho, U)`.
   - Assemble `UEqn = ddt(rho,U) [+ MRF frame acceleration] + HUEqn` and
     under-relax it with the factor `relaxationFactors/<velocityName>`.
   - If `momentumPredictor` is on, solve `UEqn == -grad(p)`.
   - Reset the cached `UEqn` to its unrelaxed form (original source and
     boundary coefficients) so that the relaxation does not affect time
     consistency in the pressure corrector.
3. **`pressureCorrector()`**: PISO loop (`while (pimple_.correct())`, i.e.
   `nCorrectors` from the region `PIMPLE` dict):
   - $r_{AU} = 1/A_U$, $\mathbf{U} = \mathbf{H}\,r_{AU}$,
     $\phi = \big(\mathbf{H}_f\,(r_{AU})_f\big) \cdot \mathbf{S}_f$ (`H` and
     `rAU` are interpolated separately).
   - Boundary flux fix-up (`correctPatchPhi.H`):
     - `closedVolume` on: set `phi = U·Sf` on all non-coupled patches except
       the patch named `interfaceShadow`.
     - Otherwise: set it on all non-coupled patches where `p` is not fixed.
   - `closedVolume` on: `correctClosedVolumePhi` rescales the flux on the
     `regionCoupledVelocityValue` / `regionCoupledMassTransferVelocityValue`
     patch so that the net flux is zero.
   - `hasSpacePatch` on (and `closedVolume` off): `correctSpaceVolumePhi`
     scales the outflow on the patch named `space` to balance the inflow.
   - Non-orthogonal loop: solve
     $\nabla\cdot\big((r_{AU})_f\nabla p\big) = \nabla\cdot\phi$. The pressure
     reference is set only when `closedVolume` is on. On the final
     non-orthogonal iteration, apply `phi -= pEqn.flux()`.
   - Relax `p` except on the final iteration, then apply
     $\mathbf{U} \mathrel{-}= r_{AU}\nabla p$, correct the boundary
     conditions and update `sigma`.
   - Print min/max/mean of `p` and `U` and the region volume.

`meshMotionCorrector()` is called by `multiRegionSystem::preSolve()` every time
step, and inside the DNA loop for `pUPimple`/`D` coupling. It calls
`mesh().controlledUpdate()` and rebuilds `Uf` from `U` and `phi` if no `Uf`
file was read. If a patch named `interface` exists and the mesh is changing,
it also prints flux and continuity diagnostics for that patch.

`correct()` (called from `preSolve()`) runs `correctPhi.H` when the mesh is
changing. `postSolve()` updates the translational MRF zones once per time step
(`correctMRF()`, `correctBoundaryVelocity(U, phi)`).

`getMinDeltaT()` returns `GREAT`: the region never limits the global time step
through `multiRegionSystem`. Courant-based time-step control happens only
inside `prePredictor()` (see above).

### Differences from `pimpleFluid`

| Aspect | `icoFluid` | `pimpleFluid` |
|--------|-----------|---------------|
| Pressure field | `p` [Pa], density-weighted momentum | `pKin` [m²/s²], `p` derived |
| Viscosity | dynamic `mu`, laminar only | `nu` via `singlePhaseTransportModel` + turbulence model |
| MRF | translational `MovingReferenceFrameZones` | `IOMRFZoneList` (`MRFProperties`) |
| `fvOptions` | not supported | supported |
| Pressure reference | `PISO { pRefPoint; pRefValue; }`, applied only if `closedVolume` | `PIMPLE { pKinRefCell; pKinRefValue; }` via `setRefCell` |
| Time-step control | inside `prePredictor`, `getMinDeltaT()` = `GREAT` | `getMinDeltaT()` |
| `HbyA` treatment | no `constrainHbyA`/`constrainPressure`/`ddtCorr`/consistent | full stock `pimpleFoam` treatment |

## Sub-Models

| Sub-model | Selection mechanism | Notes |
|-----------|---------------------|-------|
| Transport properties | `constant/<region>/transportProperties` | Uniform `rho`, `mu` (used if no `rho`/`mu` field files exist) |
| Translational MRF | `constant/<region>/translationalMRFZones` (read if present) | PID-controlled frame velocity (`KP`, `KI`, `KD`, `centreFromMesh`, ...) |
| Mesh motion | `constant/<region>/dynamicMeshDict` | Through `mesh().controlledUpdate()` |

No turbulence model and no `fvOptions`.

## Required Fields

| Field | Type | Description | Unit |
|-------|------|-------------|------|
| `U` | `volVectorField` | Velocity (name configurable through `velocityName`) | m/s |
| `p` | `volScalarField` | Pressure | Pa |

Optional / created at runtime:

| Field | Type | Description | Written |
|-------|------|-------------|---------|
| `rho` | `volScalarField` | Read if present, otherwise uniform from `transportProperties` | no |
| `mu` | `volScalarField` | Read if present, otherwise uniform from `transportProperties` | no |
| `phi` | `surfaceScalarField` | Volumetric flux, initialised as $\mathbf{U}_f\cdot\mathbf{S}_f$ | yes |
| `phiHbyA` | `surfaceScalarField` | Initialised from `phi` (not used by the algorithm) | yes |
| `pcorr` | `volScalarField` | Mesh-motion flux correction; `fixedValue` where `p` fixes value, else `zeroGradient` | yes |
| `Uf` | `surfaceVectorField` | Face velocity (read if present) | yes |
| `sigma` | `volSymmTensorField` | Cauchy stress | no |

## Case Setup

### Required dictionaries

```
constant/<region>/transportProperties   – rho, mu (dimensioned, MUST_READ)
system/<region>/fvSolution               – solvers, PIMPLE, PISO, relaxationFactors
system/<region>/fvSchemes                – schemes incl. a fluxRequired sub-dict
system/fvSolution                        – top-level PIMPLE { nOuterCorrectors; }
```

### `multiRegionProperties` entry

```
regions
(
    ( fluid (icoFluid) )
);
```

### `controlDict` (singleRegionMode)

```
singleRegionMode    yes;
regionType          icoFluid;
```

### `transportProperties`

```
rho     rho [1 -3 0 0 0 0 0] 1.173;
mu      mu  [1 -1 -1 0 0 0 0] 1.856e-5;
```

### Key region `fvSolution` settings

```
solvers
{
    p     { solver PCG; preconditioner DIC; tolerance 1e-06; relTol 0; }
    pcorr { solver PCG; preconditioner DIC; tolerance 1e-06; relTol 0; }
    U     { solver BiCGStab; preconditioner DILU; tolerance 1e-05; relTol 0; }
}

PIMPLE                          // read by pimpleControl
{
    momentumPredictor        yes;
    nCorrectors              2;   // PISO correctors
    nNonOrthogonalCorrectors 1;
}

PISO
{
    pRefValue  0;               // always required
    pRefPoint  (0 0 0);         // required if p needs a reference
}

relaxationFactors               // sub-dict must exist (may be empty)
{
    // U  0.7;                  // momentum relaxation, key = velocityName, default 1
}

closedVolume   no;              // default no; top level of fvSolution
hasSpacePatch  no;              // default no; requires a patch named "space"
// velocityName U;              // default U
```

### Key `fvSchemes` entries

```
divSchemes
{
    div(phi,U)  Gauss linear;   // name used explicitly by the convection term
}

fluxRequired                    // sub-dict must exist; each non-default entry is registered
{
    p;
    pcorr;
}
```

### Adjustable time stepping (`controlDict`)

```
adjustTimeStep  yes;      // default no
maxCo           0.5;      // default 1.0
maxDeltaT       1e-3;     // default GREAT
minDeltaT       1e-8;     // required when adjustTimeStep is on
```

## Coupling

No interface type refers to `icoFluid` by name. The partitioned p-U BCs in
`src/fvPatchFields/coupledFvPatchFields/partitionedCoupledFvPatchFields`
follow the same conventions as `icoFluid`: a `p` field in Pa, a registered
`rho` field, `rho`/`mu` in `transportProperties`, a `PISO.pRefPoint`, an
`interfaceShadow` patch and `closedVolume`. `icoFluid` also provides the
fields those BCs and the interfaces look up.

| Interface type | Fields used from the `icoFluid` region | Coupling mode |
|----------------|----------------------------------------|---------------|
| `capillaryInterface` | `U`, `p`, `rho` (field), `transportProperties` (`rho`, `mu`) via `pUPimple` | partitioned (DNA) |
| `massTransferInterface` | `U`, `p`, `phi`, `transportProperties` via the `regionCoupledMassTransfer*` BCs | partitioned (DNA) |
| `fsiInterface` | `sigma` (read by `regionCoupledTraction` on the solid side) | partitioned (DNA, field `D`) |

`regionCoupledPressureValue` / `regionCoupledMassTransferPressureValue` also
need a `g` (`uniformDimensionedVectorField`) in the fluid region, which
`icoFluid` does not create.

## Known Limitations

- **Laminar, Newtonian, incompressible only**: no turbulence model,
  constant (or user-supplied, not updated) `rho` and `mu`.
- **No `fvOptions`**: no source terms or constraints.
- **No monolithic p-U coupling**: `setCoupledEqns` is a no-op. The monolithic
  variant `pUCoupledIcoFluid` is currently disabled (commented out in
  `src/multiRegionSystem/Make/files`).
- **Time-step control bypasses `getMinDeltaT()`**: `deltaT` is changed inside
  `prePredictor()`, after the solver has already picked the time step.
  `myTimeIndex_` is only advanced in `postSolve()`, and that only runs for
  `pUPimple`/`D` DNA coupling (`postSolvePIMPLE`). Without that coupling,
  Courant evaluation and `setDeltaT` run in every outer corrector and the
  translational-MRF update in `postSolve()` never runs.
- **Hard-coded patch names**: `interfaceShadow` (excluded from the
  `closedVolume` boundary-flux fix-up), `space` (`hasSpacePatch`), and
  `interface` (diagnostics in `meshMotionCorrector`).
- **`closedVolume` needs a coupled velocity patch**: `correctClosedVolumePhi`
  looks for a `regionCoupledVelocityValue` or
  `regionCoupledMassTransferVelocityValue` patch on `U`. With none present it
  indexes patch `-1`. The same applies to `hasSpacePatch` without a `space`
  patch.
- **Pressure reference only with `closedVolume`**: the `p` equation is
  referenced only when `closedVolume` is on, while `PISO.pRefValue` is always
  required. In `setRefCell.H` the `break` sits inside the `forAll` over
  processors. If the reference point is found on more than one processor and
  not on processor 0, no processor is selected (`procWithRef = -1`).
- **Simplified `HbyA`**: no `constrainHbyA`, `constrainPressure`, `ddtCorr` or
  `adjustPhi`, so open-domain pressure BCs and flux consistency follow the
  `correctPatchPhi.H` logic only.
- **Cosmetic**: the logged velocity "min" actually prints `gMax(U)`.

## Tutorials

No tutorial currently uses `icoFluid`. The multiphase tutorials under
`tutorials/multiphaseFlow/` (e.g. `2dRisingBubble`) use the same dictionary
conventions (`PISO`, `closedVolume`, `transportProperties` with `rho`/`mu`)
but select `interTrackFluid`, which is also disabled in
`src/multiRegionSystem/Make/files`.

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
