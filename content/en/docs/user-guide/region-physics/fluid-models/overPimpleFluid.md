---
mathjax: true
title: "overPimpleFluid"
description: "Transient, incompressible, turbulent single-phase flow on moving overset (Chimera) meshes using the PIMPLE algorithm."
tags: ["region-type", "fluid", "incompressible", "turbulence", "moving-mesh", "overset"]
weight: 10
---

## Overview

`overPimpleFluid` solves the transient, incompressible Navier-Stokes equations
for a Newtonian single-phase fluid on an **overset (Chimera) mesh** using the
**PIMPLE** pressure-velocity coupling algorithm. It is the multiPhysicsFoam
port of OpenFOAM's `overPimpleDyMFoam` and shares its feature set with
[`pimpleFluid`]({{< relref "pimpleFluid" >}}): turbulence modelling (laminar, RANS, LES),
MRF and `fvOptions`. In addition, it masks hole cells (`cellMask`), tags fringe
cells (`interpolatedCells`) and can make the predicted fluxes conservative
across the overset fringe (`oversetAdjustPhi`). The region mesh must be a
`dynamicOversetFvMesh`, which makes it suitable for flows around bodies that
undergo large rigid-body motions, e.g. rotors or moving objects.

## Governing Equations

### Continuity

$$
\nabla \cdot \mathbf{U} = 0
$$

### Momentum

$$
\frac{\partial \mathbf{U}}{\partial t} + \nabla \cdot \left(\mathbf{U}\mathbf{U}\right) - \nabla \cdot \mathbf{R}
= -\chi_c\,\nabla p_\text{kin} + \mathbf{S}_U
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $\mathbf{U}$ | Velocity | m/s |
| $p_\text{kin} = p/\rho$ | Kinematic pressure | m²/s² |
| $\mathbf{R}$ | Effective stress tensor (viscous + turbulent) | m²/s² |
| $\mathbf{S}_U$ | Momentum source (from `fvOptions`) | m/s² |
| $\chi_c$ | Cell mask (`cellMask`): 0 on hole cells, 1 elsewhere | – |

Overset coupling between the mesh components is not part of the equations
above: values in fringe (interpolated) cells are supplied by the overset
cell-cell stencil (`oversetInterpolation` in `fvSchemes`), and the `overset`
patch type couples the matrices of the individual mesh zones.

### Overset masks

Both masks are derived from the cell types of the cell-cell stencil
(`Stencil::New(mesh).cellTypes()`):

| Mask | Value 0 on | Value 1 on |
|------|------------|------------|
| `cellMask` $\chi_c$ | `HOLE` and `SPECIAL` cells | all other cells |
| `interpolatedCells` $\chi_i$ | `INTERPOLATED` (fringe) cells | all other cells |

The face mask $\chi_f$ is obtained by `localMin` interpolation of $\chi_c$,
i.e. a face is masked if either neighbouring cell is a hole. The flux is
multiplied by $\chi_f$ after each pressure correction and after mesh motion.

### Stress tensor (for FSI coupling)

As in `pimpleFluid`, the Cauchy stress tensor is updated after each pressure
corrector:

$$
\boldsymbol{\sigma} = \rho \left( -p_\text{kin}\,\mathbf{I} - \mathbf{R}^\text{dev} \right)
$$

## Numerical Algorithm

`usesPIMPLE()` returns `true`, so the framework drives the region through
`prePredictor` → `momentumPredictor` → `pressureCorrector` for each outer
PIMPLE iteration (`nOuterCorrectors` from the top-level `system/fvSolution`).
`solveRegion`, `postSolve`, `correct` and `setCoupledEqns` are no-ops.

1. **Pre-predictor (mesh motion and masking)** — executed at most once per
   time step (latched on the time index) unless `moveMeshOuterCorrectors` is
   enabled; `meshMotionCorrector()` runs the same block. If the mesh changed:
   - update MRF zones;
   - recompute `cellMask`, the face mask and `interpolatedCells` from the new
     overset stencil;
   - correct `Uf` and `phi` on faces that changed from hole to active
     (`faceMaskOld`), set $\phi = \mathbf{S}_f\cdot\mathbf{U}_f$ and mask hole
     faces;
   - make $\phi$ relative to the mesh motion;
   - optionally print the mesh Courant number (`checkMeshCourantNo`).
2. **Momentum predictor** — assemble the momentum equation (cached for the
   pressure corrector), relax and constrain it, and, if `momentumPredictor`
   is on, solve it with the pressure gradient masked on hole cells:
   $A_U\mathbf{U} - \mathbf{H} = -\chi_c\nabla p_\text{kin}$.
3. **Pressure corrector (PISO loop)** — for each of `nCorrectors`:
   - $r_{AU} = 1/A_U$, $\hat{\mathbf{U}} = r_{AU}\mathbf{H}$ (via
     `constrainHbyA`), $\phi_{\hat{U}} = \hat{\mathbf{U}}_f\cdot\mathbf{S}_f$
     stored in the persistent field `phiHbyA`; MRF makes it relative;
   - if `oversetAdjustPhi` is set, `phiHbyA` is made relative to the mesh
     motion, corrected by `oversetAdjustPhi(phiHbyA, U, adjustZone)` and made
     absolute again;
   - non-orthogonal loop solving
     $\nabla\cdot\left(r_{AU,f}\nabla p_\text{kin}\right) = \nabla\cdot\phi_{\hat{U}}$;
     on the final non-orthogonal iteration
     $\phi = \phi_{\hat{U}} - \text{flux}(p\text{Eqn})$ and
     $\mathbf{U} = \chi_c\left(\hat{\mathbf{U}} - r_{AU}\,\text{reconstruct}(\text{flux}/r_{AU,f})\right)$;
     with `oversetPatchErrOutput` the overset patch flux error is reported via
     `oversetPatchPhiErr`;
   - continuity errors are evaluated with fringe and hole cells excluded
     ($\chi_i\chi_c\nabla\cdot\phi$);
   - `Uf` is reset to $\overline{\mathbf{U}}_f$ with its normal component
     replaced by $\phi/|\mathbf{S}_f|$;
   - $\phi$ is made relative to the mesh motion and multiplied by $\chi_f$.
4. **Turbulence correction** — transport model and turbulence model are
   corrected once per outer iteration (`turbCorr`).
5. **Post-processing fields** — `p = pKin·rho` and `sigma` are updated.

### Time-step control

`getMinDeltaT()` re-reads the DyM controls and returns
$\Delta t_\text{max} = \text{maxCo}/\text{Co}\cdot\Delta t$, where the
Courant number is computed from $|\phi|$ weighted by the `localMin`
interpolation of $\chi_c + \chi_i$ (as in `overPimpleDyMFoam`). The
`maxDeltaT` cap is applied by the solver's `setDeltaT.H`.

### Differences to `pimpleFluid`

| Aspect | `pimpleFluid` | `overPimpleFluid` |
|--------|---------------|-------------------|
| Mesh | any `dynamicFvMesh` (static or moving) | `dynamicOversetFvMesh` (overset stencil required) |
| Hole / fringe handling | – | `cellMask`, face mask, `interpolatedCells` |
| Pressure gradient & velocity correction | unmasked | multiplied by `cellMask` |
| Fringe flux adjustment | – | optional `oversetAdjustPhi` |
| `correctPhi` after mesh motion | `pcorr` solve if `correctPhi` | no `pcorr` solve; only face-mask correction of `Uf`/`phi` |
| `ddtCorr` | applied if `pimple.ddtCorr()` | not applied |
| `Uf` | only for moving meshes | always created |
| Continuity error | all cells | fringe and hole cells excluded |
| Non-interpolated fields | – | `HbyA`, `grad(pKin)`, `surfaceIntegrate(phi)`, `surfaceIntegrate(phiHbyA)`, `cellMask`, `cellDisplacement`, `interpolatedCells`, `cellInterpolationWeight` are added to the stencil's `nonInterpolatedFields` |

## Sub-Models

| Sub-model | Selection mechanism | Notes |
|-----------|---------------------|-------|
| Turbulence (laminar / RANS / LES) | `constant/turbulenceProperties` | Uses `incompressible::turbulenceModel` |
| Transport properties | `constant/transportProperties` | $\nu$ via `singlePhaseTransportModel`; `rho` is also read |
| Moving Reference Frame (MRF) | `constant/MRFProperties` | `IOMRFZoneList` |
| `fvOptions` | `system/fvOptions` | Explicit/implicit source terms |
| Overset stencil | `oversetInterpolation` in `system/fvSchemes` | Any `cellCellStencil`, incl. the repo's `fastOversetCellStencil` (`src/numerics/overset`) |
| Mesh motion | `constant/dynamicMeshDict` | `dynamicOversetFvMesh` with any motion solver, e.g. `multiSolidBodyMotionSolver` |

## Required Fields

| Field | Type | Description | Unit |
|-------|------|-------------|------|
| `U` | `volVectorField` | Velocity | m/s |
| `pKin` | `volScalarField` | Kinematic pressure $p/\rho$ | m²/s² |
| `zoneID` | `volScalarField` | Overset mesh-component identifier (read by the overset stencil) | – |
| `pointDisplacement` | `pointVectorField` | Required by displacement-based motion solvers | m |
| Turbulence fields | — | As required by selected model | — |

All patches of these fields on overset boundaries must be of type `overset`.

Additionally created at runtime:

| Field | Type | Description | Written |
|-------|------|-------------|---------|
| `phi` | `surfaceScalarField` | Face-normal flux | yes |
| `phiHbyA` | `surfaceScalarField` | Predicted flux (persistent) | yes |
| `Uf` | `surfaceVectorField` | Face velocity (read if present) | yes |
| `cellMask` | `volScalarField` | Hole-cell mask | yes |
| `interpolatedCells` | `volScalarField` | Fringe-cell mask | yes |
| `sigma` | `volSymmTensorField` | Cauchy stress tensor | yes |
| `p` | `volScalarField` | Pressure `pKin·rho` | no |
| `cumulativeContErr` | `uniformDimensionedScalarField` | Cumulative continuity error (`<time>/uniform`) | yes |

## Case Setup

### Required dictionaries

```
constant/transportProperties    – nu and rho
constant/turbulenceProperties   – turbulence model selection
constant/dynamicMeshDict        – dynamicOversetFvMesh + motion solver
system/fvSolution               – PIMPLE sub-dict and linear solver settings
system/fvSchemes                – discretisation schemes + oversetInterpolation
system/controlDict              – libs (overset fvMotionSolvers)
```

### `multiRegionProperties` entry

```
regions
(
    ( fluid (overPimpleFluid) )
);
```

### `controlDict` (singleRegionMode)

```
libs                (overset fvMotionSolvers);

singleRegionMode    yes;
regionType          overPimpleFluid;

adjustTimeStep      true;
maxCo               1;          // default 1.0
maxDeltaT           1e-3;       // default GREAT
```

### `transportProperties`

```
transportModel  Newtonian;
nu              1e-05;
rho             1;      // used for p = pKin*rho and sigma
```

### `dynamicMeshDict`

```
dynamicFvMesh       dynamicOversetFvMesh;

solver          multiSolidBodyMotionSolver;

multiSolidBodyMotionSolverCoeffs
{
    movingZone1
    {
        solidBodyMotionFunction rotatingMotion;
        rotatingMotionCoeffs
        {
            origin      (0.005 0.005 0.005);
            axis        (0 0 1);
            omega       100.0;
        }
    }
}
```

### Key `fvSolution` settings

```
PIMPLE
{
    momentumPredictor        false;
    nOuterCorrectors         1;
    nCorrectors              3;
    nNonOrthogonalCorrectors 0;

    // Overset controls
    oversetAdjustPhi         false;   // default false
    adjustZone               -1;      // default -1, zone used by oversetAdjustPhi
    oversetPatchErrOutput    false;   // default false

    // Moving mesh
    checkMeshCourantNo       no;      // default false
    moveMeshOuterCorrectors  no;      // default false

    // Pressure reference (pKin without fixed-value patch)
    pKinRefPoint             (0.0001 0.0001 0.001);  // or pKinRefCell (default 0)
    pKinRefValue             0.0;                    // default 0
}
```

`correctPhi` (default: `mesh.dynamic()`) is read but currently has no effect.

### Key `fvSchemes` settings

```
oversetInterpolation
{
    method                  fastOversetCellStencil;   // or inverseDistance, cellVolumeWeight, ...
    weightMethod            correctedCellVolumeWeight; // fastOversetCellStencil only
    allowInterpolatedDonors false;
    excludedReachablePatches ( hole );
}

oversetInterpolationSuppressed
{
    grad(pKin);
    surfaceIntegrate(phiHbyA);
}

fluxRequired
{
    default         no;
    pKin            ;
}
```

### Mesh preparation

The overset mesh components are combined into a single mesh and `zoneID` is
set per component, e.g. (from `twoSimpleRotors/Allrun.pre`):

```
blockMesh
topoSet
subsetMesh box -patch hole -overwrite
topoSet
restore0Dir
setFields          # writes zoneID from cellSets
```

## Coupling

No interface type in `src/regionInterfaces` refers to `overPimpleFluid`
specifically. The region registers the same fields as `pimpleFluid` (`U`,
`pKin`, `p`, `phi`, `sigma`), so interfaces and partitioned BCs that look these
up by name (`fsiInterface` via `regionCoupledTraction`/`sigma`,
`heatTransferInterface`/`defaultInterface` via additional region types such as
`transportTemperature`) can in principle be attached to it.

| Interface type | Exchanged fields | Coupling modes |
|----------------|-----------------|----------------|
| `heatTransferInterface` | `T`, heat flux (with a temperature region type stacked) | partitioned (DNA) — untested |
| `fsiInterface` | Displacement `D`, traction `sigma·n` | partitioned (DNA) — untested |
| `defaultInterface` | Any scalar (generic `regionCoupledScalarJump`/`Flux` BCs) | partitioned (DNA) — untested |

None of these combinations is covered by a tutorial.

## Known Limitations

- **Overset mesh required** — the constructor calls `Stencil::New(mesh())`, so
  the region mesh must be a `dynamicOversetFvMesh` (and the `overset` library
  must be loaded).
- **No `correctPhi` / `pcorr` solve** — after mesh motion the flux is only
  corrected on hole ↔ active faces via the face mask; the `correctPhi`
  keyword is read but not used.
- **No `ddtCorr`** — the `fvc::ddtCorr` flux correction is not added to
  `phiHbyA`.
- **Monolithic coupling not implemented** — `setCoupledEqns` is a no-op.
- **`solveRegion`, `postSolve`, `correct` are no-ops** — all work happens in
  the PIMPLE hooks.
- **Incompressible, Newtonian, single phase, no energy equation** — same
  restrictions as `pimpleFluid`; pair with a temperature region type for heat
  transfer.
- **Multi-region use untested** — only a single-region tutorial exists.

## Tutorials

| Tutorial path | Description |
|---------------|-------------|
| `tutorials/singleRegion/overPimpleFluid/twoSimpleRotors` | Two counter-rotating rotors on overset meshes in a closed box (laminar, `fastOversetCellStencil`, single-region mode) |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
