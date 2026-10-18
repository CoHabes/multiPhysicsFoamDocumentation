---
mathjax: true
title: "fsiInterface"
description: "Partitioned fluid-structure interaction interface: fluid traction drives the solid displacement, the solid displacement drives the fluid mesh motion."
tags: ["interface-type", "fsi"]
weight: 10
---

## Overview

`fsiInterface` couples a fluid region with a deforming solid region
(`solidStVK`) along a common surface, using a partitioned Dirichlet-Neumann
scheme. The solid receives the fluid traction as a Neumann condition on its
displacement `D`. The fluid interface patch is moved to the deformed solid
surface by the `interfaceTrackingFvMesh` dynamic mesh of the fluid region. The
class has very little logic of its own. It builds a patch-to-patch mapping in
the deformed solid configuration, defines a geometric interface residual for
the DNA loop, and logs the total force on the solid. The traction transfer
lives in the `regionCoupledTraction` boundary condition. The displacement
transfer lives in `movingInterfacePatches` (`src/dynamicFvMesh`).

## Interface Conditions

Let $\Gamma$ be the interface between the fluid region $f$ and the solid
region $s$. $\mathbf{x}_0$ is the undeformed position of a solid interface point.

### Kinematic condition

$$
\mathbf{x}_f = \mathbf{x}_0 + \mathbf{D}_s \quad \text{on } \Gamma
$$

The fluid interface points follow the solid point displacement `pointD`.
Velocity continuity is not imposed through `fsiInterface`. In all FSI
tutorials the fluid velocity on the interface uses `movingWallVelocity`, so the
fluid velocity at the wall equals the velocity of the moving mesh.

### Dynamic condition

$$
\mathbf{n}_s \cdot \boldsymbol{\sigma}_s = -\,\mathbf{n}_f \cdot \boldsymbol{\sigma}_f \quad \text{on } \Gamma
$$

`regionCoupledTraction` implements this as a fixed-gradient condition on `D`.
The fluid traction $\mathbf{t}_f = \boldsymbol{\sigma}_f\cdot\mathbf{n}_f$
(fluid patch normals) is mapped to the solid patch and its sign is flipped,
$\mathbf{t} = -\mathbf{t}_f$. The normal gradient is then

$$
\frac{\partial \mathbf{D}}{\partial n}
=
\frac{1}{K}\left[\,\mathbf{t} - \mathbf{n}\cdot\boldsymbol{\sigma}_s + K\,\mathbf{n}_0\cdot\nabla\mathbf{D}\,\right],
\qquad
\mathbf{n} = \frac{\mathbf{F}^{-T}\cdot\mathbf{n}_0}{\lvert\mathbf{F}^{-T}\cdot\mathbf{n}_0\rvert}
$$

At convergence of the solid Picard loop, $\mathbf{n}\cdot\boldsymbol{\sigma}_s = \mathbf{t}$
holds in the deformed configuration.

| Symbol | Description | Unit |
|--------|-------------|------|
| $\mathbf{D}_s$ | Solid total displacement (`D`, point field `pointD`) | m |
| $\boldsymbol{\sigma}_f$ | Fluid Cauchy stress (`sigma` in the fluid region) | Pa |
| $\boldsymbol{\sigma}_s$ | Solid Cauchy stress (`sigma` in the solid region) | Pa |
| $\mathbf{n}_0$ | Undeformed solid patch normal | – |
| $\mathbf{n}$ | Deformed solid patch normal | – |
| $\mathbf{F}^{-1}$ | Inverse deformation gradient (`Finv`) | – |
| $K$ | Implicit stiffness $2\mu+\lambda$ (`impK`, $1/K$ = `rImpK`) | Pa |

### Interface residual

The DNA convergence measure is the mismatch between the deformed solid
interface points and the fluid interface points. The point positions of the
side listed second in `coupledPatchPair` (side B) are mapped to side A:

$$
r_\Gamma =
\frac{\lVert \mathbf{x}_s - \mathbf{x}_f \rVert_2}
{\max\left(\lVert \mathbf{x}_s \rVert_2,\ \lVert \mathbf{x}_f \rVert_2,\ \text{SMALL}\right)},
\qquad \mathbf{x}_s = \mathbf{x}_0 + \mathbf{D}_s
$$

The norms are global sums over all interface points (`gSum`). The face-wise
field written with `outputIntResField` is $\lvert\mathbf{x}_s-\mathbf{x}_f\rvert$,
interpolated from points to faces on patch A.

## Implementation

Which side is the solid is decided at run time. Side A is assumed to be the
solid unless the mesh of side B holds a `pointVectorField` named `pointD`.

1. **Construction.** The constructor reads the `fsiInterfaceCoeffs` sub-dict,
   which must be present but is not used. It then builds the
   interface-to-interface mapping (`makeInterfaceToInterface()`):
   1. The solid mesh is moved temporarily to `points + pointD`.
   2. An `interfaceToInterfaceMapping` of type `interfaceTransferMethod`
      (`AMI` or `directMap`) is built between the two patches.
   3. The solid mesh is moved back and its `moving` flag is reset. The solid
      therefore stays total Lagrangian.
2. **Mapping updates.** `interfaces_->update()` calls
   `updateInterpolatorAndGlobalPatches()` only if mesh A is `changing()`. The
   mapping is then rebuilt every `interpolatorUpdateFrequency` time steps (0 =
   never), or when mesh A changes without moving (topology change). A rebuild
   goes through the same deformed-configuration step as in 1.
3. **DNA loop for `D`** (`multiRegionSystem::solve()`). The stand-alone PIMPLE
   call is skipped when `D` is DNA-coupled. Each DNA iteration then runs:
   1. `solvePIMPLE()` on the fluid region(s). This updates the fluid `sigma`.
   2. `assembleAndSolveEqns("D")`. The solid momentum equation is solved with
      `regionCoupledTraction` on the interface. The BC also stores the total
      force $\sum \mathbf{t}\,\lvert\mathbf{S}_0\rvert$ on the interface via
      `setTotalForce()`.
   3. `meshMotionCorrector()` on all regions. For the fluid
      `interfaceTrackingFvMesh` with `fsiInterface yes`, the solid
      `pointD` on the neighbour patch is mapped to the fluid patch. The mapping
      is built in the deformed solid configuration as in 1. The mapped
      displacement is relaxed by the patch's own acceleration model (`accType`,
      `relax` in `dynamicMeshDict`). The increment relative to the previous
      iteration is imposed on the fluid patch, and the motion solver moves the
      interior points.
   4. `interfaces_->update()`. This runs step 2 and `info()`, which prints
      `Total force on solid interface: (...)`.
   5. `dnaControl` evaluates `normIntRes("D")` against `maxIntRes`.
4. After the loop, `postSolvePIMPLE()` runs on the fluid region(s).

`correct()` is a no-op. `getMinDeltaT()` returns `GREAT`, so the interface does
not limit the time step.

If `coupled no`, `regionCoupledTraction` returns before updating its gradient.
`dnaControl` then stops after the first pass, and the interface residual
returns 0 for non-coupled interfaces. `couplingStartTime` switches coupling on
once that time is reached.

## Coupled Fields

| Field | Required boundary conditions (fluid / solid) | Description |
|-------|-----------------------------------------------|-------------|
| `D` | none (the fluid has no `D`) / `regionCoupledTraction` | Solid displacement, driven by the fluid traction |

Fields looked up, but not coupled through a patch field pair:

| Field | Region | Used by | Provided by |
|-------|--------|---------|-------------|
| `pointD` | solid | `fsiInterface` (mapping, residual, solid-side detection), `movingInterfacePatches` | `solidStVK` |
| `sigma` (`volSymmTensorField`) | fluid | `regionCoupledTraction` (fluid traction) | `pimpleFluid` (also declared by `icoFluid`, `rhoPimpleFluid`, `overPimpleFluid`, `buoyantPimpleFluid`) |
| `sigma`, `impK`, `rImpK`, `gradD`, `Finv` | solid | `regionCoupledTraction` | `solidStVK` |

Only `maxIntRes` is meaningful for `D`. The fluid has no `D` field, so jump and
flux residuals are not evaluated for this interface.
`regionCoupledTraction` uses a direct flux, so it cannot provide a flux
residual anyway.

`fsiInterface` creates no fields itself. With `outputIntResField yes`,
`dnaControl` writes `DinitialIntRes` and `DfinalIntRes` on mesh A.

## Compatible Region Types

| Region type (fluid side) | Region type (solid side) | Notes |
|--------------------------|--------------------------|-------|
| `pimpleFluid` | `solidStVK` | Used in all FSI tutorials. The fluid mesh must be `interfaceTrackingFvMesh` with `fsiInterface yes` on the interface patch |
| `pimpleFluid transportTemperature` | `solidStVK conductTemperature` | Combined with `heatTransferInterface` on the same patch pair (`heatedBeamInCrossFlow`) |
| `icoFluid`, `rhoPimpleFluid`, `overPimpleFluid`, `buoyantPimpleFluid` | `solidStVK` | These provide a `sigma` stress field. Untested: no tutorial, and the dynamic-mesh support of each region type has not been checked |

## Case Setup

### `regionInterfaceProperties` entry

From `tutorials/fluidStructureInteraction/HronTurekFsi3`:

```
partitionedCoupledPatches
(
    fluidsolid
    {
        interfaceType (fsiInterface);

        coupledPatchPair
        (
            (fluid interfaceShadow)
            (solid interface)
        );

        coupledFields
        (
            D
        );

        coupled             no;     // default yes
        couplingStartTime   2.0;    // default 0; switches coupling on at t >= 2

        fsiInterfaceCoeffs          // required, currently no entries read
        {}
    }
);

monolithicCoupledPatches
(
);

curvatureCorrectedSurfacePatches 0();   // required (no default)

interpolatorUpdateFrequency 0;          // default 1; 0 = build mapping once

interfaceTransferMethod AMI;            // default AMI; or directMap

directMapCoeffs
{}

AMICoeffs                               // <interfaceTransferMethod>Coeffs is required
{}
```

`interfaceTransferMethod`, `<method>Coeffs` and `interpolatorUpdateFrequency`
are top-level entries and apply to every interface in the file.
`beamInCrossFlow` and `heatedBeamInCrossFlow` use `directMap` and `coupled yes`.

FSI and conjugate heat transfer on the same patch pair
(`heatedBeamInCrossFlow`):

```
interfaceType (fsiInterface heatTransferInterface);
coupledFields
(
    D
    T
);
```

### `multiRegionProperties` DNA entry

```
regions
(
    ( fluid (pimpleFluid) )
    ( solid (solidStVK) )
);

DNA
{
    D
    {
        maxCoupleIter 100;

        residualControl
        {
            maxIntRes 1e-7;
            outputIntResField yes;
        }
    }
}
```

### Fluid dynamic mesh (`constant/fluid/dynamicMeshDict`)

From `HronTurekFsi3`. The `interfaceShadow` entry names the solid patch and region:

```
dynamicFvMesh interfaceTrackingFvMesh;

solver velocityLaplacian;

diffusivity quadratic inverseDistance 2(interfaceShadow cylinder);

twoDMotion yes;

motionPlaneNormal (0 0 1);

interfaceTrackingFvMeshCoeffs
    movingSurfacePatches
    (
        interfaceShadow
        {
            fsiInterface yes;               // default false

            neighbourPatch   interface;     // solid patch
            neighbourMesh    solid;         // solid region (required with fsiInterface yes)

            accType           aitken;       // fixed (default) | aitken | IQNILS
            relax             0.05;         // default 1.0

            normalMotionDir no;

            fixedSurfacePatches 0();

            interfaceDeformationLimit interfaceDeformationLimit [ 0 0 0 0 0 0 0 ] 0.0;

            interpolatorUpdateFrequency 0;  // default 0

            pointNormalsCorrectionPatches 0();

            interfaceTransferMethod AMI;    // default AMI

            directMapCoeffs{}

            AMICoeffs{}

            smoothing no;

            smoothSurfaceMeshCoeffs{}
        }
    );
```

This mapping (fluid mesh motion) is separate from the `fsiInterface` mapping
(traction and residual). Both are configured independently. The tutorials use
the same method in both places.

### Boundary conditions

Solid `0/solid/D` (`HronTurekFsi3/0/solid/orig/partitioned/D`):

```
interface
{
    type                    regionCoupledTraction;
    neighbourFieldName      sigma;
    interfaceType           (fsiInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform (0 0 0);
}
```

`interfaceType` is required. `secondOrder` and `nonOrthCorr` default to
`false`. `accType` and `relax` may also be set here to relax the traction
gradient (default `fixed`, `1.0`). None of the tutorials do this.
`neighbourFieldName` is read by the coupling manager, but
`regionCoupledTraction` always looks up the fluid field named `sigma`.

Fluid `0/fluid/U`, `0/fluid/pKin` and `0/fluid/pointMotionU`
(`HronTurekFsi3/0/fluid/orig/partitioned/`):

```
// U
interfaceShadow
{
    type            movingWallVelocity;
    value           uniform (0 0 0);
}

// pKin
interfaceShadow
{
    type zeroGradient;
    value uniform 0;
}

// pointMotionU
interfaceShadow
{
    type            fixedValue;     // overwritten each motion step by interfaceTrackingFvMesh
    value           uniform (0 0 0);
}
```

`interfaceTrackingFvMesh` requires the `pointMotionU` interface patch to be
`fixedValue`. It `refCast`s it to `fixedValuePointPatchVectorField`.

## Known Limitations

- **Partitioned only.** There is no monolithic FSI path, because monolithic
  coupling is disabled in `multiRegionSystem`.
- **Solid must be `solidStVK`.** The solid side is detected by the presence of
  `pointD`. `regionCoupledTraction` needs `impK`, `rImpK`, `gradD`, `sigma` and
  `Finv`, which only `solidStVK` provides.
- **The fluid mesh motion is not done by `fsiInterface`.** Without an
  `interfaceTrackingFvMesh` entry with `fsiInterface yes` on the fluid patch,
  the fluid mesh does not follow the solid. The interface residual then never
  drops.
- **Two independent mappings.** `fsiInterface` (traction, residual) and
  `movingInterfacePatches` (displacement) each build their own
  interface-to-interface mapping from different dictionaries. Inconsistent
  settings are not detected.
- **Mapping is built on the configuration at construction time** when
  `interpolatorUpdateFrequency 0`, as in all tutorials. With large
  deformations, the AMI/directMap weights are not refreshed. This is inferred
  from the code and has not been tested for large-deformation cases.
- **Update trigger only checks mesh A.** `changing()` and `moving()` look only
  at the region listed first in `coupledPatchPair`. If the static solid is
  listed first, the `fsiInterface` mapping is never rebuilt, whatever the
  value of `interpolatorUpdateFrequency`. The tutorials list the fluid first.
  This is inferred from the code and untested.
- **Residual normalisation depends on the coordinate origin.** $r_\Gamma$ is
  normalised by the magnitude of the absolute point positions. Moving the
  geometry away from the origin makes the same mismatch look smaller.
  Inferred from the code.
- **Total force uses the undeformed area.** The logged total force sums the
  traction times the undeformed face area $\lvert\mathbf{S}_0\rvert$. It is a
  diagnostic only.
- **`fsiInterfaceCoeffs` is required but unused**, and `correct()` is a no-op.
- **RBF transfer is unavailable** (it needs Eigen and is disabled), so only
  `AMI` and `directMap` can be used.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/fluidStructureInteraction/HronTurekFsi3` | Hron-Turek FSI3 benchmark (pimpleFluid + solidStVK), AMI transfer, coupling starts at t = 2 s |
| `tutorials/fluidStructureInteraction/beamInCrossFlow` | Elastic beam in cross-flow (pimpleFluid + solidStVK), directMap transfer |
| `tutorials/fluidStructureInteraction/heatedBeamInCrossFlow` | Heated elastic beam in cross-flow; `fsiInterface` + `heatTransferInterface` on one patch pair (pimpleFluid + transportTemperature / solidStVK + conductTemperature) |

All three are run with `./Allrun -p` (partitioned).

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
