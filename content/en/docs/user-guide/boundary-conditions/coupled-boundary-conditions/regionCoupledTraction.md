---
mathjax: true
title: "regionCoupledTraction"
description: "Partitioned traction condition for the solid displacement D on an fsiInterface: maps the fluid traction to the solid and converts it into a fixed displacement gradient."
tags: ["boundary-condition", "coupled", "fsi", "solid-mechanics"]
weight: 10
---

## Overview

`regionCoupledTraction` applies the fluid load to the solid in a partitioned
fluid-structure interaction. It acts on the solid displacement `D`
(`volVectorField`) and derives from `genericRegionCoupledFlux<vector>`, i.e.
from `fixedGradient`. It uses the "direct flux" path of the base class: the
gradient is computed entirely from the fluid traction and the solid stress
state, with no flux matching against a partner BC. There is no partner
coupled BC on the fluid side, because the fluid has no `D` field. The fluid
side follows the solid through the `interfaceTrackingFvMesh` mesh motion, and
its velocity usually uses `movingWallVelocity`. The BC belongs to
`fsiInterface` and aborts if the first entry of its `interfaceType` list is
not `fsiInterface`. The solid region must be `solidStVK`.

## Mathematical Formulation

The fluid traction is evaluated on the fluid patch, mapped to the solid patch
and its sign is flipped:

$$
\mathbf{t} = -\,\mathcal{I}_{f\to s}\!\left(\boldsymbol{\sigma}_f\cdot\mathbf{n}_f\right)
$$

The normal displacement gradient is then

$$
\frac{\partial \mathbf{D}}{\partial n} =
\frac{1}{K}\left[\,\mathbf{t} - \mathbf{n}\cdot\boldsymbol{\sigma}_s + K\,\mathbf{n}_0\cdot\nabla\mathbf{D}\,\right],
\qquad
\mathbf{n} = \frac{\mathbf{F}^{-T}\cdot\mathbf{n}_0}{\lvert\mathbf{F}^{-T}\cdot\mathbf{n}_0\rvert}
\quad \text{on } \Gamma_s
$$

This is the usual segregated traction condition: the explicit stress
$\mathbf{n}\cdot\boldsymbol{\sigma}_s$ is replaced by the implicit
$K\,\partial\mathbf{D}/\partial n$ part. At convergence of the solid Picard
loop $\mathbf{n}\cdot\boldsymbol{\sigma}_s = \mathbf{t}$ holds in the deformed
configuration. The gradient is relaxed by the acceleration model,
$g^{k} = g^{k-1} + \omega\,(g^{*} - g^{k-1})$, and the patch value is
reconstructed as $\mathbf{D}_b = \mathbf{D}_P + \mathbf{k}\cdot\nabla\mathbf{D}_P + g/\Delta_f$
(first order) or $\mathbf{D}_P + \mathbf{k}\cdot\nabla\mathbf{D}_P + (g + \mathbf{n}_0\cdot\nabla\mathbf{D}_P)/(2\Delta_f)$
(`secondOrder yes`); the $\mathbf{k}$ term only with `nonOrthCorr yes`.

| Symbol | Description | Unit |
|--------|-------------|------|
| $\mathbf{D}$ | Solid total displacement | m |
| $\boldsymbol{\sigma}_f$ | Fluid Cauchy stress, `volSymmTensorField sigma` of the fluid region | Pa |
| $\mathbf{n}_f$ | Outward unit normal of the fluid patch | – |
| $\boldsymbol{\sigma}_s$ | Solid Cauchy stress (`sigma` of the solid) | Pa |
| $\mathbf{n}_0$ | Undeformed solid patch normal | – |
| $\mathbf{n}$ | Deformed solid patch normal | – |
| $\mathbf{F}^{-1}$ | Inverse deformation gradient (`Finv`) | – |
| $K$ | Implicit stiffness (`impK`); $1/K$ = `rImpK` | Pa |
| $\nabla\mathbf{D}$ | Displacement gradient (`gradD`) | – |

## Implementation

`updateCoeffs()` is the generic one of `genericRegionCoupledFlux` with
`useDirectFlux()` returning `true`:

1. Return if updated or if the interface has `coupled no`. The gradient is
   then frozen at its last value. With `couplingStartTime`, the coupling is
   switched on once that time is reached.
2. Call `update()` on the first `interfaceType`, i.e. the `fsiInterface`.
3. Set `gradient() = flux()`, relax it, call `fixedGradient::updateCoeffs()`.

`flux()`:

- checks that the first `interfaceType` is `fsiInterface`;
- looks up the fluid field named `sigma` (hardcoded) in the neighbour region,
  takes its patch values on the neighbour patch, dots them with the
  neighbour patch normals, interpolates the result to the solid patch and
  negates it;
- looks up `impK`, `rImpK`, `gradD`, `sigma` and `Finv` on the solid patch
  (all provided by `solidStVK`);
- calls `fsiInterface::setTotalForce()` with
  $\sum \mathbf{t}\,\lvert\mathbf{S}_0\rvert$ (undeformed face areas). The
  interface logs this value in `info()`.

`evaluate()` uses the registered `grad(D)` if one exists. `solidStVK`
registers its gradient as `gradD`, so in practice `fvc::grad(D)` is computed
on every evaluation.

In the DNA loop for `D` (`multiRegionSystem::solve()`), each iteration runs
`solvePIMPLE()` on the fluid (updating the fluid `sigma`), then assembles and
solves `D` with this BC, then `meshMotionCorrector()` and
`interfaces_->update()`. Convergence is judged by the `fsiInterface`
interface residual (`maxIntRes`). The jump and flux residuals of the generic
base classes cannot be used: `rawResidual()`/`normResidual()` abort for a
direct-flux BC, and the fluid has no `D` field anyway.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `regionCoupledTraction` | yes | — |
| `interfaceType` | List of interface type names. The first entry must be `fsiInterface` | yes | — |
| `neighbourFieldName` | Read by the coupling manager, but ignored: the fluid field is always looked up as `sigma` | no | local field name (`D`) |
| `gradient` | Initial displacement gradient. Read if present; otherwise zero gradient and the patch value is set from the internal field | no | zero |
| `value` | Not read (`fixedGradient` ignores it). Written by `write()` | no | — |
| `diffusivity` | Read by the base class, not used (direct flux) | no | empty |
| `accType` | Acceleration model for the gradient: `fixed`, `aitken`, `IQNILS`. A warning is printed if missing | no | `fixed` |
| `relax` | (Initial) relaxation factor | no | `1.0` |
| `couplingReuse` | `IQNILS` only | no | `0` |
| `nonOrthCorr` | Non-orthogonal correction in `evaluate()` | no | `false` |
| `secondOrder` | Second-order reconstruction of the patch value | no | `false` |

`write()` writes `diffusivity` only if it was set, and does not write
`gradient`.

## Usage

Solid `0/solid/D` from
`tutorials/fluidStructureInteraction/HronTurekFsi3/0/solid/orig/partitioned/D`:

```
interface
{
    type                    regionCoupledTraction;
    neighbourFieldName      sigma;          // ignored, sigma is hardcoded
    interfaceType           (fsiInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform (0 0 0); // not read
}
```

The matching fluid patch has no coupled BC. In the tutorials it uses
`movingWallVelocity` for `U`, `zeroGradient` for `pKin`, and `fixedValue` for
`pointMotionU` (driven by `interfaceTrackingFvMesh` with `fsiInterface yes`).

`constant/regionInterfaceProperties`:

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

        fsiInterfaceCoeffs
        {}
    }
);
```

`constant/multiRegionProperties`:

```
DNA
{
    D
    {
        maxCoupleIter 100;

        residualControl
        {
            maxIntRes 1e-7;
        }
    }
}
```

The full setup (dynamic mesh, transfer method, coupling start time) is
described in `fsiInterface`.

## Known Limitations

- **`fsiInterface` and `solidStVK` only.** The interface type is checked;
  the solid fields `impK`, `rImpK`, `gradD`, `sigma` and `Finv` are only
  provided by `solidStVK`.
- **The fluid stress must be called `sigma`.** `neighbourFieldName` is read
  but not used.
- **No generic residuals.** Use the `fsiInterface` interface residual
  (`maxIntRes`); jump/flux residuals abort for this BC.
- **Relaxation history starts from the patch value** (*inferred*). The base
  class initialises the acceleration model with the patch values of `D` but
  relaxes the gradient. With `relax 1` (the default, used in all tutorials)
  this has no effect.
- **Total force uses the undeformed area** and is set on every call of
  `flux()`. It is a diagnostic only.
- **`gradient` is not written** to the time directories, so a restart starts
  from zero gradient until the first update.
- **Partitioned only.** There is no monolithic FSI path.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/fluidStructureInteraction/HronTurekFsi3` | Hron-Turek FSI3 benchmark (`pimpleFluid` + `solidStVK`), AMI transfer, coupling starts at t = 2 s |
| `tutorials/fluidStructureInteraction/beamInCrossFlow` | Elastic beam in cross-flow (`pimpleFluid` + `solidStVK`), directMap transfer |
| `tutorials/fluidStructureInteraction/heatedBeamInCrossFlow` | Heated elastic beam; `fsiInterface` + `heatTransferInterface` on one patch pair |

The BC is set in the `0/solid/orig/partitioned/D` variant of each tutorial
(`./Allrun -p`).

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
