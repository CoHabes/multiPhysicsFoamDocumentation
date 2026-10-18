---
mathjax: true
title: "regionCoupledVelocityValue"
description: "Partitioned Dirichlet condition for the velocity on the side-B patch of a capillaryInterface: neighbour velocity with the normal component replaced by the local mesh normal velocity."
tags: ["boundary-condition", "coupled", "multiphase", "interface-tracking"]
weight: 10
---

## Overview

`regionCoupledVelocityValue` is the Dirichlet (jump) half of the partitioned
velocity coupling across a sharp, mesh-tracked fluid-fluid interface. It acts
on a `volVectorField` (the velocity `U`) and derives from
`genericRegionCoupledJump<vector>`, i.e. from `fixedValue`. The patch value is
the neighbour velocity, interpolated onto this patch, with its normal
component replaced by the normal velocity of the local mesh. The tangential
velocity is therefore continuous. The partner condition on the other side is
`regionCoupledVelocityFlux`. Both belong to `capillaryInterface`, and this BC
aborts if the first entry of its `interfaceType` list is not
`capillaryInterface`. The pressure is coupled by `regionCoupledPressureValue`
and `regionCoupledPressureFlux` in the same `pUPimple` DNA loop.

## Mathematical Formulation

Region $A$ is the first entry of `coupledPatchPair` and carries the
finite-area mesh of the `capillaryInterface`. This BC sits on side $B$ (see
[Known Limitations](#known-limitations)). $\mathcal{I}_{A\to B}$ is the
patch-to-patch interpolation (`interfaceTransferMethod`), $\mathbf{n}_B$ the
outward unit normal of the local patch.

$$
\mathbf{U}_B^{*} = \mathcal{I}_{A\to B}(\mathbf{U}_A) - \mathbf{n}_B\left(\mathbf{n}_B\cdot\mathcal{I}_{A\to B}(\mathbf{u}_s)\right) + \mathbf{n}_B\,\frac{\phi_{\text{mesh},B}}{|\mathbf{S}_f|}
\quad \text{on } \Gamma_B
$$

$\mathbf{U}_B^{*}$ is then passed through the acceleration model:

$$
\mathbf{U}_B^{k} = \mathbf{U}_B^{k-1} + \omega\left(\mathbf{U}_B^{*} - \mathbf{U}_B^{k-1}\right)
$$

with $\omega$ = `relax` (`fixed`), the Aitken factor (`aitken`), or an IQN-ILS
update (`IQNILS`).

The surface velocity $\mathbf{u}_s$ of the `capillaryInterface` is
$\mathbf{U}_A - (\mathbf{U}_A\cdot\mathbf{n}_A)\mathbf{n}_A + w_A\mathbf{n}_A$,
where $w_A$ is the side-A mesh normal velocity. The term
$-\mathbf{n}_B(\mathbf{n}_B\cdot\mathbf{u}_s)$ therefore removes the normal part
of $\mathbf{U}_A$ up to the side-A mesh motion. The last term imposes the side-B
mesh normal velocity.

The diffusive momentum flux handed to `regionCoupledVelocityFlux` is

$$
\mathbf{F}_B = \mu_B\,(\mathbf{I}-\mathbf{n}_B\mathbf{n}_B)\cdot\frac{\partial\mathbf{U}_B}{\partial n},
\qquad
\frac{\partial\mathbf{U}_B}{\partial n} = \Delta_f\left(\mathbf{U}_B - \mathbf{U}_P - \mathbf{k}\cdot\nabla\mathbf{U}_P\right)
$$

The non-orthogonal term $\mathbf{k}\cdot\nabla\mathbf{U}_P$ is only present with
`nonOrthCorr yes`. With `secondOrder yes` the gradient is
$2\Delta_f(\mathbf{U}_B - \mathbf{U}_P - \mathbf{k}\cdot\nabla\mathbf{U}_P) - \mathbf{n}_B\cdot\nabla\mathbf{U}_P$.
In both cases the normal component is removed.

| Symbol | Description | Unit |
|--------|-------------|------|
| $\mathbf{U}$ | Velocity | m/s |
| $\mathbf{u}_s$ | Surface velocity of the `capillaryInterface` (`Us()`, on side A) | m/s |
| $\phi_{\text{mesh},B}$ | Boundary value of `fvc::meshPhi(rho, U)` if a `volScalarField rho` is registered in region B, otherwise `fvc::meshPhi(U)` | m³/s |
| $\mathbf{S}_f$ | Face area vector of the local patch | m² |
| $\mu_B$ | Diffusivity `mu` of region B (field or `transportProperties`) | Pa·s |
| $\Delta_f$ | Patch `deltaCoeffs` | 1/m |
| $\mathbf{k}$ | Non-orthogonal part of the face-to-cell vector, $\mathbf{d}-\mathbf{n}(\mathbf{n}\cdot\mathbf{d})$ | m |
| $\omega$ | Relaxation factor from the acceleration model | – |

## Implementation

`updateCoeffs()` is the generic one of `genericRegionCoupledJump`. Each call:

1. Returns immediately if the patch is already updated or the interface has
   `coupled no`. In that case the patch keeps its last value.
2. Calls `update()` on the interface named first in `interfaceType`. For
   `capillaryInterface` this refreshes the curvature and the surface velocity
   $\mathbf{u}_s$ (`correct()`).
3. Looks up the neighbour field with the **local** field name (`U`),
   interpolates its patch values onto this patch and adds `valueJump()`.
   `valueJump()` returns
   $-\mathbf{n}_B(\mathbf{n}_B\cdot\mathcal{I}_{A\to B}(\mathbf{u}_s)) + \mathbf{n}_B\,\phi_{\text{mesh},B}/|\mathbf{S}_f|$.
4. Assigns the result, relaxes it with the acceleration model and calls
   `fixedValue::updateCoeffs()`.

The matrix coefficients come from the vector specialisations in
`genericRegionCoupledJump`. `gradientInternalCoeffs()` is the `fixedValue` one
($-\Delta_f$). `gradientBoundaryCoeffs()` subtracts the normal part of the
current face gradient explicitly, using the current cell values. At
convergence only the tangential gradient enters the viscous term of the
momentum equation. `snGrad()` drops the normal part in the same way. If a
`volTensorField` named `grad(U)` is registered, it is used for the
non-orthogonal and second-order terms; otherwise `fvc::grad(U)` is computed on
every call.

`flux()` returns `snGrad()*k`. `k` comes from `diffCoeff()` with the hardcoded
diffusivity name `mu`: the patch values of a registered `volScalarField mu`,
otherwise the scalar `mu` from the region's registered `transportProperties`
dictionary.

The DNA loop for `pUPimple` (`multiRegionSystem::solve()`) runs
`solvePIMPLE()` on all regions, in which this BC is evaluated, then
`meshMotionCorrector()` and `interfaces_->update()`. The jump residual for `U`
is taken from this BC, because the side-A `U` patch carries a flux BC:

$$
r = \frac{\lVert \mathbf{U}_B - \mathcal{I}_{A\to B}(\mathbf{U}_A) - \text{valueJump}\rVert_2}
{\max\left(\min\left(\lVert\mathbf{U}_B\rVert_2,\ \lVert\mathcal{I}_{A\to B}(\mathbf{U}_A)+\text{valueJump}\rVert_2\right),\ \text{SMALL}\right)}
$$

The class also defines `updatePhi()`, which would overwrite the boundary
`phi` with the negated interpolated neighbour flux. The call in the base class
is commented out, so it is never executed.

`correctClosedVolumePhi()` (used by `icoFluid` and `interTrackFluid` when
`closedVolume` is on in `fvSolution`) looks for the patch carrying this BC (or
`regionCoupledMassTransferVelocityValue`). It sets `phi = U·Sf` there and
removes the net flux, weighted by $|\phi|$, so the net volume flux through the
interface is zero.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `regionCoupledVelocityValue` | yes | — |
| `interfaceType` | List of interface type names. The first entry is the interface this BC belongs to and must be `capillaryInterface`. The neighbour region and patch are taken from the matching interface entry in `regionInterfaceProperties` | yes | — |
| `value` | Initial patch value. Read if present. If missing, the patch values are not initialised | in practice yes | — |
| `diffusivity` | Read, but ignored: the BC always uses `mu` | no | empty |
| `neighbourFieldName` | Read by the coupling manager. Not used by this BC, which looks up the neighbour field by the local name | no | local field name |
| `jumpCoefficient` | Read by the base class, not used, not written back | no | empty |
| `accType` | Acceleration model: `fixed`, `aitken` or `IQNILS`. A warning is printed if missing | no | `fixed` |
| `relax` | (Initial) relaxation factor | no | `1.0` |
| `couplingReuse` | `IQNILS` only: number of previous time steps whose coupling data are reused | no | `0` |
| `nonOrthCorr` | Non-orthogonal correction in `snGrad()` / matrix coefficients | no | `false` |
| `secondOrder` | Second-order face-normal gradient | no | `false` |

`write()` writes `interfaceType`, the resolved `neighbourFieldName`,
`diffusivity mu`, `accType`, `relax` (and `couplingReuse` for `IQNILS`),
`nonOrthCorr`, `secondOrder` and `value`.

## Usage

`tutorials/multiphaseFlow/2dRisingBubble/0/fluidB/orig/partitioned/U`
(side B):

```
interfaceShadow
{
    type                    regionCoupledVelocityValue;
    diffusivity             "mu";                   // ignored, always mu
    interfaceType           (capillaryInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform (0 0 0);
}
```

Partner on side A, `0/fluidA/orig/partitioned/U`:

```
interface
{
    type                    regionCoupledVelocityFlux;
    diffusivity             "mu";
    interfaceType           (capillaryInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform (0 0 0);
}
```

`constant/regionInterfaceProperties`:

```
partitionedCoupledPatches
(
    fluidAfluidB
    {
        interfaceType (capillaryInterface);

        coupledPatchPair
        (
            (fluidA interface)          // side A: faMesh, VelocityFlux, PressureValue
            (fluidB interfaceShadow)    // side B: VelocityValue, PressureFlux
        );

        coupledFields
        (
            pUPimple
        );

        capillaryInterfaceCoeffs
        {
            sigma       sigma [ 1 0 -2 0 0 0 0 ] 0.07232;
        }
    }
);
```

`constant/multiRegionProperties` (`U` and `p` must exist in both regions):

```
DNA
{
    pUPimple
    {
        maxCoupleIter 10;

        residualControl
        {
            U { maxJumpRes 1e-5; maxFluxRes 1e-4; }
            p { maxJumpRes 1e-5; maxFluxRes 1e-5; }
        }
    }
}
```

## Known Limitations

- **`capillaryInterface` only.** `capInterface()` aborts if the first
  `interfaceType` entry is not `capillaryInterface`.
- **Must sit on side B** (*inferred*). $\mathbf{u}_s$ lives on the side-A
  finite-area mesh and is interpolated "from the neighbour". On side A, the
  side-A field would be mapped as if it belonged to side B.
- **`diffusivity` is ignored.** `diffusivityName()` is hardcoded to `mu`.
- **`neighbourFieldName` has no effect.** The neighbour velocity is looked up
  by the local field name.
- **`value` should be given.** The base class reads it only if present and
  does not initialise the patch otherwise.
- **`updatePhi()` is dead code**, so the boundary `phi` is not set from the
  neighbour by this BC.
- **The header example omits `interfaceType`**, which is required.
- **No runnable tutorial.** All tutorials that use this BC run on
  `interTrackFluid`, which is commented out in
  `src/multiRegionSystem/Make/files`.
- **Partitioned only.** There is no monolithic variant.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/multiphaseFlow/2dRisingBubble` | 2D rising gas bubble, `interTrackFluid` + `interTrackFluid` (region type currently disabled) |
| `tutorials/multiphaseFlow/3dDuineveldRisingBubble` | 3D rising bubble (Duineveld), same setup |
| `tutorials/multiphaseFlow/3dOscillatingBubble` | 3D oscillating bubble, same setup |

The BC is set in the `0/fluidB/orig/partitioned/U` variant of each tutorial.

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
