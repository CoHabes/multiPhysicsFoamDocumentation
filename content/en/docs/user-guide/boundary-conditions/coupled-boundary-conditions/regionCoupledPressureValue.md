---
mathjax: true
title: "regionCoupledPressureValue"
description: "Partitioned Dirichlet condition for the pressure on the side-A patch of a capillaryInterface: neighbour pressure plus the capillary, viscous and hydrostatic normal-stress jump."
tags: ["boundary-condition", "coupled", "multiphase", "surface-tension"]
weight: 10
---

## Overview

`regionCoupledPressureValue` is the Dirichlet (jump) half of the partitioned
pressure coupling across a sharp, mesh-tracked fluid-fluid interface. It acts
on a `volScalarField` (`p`, or the kinematic pressure `pKin`) and derives from
`genericRegionCoupledJump<scalar>`, i.e. from `fixedValue`. The patch value is
the interpolated neighbour pressure plus the normal-stress jump: capillary
pressure $\sigma\kappa$, the viscous normal-stress jump and a hydrostatic
correction. The partner condition on the other side is
`regionCoupledPressureFlux`. Both belong to `capillaryInterface`, and this BC
aborts if the first entry of its `interfaceType` list is not
`capillaryInterface`. The variant with phase change is
`regionCoupledMassTransferPressureValue`.

## Mathematical Formulation

Region $A$ is the first entry of `coupledPatchPair` and carries the `faMesh`
of the `capillaryInterface`. This BC sits on side $A$ (see
[Known Limitations](#known-limitations)). $\mathcal{I}_{B\to A}$ is the
patch-to-patch interpolation.

For a field named `p`:

$$
p_A^{*} = \mathcal{I}_{B\to A}(p_B) + 2(\mu_B - \mu_A)\,\nabla_s\cdot\mathbf{u}_s - \sigma\kappa + (\rho_B - \rho_A)\,(\mathbf{x}_f - \mathbf{x}_\text{ref})\cdot\mathbf{g}
\quad \text{on } \Gamma_A
$$

For a field named `pKin` the same jump is applied in kinematic form:

$$
p_{\text{kin},A}^{*} = \frac{1}{\rho_A}\Big[\rho_B\,\mathcal{I}_{B\to A}(p_{\text{kin},B}) + 2(\mu_B - \mu_A)\,\nabla_s\cdot\mathbf{u}_s - \sigma\kappa + (\rho_B - \rho_A)\,(\mathbf{x}_f - \mathbf{x}_\text{ref})\cdot\mathbf{g}\Big]
$$

The viscous term is the jump of $2\mu\,\partial u_n/\partial n$ with
$\partial u_n/\partial n = -\nabla_s\cdot\mathbf{u}_s$. The gravity term implies
that `p` excludes the hydrostatic part
$\rho\,\mathbf{g}\cdot(\mathbf{x}-\mathbf{x}_\text{ref})$ (*inferred from the
formula*). The value is then relaxed:
$p_A^{k} = p_A^{k-1} + \omega\,(p_A^{*} - p_A^{k-1})$.

| Symbol | Description | Unit |
|--------|-------------|------|
| $p$, $p_\text{kin}$ | Pressure / kinematic pressure | Pa / m²/s² |
| $\mathbf{u}_s$ | Surface velocity on the `faMesh` (`Us()`) | m/s |
| $\sigma$ | Surface tension (`sigma()` of the `capillaryInterface`) | N/m |
| $\kappa$ | Face curvature, `faMesh::faceCurvatures()` | 1/m |
| $\mu_A$, $\mu_B$ | `mu` from each region's `transportProperties` | Pa·s |
| $\rho_A$, $\rho_B$ | `rho` from each region's `transportProperties` | kg/m³ |
| $\mathbf{g}$ | `uniformDimensionedVectorField g` of region A | m/s² |
| $\mathbf{x}_f$ | Patch face centres | m |
| $\mathbf{x}_\text{ref}$ | `PISO.pRefPoint` from region A's `fvSolution` | m |

## Implementation

`updateCoeffs()` is the generic one of `genericRegionCoupledJump`: return if
updated or `coupled no`; call `update()` on the first `interfaceType`
(curvature and $\mathbf{u}_s$); interpolate the neighbour field (looked up by
the local field name) onto this patch; add `valueJump()`; relax; call
`fixedValue::updateCoeffs()`.

`valueJump()`:

- takes $\sigma$, $\kappa$ and $\mathbf{u}_s$ from the `capillaryInterface`
  and computes $\nabla_s\cdot\mathbf{u}_s$ with `fac::div`;
- reads `pRefPoint` from the `PISO` sub-dictionary of region A's
  `fvSolution` and looks up `g` in region A (both required);
- reads `mu` and `rho` of both regions from their registered
  `transportProperties` dictionaries;
- for `pKin` returns the bracket above divided by $\rho_A$, minus the
  interpolated neighbour value, so that the base class sum gives the
  kinematic formula.

It also computes the mesh flux and the patch normals, but does not use them.

The jump residual for `p` in the `pUPimple` DNA loop is taken from this BC,
because side A is checked first:
$\lVert p_A - \mathcal{I}_{B\to A}(p_B) - \text{valueJump}\rVert_2$, normalised as
described in `genericRegionCoupledJump`.

`flux()` (inherited) returns `snGrad()*k` with $k$ from `diffCoeff()` and the
hardcoded diffusivity name `v`. `diffCoeff()` maps `v` to `rho`, but the step
that should take the reciprocal is never reached, so $k = \rho_A$ instead of
$1/\rho_A$ ($k = 1$ for `pKin`). The partner `regionCoupledPressureFlux`
cancels this flux, so the error does not reach the solution.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `regionCoupledPressureValue` | yes | — |
| `interfaceType` | List of interface type names. The first entry must be `capillaryInterface` | yes | — |
| `value` | Initial patch value. Read if present. If missing, the patch values are not initialised | in practice yes | — |
| `diffusivity` | Read, but ignored: the BC always uses `v` | no | empty |
| `neighbourFieldName` | Read by the coupling manager, not used by this BC | no | local field name |
| `jumpCoefficient` | Read by the base class, not used, not written back | no | empty |
| `accType` | Acceleration model: `fixed`, `aitken`, `IQNILS`. A warning is printed if missing | no | `fixed` |
| `relax` | (Initial) relaxation factor | no | `1.0` |
| `couplingReuse` | `IQNILS` only | no | `0` |
| `nonOrthCorr` | Non-orthogonal correction in `snGrad()` / matrix coefficients | no | `false` |
| `secondOrder` | Second-order face-normal gradient | no | `false` |

Entries read from other files:

| Entry | File | Required |
|-------|------|----------|
| `PISO { pRefPoint (x y z); }` | `system/<regionA>/fvSolution` | yes |
| `g` | `constant/<regionA>/g` (registered `uniformDimensionedVectorField`) | yes |
| `mu`, `rho` | `constant/<region>/transportProperties`, both regions | yes |

## Usage

`tutorials/multiphaseFlow/2dRisingBubble/0/fluidA/orig/partitioned/p`
(side A):

```
interface
{
    type                    regionCoupledPressureValue;
    diffusivity             "v";                    // ignored, always v
    interfaceType           (capillaryInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform 0;
}
```

Partner on side B, `0/fluidB/orig/partitioned/p`:

```
interfaceShadow
{
    type                    regionCoupledPressureFlux;
    diffusivity             "v";
    interfaceType           (capillaryInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform 0;
}
```

`system/fluidA/fvSolution` must contain `PISO { ... pRefPoint (0 0 0); }`.
The interface and DNA entries are shown on the `regionCoupledVelocityValue`
page and in `capillaryInterface`.

## Known Limitations

- **`capillaryInterface` only.** `capInterface()` aborts unless the first
  `interfaceType` entry is `capillaryInterface`.
- **Must sit on side A** (*inferred*). $\sigma$, $\kappa$ and
  $\nabla_s\cdot\mathbf{u}_s$ are taken from the side-A `faMesh` index by
  index, without interpolation.
- **Constant properties.** $\mu$ and $\rho$ are single scalars from
  `transportProperties`; registered `mu` or `rho` fields are not used.
- **`PISO.pRefPoint` is required**, even if the region solver uses a
  `PIMPLE` dictionary.
- **`diffCoeff()` returns $\rho$ instead of $1/\rho$** for the diffusivity
  `v`. Without consequence for this BC (see [Implementation](#implementation)).
- **`diffusivity` and `neighbourFieldName` are ignored.**
- **The header example omits `interfaceType`**, which is required.
- **No runnable tutorial.** The tutorials use the disabled `interTrackFluid`.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/multiphaseFlow/2dRisingBubble` | 2D rising gas bubble, `interTrackFluid` + `interTrackFluid` (region type currently disabled) |
| `tutorials/multiphaseFlow/3dDuineveldRisingBubble` | 3D rising bubble (Duineveld), same setup |
| `tutorials/multiphaseFlow/3dOscillatingBubble` | 3D oscillating bubble, same setup |

The BC is set in the `0/fluidA/orig/partitioned/p` variant of each tutorial.

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
