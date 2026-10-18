---
mathjax: true
title: "regionCoupledVelocityFlux"
description: "Partitioned Neumann condition for the velocity on the side-A patch of a capillaryInterface: tangential stress balance with surface tension and the normal gradient from surface continuity."
tags: ["boundary-condition", "coupled", "multiphase", "surface-tension"]
weight: 10
---

## Overview

`regionCoupledVelocityFlux` is the Neumann (flux) half of the partitioned
velocity coupling across a sharp, mesh-tracked fluid-fluid interface. It acts
on a `volVectorField` (the velocity `U`) and derives from
`genericRegionCoupledFlux<vector>`, i.e. from `fixedGradient`. It sets the
normal velocity gradient so that the tangential viscous stress balances the
neighbour stress plus the tangential surface-tension force. The normal
component of the gradient follows from surface continuity. The partner
condition on the other side is `regionCoupledVelocityValue`. Both belong to
`capillaryInterface`, and this BC aborts if the first entry of its
`interfaceType` list is not `capillaryInterface`.

## Mathematical Formulation

Region $A$ is the first entry of `coupledPatchPair` and carries the
finite-area mesh (`faMesh`) of the `capillaryInterface`. This BC sits on side
$A$ (see [Known Limitations](#known-limitations)). $\mathbf{n}$ is the
`faMesh` face normal (`faceAreaNormals()`), $\mathcal{I}_{B\to A}$ the
patch-to-patch interpolation, and $\mathbf{F}_B$ the flux computed by
`regionCoupledVelocityValue` on side B
($\mu_B$ times the tangential part of $\partial\mathbf{U}_B/\partial n_B$).
With $\tilde{\mathbf{F}} = \mathcal{I}_{B\to A}(\mathbf{F}_B)$, the gradient is

$$
\mu_A \frac{\partial \mathbf{U}_A}{\partial n} =
-\tilde{\mathbf{F}} + \mathbf{n}(\mathbf{n}\cdot\tilde{\mathbf{F}}) + \mathbf{f}_{\sigma,t} - \mu_A\,\mathbf{n}\,(\nabla_s\cdot\mathbf{u}_s) + (\mu_B - \mu_A)\,(\nabla_s\mathbf{u}_s)\cdot\mathbf{n}
\quad \text{on } \Gamma_A
$$

The first two terms are $-(\mathbf{I}-\mathbf{n}\mathbf{n})\cdot\tilde{\mathbf{F}}$.
The tangential surface-tension force is

$$
\mathbf{f}_{\sigma,t} = \sigma\,\frac{1}{S_f}\sum_{e} \mathbf{L}_e\,c_e \;-\; \sigma\kappa\,\mathbf{n}
$$

The normal term $-\mu_A\mathbf{n}\,\nabla_s\cdot\mathbf{u}_s$ uses
$\partial u_n/\partial n = -\nabla_s\cdot\mathbf{u}_s$ (incompressibility at
the surface). The gradient is relaxed by the acceleration model:
$g^{k} = g^{k-1} + \omega\,(g^{*} - g^{k-1})$.

The patch value is reconstructed from the gradient in `evaluate()`:

$$
\mathbf{U}_A = \mathbf{U}_P + \mathbf{k}\cdot\nabla\mathbf{U}_P + \frac{1}{\Delta_f}\,\frac{\partial \mathbf{U}_A}{\partial n}
\qquad\text{or}\qquad
\mathbf{U}_A = \mathbf{U}_P + \mathbf{k}\cdot\nabla\mathbf{U}_P + \frac{1}{2\Delta_f}\left(\frac{\partial \mathbf{U}_A}{\partial n} + \mathbf{n}\cdot\nabla\mathbf{U}_P\right)
$$

(first order / `secondOrder yes`; the $\mathbf{k}$ term only with `nonOrthCorr yes`).

| Symbol | Description | Unit |
|--------|-------------|------|
| $\mathbf{U}$ | Velocity | m/s |
| $\mathbf{u}_s$ | Surface velocity on the `faMesh` (`Us()`) | m/s |
| $\sigma$ | Surface tension (`sigma()` of the `capillaryInterface`) | N/m |
| $\kappa$ | Face curvature, `faMesh::faceCurvatures()` | 1/m |
| $\mathbf{L}_e$, $c_e$ | Edge length vectors and `edgeLengthCorrection()` of the `faMesh` | m, – |
| $\mu_A$, $\mu_B$ | `mu` from each region's `transportProperties` | Pa·s |
| $\Delta_f$ | Patch `deltaCoeffs` | 1/m |
| $\mathbf{k}$ | Non-orthogonal part of the face-to-cell vector | m |

## Implementation

`updateCoeffs()` is the generic one of `genericRegionCoupledFlux`:

1. Return if already updated or if the interface has `coupled no` (the
   gradient is then frozen).
2. Call `update()` on the first `interfaceType` (refreshes curvature and
   $\mathbf{u}_s$ of the `capillaryInterface`).
3. Look up the neighbour `U` by the local field name and `refCast` its patch
   field to `genericRegionCoupledJump<vector>`. Its `flux()` is interpolated
   onto this patch and negated.
4. Add `fluxJump()` and divide by `k = diffCoeff()`. `diffCoeff()` uses the
   hardcoded name `mu`: the patch values of a registered `volScalarField mu`,
   otherwise the scalar `mu` from `transportProperties`.
5. Relax the gradient with the acceleration model and call
   `fixedGradient::updateCoeffs()`.

`fluxJump()` computes $\nabla_s\cdot\mathbf{u}_s$ (`fac::div`),
$\nabla_s\mathbf{u}_s$ (`fac::grad`) and the edge-integrated surface-tension
force on the `faMesh`. $\mu_A$ and $\mu_B$ are always read from the
`transportProperties` dictionaries, not from fields. If a region has a
`pKin` field, its `mu` must equal `rho*nu` from the same dictionary
(exact floating-point comparison), otherwise the run stops with a fatal
error. Each call prints `Us <value>` (the volume swept by $\mathbf{u}_s$ in
one time step) to the log.

`evaluate()` uses the registered `grad(U)` if it exists, otherwise
`fvc::grad(U)`.

The `U` flux residual of the `pUPimple` DNA loop is taken from this BC
(side A is checked first):

$$
r = \frac{\lVert \mu_A\,\partial_n\mathbf{U}_A + \tilde{\mathbf{F}} - \text{fluxJump}\rVert_2}
{\max\left(\min\left(\lVert\mu_A\,\partial_n\mathbf{U}_A\rVert_2,\ \lVert -\tilde{\mathbf{F}}+\text{fluxJump}\rVert_2\right),\ \text{SMALL}\right)}
$$

with $\partial_n\mathbf{U}_A$ from `snGrad()` of the patch.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `regionCoupledVelocityFlux` | yes | — |
| `interfaceType` | List of interface type names. The first entry must be `capillaryInterface` | yes | — |
| `gradient` | Initial patch-normal gradient. Read if present; otherwise zero gradient and the patch value is set from the internal field | no | zero |
| `value` | Not read (`fixedGradient` ignores it). Written by `write()` | no | — |
| `diffusivity` | Read, but ignored: the BC always uses `mu` | no | empty |
| `neighbourFieldName` | Read by the coupling manager, not used by this BC | no | local field name |
| `accType` | Acceleration model for the gradient: `fixed`, `aitken`, `IQNILS`. A warning is printed if missing | no | `fixed` |
| `relax` | (Initial) relaxation factor | no | `1.0` |
| `couplingReuse` | `IQNILS` only | no | `0` |
| `nonOrthCorr` | Non-orthogonal correction in `evaluate()` | no | `false` |
| `secondOrder` | Second-order reconstruction of the patch value | no | `false` |

`write()` does not write `gradient`. After a restart or `decomposePar` the
gradient therefore starts from zero.

## Usage

`tutorials/multiphaseFlow/2dRisingBubble/0/fluidA/orig/partitioned/U`
(side A):

```
interface
{
    type                    regionCoupledVelocityFlux;
    diffusivity             "mu";                   // ignored, always mu
    interfaceType           (capillaryInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform (0 0 0);        // not read
}
```

Partner on side B, `0/fluidB/orig/partitioned/U`:

```
interfaceShadow
{
    type                    regionCoupledVelocityValue;
    interfaceType           (capillaryInterface);
    value                   uniform (0 0 0);
}
```

The `regionInterfaceProperties` and `multiRegionProperties` entries are shown
on the `regionCoupledVelocityValue` page and in `capillaryInterface`.

## Known Limitations

- **`capillaryInterface` only.** `capInterface()` aborts unless the first
  `interfaceType` entry is `capillaryInterface`.
- **Must sit on side A** (*inferred*). The `faMesh` normals,
  $\mathbf{u}_s$, $\sigma$ and $\kappa$ are used index by index on the patch
  without interpolation, which is only valid on the side that owns the
  `faMesh`.
- **Partner must be jump-derived.** The neighbour `U` patch field is
  `refCast` to `genericRegionCoupledJump<vector>`; any other type aborts.
- **Exact `mu == rho*nu` check** for regions with `pKin`. Rounding in the
  input values can trigger the fatal error.
- **`diffusivity` and `neighbourFieldName` are ignored.**
- **Relaxation history starts from the patch value** (*inferred*). The base
  class initialises the acceleration model with the patch values but relaxes
  the gradient. With `relax 1` (the default, used in all tutorials) this has
  no effect.
- **`gradient` is not written** to the time directories.
- **No runnable tutorial.** The tutorials use the disabled `interTrackFluid`.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/multiphaseFlow/2dRisingBubble` | 2D rising gas bubble, `interTrackFluid` + `interTrackFluid` (region type currently disabled) |
| `tutorials/multiphaseFlow/3dDuineveldRisingBubble` | 3D rising bubble (Duineveld), same setup |
| `tutorials/multiphaseFlow/3dOscillatingBubble` | 3D oscillating bubble, same setup |

The BC is set in the `0/fluidA/orig/partitioned/U` variant of each tutorial.

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
