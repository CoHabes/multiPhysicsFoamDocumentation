---
mathjax: true
title: "regionCoupledPressureFlux"
description: "Partitioned Neumann condition for the pressure on the side-B patch of a capillaryInterface: the normal pressure gradient balances the local normal acceleration."
tags: ["boundary-condition", "coupled", "multiphase"]
weight: 10
---

## Overview

`regionCoupledPressureFlux` is the Neumann (flux) half of the partitioned
pressure coupling across a sharp, mesh-tracked fluid-fluid interface. It acts
on a `volScalarField` (`p` or `pKin`) and derives from
`genericRegionCoupledFlux<scalar>`, i.e. from `fixedGradient`. Although it is
built on the generic flux-matching scheme, its `fluxJump()` cancels the
neighbour flux. The normal pressure gradient therefore depends only on the
local normal acceleration $\mathbf{n}\cdot\partial\mathbf{U}/\partial t$. The
partner condition on the other side is `regionCoupledPressureValue` (or
`regionCoupledMassTransferPressureValue`, *inferred*: there is no mass-transfer
pressure-flux BC). It is used with `capillaryInterface`.

## Mathematical Formulation

The BC sits on side $B$ in all tutorials. $\mathbf{n}_B$ is the outward patch
normal and $\tilde F = \mathcal{I}_{A\to B}(F_A)$ the interpolated neighbour
flux `snGrad(p_A)*k_A`. The generic flux BC sets
$k_B\,\partial p_B/\partial n = -\tilde F + \text{fluxJump}$ with

$$
\text{fluxJump} = \tilde F - \mathbf{n}_B\cdot\left(\frac{\partial\mathbf{U}_B}{\partial t}\right)_b
$$

so that

$$
k_B\,\frac{\partial p_B}{\partial n} = -\,\mathbf{n}_B\cdot\left(\frac{\partial\mathbf{U}_B}{\partial t}\right)_b
\quad \text{on } \Gamma_B
$$

This is the normal momentum balance at the interface without the convective
and viscous terms (the code that included them is commented out). The
gradient is relaxed by the acceleration model,
$g^{k} = g^{k-1} + \omega\,(g^{*} - g^{k-1})$, and the patch value is
reconstructed as $p_B = p_P + \mathbf{k}\cdot\nabla p_P + g/\Delta_f$
(first order) or $p_P + \mathbf{k}\cdot\nabla p_P + (g + \mathbf{n}\cdot\nabla p_P)/(2\Delta_f)$
(`secondOrder yes`); the $\mathbf{k}$ term only with `nonOrthCorr yes`.

The diffusivity $k_B$ comes from `diffCoeff()` with the hardcoded name `v`:

| Field name | Intended $k_B$ | Actual $k_B$ |
|------------|----------------|--------------|
| `pKin` | $1$ | $1$ |
| any other (e.g. `p`) | $1/\rho_B$ | $\rho_B$ (registered `volScalarField rho`, else `rho` from `transportProperties`) |

For `p` the imposed gradient is therefore
$-\mathbf{n}_B\cdot\partial_t\mathbf{U}_B/\rho_B$ instead of
$-\rho_B\,\mathbf{n}_B\cdot\partial_t\mathbf{U}_B$ (see
[Known Limitations](#known-limitations)).

| Symbol | Description | Unit |
|--------|-------------|------|
| $p$ | Pressure (`p`) or kinematic pressure (`pKin`) | Pa or m²/s² |
| $(\partial\mathbf{U}/\partial t)_b$ | Patch value of `fvc::ddt(U)` in the local region | m/s² |
| $k_B$ | Diffusivity from `diffCoeff()` | see table |
| $\Delta_f$ | Patch `deltaCoeffs` | 1/m |

## Implementation

`updateCoeffs()` is the generic one of `genericRegionCoupledFlux`:

1. Return if updated or `coupled no` (gradient frozen).
2. Call `update()` on the first `interfaceType`.
3. Look up the neighbour field by the local field name and `refCast` its
   patch field to `genericRegionCoupledJump<scalar>`; interpolate its
   `flux()` and negate it.
4. Add `fluxJump()`, divide by `diffCoeff()`, relax, call
   `fixedGradient::updateCoeffs()`.

`fluxJump()` repeats step 3 to get $\tilde F$ and returns
$\tilde F - \mathbf{n}_B\cdot(\partial_t\mathbf{U})_b$. It also reads `rho` and
`mu` from the local `transportProperties` (required to exist), but does not
use them. Unlike the other capillary BCs it never calls `capInterface()`, so
it does not check the interface type.

The `p` flux residual of the `pUPimple` DNA loop comes from this BC, because
the side-A pressure patch carries a jump BC. The neighbour flux cancels in it
as well:

$$
r_\text{raw} = \left| k_B\,\partial_n p_B + \mathbf{n}_B\cdot(\partial_t\mathbf{U}_B)_b \right|
$$

normalised as described in `genericRegionCoupledFlux`.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `regionCoupledPressureFlux` | yes | — |
| `interfaceType` | List of interface type names; the first entry is updated by this BC. Tutorials use `(capillaryInterface)` | yes | — |
| `gradient` | Initial patch-normal gradient. Read if present; otherwise zero gradient and the patch value is set from the internal field | no | zero |
| `value` | Not read (`fixedGradient` ignores it). Written by `write()` | no | — |
| `diffusivity` | Read, but ignored: the BC always uses `v` | no | empty |
| `neighbourFieldName` | Read by the coupling manager, not used by this BC | no | local field name |
| `accType` | Acceleration model for the gradient: `fixed`, `aitken`, `IQNILS`. A warning is printed if missing | no | `fixed` |
| `relax` | (Initial) relaxation factor | no | `1.0` |
| `couplingReuse` | `IQNILS` only | no | `0` |
| `nonOrthCorr` | Non-orthogonal correction in `evaluate()` | no | `false` |
| `secondOrder` | Second-order reconstruction of the patch value | no | `false` |

`write()` does not write `gradient`.

## Usage

`tutorials/multiphaseFlow/2dRisingBubble/0/fluidB/orig/partitioned/p`
(side B):

```
interfaceShadow
{
    type                    regionCoupledPressureFlux;
    diffusivity             "v";                    // ignored, always v
    interfaceType           (capillaryInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform 0;              // not read
}
```

Partner on side A, `0/fluidA/orig/partitioned/p`:

```
interface
{
    type                    regionCoupledPressureValue;
    interfaceType           (capillaryInterface);
    value                   uniform 0;
}
```

The interface and DNA entries are shown on the `regionCoupledVelocityValue`
page and in `capillaryInterface`.

## Known Limitations

- **Wrong diffusivity for `p` (bug).** `diffCoeff()` renames `v` to `rho`
  before it checks for `v`, so the reciprocal is never taken and
  $k_B = \rho_B$. For a field named `p` the imposed gradient is off by a
  factor $\rho_B^2$ and has the wrong units. `pKin` is not affected. The
  tutorials use `p`.
- **Simplified normal momentum balance.** Convective and viscous
  contributions to the normal pressure gradient are not included.
- **Neighbour flux has no effect.** The interpolated neighbour flux cancels.
  It is still computed, and the neighbour patch field must derive from
  `genericRegionCoupledJump<scalar>`, otherwise the `refCast` aborts.
- **No interface-type check.** Any interface type listed first is accepted;
  it is only `update()`d.
- **`diffusivity` and `neighbourFieldName` are ignored**; `gradient` is not
  written to the time directories.
- **No runnable tutorial.** The tutorials use the disabled `interTrackFluid`.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/multiphaseFlow/2dRisingBubble` | 2D rising gas bubble, `interTrackFluid` + `interTrackFluid` (region type currently disabled) |
| `tutorials/multiphaseFlow/3dDuineveldRisingBubble` | 3D rising bubble (Duineveld), same setup |
| `tutorials/multiphaseFlow/3dOscillatingBubble` | 3D oscillating bubble, same setup |

The BC is set in the `0/fluidB/orig/partitioned/p` variant of each tutorial.

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
