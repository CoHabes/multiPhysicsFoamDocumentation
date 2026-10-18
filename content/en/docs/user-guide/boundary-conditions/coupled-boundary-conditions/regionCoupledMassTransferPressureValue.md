---
mathjax: true
title: "regionCoupledMassTransferPressureValue"
description: "Partitioned Dirichlet condition for the pressure at a phase-change interface: capillary pressure jump of regionCoupledPressureValue plus the vapour recoil term."
tags: ["boundary-condition", "coupled", "phase-change", "multiphase"]
weight: 10
---

## Overview

`regionCoupledMassTransferPressureValue` is the Dirichlet (jump) half of the
partitioned pressure coupling across a sharp liquid-vapour interface with
phase change. It acts on a `volScalarField` (`p` or `pKin`) and derives from
`genericRegionCoupledJump<scalar>`, i.e. from `fixedValue`. It applies the
same normal-stress jump as `regionCoupledPressureValue` and adds the recoil
pressure $-(1/\rho_\text{own}-1/\rho_\text{nbr})\,\dot m^2$. There is no
mass-transfer pressure-flux BC; the partner is presumably
`regionCoupledPressureFlux` (*inferred*). The BC needs a `capillaryInterface`
as the **first** `interfaceType` entry and a `massTransferInterface` on the
same patch pair.

## Mathematical Formulation

The BC sits on side A ("own"), which carries the `faMesh` of the
`capillaryInterface`; "nbr" is side B.
$\tilde{\dot m} = \mathcal{I}_{\text{nbr}\to\text{own}}(\dot m)$.

For a field named `p`:

$$
p_\text{own}^{*} = \mathcal{I}(p_\text{nbr}) + 2(\mu_\text{nbr}-\mu_\text{own})\,\nabla_s\cdot\mathbf{u}_s - \sigma\kappa + (\rho_\text{nbr}-\rho_\text{own})\,(\mathbf{x}_f-\mathbf{x}_\text{ref})\cdot\mathbf{g} - \left(\frac{1}{\rho_\text{own}}-\frac{1}{\rho_\text{nbr}}\right)\tilde{\dot m}^2
\quad \text{on } \Gamma
$$

For a field named `pKin`:

$$
p_{\text{kin,own}}^{*} = \frac{1}{\rho_\text{own}}\Big[\rho_\text{nbr}\,\mathcal{I}(p_\text{kin,nbr}) + 2(\mu_\text{nbr}-\mu_\text{own})\,\nabla_s\cdot\mathbf{u}_s - \sigma\kappa + (\rho_\text{nbr}-\rho_\text{own})\,(\mathbf{x}_f-\mathbf{x}_\text{ref})\cdot\mathbf{g} - \left(\frac{1}{\rho_\text{own}}-\frac{1}{\rho_\text{nbr}}\right)\tilde{\dot m}^2\Big]
$$

followed by relaxation, $p^{k} = p^{k-1} + \omega\,(p^{*} - p^{k-1})$.

| Symbol | Description | Unit |
|--------|-------------|------|
| $p$, $p_\text{kin}$ | Pressure / kinematic pressure | Pa / m²/s² |
| $\dot m$ | Interfacial mass flux, `mDotS()` of the `massTransferInterface` | kg/(m² s) |
| $\mathbf{u}_s$ | Surface velocity on the `faMesh` (`Us()`) | m/s |
| $\sigma$, $\kappa$ | Surface tension and face curvature from the `capillaryInterface` | N/m, 1/m |
| $\mu$, $\rho$ | `mu`, `rho` from each region's `transportProperties` | Pa·s, kg/m³ |
| $\mathbf{g}$ | `uniformDimensionedVectorField g` of the local region | m/s² |
| $\mathbf{x}_\text{ref}$ | `PISO.pRefPoint` from the local `fvSolution` | m |

## Implementation

`updateCoeffs()` is the generic one of `genericRegionCoupledJump`: return if
updated or `coupled no`; call `update()` on the first `interfaceType` (the
`capillaryInterface`); interpolate the neighbour field (looked up by the
local field name); add `valueJump()`; relax; call
`fixedValue::updateCoeffs()`.

`valueJump()` is that of `regionCoupledPressureValue` with the extra recoil
term. It aborts unless the first `interfaceType` is `capillaryInterface`,
finds the `massTransferInterface` by type name and interpolates `mDotS()`
from the neighbour side. `PISO.pRefPoint`, `g`, and `mu`/`rho` of both
regions are required. The mesh flux and patch normals are computed but not
used.

Since only the first interface is `update()`d, $\dot m$ is not recomputed by
this BC. It changes only when a BC that lists `massTransferInterface` first
(e.g. `regionCoupledPhaseChangeHeatFlux` on `T`) is updated. Until then the
recoil term is zero and the BC behaves like `regionCoupledPressureValue`.

The `p` jump residual of the `pUPimple` DNA loop is taken from this BC
(side A is checked first). `flux()` uses the hardcoded diffusivity `v`, for
which `diffCoeff()` returns $\rho$ instead of $1/\rho$ (1 for `pKin`); the
partner `regionCoupledPressureFlux` cancels this flux.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `regionCoupledMassTransferPressureValue` | yes | — |
| `interfaceType` | List of interface type names. First entry must be `capillaryInterface`; `massTransferInterface` must also be listed | yes | — |
| `value` | Initial patch value. Read if present. If missing, the patch values are not initialised | in practice yes | — |
| `diffusivity` | Read, but ignored: the BC always uses `v` | no | empty |
| `neighbourFieldName` | Read by the coupling manager, not used by this BC | no | local field name |
| `jumpCoefficient` | Read by the base class, not used, not written back | no | empty |
| `accType` | Acceleration model: `fixed`, `aitken`, `IQNILS`. A warning is printed if missing | no | `fixed` |
| `relax` | (Initial) relaxation factor | no | `1.0` |
| `couplingReuse` | `IQNILS` only | no | `0` |
| `nonOrthCorr` | Non-orthogonal correction in `snGrad()` / matrix coefficients | no | `false` |
| `secondOrder` | Second-order face-normal gradient | no | `false` |

Entries read from other files: `PISO { pRefPoint (...); }` in the local
`fvSolution`, `g` in the local `constant/`, and `mu`, `rho` in
`transportProperties` of both regions.

## Usage

No tutorial uses this BC. Untested entry on side A, following
`massTransferInterface`:

```
// 0/fluidA/p
interface
{
    type                    regionCoupledMassTransferPressureValue;
    interfaceType           (capillaryInterface massTransferInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform 0;
}
```

Partner on side B (*inferred*):

```
// 0/fluidB/p
interfaceShadow
{
    type                    regionCoupledPressureFlux;
    interfaceType           (capillaryInterface massTransferInterface);
    value                   uniform 0;
}
```

The `regionInterfaceProperties` entry is shown on the
`regionCoupledMassTransferVelocityValue` page and in `massTransferInterface`.

## Known Limitations

- **Untested, no tutorial**; the intended region type `interTrackFluid` is
  disabled.
- **$\dot m$ is only updated through another BC.** See
  [Implementation](#implementation).
- **`mDotS()` is interpolated from the neighbour although it is defined on
  side A** (*inferred*), where this BC sits. Harmless only for identical
  conformal patches.
- **Must sit on side A** (*inferred*): $\sigma$, $\kappa$ and
  $\nabla_s\cdot\mathbf{u}_s$ are used without interpolation.
- **No dedicated flux partner.** `regionCoupledPressureFlux` does not
  account for phase change.
- **Constant properties**, `PISO.pRefPoint` required, `diffusivity` and
  `neighbourFieldName` ignored (as for `regionCoupledPressureValue`).

## Tutorials

| Tutorial | Description |
|----------|-------------|
| — | No tutorial uses `regionCoupledMassTransferPressureValue`. |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
