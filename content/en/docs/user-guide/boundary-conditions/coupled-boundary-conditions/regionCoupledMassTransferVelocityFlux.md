---
mathjax: true
title: "regionCoupledMassTransferVelocityFlux"
description: "Partitioned Neumann condition for the velocity at a phase-change interface: capillary tangential stress balance plus the momentum flux carried by the interfacial mass flux."
tags: ["boundary-condition", "coupled", "phase-change", "multiphase"]
weight: 10
---

## Overview

`regionCoupledMassTransferVelocityFlux` is the Neumann (flux) half of the
partitioned velocity coupling across a sharp liquid-vapour interface with
phase change. It acts on a `volVectorField` (`U`) and derives from
`genericRegionCoupledFlux<vector>`, i.e. from `fixedGradient`. It uses the
same stress balance as `regionCoupledVelocityFlux` and adds the momentum
flux $\dot m\,[\![\mathbf{u}_\text{rel}]\!]$ that the mass flux carries
across the interface. The partner condition on the other side is
`regionCoupledMassTransferVelocityValue`. The BC needs a `capillaryInterface`
as the **first** `interfaceType` entry and a `massTransferInterface` on the
same patch pair.

## Mathematical Formulation

"own" is the side the BC sits on (side A, which carries the `faMesh`) and
"nbr" the other side. $\mathbf{n}$ is the `faMesh` face normal,
$\mathbf{n}_\varepsilon = \mathbf{n}$ if the local region is `meshA` of the
`capillaryInterface` (otherwise $-\mathbf{n}$). $\tilde{\mathbf{F}}$ is the
interpolated diffusive flux of the partner BC ($\mu_\text{nbr}$ times the
tangential face-normal gradient) and
$\tilde{\dot m} = \mathcal{I}_{\text{nbr}\to\text{own}}(\dot m)$.

$$
\mu_\text{own}\frac{\partial \mathbf{U}_\text{own}}{\partial n} =
-(\mathbf{I}-\mathbf{n}\mathbf{n})\cdot\tilde{\mathbf{F}} + \mathbf{f}_{\sigma,t} - \mu_\text{own}\,\mathbf{n}\,(\nabla_s\cdot\mathbf{u}_s) + (\mu_\text{nbr}-\mu_\text{own})\,(\nabla_s\mathbf{u}_s)\cdot\mathbf{n} + \tilde{\dot m}\left[\left(\frac{\tilde{\dot m}}{\rho_\text{own}}\mathbf{n}_\varepsilon - \mathbf{U}_\text{own}\right) - \left(\frac{\tilde{\dot m}}{\rho_\text{nbr}}\mathbf{n}_\varepsilon - \mathbf{U}_\text{nbr}\right)\right]
$$

with the tangential surface-tension force

$$
\mathbf{f}_{\sigma,t} = \sigma\,\frac{1}{S_f}\sum_{e} \mathbf{L}_e\,c_e \;-\; \sigma\kappa\,\mathbf{n}
$$

$\mathbf{U}_\text{own}$ are the current patch values (previous iteration) and
$\mathbf{U}_\text{nbr}$ the neighbour patch values **without interpolation**
(see [Known Limitations](#known-limitations)). The gradient is relaxed,
$g^{k} = g^{k-1} + \omega\,(g^{*} - g^{k-1})$, and the patch value is
reconstructed from it in `evaluate()` as described for
`regionCoupledVelocityFlux`.

| Symbol | Description | Unit |
|--------|-------------|------|
| $\mathbf{U}$ | Velocity | m/s |
| $\dot m$ | Interfacial mass flux, `mDotS()` of the `massTransferInterface` | kg/(m² s) |
| $\mathbf{u}_s$ | Surface velocity on the `faMesh` (`Us()`) | m/s |
| $\sigma$, $\kappa$ | Surface tension and face curvature from the `capillaryInterface` | N/m, 1/m |
| $\mathbf{L}_e$, $c_e$ | `faMesh` edge length vectors and `edgeLengthCorrection()` | m, – |
| $\mu$, $\rho$ | `mu`, `rho` from each region's `transportProperties` | Pa·s, kg/m³ |

## Implementation

`updateCoeffs()` is the generic one of `genericRegionCoupledFlux`: return if
updated or `coupled no`; call `update()` on the first `interfaceType` (the
`capillaryInterface`); interpolate and negate the neighbour `flux()` (the
neighbour patch field is `refCast` to `genericRegionCoupledJump<vector>`);
add `fluxJump()`; divide by `diffCoeff()` (hardcoded name `mu`: registered
`volScalarField mu`, else `transportProperties`); relax; call
`fixedGradient::updateCoeffs()`.

`fluxJump()`:

- aborts unless the first `interfaceType` is `capillaryInterface`; finds the
  `massTransferInterface` by type name;
- computes $\nabla_s\cdot\mathbf{u}_s$, $\nabla_s\mathbf{u}_s$ and the
  edge-integrated surface-tension force on the `faMesh`;
- reads `mu` and `rho` of both regions from `transportProperties`; if a
  region has a `pKin` field, its `mu` must equal `rho*nu` exactly, otherwise
  the run stops;
- interpolates `mDotS()` from the neighbour side.

Because only the first interface is `update()`d, $\dot m$ is not recomputed by
this BC. It changes only when a BC that lists `massTransferInterface` first
(e.g. `regionCoupledPhaseChangeHeatFlux` on `T`) is updated.

The `U` flux residual of the `pUPimple` DNA loop is taken from this BC if it
sits on side A:
$\lVert\mu_\text{own}\partial_n\mathbf{U}_\text{own} + \tilde{\mathbf{F}} - \text{fluxJump}\rVert_2$,
normalised as in `genericRegionCoupledFlux`.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `regionCoupledMassTransferVelocityFlux` | yes | — |
| `interfaceType` | List of interface type names. First entry must be `capillaryInterface`; `massTransferInterface` must also be listed | yes | — |
| `gradient` | Initial patch-normal gradient. Read if present; otherwise zero gradient and the patch value is set from the internal field | no | zero |
| `value` | Not read (`fixedGradient` ignores it). Written by `write()` | no | — |
| `diffusivity` | Read, but ignored: the BC always uses `mu` | no | empty |
| `neighbourFieldName` | Read by the coupling manager, not used by this BC | no | local field name |
| `accType` | Acceleration model for the gradient: `fixed`, `aitken`, `IQNILS`. A warning is printed if missing | no | `fixed` |
| `relax` | (Initial) relaxation factor | no | `1.0` |
| `couplingReuse` | `IQNILS` only | no | `0` |
| `nonOrthCorr` | Non-orthogonal correction in `evaluate()` | no | `false` |
| `secondOrder` | Second-order reconstruction of the patch value | no | `false` |

`write()` does not write `gradient`.

## Usage

No tutorial uses this BC. Untested entry on side A, following
`massTransferInterface`:

```
// 0/fluidA/U
interface
{
    type                    regionCoupledMassTransferVelocityFlux;
    interfaceType           (capillaryInterface massTransferInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform (0 0 0);
}
```

Partner on side B:

```
// 0/fluidB/U
interfaceShadow
{
    type                    regionCoupledMassTransferVelocityValue;
    interfaceType           (capillaryInterface massTransferInterface);
    value                   uniform (0 0 0);
}
```

The `regionInterfaceProperties` entry is shown on the
`regionCoupledMassTransferVelocityValue` page and in `massTransferInterface`.

## Known Limitations

- **Neighbour velocity is not interpolated (bug).** $\mathbf{U}_\text{nbr}$
  in the mass-flux term is the raw neighbour patch field, combined index by
  index with local fields. This is only correct for conformal patches with
  identical face order; otherwise the sizes may not even match.
- **`mDotS()` is interpolated from the neighbour although it is defined on
  side A** (*inferred*), where this BC sits. Harmless only for identical
  conformal patches.
- **Must sit on side A** (*inferred*): `faMesh` normals, $\mathbf{u}_s$,
  $\sigma$, $\kappa$ are used without interpolation.
- **$\dot m$ is only updated through another BC.** See
  [Implementation](#implementation).
- **Exact `mu == rho*nu` check** for regions with `pKin`.
- **Relaxation history starts from the patch value** (*inferred*). The base
  class initialises the acceleration model with the patch values but relaxes
  the gradient; harmless with `relax 1`.
- **Untested, no tutorial**; `interTrackFluid` is disabled.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| — | No tutorial uses `regionCoupledMassTransferVelocityFlux`. |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
