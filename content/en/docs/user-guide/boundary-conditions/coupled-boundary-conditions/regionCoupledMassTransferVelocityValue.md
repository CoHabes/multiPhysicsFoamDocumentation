---
mathjax: true
title: "regionCoupledMassTransferVelocityValue"
description: "Partitioned Dirichlet condition for the velocity at a phase-change interface: neighbour velocity plus the normal velocity jump caused by the interfacial mass flux."
tags: ["boundary-condition", "coupled", "phase-change", "multiphase"]
weight: 10
---

## Overview

`regionCoupledMassTransferVelocityValue` is the Dirichlet (jump) half of the
partitioned velocity coupling across a sharp liquid-vapour interface with
phase change. It acts on a `volVectorField` (`U`) and derives from
`genericRegionCoupledJump<vector>`, i.e. from `fixedValue`. The patch value is
the interpolated neighbour velocity plus the jump in normal velocity that a
mass flux $\dot m$ produces between two phases of different density. The
partner condition on the other side is
`regionCoupledMassTransferVelocityFlux`. The BC needs both a
`capillaryInterface` and a `massTransferInterface` on the same patch pair,
and the pressure is coupled by `regionCoupledMassTransferPressureValue` and
(*inferred*) `regionCoupledPressureFlux`. It is the phase-change counterpart
of `regionCoupledVelocityValue`.

## Mathematical Formulation

"own" is the side the BC sits on and "nbr" the other side.
$\mathcal{I}_{\text{nbr}\to\text{own}}$ is the patch-to-patch interpolation.
$\mathbf{n}_\varepsilon$ is the local patch normal if the local region is
`meshA` of the `capillaryInterface`, and its negative otherwise. In both
cases it points out of region A.

$$
\mathbf{U}_\text{own}^{*} = \mathcal{I}_{\text{nbr}\to\text{own}}(\mathbf{U}_\text{nbr}) + \left(\frac{1}{\rho_\text{own}} - \frac{1}{\rho_\text{nbr}}\right)\tilde{\dot m}\,\mathbf{n}_\varepsilon,
\qquad
\tilde{\dot m} = \mathcal{I}_{\text{nbr}\to\text{own}}(\dot m)
\quad \text{on } \Gamma
$$

followed by relaxation,
$\mathbf{U}^{k} = \mathbf{U}^{k-1} + \omega\,(\mathbf{U}^{*} - \mathbf{U}^{k-1})$.

Unlike `regionCoupledVelocityValue`, no surface-velocity or mesh-velocity
correction of the normal component is applied. That version is still in the
source but commented out.

The diffusive momentum flux handed to `regionCoupledMassTransferVelocityFlux`
is the same as for `regionCoupledVelocityValue`: $\mu_\text{own}$ times the
tangential part of the face-normal gradient (normal component removed in the
vector `snGrad()` of `genericRegionCoupledJump`).

| Symbol | Description | Unit |
|--------|-------------|------|
| $\mathbf{U}$ | Velocity | m/s |
| $\dot m$ | Interfacial mass flux, `mDotS()` of the `massTransferInterface` (area field on the side-A `faMesh`) | kg/(m² s) |
| $\rho_\text{own}$, $\rho_\text{nbr}$ | `rho` from each region's `transportProperties` | kg/m³ |
| $\mathbf{n}_\varepsilon$ | Interface normal pointing out of region A | – |
| $\mu_\text{own}$ | Diffusivity `mu` (field or `transportProperties`) | Pa·s |

## Implementation

`updateCoeffs()` is the generic one of `genericRegionCoupledJump`: return if
updated or `coupled no`; call `update()` on the **first** `interfaceType`;
interpolate the neighbour `U` (looked up by the local field name); add
`valueJump()`; relax; call `fixedValue::updateCoeffs()`.

`valueJump()`:

- finds the `capillaryInterface` and the `massTransferInterface` by type name
  on this patch (`rgInterface(<type>)`), independent of their order in
  `interfaceType`;
- interpolates `mDotS()` from the neighbour side;
- reads `rho` of both regions from `transportProperties`;
- `refCast`s the neighbour `U` patch field to `genericRegionCoupledFlux<vector>`
  and uses its values only for a diagnostic line;
- prints `U Dirichlet boundary contributions: ...` (integrated volume fluxes
  of the neighbour velocity, of the jump term and of mesh motion plus
  $\dot m/\rho$) on every call.

Only the interface listed first in `interfaceType` is `update()`d by this
BC. With `(capillaryInterface massTransferInterface)` the
`massTransferInterface` is not corrected here. $\dot m$ is only recomputed
when a BC that lists `massTransferInterface` first (e.g.
`regionCoupledPhaseChangeHeatFlux` on `T`) is updated. Otherwise it stays at
its initial value of zero.

The class has two more functions that are not used in the solution:

- `updatePhi()` would set the boundary `phi` to the negated interpolated
  neighbour flux plus $\pm(1/\rho_\text{own}-1/\rho_\text{nbr})\dot m\,|\mathbf{S}_f|$.
  The call in the base class is commented out.
- The member `correctClosedVolumePhi(phi, U, p, rAU)` would correct the
  interface flux to the net volume flux $-\dot m/\rho\,|\mathbf{S}_f|$. It is
  never called. The global `correctClosedVolumePhi()` used by `icoFluid` and
  `interTrackFluid` (`closedVolume` on) recognises this BC, but forces the net
  volume flux through the patch to **zero**.

The `U` jump residual is evaluated on this BC if the side-A `U` patch is not
jump-derived (as in the intended setup):
$\lVert\mathbf{U}_\text{own} - \mathcal{I}(\mathbf{U}_\text{nbr}) - \text{valueJump}\rVert_2$,
normalised as in `genericRegionCoupledJump`.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `regionCoupledMassTransferVelocityValue` | yes | — |
| `interfaceType` | List of interface type names. Must contain `capillaryInterface` and `massTransferInterface`; the first entry is the one this BC `update()`s | yes | — |
| `value` | Initial patch value. Read if present. If missing, the patch values are not initialised | in practice yes | — |
| `diffusivity` | Read, but ignored: the BC always uses `mu` | no | empty |
| `neighbourFieldName` | Read by the coupling manager, not used by this BC | no | local field name |
| `jumpCoefficient` | Read by the base class, not used, not written back | no | empty |
| `accType` | Acceleration model: `fixed`, `aitken`, `IQNILS`. A warning is printed if missing | no | `fixed` |
| `relax` | (Initial) relaxation factor | no | `1.0` |
| `couplingReuse` | `IQNILS` only | no | `0` |
| `nonOrthCorr` | Non-orthogonal correction in `snGrad()` / matrix coefficients | no | `false` |
| `secondOrder` | Second-order face-normal gradient | no | `false` |

`rho` (and `mu`) must be present in `transportProperties` of both regions.

## Usage

No tutorial uses this BC. The entry below follows `massTransferInterface`
and has not been run. Side B (second in `coupledPatchPair`):

```
// 0/fluidB/U
interfaceShadow
{
    type                    regionCoupledMassTransferVelocityValue;
    interfaceType           (capillaryInterface massTransferInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform (0 0 0);
}
```

Partner on side A:

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

`regionInterfaceProperties` (excerpt):

```
fluidAfluidB
{
    interfaceType (capillaryInterface massTransferInterface);

    coupledPatchPair
    (
        (fluidA interface)
        (fluidB interfaceShadow)
    );

    coupledFields
    (
        pUPimple
        T
    );

    capillaryInterfaceCoeffs    { sigma sigma [1 0 -2 0 0 0 0] 0.07232; }
    massTransferInterfaceCoeffs { hlv hlv [0 2 -2 0 0 0 0] 2.257e6; TSat TSat [0 0 0 1 0 0 0] 373.15; }
}
```

The full case setup, including the `T` coupling that drives $\dot m$, is in
`massTransferInterface`.

## Known Limitations

- **Untested, no tutorial.** The interface-tracking region type it was
  written for (`interTrackFluid`) is disabled in
  `src/multiRegionSystem/Make/files`.
- **$\dot m$ is only updated through another BC.** See
  [Implementation](#implementation).
- **Must sit on side B** (*inferred*). `mDotS()` lives on the side-A `faMesh`
  and is interpolated "from the neighbour".
- **Neighbour must be flux-derived.** The neighbour `U` patch field is
  `refCast` to `genericRegionCoupledFlux<vector>` (for a log line only); any
  other type aborts.
- **Closed-volume correction ignores phase change.** With `closedVolume` on,
  the global `correctClosedVolumePhi()` forces zero net volume flux through
  this patch, which contradicts a non-zero $\dot m$ (*inferred*).
- **No normal mesh-velocity term.** The normal velocity is not tied to the
  mesh motion as in `regionCoupledVelocityValue` (*inferred from the
  commented-out code; intent unclear*).
- **Log output on every call**; `diffusivity` and `neighbourFieldName` are
  ignored.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| — | No tutorial uses `regionCoupledMassTransferVelocityValue`. |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
