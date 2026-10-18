---
mathjax: true
title: "regionCoupledPhaseChangeHeatFlux"
description: "Partitioned flux-side temperature BC of massTransferInterface: heat-flux balance across a phase-change interface with the latent-heat flux jump computed by the interface."
tags: ["boundary-condition", "coupled", "heat-transfer", "phase-change"]
weight: 10
---

## Overview

`regionCoupledPhaseChangeHeatFlux` is the flux side (Neumann side) of the
temperature coupling of `massTransferInterface` (evaporation/condensation
between two fluid regions). It derives from
`genericRegionCoupledFlux<scalar>` (and so from `fixedGradient`). It fixes the
normal temperature gradient so that the conductive heat fluxes of both sides
add up to the interfacial heat-flux jump $q_{sat}$, which
`massTransferInterface` computes and consumes as latent heat. The interface
found through the first `interfaceType` entry must be a
`massTransferInterface`. The partner on the neighbour region must be a
`genericRegionCoupledJump<scalar>`-derived BC. The class does not prescribe
which one: `regionCoupledScalarJump` with `interfaceType
(massTransferInterface)` and `diffusivity k` fits (*inferred*).
`regionCoupledTemperatureJump` only works if a `heatTransferInterface` on the
same patch pair is listed first in its own `interfaceType`.

## Mathematical Formulation

Let $B$ be this patch and $A$ the neighbour (jump-side) patch:

$$
k_B \frac{\partial T_B}{\partial n_B} + k_A \frac{\partial T_A}{\partial n_A} = q_{sat}
\quad \text{on } \Gamma
$$

imposed as the gradient

$$
\left.\frac{\partial T_B}{\partial n}\right|_b
= \frac{-\,\mathcal{I}_{A\to B}\!\left(k_A\, \partial_n T_A\right) + q_{sat}}{k_B}
$$

followed by relaxation with the acceleration model. The flux jump is the
heat-flux imbalance at saturation computed by `massTransferInterface`:

$$
q_{sat} = k_A \left.\frac{\partial T_A}{\partial n_A}\right|_{T_\Gamma = T_{sat}} + k_B \left.\frac{\partial T_B}{\partial n_B}\right|_{T_\Gamma = T_{sat}},
\qquad
\dot m = -\frac{q_{sat}}{h_{lv}}
$$

The value continuity $T_A = T_B$ comes from the jump-side partner. This BC
does not fix $T_\Gamma$ to $T_{sat}$.

| Symbol | Description | Unit |
|--------|-------------|------|
| $T$ | Temperature (field `T` in both regions) | K |
| $k_B$ | Thermal conductivity, diffusivity name fixed to `k` | W/(m K) |
| $k_A\,\partial_n T_A$ | `flux()` of the jump-side partner (its own `diffusivity`) | W/m² |
| $q_{sat}$ | `massTransferInterface::fluxMTJump()` (area field `<patchA>fluxMTJump`) | W/m² |
| $T_{sat}$, $h_{lv}$ | Saturation temperature and latent heat of the interface | K, J/kg |
| $\dot m$ | Interfacial mass flux of the interface | kg/(m² s) |

## Implementation

Run-time work is inherited from `genericRegionCoupledFlux<scalar>` (see
there). The class adds:

- **`fluxJump()`**: returns the internal field of
  `massTransInterface().fluxMTJump()`, an `areaScalarField` on the
  interface's finite-area mesh. It is used **as is**, without interpolation
  between the two sides.
- **`massTransInterface()`**: checks that the interface matched through the
  first `interfaceType` entry is a `massTransferInterface`, otherwise
  `FatalError`. Called from `fluxJump()`, i.e. in every `updateCoeffs()` and
  every residual evaluation.
- **`diffusivityName()`**: hard-coded to `k`. The `diffusivity` keyword is
  read by the base class but ignored. `diffCoeff()` takes a `volScalarField k`
  of the local region if one is registered, otherwise the uniform `k` from
  the region's `transportProperties`.

Per DNA iteration of `T`, `updateCoeffs()` first calls the interface's
`update()`. This runs `massTransferInterface::correct()` → `updateMDotS()`,
which recomputes $q_{sat}$ from copies of both `T` fields with the interface
patches set to $T_{sat}$ (and sets $\dot m$ to zero at time index 1). The
jump-side partner triggers the same update. The BC then `refCast`s the
neighbour `T` patch field to `genericRegionCoupledJumpFvPatchField<scalar>`,
interpolates and negates its `flux()`, adds $q_{sat}$, divides by $k_B$,
relaxes and fixes the gradient. The flux residual
$|k_B g_B + \mathcal{I}(q_A) - q_{sat}|$ is checked against `maxFluxRes`.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `regionCoupledPhaseChangeHeatFlux` | yes | — |
| `interfaceType` | Interface type list. The **first** entry must be `massTransferInterface`. | yes | — |
| `diffusivity` | Read by the base class but **ignored** (`k` is always used) | no | — |
| `neighbourFieldName` | Read and written, not used | no | local field name (`T`) |
| `accType` | `fixed`, `aitken` or `IQNILS` | no | `fixed` (with warning) |
| `relax` | Initial / fixed relaxation factor | no | `1.0` |
| `couplingReuse` | `IQNILS` only: reused time steps | no | `0` |
| `nonOrthCorr` | Non-orthogonal correction in `evaluate()` | no | `false` |
| `secondOrder` | Second-order value evaluation | no | `false` |
| `gradient` | Initial gradient; if given, the initial patch value is evaluated from it | no | zero |
| `value` | Not read (only written) | no | — |

The interface parameters ($T_{sat}$, $h_{lv}$) are read from
`massTransferInterfaceCoeffs` in `constant/regionInterfaceProperties`; see
`docs/interfaces/massTransferInterface.md`.

## Usage

There is no tutorial. Entry pattern (from
`docs/interfaces/massTransferInterface.md`):

```
// 0/fluidB/T
interfaceShadow
{
    type                    regionCoupledPhaseChangeHeatFlux;
    interfaceType           (massTransferInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform 373.15;
}
```

Partner on the other region (*inferred*: the class only requires a
`genericRegionCoupledJump<scalar>`-derived BC):

```
// 0/fluidA/T
interface
{
    type                    regionCoupledScalarJump;
    interfaceType           (massTransferInterface);
    diffusivity             "k";
    accType                 aitken;
    relax                   0.75;
    value                   uniform 373.15;
}
```

Both regions need a uniform `k` in their `transportProperties` (read by the
interface and by this BC). The patch pair is declared in
`partitionedCoupledPatches` with `interfaceType (... massTransferInterface)`,
`coupledFields (T ...)` and a `massTransferInterfaceCoeffs` sub-dictionary.

## Known Limitations

- **No tutorial; untested.** No case uses this BC or `massTransferInterface`.
  The `interTrackFluid` region type it was written for is commented out in
  `src/multiRegionSystem/Make/files`.
- **Needs `massTransferInterface` first** in `interfaceType`; any other type
  aborts in `fluxJump()`.
- **Conformal patches assumed.** $q_{sat}$ lives on the finite-area mesh of
  the interface (built on side A) and is returned without mapping. If this BC
  is on side B, the result is only correct when both patches have the same
  faces in the same order (*inferred*). Otherwise the field sizes do not
  match.
- **Diffusivity name fixed to `k`.** The `diffusivity` keyword is ignored,
  yet `write()` writes `diffusivity k;`. A registered `volScalarField k` takes
  precedence over the `transportProperties` entry. In a turbulent region that
  field is the turbulent kinetic energy, not the conductivity (*inferred*).
- **Interface update on every call.** The interface's `correct()`
  (recomputing $q_{sat}$ and $\dot m$, with an `Info` line) runs in every
  `updateCoeffs()` of every BC coupled through this interface.
- **$T_{sat}$ is not imposed** on the interface temperature; it only enters
  $q_{sat}$.
- Inherited from `genericRegionCoupledFlux`: partner must be a jump BC,
  `value` ignored on input, `gradient` not written, first relaxation uses
  the patch values as reference, `neighbourFieldName` unused, `IQNILS` not
  parallel-consistent, partitioned only.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| — | No tutorial uses `regionCoupledPhaseChangeHeatFlux`. |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
