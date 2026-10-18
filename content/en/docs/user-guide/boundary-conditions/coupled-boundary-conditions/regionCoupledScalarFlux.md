---
mathjax: true
title: "regionCoupledScalarFlux"
description: "Partitioned flux-side coupling BC for any scalar field: fixes the normal gradient so that the diffusive flux balances the flux of the value-side partner (zero flux jump), for use with any interface type."
tags: ["boundary-condition", "coupled", "scalar-transport", "partitioned"]
weight: 10
---

## Overview

`regionCoupledScalarFlux` is the concrete scalar flux side (Neumann side) of
the partitioned DNA coupling. It derives from
`genericRegionCoupledFlux<scalar>` (and so from `fixedGradient`). It fixes the
patch-normal gradient so that the local diffusive flux is equal and opposite
to the flux of the value-side partner, i.e. flux continuity without a source
at the interface. It does not check the interface type and is normally used
with `defaultInterface`. The partner on the neighbour patch must be a
`genericRegionCoupledJump<scalar>`-derived BC, usually
`regionCoupledScalarJump`.

## Mathematical Formulation

Let $F$ be this patch and $J$ the neighbour (jump-side) patch:

$$
\Gamma_F \frac{\partial \phi_F}{\partial n_F} + \Gamma_J \frac{\partial \phi_J}{\partial n_J} = 0
\quad \text{on } \Gamma
$$

imposed as the gradient

$$
\left.\frac{\partial \phi_F}{\partial n}\right|_b
= -\frac{\mathcal{I}_{J\to F}\!\left(\Gamma_J\, \partial_n \phi_J\right)}{\Gamma_F}
$$

followed by relaxation with the acceleration model. `fluxJump()` returns zero.

| Symbol | Description | Unit |
|--------|-------------|------|
| $\phi$ | Coupled scalar field (same name in both regions) | field-dependent |
| $\Gamma_F$ | Diffusivity named by this BC's `diffusivity` | field-dependent |
| $\Gamma_J\,\partial_n\phi_J$ | `flux()` of the jump-side partner, with the partner's `diffusivity` | field-dependent |
| $\mathcal{I}_{J\to F}$ | Interface-to-interface interpolation | – |

## Implementation

Everything except `fluxJump()` is inherited from
`genericRegionCoupledFlux<scalar>`. In each DNA iteration `updateCoeffs()`
updates the interface, `refCast`s the neighbour patch field to
`genericRegionCoupledJumpFvPatchField<scalar>`, takes its `flux()`,
interpolates and negates it, adds `fluxJump()` (= 0), divides by the local
`diffCoeff()`, relaxes the gradient and fixes it. `evaluate()` sets the
patch value from the gradient with optional non-orthogonal and second-order
corrections. The flux residual
$|\Gamma_F\,g_F + \mathcal{I}(q_J)|$ is checked against `maxFluxRes`.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `regionCoupledScalarFlux` | yes | — |
| `interfaceType` | Interface type list, e.g. `(defaultInterface)`. Only the first entry is used; any type is accepted. | yes | — |
| `diffusivity` | Diffusivity name (`volScalarField` or `transportProperties` entry) | effectively yes | `word::null` (lookup fails) |
| `neighbourFieldName` | Read and written, not used | no | local field name |
| `accType` | `fixed`, `aitken` or `IQNILS` | no | `fixed` (with warning) |
| `relax` | Initial / fixed relaxation factor | no | `1.0` |
| `couplingReuse` | `IQNILS` only: reused time steps | no | `0` |
| `nonOrthCorr` | Non-orthogonal correction in `evaluate()` | no | `false` |
| `secondOrder` | Second-order value evaluation | no | `false` |
| `gradient` | Initial gradient; if given, the initial patch value is evaluated from it | no | zero |
| `value` | Not read (only written) | no | — |

## Usage

From `tutorials/electrochemistry/NTGKBatteryWithTabs/0/negTabRegion/orig/partitioned/T`
(flux side on the tab):

```
toCellFromNeg
{
    type                    regionCoupledScalarFlux;
    interfaceType           (defaultInterface);
    diffusivity             "k";
    secondOrder             false;
    nonOrthCorr             false;
    value                   $internalField;
}
```

Partner in `0/batteryRegion/orig/partitioned/T`:

```
toNegTab
{
    type                    regionCoupledScalarJump;
    interfaceType           (defaultInterface);
    diffusivity             "k";
    accType                 aitken;
    relax                   0.75;
    secondOrder             false;
    nonOrthCorr             false;
    value                   $internalField;
}
```

The `regionInterfaceProperties` (`cellToNegTab`, `defaultInterface`,
`coupledFields (phiPos phiNeg T)`) and the `DNA` entries are shown in
`regionCoupledScalarJump`. In this tutorial `k` is a uniform entry in the
`transportProperties` of the `NTGK` and `conductPotentialTemperature`
regions.

## Known Limitations

- **Partner must be a jump BC** of the same field (bad-cast error
  otherwise).
- **`diffusivity` must be set**; the default `word::null` lookup fails.
- **No interfacial source.** `fluxJump()` is always zero. Use a derived BC
  with its own `fluxJump()` (e.g. `regionCoupledPhaseChangeHeatFlux`) for a
  flux jump.
- **No consistency check** between the diffusivities named on the two
  sides.
- **`value` is ignored on input** and **`gradient` is not written**. On
  restart the coupled gradient starts from zero unless `gradient` is given.
- **First relaxation uses the patch values as reference** (inherited; only
  matters with `relax` ≠ 1 or `aitken`/`IQNILS` on this side).
- Inherited from `genericRegionCoupledFlux`: `neighbourFieldName` unused,
  `nonOrthCorr`/`secondOrder` only in the evaluated value, `IQNILS` not
  parallel-consistent, partitioned only.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/electrochemistry/NTGKBatteryWithTabs` | `phiPos`: flux side on the battery (`diffusivity sigmaPos`, initial `gradient uniform 0.1266`). `phiNeg`: flux side on the negative tab (`diffusivity sigma`). `T`: flux side on both tabs (`diffusivity k`). All with `defaultInterface` and `directMap`. |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
