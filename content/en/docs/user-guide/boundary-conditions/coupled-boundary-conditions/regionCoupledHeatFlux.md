---
mathjax: true
title: "regionCoupledHeatFlux"
description: "Partitioned flux-side temperature BC of heatTransferInterface: fixes the normal temperature gradient so that the conductive heat flux balances the heat flux of the value-side partner."
tags: ["boundary-condition", "coupled", "heat-transfer", "conjugate-heat-transfer"]
weight: 10
---

## Overview

`regionCoupledHeatFlux` is the flux side (Neumann side) of the temperature
coupling used with `heatTransferInterface` for conjugate heat transfer and
thermo-FSI. It derives from `genericRegionCoupledFlux<scalar>` (and so from
`fixedGradient`). It fixes the normal temperature gradient so that the local
conductive heat flux is equal and opposite to the heat flux of the
value-side partner on the neighbour region, without an interfacial source.
The partner must be a `genericRegionCoupledJump<scalar>`-derived BC, in
practice `regionCoupledTemperatureJump`. Unlike its partner, this BC does
**not** check the interface type (see Implementation).

## Mathematical Formulation

Let $B$ be this patch and $A$ the neighbour (jump-side) patch:

$$
\kappa_A \frac{\partial T_A}{\partial n_A} + \kappa_B \frac{\partial T_B}{\partial n_B} = 0
\quad \text{on } \Gamma
$$

imposed as the gradient

$$
\left.\frac{\partial T_B}{\partial n}\right|_b
= -\frac{\mathcal{I}_{A\to B}\!\left(\kappa_A\, \partial_n T_A\right)}{\kappa_B}
$$

followed by relaxation with the acceleration model. `fluxJump()` returns zero.

| Symbol | Description | Unit |
|--------|-------------|------|
| $T$ | Temperature (field `T` in both regions) | K |
| $\kappa_B$ | Diffusivity named by `diffusivity`, default `kappa` | W/(m K) |
| $\kappa_A\,\partial_n T_A$ | `flux()` of the jump-side partner | W/m² |
| $\mathcal{I}_{A\to B}$ | Interface-to-interface interpolation | – |

## Implementation

Run-time work is inherited from `genericRegionCoupledFlux<scalar>` (see
there). The class adds:

- **`fluxJump()`**: returns a zero field.
- **`diffusivityName()`**: the `diffusivity` keyword if set, otherwise
  `kappa`. `diffCoeff()` uses a `volScalarField` with that name on the local
  region if one is registered, otherwise a uniform `dimensionedScalar` from
  the region's `transportProperties`.
- **`heatTransInterface()`**: a type check plus cast to
  `heatTransferInterface`, as in `regionCoupledTemperatureJump`. It is
  declared and defined but **never called**, so the BC accepts any interface
  type listed first in `interfaceType`.

Per DNA iteration of `T`, `updateCoeffs()` updates the interface,
`refCast`s the neighbour `T` patch field to
`genericRegionCoupledJumpFvPatchField<scalar>`, interpolates and negates its
`flux()`, divides by the local $\kappa$, relaxes the gradient and fixes it.
`evaluate()` sets the patch temperature from the gradient (with the optional
non-orthogonal and second-order corrections). The flux residual
$|\kappa_B g_B + \mathcal{I}(q_A)|$ is checked against `maxFluxRes`.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `regionCoupledHeatFlux` | yes | — |
| `interfaceType` | Interface type list; the first entry selects the interface (normally `heatTransferInterface`; not checked) | yes | — |
| `diffusivity` | Name of the thermal conductivity field / `transportProperties` entry | no | `kappa` |
| `neighbourFieldName` | Read and written, not used | no | local field name (`T`) |
| `accType` | `fixed`, `aitken` or `IQNILS` | no | `fixed` (with warning) |
| `relax` | Initial / fixed relaxation factor | no | `1.0` |
| `couplingReuse` | `IQNILS` only: reused time steps | no | `0` |
| `nonOrthCorr` | Non-orthogonal correction in `evaluate()` | no | `false` |
| `secondOrder` | Second-order value evaluation | no | `false` |
| `gradient` | Initial gradient; if given, the initial patch value is evaluated from it | no | zero |
| `value` | Not read (only written) | no | — |

## Usage

From `tutorials/fluidStructureInteraction/heatedBeamInCrossFlow/0/solid/orig/partitioned/T`:

```
interface
{
    type                    regionCoupledHeatFlux;
    interfaceType           (heatTransferInterface);
    diffusivity             "kappa";
    secondOrder             false;
    nonOrthCorr             false;
    value                   $internalField;
}
```

Partner in `0/fluid/orig/partitioned/T`:

```
interfaceShadow
{
    type                    regionCoupledTemperatureJump;
    interfaceType           (heatTransferInterface);
    diffusivity             "kappa";
    accType                 aitken;
    relax                   0.75;
    secondOrder             false;
    nonOrthCorr             false;
    value                   $internalField;
}
```

Interface in `constant/regionInterfaceProperties` (heat transfer stacked on
FSI):

```
fluidsolid
{
    interfaceType (fsiInterface heatTransferInterface);
    coupledPatchPair
    (
        (fluid interfaceShadow)
        (solid interface)
    );
    coupledFields
    (
        D
        T
    );
    coupled             yes;
    fsiInterfaceCoeffs
    {}
    heatTransferInterfaceCoeffs
    {
        jumpCoefficient               jumpCoefficient [ 0 0 0 0 0 0 0 ] 1;
    }
}
```

## Known Limitations

- **No interface-type check.** `heatTransInterface()` is never called, so a
  wrong first `interfaceType` entry is not caught here. It is usually caught
  by the partner `regionCoupledTemperatureJump`, which does check.
- **Partner must be a jump BC** of `T` (bad-cast error otherwise).
- **No interfacial heat source or contact resistance**: `fluxJump()` is zero.
  For phase change use `regionCoupledPhaseChangeHeatFlux`.
- **No unit check on `kappa`** across the two regions.
- **`value` is ignored on input** and **`gradient` is not written** (restart
  starts from a zero gradient).
- **First relaxation uses the patch values as reference** (inherited; the
  tutorials use the default `fixed`/`relax 1` on this side, where it has no
  effect).
- Inherited from `genericRegionCoupledFlux`: `neighbourFieldName` unused,
  corrections only in the evaluated value, `IQNILS` not parallel-consistent,
  partitioned only.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/conjugateHeatTransfer/flowOverHeatedPlate` | Solid side (`conductTemperature`, patch `top`), default `kappa`, no relaxation. |
| `tutorials/conjugateHeatTransfer/ShellAndTubeHeatExchanger` | Solid side towards both fluids (`solid_to_inner`, `solid_to_outer`), `nonOrthCorr true`, `AMI`. |
| `tutorials/fluidStructureInteraction/heatedBeamInCrossFlow` | Solid side (`solidStVK conductTemperature`) of the thermo-FSI interface. |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
