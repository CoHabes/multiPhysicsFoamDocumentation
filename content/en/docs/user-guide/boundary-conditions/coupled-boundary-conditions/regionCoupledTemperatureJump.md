---
mathjax: true
title: "regionCoupledTemperatureJump"
description: "Partitioned value-side temperature BC of heatTransferInterface: fixes the patch temperature to the interpolated neighbour temperature times the interface's constant jump coefficient."
tags: ["boundary-condition", "coupled", "heat-transfer", "conjugate-heat-transfer"]
weight: 10
---

## Overview

`regionCoupledTemperatureJump` is the value side (Dirichlet side) of the
temperature coupling of `heatTransferInterface`, used for conjugate heat
transfer and thermo-FSI. It derives from `genericRegionCoupledJump<scalar>`
(and so from `fixedValue`). It fixes the patch temperature to the neighbour
temperature, interpolated onto this patch and scaled by the interface's
`jumpCoefficient` $c$ ($c = 1$: temperature continuity). It also supplies its
heat flux $\kappa\,\partial T/\partial n$ to the partner
`regionCoupledHeatFlux` on the neighbour region. The interface found through
the first `interfaceType` entry must be a `heatTransferInterface`.

## Mathematical Formulation

Let $A$ be this patch and $B$ the neighbour (flux-side) patch:

$$
T_A = c\;\mathcal{I}_{B\to A}(T_B) \quad \text{on } \Gamma
$$

implemented as $T_A = \mathcal{I}_{B\to A}(T_B) + j$ with

$$
j = \mathcal{I}_{B\to A}(T_B)\,(c - 1),
$$

followed by relaxation with the acceleration model. The heat flux handed to
`regionCoupledHeatFlux` is

$$
q_A = \kappa_A \left.\frac{\partial T_A}{\partial n}\right|_b
$$

with the coupled `snGrad` of `genericRegionCoupledJump` (optional
non-orthogonal and second-order corrections).

| Symbol | Description | Unit |
|--------|-------------|------|
| $T$ | Temperature (field `T` in both regions) | K |
| $c$ | `heatTransferInterfaceCoeffs/jumpCoefficient` (uniform, constant) | – |
| $\kappa_A$ | Diffusivity named by `diffusivity`, default `kappa` | W/(m K) |
| $\mathcal{I}_{B\to A}$ | Interface-to-interface interpolation | – |

## Implementation

Run-time work is inherited from `genericRegionCoupledJump<scalar>` (see
there). The class adds:

- **`valueJump()`**: looks up the neighbour `T` under the local field name,
  interpolates its neighbour patch values, and returns
  $\mathcal{I}(T_B)\,(c - 1)$. Here $c$ is
  `heatTransInterface().jumpCoefficient().value()`, a single scalar.
- **`heatTransInterface()`**: checks that the interface matched through the
  first `interfaceType` entry has type `heatTransferInterface`, otherwise
  `FatalError` ("BC can only be used in combination with a
  heatTransferInterface"). It is called from `valueJump()`, so the check
  runs on the first `updateCoeffs()` and on every residual evaluation.
- **`diffusivityName()`**: returns the `diffusivity` keyword if set,
  otherwise `kappa`. `diffCoeff()` uses the `volScalarField` with this name
  on the local region if one is registered, else a uniform
  `dimensionedScalar` of that name from the region's `transportProperties`.

Per DNA iteration of `T`: the interface is updated (its `correct()` is a
no-op), the value $\mathcal{I}(T_B) + j$ is set, relaxed (`accType`, `relax`)
and fixed. The jump residual $|T_A - \mathcal{I}(T_B) - j|$ is checked against
`maxJumpRes`.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `regionCoupledTemperatureJump` | yes | — |
| `interfaceType` | Interface type list. The **first** entry must be `heatTransferInterface`. | yes | — |
| `diffusivity` | Name of the thermal conductivity field / `transportProperties` entry | no | `kappa` |
| `neighbourFieldName` | Read and written, not used | no | local field name (`T`) |
| `jumpCoefficient` | Read by the base class but **ignored** (the coefficient comes from the interface) | no | `word::null` |
| `accType` | `fixed`, `aitken` or `IQNILS` | no | `fixed` (with warning) |
| `relax` | Initial / fixed relaxation factor | no | `1.0` |
| `couplingReuse` | `IQNILS` only: reused time steps | no | `0` |
| `nonOrthCorr` | Non-orthogonal correction of `snGrad` | no | `false` |
| `secondOrder` | Second-order `snGrad` | no | `false` |
| `value` | Initial patch value (read if present, otherwise uninitialised) | effectively yes | — |

## Usage

From `tutorials/conjugateHeatTransfer/flowOverHeatedPlate/0/fluid/orig/partitioned/T`:

```
interface
{
    type                    regionCoupledTemperatureJump;
    interfaceType           (heatTransferInterface);

    accType                 aitken;
    relax                   0.75;

    secondOrder             false;
    nonOrthCorr             false;
    value                   $internalField;
}
```

Partner in `0/solid/orig/partitioned/T`:

```
top
{
    type                    regionCoupledHeatFlux;
    interfaceType           (heatTransferInterface);

    secondOrder             false;
    nonOrthCorr             false;
    value                   $internalField;
}
```

Interface in `constant/regionInterfacePropertiesPartitioned`:

```
partitionedCoupledPatches
(
    fluidsolid
    {
        interfaceType (heatTransferInterface);
        coupledPatchPair
        (
            (fluid interface)
            (solid top)
        );
        coupledFields
        (
            T
        );
        heatTransferInterfaceCoeffs
        {
            jumpCoefficient               jumpCoefficient [ 0 0 0 0 0 0 0 ] 1;
        }
    }
);
```

and the `T` entry in the `DNA` block of `constant/multiRegionProperties`
(`maxCoupleIter`, `maxJumpRes`, `maxFluxRes`). For thermo-FSI the interface
entry lists `interfaceType (fsiInterface heatTransferInterface);`, while the
`T` patches keep `interfaceType (heatTransferInterface);`
(`tutorials/fluidStructureInteraction/heatedBeamInCrossFlow`).

## Known Limitations

- **Needs `heatTransferInterface` first** in the BC's `interfaceType` list;
  any other type aborts in `valueJump()`.
- **Constant, uniform jump coefficient.** $c$ is a single
  `dimensionedScalar` read once by the interface. The `jumpCoefficient` patch
  keyword (field lookup) of `regionCoupledScalarJump` is not used here.
- **No unit check on `kappa`.** The flux balance with the partner only
  conserves energy if both sides name conductivities in the same units
  (e.g. not `DT` in m²/s on one side).
- **Missing diffusivity is fatal** when the partner first asks for `flux()`
  (no `kappa` field and no `kappa` in `transportProperties`).
- Inherited from `genericRegionCoupledJump`: `neighbourFieldName` unused,
  side detection by patch name only, acceleration state carried over time
  steps, `IQNILS` not parallel-consistent, partitioned only (the
  `*Monolithic` tutorial variants use other BCs).

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/conjugateHeatTransfer/flowOverHeatedPlate` | Fluid side of `fluid (pimpleFluid transportTemperature)` / `solid (conductTemperature)`, Aitken 0.75, `directMap`. |
| `tutorials/conjugateHeatTransfer/ShellAndTubeHeatExchanger` | Fluid side on `inner_to_solid` and `outer_to_solid`, `nonOrthCorr true`, default `fixed` relaxation, `AMI`. |
| `tutorials/fluidStructureInteraction/heatedBeamInCrossFlow` | Fluid side of the thermo-FSI interface (`diffusivity "kappa"`, Aitken 0.75), stacked with `fsiInterface`. |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
