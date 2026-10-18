---
mathjax: true
title: "regionCoupledScalarJump"
description: "Partitioned value-side coupling BC for any scalar field: fixes the patch value to the interpolated neighbour value times an optional partition coefficient, for use with any interface type."
tags: ["boundary-condition", "coupled", "scalar-transport", "partitioned"]
weight: 10
---

## Overview

`regionCoupledScalarJump` is the concrete scalar value side (Dirichlet side)
of the partitioned DNA coupling. It derives from
`genericRegionCoupledJump<scalar>` (and so from `fixedValue`). It fixes the
patch value to the neighbour value interpolated onto this patch, scaled by a
jump (partition) coefficient $c$, which defaults to 1, i.e. continuity. It does
not check the interface type, so it works with any interface, in practice
`defaultInterface`. It is the generic choice for potentials, temperatures
and concentrations that need no interface physics. The partner on the
neighbour patch is `regionCoupledScalarFlux` (or any other
`genericRegionCoupledFlux<scalar>`-derived BC, e.g.
`regionCoupledPhaseChangeHeatFlux`).

## Mathematical Formulation

Let $J$ be this patch and $F$ the neighbour (flux-side) patch:

$$
\phi_J = c\;\mathcal{I}_{F\to J}(\phi_F) \quad \text{on } \Gamma
$$

implemented as $\phi_J = \mathcal{I}_{F\to J}(\phi_F) + j$ with the value jump

$$
j = \mathcal{I}_{F\to J}(\phi_F)\,(c - 1).
$$

The value is then relaxed by the acceleration model. The flux handed to the
partner is

$$
q_J = \Gamma_J \left.\frac{\partial \phi_J}{\partial n}\right|_b
$$

with the coupled `snGrad` of `genericRegionCoupledJump` (optional
non-orthogonal and second-order corrections).

| Symbol | Description | Unit |
|--------|-------------|------|
| $\phi$ | Coupled scalar field (same name in both regions) | field-dependent |
| $c$ | Jump (partition) coefficient, keyword `jumpCoefficient` | – |
| $\mathcal{I}_{F\to J}$ | Interface-to-interface interpolation | – |
| $\Gamma_J$ | Diffusivity named by `diffusivity` | field-dependent |

## Implementation

All run-time work (interface lookup, `updateCoeffs()`, relaxation, `snGrad`,
`flux()`, `diffCoeff()`, residuals) is inherited from
`genericRegionCoupledJump<scalar>`. The class only implements
`valueJump()`:

1. Looks up the neighbour field under the local field name and interpolates
   its neighbour patch values onto this patch.
2. Determines $c$ per face:
   - if a `volScalarField` named `jumpCoefficient` is registered on the
     **local** region, its local patch values are used;
   - else, if `jumpCoefficient` is set, a uniform `dimensionedScalar` of that
     name is read from the local region's `transportProperties`;
   - else $c = 1$.
3. Returns $\mathcal{I}_{F\to J}(\phi_F)\,(c - 1)$.

`valueJump()` is evaluated in `updateCoeffs()` and again in each residual
evaluation (jump residual $|\phi_J - \mathcal{I}(\phi_F) - j|$, checked
against `maxJumpRes`).

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `regionCoupledScalarJump` | yes | — |
| `interfaceType` | Interface type list, e.g. `(defaultInterface)`. Only the first entry is used to find the interface; any type is accepted. | yes | — |
| `diffusivity` | Diffusivity name (`volScalarField` or `transportProperties` entry) used by `flux()` for the partner | effectively yes | `word::null` (lookup fails) |
| `jumpCoefficient` | Name of the partition-coefficient field / `transportProperties` entry | no | `word::null` → $c = 1$ |
| `neighbourFieldName` | Read and written, not used | no | local field name |
| `accType` | `fixed`, `aitken` or `IQNILS` | no | `fixed` (with warning) |
| `relax` | Initial / fixed relaxation factor | no | `1.0` |
| `couplingReuse` | `IQNILS` only: reused time steps | no | `0` |
| `nonOrthCorr` | Non-orthogonal correction of `snGrad` | no | `false` |
| `secondOrder` | Second-order `snGrad` | no | `false` |
| `value` | Initial patch value (read if present, otherwise uninitialised) | effectively yes | — |

## Usage

From `tutorials/electrochemistry/NTGKBatteryWithTabs/0/posTabRegion/orig/partitioned/phiPos`
(jump side on the tab):

```
toCellFromPos
{
    type                    regionCoupledScalarJump;
    interfaceType           (defaultInterface);
    diffusivity             "sigma";
    accType                 aitken;
    relax                   0.1;
    secondOrder             false;
    nonOrthCorr             false;
    value                   $internalField;
}
```

Partner in `0/batteryRegion/orig/partitioned/phiPos`:

```
toPosTab
{
    type                    regionCoupledScalarFlux;
    interfaceType           (defaultInterface);
    diffusivity             "sigmaPos";
    secondOrder             false;
    nonOrthCorr             false;
    value                   $internalField;
    gradient                uniform 0.1266;
}
```

Interface entry in `constant/regionInterfaceProperties`:

```
cellToPosTab
{
    interfaceType (defaultInterface);
    coupledPatchPair
    (
        (batteryRegion toPosTab)
        (posTabRegion  toCellFromPos)
    );
    coupledFields
    (
        phiPos
        phiNeg
        T
    );
    defaultInterfaceCoeffs
    {}
}
```

and a DNA entry per coupled field in `constant/multiRegionProperties`:

```
DNA
{
    phiPos
    {
        maxCoupleIter 50;
        residualControl
        {
            maxJumpRes 1e-05;
            outputJumpResField no;
            maxFluxRes 1e-05;
            outputFluxResField no;
        }
    }
}
```

For a partition coefficient (e.g. Henry's law $C_A = H\,C_B$), add
`jumpCoefficient H;` and provide `H` either as a `volScalarField` on the
jump-side region or as `H H [0 0 0 0 0 0 0] 0.5;` in its
`transportProperties`.

## Known Limitations

- **No interface-type check.** Any interface type listed first in
  `interfaceType` is accepted. The interface's own physics, if any, is not
  used, except that its `correct()` is called through `update()`.
- **`diffusivity` must be set.** The default `word::null` makes the
  `diffCoeff()` lookup in `transportProperties` fail as soon as the partner
  asks for `flux()`.
- **Diffusivities are not checked for consistency.** Each side uses whatever
  its `diffusivity` names (in NTGK, `sigmaPos` on the battery and `sigma` on
  the tab). Flux conservation requires the same units on both sides.
- **Jump coefficient is taken from the jump side only** (local field or local
  `transportProperties`). With a `volScalarField`, the local patch values
  are used and not interpolated.
- **`jumpCoefficient` is not written** by `write()`, so it disappears from
  written time directories and after `decomposePar`.
- Inherited limitations of `genericRegionCoupledJump`: `neighbourFieldName`
  unused, side detection by patch name only, `IQNILS` not parallel-consistent,
  partitioned only.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/electrochemistry/NTGKBatteryWithTabs` | `batteryRegion (NTGK)` coupled to two `conductPotentialTemperature` tabs through `defaultInterface`. `phiPos`: jump side on the positive tab (`diffusivity sigma`, Aitken 0.1). `phiNeg`: jump side on the battery (`diffusivity sigmaNeg`, Aitken 0.75). `T`: jump side on the battery at both tabs (`diffusivity k`, Aitken 0.75). `directMap` transfer. |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
