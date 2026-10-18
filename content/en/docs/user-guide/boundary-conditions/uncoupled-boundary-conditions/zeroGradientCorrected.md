---
mathjax: true
title: "zeroGradientCorrected"
description: "Zero-gradient condition whose patch value includes a non-orthogonal correction evaluated on a patch sub-mesh; foam-extend interface-tracking legacy, compiled into libinterfaceTrackingFvMesh."
tags: ["boundary-condition", "uncoupled", "interface-tracking", "non-orthogonal-correction"]
weight: 10
---

## Overview

`zeroGradientCorrected` is a zero-gradient condition for fields of any type
(scalar, vector, sphericalTensor, symmTensor, tensor). Instead of copying the
near-wall cell value, it sets the patch value to the cell value plus a
non-orthogonal correction $\mathbf{c}\cdot\nabla\phi_P$, evaluated on a small
sub-mesh around the patch. The same correction enters the matrix through
`valueBoundaryCoeffs`. It comes from the foam-extend interface-tracking
solver (`interTrackFoam`, Tuković). The BC is compiled into
`libinterfaceTrackingFvMesh`. It derives from the common base
`correctedFvPatchField` (type `corrected`). The base class, the sub-mesh
machinery, the `corrected` type and the porting status (`interTrackFluid`
disabled, `patchCorrected` snGrad not ported) are described on the
[fixedValueCorrected]({{< relref "fixedValueCorrected" >}}) page.

## Mathematical Formulation

$$
\left.\frac{\partial\phi}{\partial n}\right|_b = 0 \quad \text{on } \partial\Omega,
\qquad
\phi_b = \phi_P + \phi_c,
\qquad
\phi_c = \mathbf{c}\cdot(\nabla\phi)_P,
\quad
\mathbf{c} = \boldsymbol{\Delta} - \mathbf{n}\,(\mathbf{n}\cdot\boldsymbol{\Delta})
$$

Matrix coefficients:

| Coefficient | Value |
|-------------|-------|
| `valueInternalCoeffs` | $1$ |
| `valueBoundaryCoeffs` | $\phi_c$ |
| `gradientInternalCoeffs` | $0$ |
| `gradientBoundaryCoeffs` | $0$ |
| `snGrad()` | $0$ |

| Symbol | Description | Unit |
|--------|-------------|------|
| $\phi$ | Patch field | field unit |
| $\mathbf{n}$ | Patch unit normal | – |
| $\boldsymbol{\Delta}$ | Patch delta vector (`patch().delta()`) | m |
| $\mathbf{c}$ | Correction vector | m |
| $\phi_c$ | Stored correction (`corrVecGrad`), updated in `updateCoeffs()` | field unit |

**Note on OpenFOAM.com.** `fvPatch::delta()` is patch-normal for every
non-coupled patch in OpenFOAM.com, so $\mathbf{c} \equiv \mathbf{0}$ and
$\phi_c \equiv 0$. The coefficients reduce to those of stock `zeroGradient`.

## Implementation

- **Construction:** the base class reads `volFieldName` (required) and the
  optional `corrVecGrad`. The patch value is set to $\phi_P + \phi_c$. `value`
  is not read.
- **`updateCoeffs()`** (base class): recomputes $\mathbf{c}$ if the mesh moves,
  rebuilds $\phi_c$ on the patch sub-mesh (gradient scheme
  `snGradCorr(<volFieldName>)`). It is called whenever an `fvMatrix` for the
  field is constructed, so `valueBoundaryCoeffs` uses the current $\phi_c$.
- **`evaluate()`:** declared **without** the `Pstream::commsTypes` argument. It
  does not override `fvPatchField::evaluate(const Pstream::commsTypes)`, which
  `correctBoundaryConditions()` calls. The base version runs, and it only
  resets the `updated` flag. **The patch value is not recomputed** after a
  solve. It stays at the construction value $\phi_P + \phi_c$ of the initial
  field, so it does not behave like a zero-gradient value at all.
- **Mapping:** `autoMap`/`rmap` map the value only.
- **Sub-mesh motion:** `interfaceTrackingFvMesh::update()` moves the
  sub-mesh only for `U` and `p`.
- **`write()`:** `type`, `volFieldName`, `corrVecGrad`. No `value` is written.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `zeroGradientCorrected` | yes | – |
| `volFieldName` | Name of the volume field this patch field belongs to | yes | – (no default) |
| `corrVecGrad` | Initial stored correction $\phi_c$ | no | zero |

The library must be loaded with
`libs ("libinterfaceTrackingFvMesh.so");` in `system/controlDict`.

## Usage

```
outlet
{
    type            zeroGradientCorrected;
    volFieldName    U;
}
```

## Known Limitations

- **Patch value frozen at its initial value**, because the `evaluate()`
  signature does not match (see Implementation). Explicit operators that use
  patch values (`fvc::grad`, `fvc::interpolate`, face fluxes) see the initial
  value.
- **No correction in OpenFOAM.com** (patch-normal `delta()`). The cost of the
  sub-mesh gradient remains on every `updateCoeffs()`.
- **Interface-tracking context not available:** `interTrackFluid` is disabled
  and `patchCorrected` is not ported (see
  [fixedValueCorrected]({{< relref "fixedValueCorrected#known-limitations" >}})).
- **`volFieldName` has no default** and must match the field name.
- **Sub-mesh motion only for `U` and `p`.**

In practice, use the stock `zeroGradient` instead.

## Tutorials

No tutorial uses `zeroGradientCorrected`.

## References

- Tuković, Ž., Jasak, H., *A moving mesh finite volume interface tracking method for surface tension dominated interfacial fluid flow*, Computers & Fluids 55, 2012. DOI: [10.1016/j.compfluid.2011.11.003](https://doi.org/10.1016/j.compfluid.2011.11.003)
- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
