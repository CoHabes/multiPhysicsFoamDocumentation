---
mathjax: true
title: "fixedGradientCorrected"
description: "Fixed-gradient condition whose patch value includes a non-orthogonal correction evaluated on a patch sub-mesh; foam-extend interface-tracking legacy, compiled into libinterfaceTrackingFvMesh."
tags: ["boundary-condition", "uncoupled", "interface-tracking", "non-orthogonal-correction"]
weight: 10
---

## Overview

`fixedGradientCorrected` is a Neumann condition for fields of any type
(scalar, vector, sphericalTensor, symmTensor, tensor). It prescribes the
patch-normal gradient. Its patch value and value coefficients include a
non-orthogonal correction $\mathbf{c}\cdot\nabla\phi_P$, evaluated on a small
sub-mesh around the patch. Unlike [fixedValueCorrected]({{< relref "fixedValueCorrected" >}}),
the correction enters the matrix directly, through `valueBoundaryCoeffs`. It
comes from the foam-extend interface-tracking solver (`interTrackFoam`,
Tuković). The BC is compiled into `libinterfaceTrackingFvMesh`. It derives
from the common base `correctedFvPatchField` (type `corrected`). The base
class, the sub-mesh machinery, the `corrected` type and the porting status
(`interTrackFluid` disabled, `patchCorrected` snGrad not ported) are
described on the [fixedValueCorrected]({{< relref "fixedValueCorrected" >}}) page.

## Mathematical Formulation

$$
\left.\frac{\partial\phi}{\partial n}\right|_b = g \quad \text{on } \partial\Omega,
\qquad
\phi_b = \phi_P + \phi_c + \delta\, g
$$

with the correction from the base class

$$
\mathbf{c} = \boldsymbol{\Delta} - \mathbf{n}\,(\mathbf{n}\cdot\boldsymbol{\Delta}),
\qquad
\phi_c = \mathbf{c}\cdot(\nabla\phi)_P
$$

Matrix coefficients:

| Coefficient | Value |
|-------------|-------|
| `valueInternalCoeffs` | $1$ |
| `valueBoundaryCoeffs` | $\phi_c + \delta\,g$ |
| `gradientInternalCoeffs` | $0$ |
| `gradientBoundaryCoeffs` | $g$ |
| `snGrad()` | $g$ |

| Symbol | Description | Unit |
|--------|-------------|------|
| $\phi$ | Patch field | field unit |
| $g$ | Prescribed normal gradient (`gradient`) | field unit/m |
| $\mathbf{n}$ | Patch unit normal | – |
| $\boldsymbol{\Delta}$ | Patch delta vector (`patch().delta()`) | m |
| $\mathbf{c}$ | Correction vector | m |
| $\phi_c$ | Stored correction (`corrVecGrad`), updated in `updateCoeffs()` | field unit |
| $\delta$ | $1/$`patch().deltaCoeffs()` | m |

**Note on OpenFOAM.com.** `fvPatch::delta()` is patch-normal for every
non-coupled patch in OpenFOAM.com, so $\mathbf{c} \equiv \mathbf{0}$ and
$\phi_c \equiv 0$. The coefficients reduce to those of stock `fixedGradient`.

## Implementation

- **Construction:** the base class reads `volFieldName` (required) and the
  optional `corrVecGrad`. This class then reads `gradient` (required). If
  `value` is given, it is used as the initial patch value. Otherwise
  $\phi_P + \delta g + \phi_c$ is used, with $\phi_c = 0$ at that point unless
  `corrVecGrad` was given. `nGradInternal` is set to $g$.
- **`updateCoeffs()`** (base class): recomputes $\mathbf{c}$ if the mesh moves,
  rebuilds $\phi_c$ on the patch sub-mesh (gradient scheme
  `snGradCorr(<volFieldName>)`). It is called whenever an `fvMatrix` for the
  field is constructed, so `valueBoundaryCoeffs` always uses the current
  $\phi_c$.
- **`evaluate()`:** the class declares `evaluate()` **without** the
  `Pstream::commsTypes` argument. It therefore does not override
  `fvPatchField::evaluate(const Pstream::commsTypes)`, which
  `correctBoundaryConditions()` calls. In practice the base `evaluate` runs,
  and that only resets the `updated` flag. **The patch value is not
  recomputed** from $\phi_P + \phi_c + \delta g$ after a solve. It keeps its
  initial value unless code assigns it explicitly. Only the matrix
  coefficients follow the formula above.
- **Mapping:** `autoMap`/`rmap` map the value and `gradient`.
- **Sub-mesh motion:** `interfaceTrackingFvMesh::update()` moves the
  sub-mesh only for `U` and `p`.
- **`write()`:** `type`, `volFieldName`, `corrVecGrad`, `gradient`, `value`.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `fixedGradientCorrected` | yes | – |
| `volFieldName` | Name of the volume field this patch field belongs to | yes | – (no default) |
| `gradient` | Prescribed normal gradient $g$ | yes | – |
| `value` | Initial patch value | no | $\phi_P + \delta g + \phi_c$ |
| `corrVecGrad` | Initial stored correction $\phi_c$ | no | zero |

The library must be loaded with
`libs ("libinterfaceTrackingFvMesh.so");` in `system/controlDict`.

## Usage

```
outlet
{
    type            fixedGradientCorrected;
    volFieldName    p;
    gradient        uniform 0;
    value           uniform 0;
}
```

## Known Limitations

- **Patch value not updated by `correctBoundaryConditions()`** because the
  `evaluate()` signature does not match (see Implementation). Explicit
  operators that use patch values (`fvc::grad`, `fvc::interpolate`, face
  fluxes) see a stale value.
- **No correction in OpenFOAM.com** (patch-normal `delta()`). The cost of the
  sub-mesh gradient remains on every `updateCoeffs()`.
- **Interface-tracking context not available:** `interTrackFluid` is disabled
  and `patchCorrected` is not ported (see
  [fixedValueCorrected]({{< relref "fixedValueCorrected#known-limitations" >}})).
- **`volFieldName` has no default** and must match the field name.
- **Sub-mesh motion only for `U` and `p`.**

In practice, use the stock `fixedGradient` instead.

## Tutorials

No tutorial uses `fixedGradientCorrected`.

## References

- Tuković, Ž., Jasak, H., *A moving mesh finite volume interface tracking method for surface tension dominated interfacial fluid flow*, Computers & Fluids 55, 2012. DOI: [10.1016/j.compfluid.2011.11.003](https://doi.org/10.1016/j.compfluid.2011.11.003)
- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
