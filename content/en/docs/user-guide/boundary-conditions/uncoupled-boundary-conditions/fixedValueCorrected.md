---
mathjax: true
title: "fixedValueCorrected"
description: "Fixed-value condition with a stored non-orthogonal correction of the patch-normal gradient, meant for interface-tracking meshes together with the patchCorrected snGrad scheme (not ported)."
tags: ["boundary-condition", "uncoupled", "interface-tracking", "non-orthogonal-correction"]
weight: 10
---

## Overview

`fixedValueCorrected` is a Dirichlet condition for fields of any type
(scalar, vector, sphericalTensor, symmTensor, tensor). On top of the fixed
value it stores a non-orthogonal correction $\mathbf{c}\cdot\nabla\phi_P$ of the
near-wall value. The correction is evaluated on a small sub-mesh around the
patch. It comes from the foam-extend interface-tracking solver
(`interTrackFoam`, Tuković), where it is used on patches whose near-wall
cells become skewed as the mesh follows a moving interface. The fixed-value
matrix coefficients themselves are not corrected. The correction only reaches
the discretisation through `snGradCorrection()`, which is called by the
`patchCorrected` snGrad scheme. That scheme is **not ported** (see
[Known Limitations](#known-limitations)).

The BC is compiled into `libinterfaceTrackingFvMesh`, not
`libmultiRegionSystem`, so this library must be loaded. It is one of three
types derived from the common base `correctedFvPatchField`:

| Type | Page |
|------|------|
| `fixedValueCorrected` | this page |
| `fixedGradientCorrected` | [fixedGradientCorrected]({{< relref "fixedGradientCorrected" >}}) |
| `zeroGradientCorrected` | [zeroGradientCorrected]({{< relref "zeroGradientCorrected" >}}) |
| `corrected` | the base class itself, registered as a type but not usable on its own (see below) |

## Mathematical Formulation

$$
\phi_b = \phi_\text{fixed} \quad \text{on } \partial\Omega
$$

With $\mathbf{n}$ the patch unit normal, $\boldsymbol{\Delta}$ the patch delta
vector and $\delta = 1/\texttt{deltaCoeffs}$, the base class defines the
correction vector and the stored correction

$$
\mathbf{c} = \boldsymbol{\Delta} - \mathbf{n}\,(\mathbf{n}\cdot\boldsymbol{\Delta}),
\qquad
\phi_c = \mathbf{c}\cdot(\nabla\phi)_P
$$

and the corrected normal gradient

$$
\left.\frac{\partial\phi}{\partial n}\right|_b
 = \frac{\phi_b - \phi_P}{\delta} - \frac{\phi_c}{\delta}
 = \frac{\phi_b - (\phi_P + \phi_c)}{\delta}
$$

The matrix coefficients are those of plain `fixedValue`:

| Coefficient | Value |
|-------------|-------|
| `valueInternalCoeffs` | $0$ |
| `valueBoundaryCoeffs` | $\phi_b$ |
| `gradientInternalCoeffs` | $-1/\delta$ |
| `gradientBoundaryCoeffs` | $\phi_b/\delta$ |
| `snGradCorrection()` | $-\phi_c/\delta$ |

| Symbol | Description | Unit |
|--------|-------------|------|
| $\phi$ | Patch field | field unit |
| $\mathbf{n}$ | Patch unit normal | – |
| $\boldsymbol{\Delta}$ | Patch delta vector (`patch().delta()`) | m |
| $\mathbf{c}$ | Correction vector (`correctionVectors()`) | m |
| $(\nabla\phi)_P$ | Cell gradient next to the patch, computed on the patch sub-mesh | field unit/m |
| $\phi_c$ | Stored correction (`corrVecGrad`) | field unit |
| $\delta$ | $1/$`patch().deltaCoeffs()` | m |

**Note on OpenFOAM.com.** In OpenFOAM.com (v2512 and v2606),
`fvPatch::delta()` returns the patch-normal delta for every non-coupled patch.
So $\mathbf{c} \equiv \mathbf{0}$ and $\phi_c \equiv 0$ identically. The BC then
behaves exactly like stock `fixedValue`, but still does the sub-mesh work
described below on every update. The correction is only non-zero in
foam-extend, where `delta()` is the full $\mathbf{C}_f - \mathbf{C}_P$.

## Implementation

**Base class `correctedFvPatchField` (type `corrected`).**

- *Construction:* reads `volFieldName` (required) and, if present,
  `corrVecGrad`. It computes the correction vectors $\mathbf{c}$. It does
  **not** read `value` or `patchType`.
- *`updateCoeffs()`:* if the mesh is moving, recomputes $\mathbf{c}$. Then it
  always calls `updateCorrVecGrad()`:
  1. On first use, builds a patch sub-mesh (`fvMeshSubset`) from the
     patch face cells plus all their neighbouring cells. The sub-mesh takes
     its `fvSchemes`/`fvSolution` from the region mesh.
  2. Interpolates the volume field named `volFieldName` onto the sub-mesh.
  3. For each component, computes `fvc::grad` with the gradient scheme
     named `snGradCorr(<volFieldName>)` (normally covered by
     `gradSchemes { default ...; }`).
  4. Stores $\phi_c = \mathbf{c}\cdot(\nabla\phi)_P$ and
     $\mathbf{n}\cdot(\nabla\phi)_P$ (`nGradInternal`, not used by any active
     code path).
- *`movePatchSubMesh()`:* moves the sub-mesh points to the current mesh
  points. `interfaceTrackingFvMesh::update()` calls it after each mesh motion
  for every `fixedValueCorrected`, `fixedGradientCorrected` or
  `zeroGradientCorrected` patch of the fields **`U`** and **`p`** only. The
  sub-meshes of other fields are not moved.
- *`write()`:* writes `type`, `volFieldName` and `corrVecGrad`.

**`fixedValueCorrected`.**

- *Construction:* reads `value` (required) in addition to the base keywords.
- `fixesValue()` returns `true`, so no pressure reference is needed if the BC
  is used on `p`.
- All assignment operators (`=`, `+=`, `-=`, `*=`, `/=` with fields, patch
  fields or constants) are no-ops, as in stock `fixedValue`. The value can
  only be changed with `==`.
- `evaluate()` is not overridden, so the value stays fixed.
- *`write()`:* base entries plus `value`.

**`patchCorrectedSnGrad` (not compiled).** The source is in
`src/dynamicFvMesh/interfaceTrackingFvMesh/patchCorrectedSnGrad`. It is a
central snGrad scheme with explicit non-orthogonal correction, using the
gradient scheme `snGradCorr(<field>)`. On patches of type
`fixedValueCorrected` it replaces the boundary correction by
`snGradCorrection()`. This is the only code path through which the stored
correction would enter a discretised equation. The entry in
`src/dynamicFvMesh/Make/files` is commented out ("NEEDS TO BE PORTED STILL"),
so the scheme cannot be selected.

### The `corrected` type

`makePatchFields(corrected)` registers the base class as type `corrected` for
scalar, vector, sphericalTensor, symmTensor and tensor fields. The base class
does not implement any matrix coefficients (`valueInternalCoeffs` etc. fall
back to `fvPatchField`, which raises `NotImplemented`). It does not read or
write `value`, so its patch values are uninitialised. It also does not
override `evaluate()`. `type corrected;` is therefore not a usable boundary
condition. Use one of the three derived types.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `fixedValueCorrected` | yes | – |
| `volFieldName` | Name of the volume field this patch field belongs to (e.g. `U`, `p`). Used to look up the field for the sub-mesh gradient. | yes | – (no default) |
| `value` | Fixed patch value | yes | – |
| `corrVecGrad` | Initial stored correction $\phi_c$ (written on output; overwritten at the first update) | no | zero |

Library loading (`system/controlDict` of the case):

```
libs
(
    "libinterfaceTrackingFvMesh.so"
);
```

## Usage

```
movingWall
{
    type            fixedValueCorrected;
    volFieldName    U;
    value           uniform (0 0 0);
}
```

The field must be the one named in `volFieldName`. The sub-mesh is only kept
in step with mesh motion if the field is `U` or `p` and the region mesh is an
`interfaceTrackingFvMesh`.

## Known Limitations

- **No effect in OpenFOAM.com.** The patch delta is patch-normal, so the
  correction vectors and $\phi_c$ are identically zero. The BC is equivalent
  to `fixedValue`, at the cost of a sub-mesh construction and a gradient
  evaluation per component on every `updateCoeffs()`.
- **`patchCorrected` snGrad scheme not ported.** Even with a non-zero
  correction, nothing in the compiled code would use `snGradCorrection()`.
  The matrix coefficients are those of plain `fixedValue`.
- **Interface-tracking region type disabled.** The BCs were meant for
  `interTrackFluid`, which is commented out in
  `src/multiRegionSystem/Make/files`. The `multiphaseFlow` tutorials that need
  it (`2dRisingBubble`, `3dDuineveldRisingBubble`, `3dOscillatingBubble`)
  cannot currently be run. Only `3dDuineveldRisingBubble` mentions
  `patchCorrected`, and only in comments in its base-mesh `fvSchemes`.
- **`volFieldName` has no default** and must match the field name. Otherwise
  the field lookup in `updateCoeffs()` fails.
- **Sub-mesh motion only for `U` and `p`.** For any other field, the sub-mesh
  keeps its initial geometry while the mesh moves.
- **`corrected` is not a standalone BC** (see above).

In practice, use the stock `fixedValue` instead.

## Tutorials

No tutorial uses `fixedValueCorrected` (or any `*Corrected` type). The FSI
tutorials `beamInCrossFlow`, `heatedBeamInCrossFlow` and `HronTurekFsi3` load
`libinterfaceTrackingFvMesh.so` for their `interfaceTrackingFvMesh` fluid mesh,
but their fields use stock and `regionCoupled*` boundary conditions.

## References

- Tuković, Ž., Jasak, H., *A moving mesh finite volume interface tracking method for surface tension dominated interfacial fluid flow*, Computers & Fluids 55, 2012. DOI: [10.1016/j.compfluid.2011.11.003](https://doi.org/10.1016/j.compfluid.2011.11.003)
- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
