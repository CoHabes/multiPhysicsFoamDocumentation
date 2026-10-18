---
mathjax: true
title: "solidSymmetry"
description: "Symmetry-plane condition for solid displacement (vector) and scalar fields, with a non-orthogonal correction and an optional second-order correction (vector only)."
tags: ["boundary-condition", "uncoupled", "solid-mechanics", "symmetry"]
weight: 10
---

## Overview

`solidSymmetry` is a symmetry-plane condition ported from the foam-extend
solid-mechanics toolbox (Cardiff, Tuković), where it is also known as
`symmetryDisplacement`. It exists for two field types, both registered under
the same type name:

- `solidSymmetryFvPatchVectorField`, intended for the total displacement `D`
  of a solid region (`solidStVK`). It mirrors the near-wall displacement across
  the plane, so the normal displacement on the patch is zero and the
  tangential displacement is free. An optional second-order correction is
  available.
- `solidSymmetryFvPatchScalarField`, the scalar counterpart. It sets the patch
  value to the corrected near-wall value (zero normal gradient).

Both classes derive from the stock `symmetry` constraint
(`symmetryFvPatchField`). They only replace `snGrad()` and `evaluate()`, so the
patch must be of mesh type `symmetry`. The class keeps the matrix
coefficients (implicit part) of the stock condition.

## Mathematical Formulation

Let $\mathbf{n}$ be the patch unit normal, $\boldsymbol{\Delta}$ the patch delta
vector (face centre minus cell centre), $\delta = 1/\texttt{deltaCoeffs}$ the
cell-to-face distance, and $\mathbf{R} = \mathbf{I} - 2\,\mathbf{n}\mathbf{n}$
the reflection tensor. The non-orthogonal correction vector and the corrected
near-wall value are

$$
\mathbf{k} = (\mathbf{I} - \mathbf{n}\mathbf{n}) \cdot \boldsymbol{\Delta},
\qquad
\mathbf{D}_P^{*} = \mathbf{D}_P + \mathbf{k} \cdot (\nabla \mathbf{D})_P
$$

### Vector field, first order (default)

$$
\mathbf{D}_b = \tfrac{1}{2}\left(\mathbf{D}_P^{*} + \mathbf{R}\,\mathbf{D}_P^{*}\right)
             = (\mathbf{I} - \mathbf{n}\mathbf{n})\,\mathbf{D}_P^{*},
\qquad
\left.\frac{\partial \mathbf{D}}{\partial n}\right|_b
  = \frac{\mathbf{R}\,\mathbf{D}_P^{*} - \mathbf{D}_P^{*}}{2\,\delta}
  = -\frac{(\mathbf{n}\cdot\mathbf{D}_P^{*})\,\mathbf{n}}{\delta}
$$

### Vector field, second order (`secondOrder yes`)

$$
\mathbf{D}_b = (\mathbf{I} - \mathbf{n}\mathbf{n}) \cdot
  \left( \mathbf{D}_P^{*} + \tfrac{\delta}{2}\, \mathbf{n}\cdot(\nabla\mathbf{D})_P \right)
$$

$$
\left.\frac{\partial \mathbf{D}}{\partial n}\right|_b
  = -\frac{2\,(\mathbf{n}\cdot\mathbf{D}_P^{*})\,\mathbf{n}}{\delta} - \mathbf{n}\left(\mathbf{n} \cdot \left[\mathbf{n}\cdot(\nabla\mathbf{D})_P\right]\right)
$$

The tangential part of the value is extrapolated over half a cell with the
mean of the zero face gradient and the cell gradient. The normal part of the
gradient is the one-sided second-order estimate
$2(\mathbf{D}_b - \mathbf{D}_P)\cdot\mathbf{n}/\delta - (\partial_n D_n)_P$.

### Scalar field

A scalar is invariant under reflection, so $\mathbf{R}\,T = T$ and the
formulas above reduce to

$$
T_b = T_P + \mathbf{k}\cdot(\nabla T)_P,
\qquad
\left.\frac{\partial T}{\partial n}\right|_b = 0
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $\mathbf{D}$, $T$ | Patch field (vector / scalar) | field unit |
| $\mathbf{n}$ | Patch unit normal (`patch().nf()`) | – |
| $\boldsymbol{\Delta}$ | Patch delta vector (`patch().delta()`) | m |
| $\delta$ | $1/$`patch().deltaCoeffs()` | m |
| $\mathbf{k}$ | Non-orthogonal correction vector | m |
| $(\nabla\mathbf{D})_P$ | Cell-centre gradient from the registered field `gradD` | field unit/m |
| $(\nabla T)_P$ | Cell-centre gradient from the registered field `grad(<field>)` | field unit/m |

**Note on OpenFOAM.com.** In OpenFOAM.com (v2512 and v2606),
`fvPatch::delta()` returns the patch-normal delta
$\mathbf{n}(\mathbf{n}\cdot(\mathbf{C}_f-\mathbf{C}_P))$ for all non-coupled
patches. Then $\mathbf{k} \equiv \mathbf{0}$, $\mathbf{D}_P^{*} = \mathbf{D}_P$,
and the non-orthogonal correction does nothing. It was written for
foam-extend, where `delta()` is the full $\mathbf{C}_f - \mathbf{C}_P$. With
the default first-order setting, the vector condition therefore gives exactly
the values and normal gradient of the stock `symmetry` condition. The scalar
condition is equivalent to stock `symmetry` / `zeroGradient` for scalars.

## Implementation

**Construction.** The dictionary constructor calls the stock `symmetry`
constructor. That constructor ignores `value` and evaluates the field once
with the stock symmetry formula. The constructor then checks that the patch
is exactly a `symmetryFvPatch` (`isType`). Any other patch type, including
`symmetryPlane`, stops the run with a `FatalIOError`. The vector version
reads the optional `secondOrder` switch and prints it. Both versions print
"Symmetry boundary condition with non-orthogonal correction".

**`evaluate()`** (called from `correctBoundaryConditions()`):

1. Calls `updateCoeffs()` if needed. This is the base no-op.
2. Looks up the cell gradient on the patch.
   - Vector: the `volTensorField` named **`gradD`**, whatever the name of the
     field the BC is attached to.
   - Scalar: the `volVectorField` named **`grad(<fieldName>)`**, e.g.
     `grad(T)`.

   The lookup is mandatory. If the field is not registered on the mesh, the
   run stops with a fatal error.
3. Builds $\mathbf{D}_P^{*}$ and assigns the patch value with the formula
   above.
4. Vector only: calls the base `evaluate()`, which resets the `updated` flag.
   The scalar version does not call it.

**`snGrad()`** returns the normal gradient above. It is used wherever the
patch gradient is requested explicitly. The matrix coefficients are inherited
from the stock `transform`/`basicSymmetry` implementation: for vectors the
implicit diagonal is $|n_i|$ per component. Because
`gradientBoundaryCoeffs()` is `snGrad() - gradientInternalCoeffs()·D_P`, the
overridden `snGrad()` (including the second-order term) enters
`fvm::laplacian` as an explicit source. For scalars, the inherited
coefficients are those of a zero-gradient condition.

**In `solidStVK`.** The momentum equation uses
`fvm::laplacian(impKf, D, "laplacian(DD,D)")`, which picks up the coefficients
above. `gradD = fvc::grad(D)` is recomputed in `postSolve()` from the patch
values set in `evaluate()`. `solidStVK` registers `gradD`, so the vector BC
works on its `D` field. It also works on any other vector field of that
region, but it always uses `gradD`.

## Dictionary Entries

Vector version:

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `solidSymmetry` | yes | – |
| `secondOrder` | Second-order correction of value and normal gradient (`Switch`) | no | `false` |
| `patchType` | Generic `fvPatchField` keyword, read by the stock base class | no | – |
| `value` | Not read. The initial value comes from the stock symmetry evaluation. Written on output. | no | – |

Scalar version:

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `solidSymmetry` | yes | – |
| `patchType` | Generic `fvPatchField` keyword | no | – |
| `value` | Not read (see above). Written on output. | no | – |

The scalar version does not read `secondOrder`.

Required registered fields:

| Field | Type | Needed by | Provided by |
|-------|------|-----------|-------------|
| `gradD` | `volTensorField` | vector version | `solidStVK` |
| `grad(<fieldName>)` | `volVectorField` | scalar version | no shipped region type |

## Usage

The mesh patch must be of type `symmetry` (e.g. `symmetry symmetry (...)` in
`blockMeshDict`). From
`tutorials/fluidStructureInteraction/beamInCrossFlow/0/solid/orig/partitioned/D`:

```
symmetry
{
    type            solidSymmetry;
    patchType       symmetry;
    value           uniform (0 0 0);
}
```

With the second-order correction:

```
symmetry
{
    type            solidSymmetry;
    secondOrder     yes;
    value           uniform (0 0 0);
}
```

For other fields of the solid region on the same patch, the tutorials use the
stock condition (e.g. `T` in `heatedBeamInCrossFlow`: `type symmetry;`).

## Known Limitations

- **Non-orthogonal correction has no effect in OpenFOAM.com.** The patch delta
  is patch-normal, so $\mathbf{k} = \mathbf{0}$ (see the note above). With
  `secondOrder no`, the vector BC is numerically identical to stock `symmetry`.
  It still needs `gradD`.
- **Hard-coded gradient name (vector).** The vector version always reads
  `gradD`. It is only usable in a region that registers a `volTensorField`
  named `gradD`, i.e. `solidStVK`.
- **Scalar version needs `grad(<field>)`.** No compiled region type registers
  a `volVectorField` with that name. Unless something else in the case
  provides it, the scalar BC fails on its first evaluation. Because it is
  equivalent to the stock condition in OpenFOAM.com, use `symmetry` for
  scalar fields instead.
- **Patch type `symmetry` only.** `symmetryPlane` and other patch types are
  rejected.
- **`value` is ignored on read.** The initial patch value is always the stock
  symmetry evaluation of the internal field.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/fluidStructureInteraction/beamInCrossFlow` | `D` of the `solidStVK` solid on the `symmetry` patch (first order). |
| `tutorials/fluidStructureInteraction/heatedBeamInCrossFlow` | Same `D` setup in a thermo-FSI case. `T` on the same patch uses stock `symmetry`. |

No tutorial uses the scalar version or `secondOrder yes`.

## References

- Cardiff, P., Karač, A., De Jaeger, P., Jasak, H., Nagy, J., Ivanković, A., Tuković, Ž., *An open-source finite volume toolbox for solid mechanics and fluid-solid interaction simulations*, arXiv:1808.10736, 2018.
- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
