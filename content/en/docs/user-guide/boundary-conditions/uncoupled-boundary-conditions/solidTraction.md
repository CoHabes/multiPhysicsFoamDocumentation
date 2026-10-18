---
mathjax: true
title: "solidTraction"
description: "Prescribed traction and pressure on a solid displacement field, imposed as a fixed normal gradient derived from the current stress state (total Lagrangian, solidStVK)."
tags: ["boundary-condition", "uncoupled", "solid-mechanics", "traction"]
weight: 10
---

## Overview

`solidTraction` imposes a known surface traction $\mathbf{t}$ and/or normal
pressure $p$ on the displacement field `D` of a solid region. It is the
uncoupled counterpart of the FSI traction BC: use it for loaded or
traction-free solid surfaces that are not coupled to another region. The class
derives from `fixedGradient` (vector only). In each update it recomputes the
displacement normal gradient so that the stress on the patch matches the
prescribed load. The traction and the pressure can each be given as a
constant field, as a time series, or as a separate volume field. It was
ported from the foam-extend solid-mechanics toolbox (Cardiff, Tuković). It
reads the fields of `solidStVK` (`impK`, `rImpK`, `gradD`, `sigma`, `Finv`).

## Mathematical Formulation

The target condition is a traction balance in the current (deformed)
configuration:

$$
\mathbf{n}_c \cdot \boldsymbol{\sigma} = \mathbf{t} - p\,\mathbf{n}_c
\quad \text{on } \partial\Omega_t,
\qquad
\mathbf{n}_c = \frac{\mathbf{F}^{-T}\cdot\mathbf{n}_0}{\lvert \mathbf{F}^{-T}\cdot\mathbf{n}_0 \rvert}
$$

It is imposed through the normal gradient of `D`, in the segregated
"implicit stiffness" form used by `solidStVK`:

$$
\mathbf{g}^{\,new} = \frac{1}{K}\left[ \left(\mathbf{t} - p\,\mathbf{n}_c\right) - \mathbf{n}_c\cdot\boldsymbol{\sigma} + K\, \mathbf{n}_0\cdot\nabla\mathbf{D} \right],
\qquad
\mathbf{g} \leftarrow \omega\,\mathbf{g}^{\,new} + (1-\omega)\,\mathbf{g}^{\,old}
$$

$\boldsymbol{\sigma}$ and $\nabla\mathbf{D}$ are the current patch values.
When the iteration has converged, $\mathbf{g} = \mathbf{n}_0\cdot\nabla\mathbf{D}$,
so the bracket requires $\mathbf{n}_c\cdot\boldsymbol{\sigma} = \mathbf{t} - p\,\mathbf{n}_c$.
The $K\,\mathbf{n}_0\cdot\nabla\mathbf{D}$ term mirrors the implicit/explicit
Laplacian pair in the `solidStVK` momentum equation.

The patch value is then reconstructed from the gradient:

$$
\mathbf{D}_b =
\begin{cases}
\mathbf{D}_P + \mathbf{k}\cdot(\nabla\mathbf{D})_P + \delta\,\mathbf{g} & \text{first order} \\[1ex]
\mathbf{D}_P + \mathbf{k}\cdot(\nabla\mathbf{D})_P + \tfrac{\delta}{2}\left(\mathbf{g} + \mathbf{n}_0\cdot(\nabla\mathbf{D})_P\right) & \texttt{secondOrder yes}
\end{cases}
$$

with $\mathbf{k} = (\mathbf{I}-\mathbf{n}_0\mathbf{n}_0)\cdot\boldsymbol{\Delta}$.

| Symbol | Description | Unit |
|--------|-------------|------|
| $\mathbf{D}$ | Total displacement (patch field) | m |
| $\mathbf{t}$ | Prescribed traction (`traction`, `tractionSeries` or `tractionField`) | Pa |
| $p$ | Prescribed pressure, positive = compressive (`pressure`, `pressureSeries` or `pressureField`) | Pa |
| $\mathbf{n}_0$ | Patch unit normal of the undeformed mesh | – |
| $\mathbf{n}_c$ | Deformed unit normal from `Finv` | – |
| $\boldsymbol{\sigma}$ | Cauchy stress, patch values of `sigma` | Pa |
| $\mathbf{F}^{-1}$ | Inverse deformation gradient, patch values of `Finv` | – |
| $K$ | Implicit stiffness, patch values of `impK` ($1/K$ from `rImpK`) | Pa |
| $\mathbf{g}$ | Patch normal gradient (`gradient`) | – |
| $\omega$ | `relaxationFactor` | – |
| $\delta$ | $1/$`patch().deltaCoeffs()` | m |
| $\boldsymbol{\Delta}$, $\mathbf{k}$ | Patch delta vector and non-orthogonal correction vector | m |

**Note on OpenFOAM.com.** In OpenFOAM.com, `fvPatch::delta()` is
patch-normal for non-coupled patches, so $\mathbf{k} \equiv \mathbf{0}$ and the
non-orthogonal term in $\mathbf{D}_b$ has no effect.

## Implementation

**Construction.** The dictionary constructor reads `secondOrder`,
`setEffectiveTraction` and `relaxationFactor`. It reads `gradient` and `value`
if present, otherwise it starts from zero gradient and the patch-internal
value. It then sets up the traction and the pressure. For each, exactly one
of three forms is allowed. Giving two forms of the same quantity stops the
run (`FatalError`):

| Traction | Pressure | Behaviour |
|----------|----------|-----------|
| `traction` | `pressure` | Patch field (`uniform` or `nonuniform`), constant in time |
| `tractionSeries {...}` | `pressureSeries {...}` | `interpolationTable` sub-dictionary, interpolated at `time().timeOutputValue()`. The value is uniform over the patch. |
| `tractionField <name>` | `pressureField <name>` | A `volVectorField` / `volScalarField` with that name is constructed by the BC (`MUST_READ` from the time directory, `AUTO_WRITE`), and its boundary values on this patch are used |

If neither a series nor a field is given, `traction` and `pressure` are read
as plain patch fields, so both are required.

**`updateCoeffs()`** (skipped if already updated in this evaluation):

1. Once per time step (tracked through `curTimeIndex_`), calls
   `correctBoundaryConditions()` on the `tractionField` / `pressureField`, if
   one is used.
2. Updates $\mathbf{t}$ and $p$ from the field or the series.
3. If `setEffectiveTraction` is on, subtracts the patch values of a
   `volScalarField` named `p` from the pressure. This applies the effective
   rather than the total traction, which is meant for poro-elastic solids.
4. Looks up `impK`, `rImpK`, `gradD`, `sigma` and `Finv` on the patch. All are
   mandatory.
5. Computes $\mathbf{n}_c$ and the relaxed gradient $\mathbf{g}$ above, then
   calls `fixedGradient::updateCoeffs()`.

**`evaluate()`** calls `updateCoeffs()` if needed. It looks up the cell
gradient field named `grad<fieldName>` (i.e. `gradD` for `D`), sets
$\mathbf{D}_b$ by the formula above and resets the `updated` flag. The matrix
coefficients are those of the stock `fixedGradient` condition: the gradient
$\mathbf{g}$ enters the Laplacian as an explicit flux.

**Interaction with `solidStVK`.** `D` boundary conditions are updated when the
momentum equation for `D` is assembled. `sigma`, `gradD` and `Finv` come from
the previous Picard / DNA iteration (updated in `solidStVK::postSolve()`), so
the traction is enforced iteratively. It converges together with the `PICARD`
loop on `D`. The `relaxationFactor` gives extra under-relaxation of the
gradient if this loop oscillates.

**`write()`** writes the form actually used for traction and pressure
(`tractionField`/`tractionSeries`/`traction`, and likewise for pressure), then
`secondOrder`, `setEffectiveTraction`, `relaxationFactor`, `value` and
`gradient`.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `solidTraction` | yes | – |
| `traction` | Traction patch field (vector) | yes, unless `tractionSeries` or `tractionField` is given | – |
| `tractionSeries` | `interpolationTable<vector>` sub-dictionary (`file`, optional `outOfBounds`, `readerType`) | alternative to `traction` | – |
| `tractionField` | Name of a `volVectorField` to read from the time directory | alternative to `traction` | – |
| `pressure` | Pressure patch field (scalar, positive = compressive) | yes, unless `pressureSeries` or `pressureField` is given | – |
| `pressureSeries` | `interpolationTable<scalar>` sub-dictionary | alternative to `pressure` | – |
| `pressureField` | Name of a `volScalarField` to read from the time directory | alternative to `pressure` | – |
| `secondOrder` | Second-order reconstruction of the patch value | no | `false` |
| `setEffectiveTraction` | Subtract the region's `p` field from the pressure | no | `false` |
| `relaxationFactor` | Under-relaxation $\omega$ of the gradient | no | `1.0` |
| `gradient` | Initial normal gradient | no | zero |
| `value` | Initial patch value | no | patch-internal value |

In the series sub-dictionary, `outOfBounds` defaults to `warn` and
`readerType` to `openFoam`. Both are standard OpenFOAM `interpolationTable`
keywords.

Required registered fields on the solid mesh:

| Field | Type | Provided by |
|-------|------|-------------|
| `impK`, `rImpK` | `volScalarField` | `solidStVK` |
| `gradD` | `volTensorField` | `solidStVK` |
| `sigma` | `volSymmTensorField` | `solidStVK` |
| `Finv` | `volTensorField` | `solidStVK` |
| `grad<fieldName>` | `volTensorField` | `solidStVK` (`gradD`, when the BC is on `D`) |
| `p` | `volScalarField` | only with `setEffectiveTraction yes`; no shipped solid region type provides it |

## Usage

Traction-free surface (`0/<solid>/D`):

```
freeSurface
{
    type            solidTraction;
    traction        uniform (0 0 0);
    pressure        uniform 0;
    value           uniform (0 0 0);
}
```

Uniform pressure ramped in time from a table file:

```
loadedSurface
{
    type            solidTraction;
    traction        uniform (0 0 0);
    pressureSeries
    {
        file        "$FOAM_CASE/constant/solid/pressureSeries";
        outOfBounds clamp;
    }
    relaxationFactor 0.8;
    value           uniform (0 0 0);
}
```

with `constant/solid/pressureSeries` in the `openFoam` table format, e.g.
`( (0 0) (1 1e4) )`.

## Known Limitations

- **Only for `solidStVK`.** The fields `impK`, `rImpK`, `gradD`, `sigma` and
  `Finv` are looked up unconditionally. No other compiled region type
  provides all of them.
- **Field and series data are uniform or per face only.** `tractionSeries` /
  `pressureSeries` give a value that is uniform over the patch.
  `traction` / `pressure` are constant in time. For space- and time-varying
  loads, use `tractionField` / `pressureField` and drive that field's own
  boundary conditions (they are corrected once per time step).
- **`setEffectiveTraction` needs a `p` field** on the solid mesh. No shipped
  solid region type registers one.
- **Non-orthogonal correction has no effect in OpenFOAM.com** (patch-normal
  `delta()`; see the note above).
- **Copies lose the field link.** The mapping and copy constructors do not
  copy `tractionFieldPtr_` / `pressureFieldPtr_`. A cloned patch field (e.g. a
  stored old-time or previous-iteration copy) keeps the last
  traction/pressure values and stops following the field. The original
  boundary condition of `D` is not affected.
- **Inconsistent gradient names.** `updateCoeffs()` always reads `gradD`,
  while `evaluate()` reads `grad<fieldName>`. On a field other than `D` the two
  differ.

## Tutorials

No tutorial uses `solidTraction`. The FSI tutorials use
`regionCoupledTraction` on the coupled solid surface and `fixedValue` /
`solidSymmetry` elsewhere.

## References

- Cardiff, P., Karač, A., De Jaeger, P., Jasak, H., Nagy, J., Ivanković, A., Tuković, Ž., *An open-source finite volume toolbox for solid mechanics and fluid-solid interaction simulations*, arXiv:1808.10736, 2018.
- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
