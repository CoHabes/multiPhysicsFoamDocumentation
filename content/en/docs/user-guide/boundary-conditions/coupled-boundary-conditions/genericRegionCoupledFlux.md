---
mathjax: true
title: "genericRegionCoupledFlux"
description: "Base class of all partitioned flux-side (Neumann) coupling BCs: sets the patch-normal gradient so that the diffusive flux balances the interpolated flux of the value-side partner plus a flux jump, and relaxes it with an acceleration model."
tags: ["boundary-condition", "coupled", "generic", "partitioned"]
weight: 10
---

## Overview

`genericRegionCoupledFlux<Type>` is the flux side (Neumann side) of the
partitioned Dirichlet–Neumann (DNA) coupling. It derives from `fixedGradient`
and from `interfaceToInterfaceCoupleManager`. In each coupling iteration it
reads the diffusive flux of the value-side partner (a
`genericRegionCoupledJump`-derived BC on the neighbour patch), interpolates it
onto this patch, negates it, adds a *flux jump* $j_q$ and divides by the
local diffusivity. The result is relaxed with an acceleration model and fixed
as the patch-normal gradient. Every partitioned flux-side BC derives from it:
`regionCoupledScalarFlux`, `regionCoupledHeatFlux`,
`regionCoupledPhaseChangeHeatFlux`, `regionCoupledVelocityFlux`,
`regionCoupledPressureFlux`, `regionCoupledMassTransferVelocityFlux` and
`regionCoupledTraction`.

The template is also registered as a selectable patch field type
`genericRegionCoupledFlux` for `scalar` and `vector` fields. Both `fluxJump()`
and `flux()` are `notImplemented`, so the bare type aborts on the first
`updateCoeffs()` (see Known Limitations). In a case, use a derived BC such as
`regionCoupledScalarFlux`.

## Mathematical Formulation

Let $F$ be the patch with this BC and $J$ the neighbour patch with the
value-side partner. The unrelaxed gradient in coupling iteration $k$ is

$$
\left.\frac{\partial \tilde\phi_F}{\partial n}\right|^{k}
= \frac{-\,\mathcal{I}_{J\to F}\!\left(q_J\right) + j_q}{\Gamma_F}
\quad \text{on } \Gamma,
\qquad q_J = \Gamma_J \left.\frac{\partial \phi_J}{\partial n}\right|_b
$$

so that, at convergence,

$$
\Gamma_F \frac{\partial \phi_F}{\partial n_F} + \Gamma_J \frac{\partial \phi_J}{\partial n_J} = j_q .
$$

$q_J$ is the `flux()` of the jump-side partner, computed with its own
non-orthogonal/second-order settings (see `genericRegionCoupledJump`). The
gradient is then relaxed (fixed, Aitken or IQN-ILS, as for
`genericRegionCoupledJump`):

$$
g_F^{\,k} = g_F^{\,k-1} + \omega^k \left(\tilde g_F^{\,k} - g_F^{\,k-1}\right).
$$

The patch value is evaluated from the gradient:

$$
\phi_b =
\begin{cases}
\phi_P + \mathbf{k}\cdot\nabla\phi_P + \dfrac{g_F}{\Delta_f} & \texttt{secondOrder false} \\[8pt]
\phi_P + \mathbf{k}\cdot\nabla\phi_P + \dfrac{g_F + \mathbf{n}\cdot\nabla\phi_P}{2\,\Delta_f} & \texttt{secondOrder true}
\end{cases}
$$

The term $\mathbf{k}\cdot\nabla\phi_P$ ($\mathbf{k} = \mathbf{d} - \mathbf{n}(\mathbf{n}\cdot\mathbf{d})$)
is only present with `nonOrthCorr true`.

| Symbol | Description | Unit |
|--------|-------------|------|
| $\phi$ | Coupled field (`scalar` or `vector`), same name in both regions | field-dependent |
| $g_F$ | Fixed patch-normal gradient (`gradient()`) | field-dependent / m |
| $q_J$ | Diffusive flux of the jump-side partner | field-dependent |
| $\mathcal{I}_{J\to F}$ | Interface-to-interface interpolation of the owning interface | – |
| $j_q$ | Flux jump, `fluxJump()` of the derived class | as $q$ |
| $\Gamma_F$, $\Gamma_J$ | Diffusivities from `diffCoeff()` on each side | field-dependent |
| $\Delta_f$ | `patch().deltaCoeffs()` | 1/m |
| $\mathbf{n}$, $\mathbf{d}$ | Outward unit face normal, face-to-cell-centre vector | –, m |

The interface residual checked by the DNA loop is

$$
r = \left|\Gamma_F\, g_F + \mathcal{I}_{J\to F}(q_J) - j_q\right|,
\qquad
\|r\| = \frac{\|r\|_2}{\max\!\left(\min\!\left(\|\Gamma_F\, g_F\|_2,\ \|-\mathcal{I}_{J\to F}(q_J) + j_q\|_2\right),\ \texttt{SMALL}\right)}
$$

with the global $L_2$ norm over the patch faces. The own flux uses the
stored (relaxed) gradient $g_F$, not a gradient recomputed from the
solution.

## Implementation

**Finding the interface.** As for `genericRegionCoupledJump`: the interface
is the registered `regionInterfaceType` whose type equals the **first**
`interfaceType` entry and that contains this region/patch on side A or B.
No match, or more than one match of that type, is a `FatalError`. The
neighbour region/patch is the other side of that interface.

**`updateCoeffs()`**:

1. Returns at once if already updated or if the interface is not `coupled`
   (the gradient keeps its last value).
2. Calls `update()` on the matched interface (topology update, curvature
   update if a finite-area mesh exists, then the interface's `correct()`).
3. If the derived class returns `useDirectFlux() == true`, the gradient is
   set to `flux()` directly. In this repository only `regionCoupledTraction`
   does that.
4. Otherwise it looks up the neighbour field under the **local field name**,
   `refCast`s its patch field on the neighbour patch to
   `genericRegionCoupledJumpFvPatchField<Type>` and calls its `flux()`. The
   flux is interpolated with `interpolateFromNbrField`, multiplied by $-1$,
   `fluxJump()` is added, and the result is divided by `diffCoeff()`.
5. Relaxes the gradient in place with the acceleration model
   (`accModel_->relax(this->gradient())`).
6. Calls `fixedGradientFvPatchField::updateCoeffs()`.

**`evaluate()`** is specialised for `scalar` and `vector`. It calls
`updateCoeffs()` if needed and sets the patch value with the formula above,
using `grad(<field>)` if registered, otherwise `fvc::grad(<field>)`. The
`vector` version first calls `fixedGradientFvPatchVectorField::evaluate()` and
then overwrites the value. The matrix coefficients (`valueInternalCoeffs`,
`valueBoundaryCoeffs`, `gradientInternalCoeffs`, `gradientBoundaryCoeffs`) are
not overridden. The diffusion term therefore receives exactly $g_F$, while
the corrections only enter through the evaluated patch value.

**`diffCoeff()`** is identical to the jump side: name from
`diffusivityName()` (the `diffusivity` keyword unless overridden), `v` →
`rho`, `volScalarField` of that name if registered, else a uniform
`dimensionedScalar` from the region's `transportProperties`, and 1 for a
field named `pKin`.

**Residuals.** `rawResidual()` and `normResidual()` are called through
`regionInterfaceType::rawFluxRes/normFluxRes`. `dnaControl` compares the
normalised value against `maxFluxRes` and writes the raw value to
`<field>initialFluxRes` / `<field>finalFluxRes` when `outputFluxResField` is
on. `ofNormResidual()` exists but is not called. With `useDirectFlux()` all
residual functions stop with a `FatalError`, so a derived class that uses
the direct flux must not be checked through the generic flux residual.

**Acceleration.** The same `accelerationModel` family as on the jump side
(`fixed`, `aitken`, `IQNILS`), applied to the gradient field. See
`genericRegionCoupledJump` for the update formulas.

## Dictionary Entries

All keywords read by the dictionary constructor (including those of
`interfaceToInterfaceCoupleManager`, `fixedGradient` and the acceleration
model):

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `genericRegionCoupledFlux` (or a derived type) | yes | — |
| `interfaceType` | List of interface type names; only the **first** entry is used | yes | — |
| `neighbourFieldName` | Read and written back, but not used (the neighbour field is looked up under the local name) | no | local field name |
| `diffusivity` | Name of the diffusivity (`volScalarField` or `transportProperties` entry) | no (but see Known Limitations) | `word::null` |
| `accType` | Acceleration model: `fixed`, `aitken`, `IQNILS` | no | `fixed` (a warning is printed when missing) |
| `relax` | Initial / fixed relaxation factor $\omega_0$ | no | `1.0` |
| `couplingReuse` | `IQNILS` only: number of past time steps whose modes are kept | no | `0` |
| `nonOrthCorr` | Non-orthogonal correction in `evaluate()` | no | `false` |
| `secondOrder` | Second-order patch-value evaluation | no | `false` |
| `gradient` | Initial gradient (`fixedGradient`). If present, the patch value is evaluated from it at construction. If absent, the gradient is zero and the patch value is the patch-internal value. | no | zero |
| `value` | Written by `write()`, but **not read** by the constructor (`fixedGradient` constructs the value from `gradient` or by extrapolation) | no | — |

`write()` writes `type`, `interfaceType`, `neighbourFieldName`, `diffusivity`
(if non-empty), `accType`, `relax` (`couplingReuse` for `IQNILS`),
`nonOrthCorr`, `secondOrder` and `value`. The `gradient` is **not** written.

## Usage

As with the jump side, only derived types are usable. Generic layout (shown
with the concrete scalar BC):

```
// 0/regionB/C
interfaceB
{
    type            regionCoupledScalarFlux;   // genericRegionCoupledFlux-derived
    interfaceType   (defaultInterface);
    diffusivity     "D";
    nonOrthCorr     false;
    secondOrder     false;
    value           uniform 0;
}
```

Partner on the neighbour region (must be `genericRegionCoupledJump`-derived,
same field name and type):

```
// 0/regionA/C
interfaceA
{
    type            regionCoupledScalarJump;
    interfaceType   (defaultInterface);
    diffusivity     "D";
    accType         aitken;
    relax           0.75;
    value           uniform 0;
}
```

The interface entry in `constant/regionInterfaceProperties` is the same as
shown for `genericRegionCoupledJump`.

## Known Limitations

- **Not usable on its own.** `fluxJump()` and `flux()` call
  `notImplemented` (a `FatalError`). `updateCoeffs()` calls `fluxJump()`, so
  the bare type aborts on the first coupled solve.
- **The partner must be a jump BC.** The neighbour patch field is
  `refCast` to `genericRegionCoupledJumpFvPatchField<Type>`. Any other
  neighbour BC (including another flux BC) stops the run with a bad-cast
  error.
- **`diffusivity` effectively required.** With the default `word::null` the
  `diffCoeff()` lookup fails unless the derived class overrides
  `diffusivityName()`.
- **First relaxation uses the wrong reference.** The acceleration model is
  initialised with the patch **values** (`initialize(*this)`), but it relaxes
  the **gradient**. In the very first coupled iteration
  $g = \phi_b + \omega(\tilde g - \phi_b)$ for $\omega \ne 1$. With the
  default `fixed`/`relax 1` (all tutorials) this has no effect.
- **`gradient` is not written.** On restart the flux BC starts from a zero
  gradient (or from `gradient` if added by hand), not from the last
  coupled gradient.
- **`value` is ignored on input.** The initial patch value comes from
  `gradient` or zero-gradient extrapolation.
- **Corrections are explicit.** `nonOrthCorr`/`secondOrder` only change the
  evaluated patch value. The fixed-gradient matrix coefficients are unchanged.
- **`diffusivity "v"` is not inverted** (same dead check as on the jump side:
  $\Gamma = \rho$ instead of $1/\rho$).
- **`neighbourFieldName` has no effect**; the same field name is required on
  both sides.
- **`IQNILS` is not parallel-consistent** (processor-local `sum`).
- **Partitioned only.**

## Tutorials

| Tutorial | Description |
|----------|-------------|
| — | No tutorial uses the bare `genericRegionCoupledFlux` type. Derived BCs are used in `tutorials/conjugateHeatTransfer/*`, `tutorials/electrochemistry/NTGKBatteryWithTabs` and `tutorials/fluidStructureInteraction/heatedBeamInCrossFlow` (see `regionCoupledScalarFlux`, `regionCoupledHeatFlux`). |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
