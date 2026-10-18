---
mathjax: true
title: "genericRegionCoupledJump"
description: "Base class of all partitioned value-side (Dirichlet) coupling BCs: sets the patch value to the interpolated neighbour value plus a value jump, relaxes it with an acceleration model and supplies its diffusive flux to the flux-side partner."
tags: ["boundary-condition", "coupled", "generic", "partitioned"]
weight: 10
---

## Overview

`genericRegionCoupledJump<Type>` is the value side (Dirichlet side) of the
partitioned Dirichlet–Neumann (DNA) coupling. It derives from `fixedValue` and
from `interfaceToInterfaceCoupleManager`. In each coupling iteration it sets
the patch value to the neighbour patch value, interpolated onto this patch,
plus a *value jump* $j$. It then relaxes this value with an acceleration
model. It also provides the diffusive flux $\Gamma\,\partial\phi/\partial n$
that the flux-side partner (`genericRegionCoupledFlux` or a class derived
from it) reads on the other side of the interface. Every partitioned value-side
BC derives from it: `regionCoupledScalarJump`, `regionCoupledTemperatureJump`,
`regionCoupledVelocityValue`, `regionCoupledPressureValue`,
`regionCoupledMassTransferVelocityValue` and
`regionCoupledMassTransferPressureValue`.

The template is also registered as a selectable patch field type
`genericRegionCoupledJump` for `scalar` and `vector` fields. It cannot be used
on its own, though: `valueJump()` is `notImplemented` and aborts the run on
the first `updateCoeffs()` (see Known Limitations). In a case, use a derived
BC such as `regionCoupledScalarJump` instead.

## Mathematical Formulation

Let $J$ be the patch with this BC and $F$ the neighbour patch with the
flux-side partner. In coupling iteration $k$ the BC computes

$$
\tilde\phi_J^{\,k} = \mathcal{I}_{F\to J}\!\left(\phi_F\right) + j
\quad \text{on } \Gamma
$$

and fixes the relaxed value

$$
\phi_J^{\,k} = \phi_J^{\,k-1} + \omega^k \left(\tilde\phi_J^{\,k} - \phi_J^{\,k-1}\right)
$$

(or the IQN-ILS update; see Implementation). The flux handed to the partner is

$$
q_J = \Gamma_J \left.\frac{\partial \phi_J}{\partial n}\right|_b
$$

with the coupled patch-normal gradient (scalar fields)

$$
\left.\frac{\partial \phi}{\partial n}\right|_b =
\begin{cases}
\Delta_f \left(\phi_b - \phi_P - \mathbf{k}\cdot\nabla\phi_P\right) & \texttt{secondOrder false} \\[4pt]
2\,\Delta_f \left(\phi_b - \phi_P - \mathbf{k}\cdot\nabla\phi_P\right) - \mathbf{n}\cdot\nabla\phi_P & \texttt{secondOrder true}
\end{cases}
\qquad
\mathbf{k} = \mathbf{d} - \mathbf{n}\,(\mathbf{n}\cdot\mathbf{d})
$$

The term $\mathbf{k}\cdot\nabla\phi_P$ is only present with `nonOrthCorr true`.
For `vector` fields the same expressions are used, and then the normal component is
removed: $\partial_n\boldsymbol{\phi} \leftarrow (\mathbf{I} - \mathbf{n}\mathbf{n})\cdot\partial_n\boldsymbol{\phi}$.

| Symbol | Description | Unit |
|--------|-------------|------|
| $\phi$ | Coupled field (`scalar` or `vector`), same name in both regions | field-dependent |
| $\mathcal{I}_{F\to J}$ | Interface-to-interface interpolation (`AMI` or `directMap`) of the owning interface | – |
| $j$ | Value jump, `valueJump()` of the derived class | as $\phi$ |
| $\omega^k$ | Relaxation factor of the acceleration model | – |
| $\Gamma_J$ | Diffusivity of side $J$, from `diffCoeff()` | field-dependent |
| $\phi_P$, $\nabla\phi_P$ | Patch-internal cell value and cell gradient | – |
| $\Delta_f$ | `patch().deltaCoeffs()` | 1/m |
| $\mathbf{n}$, $\mathbf{d}$ | Outward unit face normal, face-to-cell-centre vector (`patch().delta()`) | –, m |

The interface residual checked by the DNA loop is

$$
r = \left|\phi_J - \mathcal{I}_{F\to J}(\phi_F) - j\right|,
\qquad
\|r\| = \frac{\|r\|_2}{\max\!\left(\min\!\left(\|\phi_J\|_2,\ \|\mathcal{I}_{F\to J}(\phi_F) + j\|_2\right),\ \texttt{SMALL}\right)}
$$

where $\|\cdot\|_2$ is the global (all-processor) $L_2$ norm over the patch faces.

## Implementation

**Finding the interface.** The coupling data are taken from a
`regionInterfaceType` object. `interfaceToInterfaceCoupleManager` looks up all
registered interfaces and keeps the one whose type equals the **first** entry
of the BC's `interfaceType` list and whose side A or side B is this
region/patch. If there is no match, or more than one interface of that type
contains this patch, the run stops with a `FatalError`. The neighbour region
and patch are taken from the matched interface: if this patch's name equals
`patchA` of the interface, the neighbour is side B, otherwise side A.

**`updateCoeffs()`**, called when the coupled equation is assembled in the
DNA loop:

1. Returns at once if the BC is already updated or the interface is not
   `coupled`. In that case the patch keeps its last value.
2. Calls `update()` on the matched interface. This rebuilds the global
   patches if the mesh topology changed, updates the curvature if the
   interface has a finite-area mesh, and calls the interface's `correct()`.
3. If the derived class returns `useDirectValue() == true`, the value is set
   to `value()` directly. No derived class in this repository does that; the
   base class returns `false`.
4. Otherwise it looks up the neighbour field **under the local field name**,
   interpolates the neighbour patch values with
   `interpolateFromNbrField`, adds `valueJump()` and assigns the result.
5. Relaxes the patch values in place with the acceleration model
   (`accModel_->relax(*this)`).
6. Calls `fixedValueFvPatchField::updateCoeffs()`.

**Matrix coefficients.** `snGrad()` and `gradientBoundaryCoeffs()` are
specialised for `scalar` and `vector`. They use the cell gradient
`grad(<field>)` if a field of that name is registered, otherwise they compute
`fvc::grad(<field>)` on the fly. They add the non-orthogonal and
second-order corrections of the formula above, so the implicit
`fvm::laplacian` of the coupled equation sees the same face gradient that
`flux()` reports. `gradientInternalCoeffs()` is inherited from `fixedValue`
($-\Delta_f$).

**`flux()`** returns `snGrad() * diffCoeff()`. The flux-side partner calls it
through a `refCast` to `genericRegionCoupledJumpFvPatchField<Type>`.

**`diffCoeff()`** resolves the diffusivity name returned by the virtual
`diffusivityName()`. The base class returns the `diffusivity` keyword; derived
classes may override it (`regionCoupledTemperatureJump` falls back to
`kappa`):

1. The name `v` is replaced by `rho`.
2. If a `volScalarField` with that name is registered on the local region,
   its patch values are used.
3. Otherwise a uniform `dimensionedScalar` with that name is read from the
   region's registered `transportProperties` dictionary. A missing entry or
   dictionary is a fatal error.
4. For a field named `pKin` the diffusivity is set to 1.

**Residuals.** `rawResidual()` (per face) and `normResidual()` (normalised,
see above) are called by `regionInterfaceType::rawJumpRes/normJumpRes`. The
DNA loop (`dnaControl`) compares the normalised residual against
`maxJumpRes` and writes the raw residual to `<field>initialJumpRes` /
`<field>finalJumpRes` if `outputJumpResField` is on. `ofNormResidual()`
(an OpenFOAM-style normalisation) is implemented but not called anywhere.
All three residual functions stop with a `FatalError` if
`useDirectValue()` is true.

### Acceleration models

The acceleration model is built from the patch dictionary by
`accelerationModel<Type>::New`. It is initialised with the initial patch
values and stores the last relaxed value $x^{k-1}$. In each call it forms the
residual $r^k = \tilde x^k - x^{k-1}$. The iteration counter is reset to 1 on
the first call in a new time step; $x^{k-1}$ is **not** reset, so the first
iteration of a time step is relaxed against the final value of the previous
one.

| `accType` | Update |
|-----------|--------|
| `fixed` (default) | $x^k = x^{k-1} + \omega_0\, r^k$ with $\omega_0 =$ `relax`. With `relax 1` this is no relaxation. |
| `aitken` | $x^k = x^{k-1} + \omega^k r^k$ with $\omega^k = \left| -\omega^{k-1} \dfrac{\sum r^{k-1}\cdot(r^k - r^{k-1})}{\sum |r^k - r^{k-1}|^2 + \texttt{SMALL}} \right|$, sums over all processors. In the first two iterations of a time step $\omega^k = \max(\omega^{k-1}, \omega_0)$. $\omega^k$ is clipped to 1 from above, and reset to $\omega_0$ if it falls below `SMALL`. |
| `IQNILS` | Interface quasi-Newton with inverse Jacobian from a least-squares model. Iteration 1 of a time step drops stored modes older than `couplingReuse` time steps. Iteration 2 stores the reference $(\tilde x_{ref}, r_{ref})$. From iteration 3 the columns $V_i = r^k - r_{ref}$, $W_i = \tilde x^k - \tilde x_{ref}$ are added. With more than one column, $c = \arg\min\|V c + r^k\|$ is solved by modified Gram–Schmidt QR and $x^k = \tilde x^k + W c$. Otherwise the fixed update with $\omega_0$ is used. |

## Dictionary Entries

All keywords read by the dictionary constructor (including those of
`interfaceToInterfaceCoupleManager` and the acceleration model):

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `genericRegionCoupledJump` (or a derived type) | yes | — |
| `interfaceType` | List of interface type names, e.g. `(defaultInterface)`. Only the **first** entry is used to find the interface. Must be a list. | yes | — |
| `neighbourFieldName` | Name of the neighbour field. Read and written back, but not used: the neighbour field is always looked up under the local field name. | no | local field name |
| `diffusivity` | Name of the diffusivity (`volScalarField` or `transportProperties` entry) used by `flux()` | no (but see Known Limitations) | `word::null` |
| `jumpCoefficient` | Name of a jump-coefficient field/entry. Stored by the base class, used only by `regionCoupledScalarJump`. | no | `word::null` |
| `accType` | Acceleration model: `fixed`, `aitken`, `IQNILS` | no | `fixed` (a warning is printed when missing) |
| `relax` | Initial / fixed relaxation factor $\omega_0$ | no | `1.0` |
| `couplingReuse` | `IQNILS` only: number of past time steps whose modes are kept | no | `0` |
| `nonOrthCorr` | Non-orthogonal correction $\mathbf{k}\cdot\nabla\phi_P$ in `snGrad` | no | `false` |
| `secondOrder` | Second-order one-sided `snGrad` | no | `false` |
| `value` | Initial patch value. Read if present; otherwise the patch values are left uninitialised, so always give it. | effectively yes | — |

`write()` writes `type`, `interfaceType`, `neighbourFieldName`, `diffusivity`
(if the resolved name is non-empty), `accType`, `relax` (`couplingReuse` for
`IQNILS`), `nonOrthCorr`, `secondOrder` and `value`. `jumpCoefficient` is
**not** written.

## Usage

Direct use only makes sense with a derived class; the base type aborts in
`valueJump()`. The generic entry layout (shown with the concrete scalar BC)
is:

```
// 0/regionA/C
interfaceA
{
    type            regionCoupledScalarJump;   // genericRegionCoupledJump-derived
    interfaceType   (defaultInterface);
    diffusivity     "D";
    accType         aitken;
    relax           0.75;
    nonOrthCorr     false;
    secondOrder     false;
    value           uniform 0;
}
```

The partner on the neighbour region must be a
`genericRegionCoupledFlux`-derived BC of the same field type, e.g.
`regionCoupledScalarFlux`:

```
// 0/regionB/C
interfaceB
{
    type            regionCoupledScalarFlux;
    interfaceType   (defaultInterface);
    diffusivity     "D";
    nonOrthCorr     false;
    secondOrder     false;
    value           uniform 0;
}
```

and the patch pair must be declared in `constant/regionInterfaceProperties`
under `partitionedCoupledPatches`, with the field in `coupledFields`:

```
partitionedCoupledPatches
(
    regionAregionB
    {
        interfaceType (defaultInterface);
        coupledPatchPair ((regionA interfaceA) (regionB interfaceB));
        coupledFields (C);
        defaultInterfaceCoeffs {}
    }
);
```

## Known Limitations

- **Not usable on its own.** `genericRegionCoupledJump` is registered for
  `scalar` and `vector`, but `valueJump()` calls `notImplemented`, which is a
  `FatalError`. `updateCoeffs()` calls it in every coupled iteration, so a case
  using the bare type aborts on the first coupled solve.
- **`diffusivity` effectively required.** The default is `word::null`. Unless
  a derived class overrides `diffusivityName()`, `flux()` (called by the
  partner) then looks up an empty name in `transportProperties` and stops.
- **`neighbourFieldName` has no effect.** Both coupled patches must carry a
  field of the same name.
- **`diffusivity "v"` is not inverted.** The name `v` is mapped to `rho`. The
  following check for `v`, meant to return $1/\rho$, compares the already
  replaced name and is never true, so $\Gamma = \rho$ is used.
- **Vector `snGrad` is tangential only.** For `vector` fields the normal
  component of $\partial_n\boldsymbol\phi$ is removed in `snGrad()` and
  `gradientBoundaryCoeffs()`. This suits the velocity BCs built on it, but not
  a general vector diffusion problem.
- **Acceleration state carries over time steps.** $x^{k-1}$ is not reset at a
  new time step, so `relax < 1` damps the transient change of the interface
  value in the first iteration of every step.
- **`IQNILS` is not parallel-consistent.** The QR products use processor-local
  `sum` instead of `gSum`, so in a decomposed run each processor computes its
  own coefficients.
- **Side detection by patch name only.** The neighbour side is decided by
  comparing the local patch name with `patchA` of the interface. Two coupled
  patches with the same name in both regions would resolve to the wrong
  neighbour (inferred, untested).
- **`jumpCoefficient` is not written** by `write()`, so it is lost after
  `decomposePar`/`reconstructPar` round trips or in written time directories.
- **Jump–jump pairing.** Nothing prevents a jump BC on both sides, but then no
  flux BC exists and `normFluxRes` returns `GREAT` with a warning, so the DNA
  loop always runs `maxCoupleIter` iterations (inferred from the code).
- **Partitioned only.** The monolithic solve is commented out in
  `multiRegionSystem::solve()`.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| — | No tutorial uses the bare `genericRegionCoupledJump` type. Derived BCs are used in `tutorials/conjugateHeatTransfer/*`, `tutorials/electrochemistry/NTGKBatteryWithTabs` and `tutorials/fluidStructureInteraction/heatedBeamInCrossFlow` (see `regionCoupledScalarJump`, `regionCoupledTemperatureJump`). |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
- Küttler, U., Wall, W. A., *Fixed-point fluid–structure interaction solvers with dynamic relaxation*, Computational Mechanics 43, 2008. DOI: [10.1007/s00466-008-0255-5](https://doi.org/10.1007/s00466-008-0255-5)
- Degroote, J., Bathe, K.-J., Vierendeels, J., *Performance of a new partitioned procedure versus a monolithic procedure in fluid–structure interaction*, Computers & Structures 87, 2009. DOI: [10.1016/j.compstruc.2008.11.013](https://doi.org/10.1016/j.compstruc.2008.11.013)
