---
mathjax: true
title: "paraboloidInletVelocity"
description: "Fixed-value inlet velocity with a doubly parabolic (paraboloid) x-velocity profile over the y-z inlet plane and an optional cosine ramp in time."
tags: ["boundary-condition", "uncoupled", "fluid-flow", "inlet"]
weight: 10
---

## Overview

`paraboloidInletVelocity` is a `fixedValue` condition for a velocity field. It
prescribes a velocity in the $+x$ direction whose magnitude is parabolic in
both $y$ and $z$, as found in the developed laminar inflow of a rectangular
channel. The peak velocity can be ramped up smoothly from zero with a
half-cosine. The BC is based on the `setInletVelocity` function object of
Ž. Tuković and was taken over from the foam-extend solid-mechanics toolbox. It
is used for the 3-D FSI benchmark "beam in cross-flow". The profile is written
in global coordinates: the inlet must be normal to $x$, and the channel must be
placed as described below.

## Mathematical Formulation

$$
\mathbf{U}_b = U(t)\;
\frac{y\,(2y_m - y)\,\left(z_m^2 - z^2\right)}{y_m^2\, z_m^2}\;
\mathbf{e}_x
$$

$(y, z)$ are the global coordinates of each face centre. The profile reaches
$U(t)$ at $y = y_m$, $z = 0$. It is zero at $y = 0$, $y = 2y_m$ and
$z = \pm z_m$. The time function is

$$
U(t) =
\begin{cases}
\dfrac{U_\text{max}}{2}\left[1 - \cos\left(\dfrac{\pi\, t}{t_\text{max}}\right)\right] & t < t_\text{end} \\[2ex]
U_\text{max} & t \ge t_\text{end}
\end{cases}
$$

With $t_\text{end} = t_\text{max}$ (as in the tutorials), this is a smooth ramp
from $0$ to $U_\text{max}$, after which the velocity is constant. With
$t_\text{end} > t_\text{max}$, the velocity follows the cosine further (it
oscillates between $0$ and $U_\text{max}$) until $t_\text{end}$, and then jumps
to $U_\text{max}$. With $t_\text{end} < t_\text{max}$, the ramp is cut off and
the velocity jumps to $U_\text{max}$ at $t_\text{end}$. With
$t_\text{end} \le 0$, $U(t) = U_\text{max}$ from the start.

| Symbol | Description | Unit |
|--------|-------------|------|
| $\mathbf{U}_b$ | Patch velocity | m/s |
| $U_\text{max}$ | `maxVelocity`, peak velocity at full load | m/s |
| $y_m$ | `yWidth`: $y$ position of the peak, i.e. half the channel height for a channel $0 \le y \le 2y_m$ | m |
| $z_m$ | `zWidth`: half-width in $z$ for a channel $-z_m \le z \le z_m$ (or the full width with a symmetry plane at $z=0$) | m |
| $t_\text{end}$ | `timeVaryingEndTime` | s |
| $t_\text{max}$ | `timeAtMaxVelocity` | s |
| $t$ | `time().value()` | s |

## Implementation

**`updateCoeffs()`** computes $U(t)$ from the current time value and
evaluates the formula for every face centre. It assigns the result with
`operator==` and calls `fixedValue::updateCoeffs()`. The profile does not
depend on the patch normal or the internal field.

**`snGrad()`** returns $(\mathbf{U}_b - \mathbf{U}_P^{*})\,/\,\delta$, where
$\mathbf{U}_P^{*} = \mathbf{U}_P + \mathbf{k}\cdot(\nabla\mathbf{U})_P$ and
$\mathbf{k} = (\mathbf{I}-\mathbf{n}\mathbf{n})\cdot\boldsymbol{\Delta}$. The
correction is only applied if a `volTensorField` named `grad(<U>)` is
registered. Otherwise the plain $(\mathbf{U}_b-\mathbf{U}_P)/\delta$ is used.

**`gradientBoundaryCoeffs()`** returns
$(\mathbf{U}_b - \mathbf{k}\cdot(\nabla\mathbf{U})_P)/\delta$. If no
`grad(<U>)` field is registered, it computes `fvc::grad(U)` on the fly for
this. `gradientInternalCoeffs()` and the value coefficients are inherited
from `fixedValue`.

A second-order variant exists in both functions, but it is hard-coded off
(`secondOrder_ = false`) and cannot be selected from the dictionary.

**Note on OpenFOAM.com.** `fvPatch::delta()` is patch-normal for non-coupled
patches, so $\mathbf{k} \equiv \mathbf{0}$. The non-orthogonal correction does
nothing, and the gradient coefficients are those of a plain `fixedValue`
condition. The on-the-fly `fvc::grad(U)` in `gradientBoundaryCoeffs()` still
costs time on every matrix assembly.

**`write()`** writes `maxVelocity`, `yWidth`, `zWidth`, `timeVaryingEndTime`,
`timeAtMaxVelocity`, then `type` and `value`.

## Dictionary Entries

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `paraboloidInletVelocity` | yes | – |
| `maxVelocity` | Peak velocity $U_\text{max}$ | yes | – |
| `yWidth` | $y_m$, $y$ coordinate of the peak (half channel height) | yes | – |
| `zWidth` | $z_m$, $z$ half-width of the profile | yes | – |
| `timeVaryingEndTime` | $t_\text{end}$. `0` gives a constant velocity from the start. | yes | – |
| `timeAtMaxVelocity` | $t_\text{max}$, end of the cosine ramp | yes if `timeVaryingEndTime > 0` | `0` (not read otherwise) |
| `value` | Initial patch value (required by `fixedValue`, overwritten at the first update) | yes | – |

## Usage

From `tutorials/fluidStructureInteraction/beamInCrossFlow/0/fluid/orig/partitioned/U`.
The fluid domain is $0 \le y \le 0.4$, $-0.4 \le z \le 0$ with a symmetry
plane at $z = 0$, so `yWidth 0.2` puts the peak at mid-height and `zWidth 0.4`
puts it on the symmetry plane:

```
inlet
{
    type                paraboloidInletVelocity;
    maxVelocity         0.3;
    yWidth              0.2;
    zWidth              0.4;
    timeVaryingEndTime  1.0;
    timeAtMaxVelocity   1.0;
    value               uniform (0 0 0);
}
```

## Known Limitations

- **x-direction and global coordinates only.** The velocity always points
  along $+\mathbf{e}_x$. $y$ and $z$ are absolute coordinates, so the channel
  must start at $y = 0$ and be centred on $z = 0$. The profile cannot be
  translated or rotated.
- **No clipping.** Faces outside $0 \le y \le 2y_m$, $|z| \le z_m$ get a
  negative (reversed) velocity.
- **`timeAtMaxVelocity` must be positive** when `timeVaryingEndTime > 0`.
  Otherwise $U(t)$ divides by zero.
- **Uses `time().value()`**, not the user time (`timeOutputValue()`), unlike
  `solidTraction`'s series.
- **Second-order correction not selectable** (hard-coded off).
- **Non-orthogonal correction has no effect in OpenFOAM.com**, but
  `gradientBoundaryCoeffs()` still computes `fvc::grad(U)` on each call unless
  `grad(U)` is registered (see the note above).

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/fluidStructureInteraction/beamInCrossFlow` | Inlet of the `pimpleFluid` region, ramped to 0.3 m/s over 1 s. |
| `tutorials/fluidStructureInteraction/heatedBeamInCrossFlow` | Same inlet in the thermo-FSI variant (`pimpleFluid transportTemperature`). |

## References

- Cardiff, P., Karač, A., De Jaeger, P., Jasak, H., Nagy, J., Ivanković, A., Tuković, Ž., *An open-source finite volume toolbox for solid mechanics and fluid-solid interaction simulations*, arXiv:1808.10736, 2018.
- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
