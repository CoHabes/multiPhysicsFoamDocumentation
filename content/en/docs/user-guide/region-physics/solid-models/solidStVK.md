---
mathjax: true
title: "solidStVK"
description: "Transient, geometrically nonlinear St. Venant-Kirchhoff solid in a total Lagrangian displacement formulation, intended as the structural side of partitioned FSI."
tags: ["region-type", "solid-mechanics", "fsi", "nonlinear-elasticity", "static-mesh"]
weight: 10
---

## Overview

`solidStVK` solves the transient momentum balance for an elastic solid
undergoing large deformations. It uses the **St. Venant-Kirchhoff**
hyperelastic law and a **total Lagrangian** formulation: the unknown is the
total displacement `D`, and every equation is assembled on the undeformed
(reference) mesh, which is never moved. The discretisation follows the
segregated finite-volume solid-mechanics approach of foam-extend
(Cardiff et al.). It uses an implicit "stiffness" Laplacian with explicit
deferred correction, and a Rhie-Chow-type momentum stabilisation. The region
type is built to act as the structural partner of a fluid region, coupled
through `fsiInterface`. On its own it does not solve anything (see
[Known Limitations](#known-limitations)).

## Governing Equations

### Momentum (total Lagrangian)

The equation assembled in `setCoupledEqns()` is

$$
\rho \frac{\partial^2 \mathbf{D}}{\partial t^2} =
\underbrace{\nabla \cdot \left( K \nabla \mathbf{D} \right)}_{\text{implicit}} - \underbrace{\nabla \cdot \left( K \nabla \mathbf{D} \right)}_{\text{explicit}} + \nabla \cdot \left( J\, \mathbf{F}^{-1} \cdot \boldsymbol{\sigma} \right) + \mathbf{S}_\text{stab}
$$

with the implicit stiffness $K = 2\mu + \lambda$. The implicit and explicit
Laplacians cancel at convergence, so the converged solution satisfies
$\rho\,\ddot{\mathbf{D}} = \nabla_0 \cdot \mathbf{P}^T$ with the first
Piola-Kirchhoff stress $\mathbf{P} = J \boldsymbol{\sigma} \mathbf{F}^{-T}$.
No body force (gravity) term is included.

### Kinematics

$$
\mathbf{F} = \mathbf{I} + (\nabla \mathbf{D})^T, \qquad
J = \det \mathbf{F}
$$

### St. Venant-Kirchhoff constitutive law

$$
\mathbf{C} = \mathbf{F}^T \mathbf{F}, \qquad
\mathbf{E} = \tfrac{1}{2}\left( \mathbf{C} - \mathbf{I} \right), \qquad
\mathbf{S} = 2\mu\,\mathbf{E} + \lambda\,\mathrm{tr}(\mathbf{E})\,\mathbf{I}
$$

$$
\boldsymbol{\sigma} = \frac{1}{J}\, \mathbf{F}\, \mathbf{S}\, \mathbf{F}^T
$$

### Lamé parameters

$$
\mu = \frac{E}{2(1+\nu)}, \qquad
\lambda =
\begin{cases}
\dfrac{\nu E}{(1+\nu)(1-2\nu)} & \text{plane strain / 3-D} \\[2ex]
\dfrac{\nu E}{(1+\nu)(1-\nu)} & \text{plane stress (\texttt{planeStress yes})}
\end{cases}
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $\mathbf{D}$ | Total displacement | m |
| $\rho$ | Solid density | kg/m³ |
| $E$ | Young's modulus | Pa |
| $\nu$ | Poisson's ratio | – |
| $\mu,\ \lambda$ | Lamé parameters | Pa |
| $K$ | Implicit stiffness $2\mu+\lambda$ (`impK`, `impKf`) | Pa |
| $\mathbf{F}$ | Deformation gradient | – |
| $J$ | Jacobian $\det\mathbf{F}$ | – |
| $\mathbf{E}$ | Green-Lagrange strain | – |
| $\mathbf{S}$ | 2nd Piola-Kirchhoff stress | Pa |
| $\boldsymbol{\sigma}$ | Cauchy stress (`sigma`) | Pa |
| $\mathbf{S}_\text{stab}$ | Momentum stabilisation term | Pa/m |

### Momentum stabilisation

$\mathbf{S}_\text{stab}$ is evaluated explicitly on the displacement
**increment** `DD` $= \mathbf{D} - \mathbf{D}^\text{old}$ and its gradient
`gradDD`. Its diffusivity is $\gamma = s\,K$, where $s$ is the scale factor,
interpolated to faces as $\gamma_f$:

| `type` | Term |
|--------|------|
| `RhieChow` (default) | $\nabla\cdot(\gamma_f \nabla \mathbf{DD}) - \nabla\cdot\left(\gamma_f\, \mathbf{S}_f \cdot (\nabla\mathbf{DD})_f\right)$, i.e. compact minus non-compact Laplacian |
| `JamesonSchmidtTurkel` | $-\nabla\cdot\left(\lvert\mathbf{S}_f\rvert\, \nabla \left[\nabla\cdot(\gamma_f \nabla \mathbf{DD})\right]\right)$ (fourth-order damping) |
| `Laplacian` | $\nabla\cdot(\gamma_f \nabla \mathbf{DD})$ |
| `none` | $\mathbf{0}$ |

On faces where $K$ differs between owner and neighbour (material
interfaces), $\gamma_f$ is replaced by
`interfaceScaleFactor`·$\tfrac{1}{2}(K_\text{own}+K_\text{nei})$. On coupled
non-processor patches this replacement is applied to every face.

## Numerical Algorithm

`usesPIMPLE()` returns `false`. `solveRegion()`, `prePredictor()`,
`momentumPredictor()`, `pressureCorrector()` and `meshMotionCorrector()` are
all no-ops. `selfCoupledFieldNames()` is not overridden, so it returns an empty
list. The displacement equation is solved **only** when `D` appears in the
DNA (partitioned) coupled field list. In that case each DNA iteration of
`multiRegionSystem::solve()` runs:

1. `solvePIMPLE()` on the fluid region(s). `D` is a special FSI field, and
   stand-alone PIMPLE is skipped whenever `D` is DNA-coupled.
2. `assembleAndSolveEqns<fvMatrix, vector>("D")`. For `solidStVK` this is a
   Picard loop:
   1. `setCoupledEqns("D")` stores `D.prevIter()` and assembles the momentum
      equation above. It registers the equation as
      `D<mesh>MeshsolidStVKTypeEqn`.
   2. The equation is relaxed (`relaxationFactors.equations`) and solved.
   3. `postSolve()` relaxes the field `D` (`relaxationFactors.fields`). It then
      updates `DD`, `gradD = fvc::grad(D)`, `gradDD`, `F`, `Finv`, `J`,
      `sigma`, interpolates `D` to `pointD` (`volPointInterpolation`), updates
      `pointDD` and sets the solid velocity `U = fvc::ddt(D)`.
   4. The loop repeats until the region's `PICARD.D` criteria are met
      (`maxCorr`, `tolerance`, `relTol`). If no `PICARD` entry is given, it
      runs one pass per DNA iteration.
3. `meshMotionCorrector()` on all regions. The fluid mesh follows `pointD` of
   the solid (via `interfaceTrackingFvMesh`). Then `interfaces_->update()`.
4. DNA convergence is checked with the `fsiInterface` residual: the
   mismatch between the deformed solid interface points
   (`points + pointD`) and the fluid interface points.

`correct()`, which runs from `regionTypeList::preSolve()` once per time step,
only calls `D.correctBoundaryConditions()`.

`getMinDeltaT()` returns `GREAT`, so the solid never limits the time step.

## Sub-Models

| Sub-model | Selection mechanism | Notes |
|-----------|---------------------|-------|
| Momentum stabilisation | `stabilisation` sub-dict in `mechanicalProperties` | `none`, `RhieChow`, `JamesonSchmidtTurkel`, `Laplacian` |

No `fvOptions`, thermal-stress, plasticity or other constitutive models are
available. The material law is fixed to St. Venant-Kirchhoff.

## Required Fields

| Field | Type | Description | Unit |
|-------|------|-------------|------|
| `D` | `volVectorField` | Total displacement (**must be read**) | m |

Created at runtime. All of these are read if present and auto-written unless
stated otherwise:

| Field | Type | Description | Unit | Default |
|-------|------|-------------|------|---------|
| `DD` | `volVectorField` | Displacement increment over the time step | m | 0 |
| `gradD` | `volTensorField` | Displacement gradient | – | 0 |
| `gradDD` | `volTensorField` | Gradient of displacement increment | – | 0 |
| `F` | `volTensorField` | Deformation gradient | – | `I` |
| `Finv` | `volTensorField` | Inverse deformation gradient | – | `I` |
| `J` | `volScalarField` | Jacobian | – | 1 |
| `sigma` | `volSymmTensorField` | Cauchy stress | Pa | 0 |
| `U` | `volVectorField` | Solid velocity $\partial\mathbf{D}/\partial t$ | m/s | 0 |
| `impK` | `volScalarField` | Implicit stiffness $2\mu+\lambda$ | Pa | uniform |
| `impKf` | `surfaceScalarField` | Face implicit stiffness | Pa | uniform |
| `rImpK` | `volScalarField` | $1/K$ (read if present, **not written**) | 1/Pa | $1/K$ |
| `pointD` | `pointVectorField` | Displacement interpolated to points | m | from `D` |
| `pointDD` | `pointVectorField` | Point displacement increment | m | from `DD` |

`pointD` is the field that `fsiInterface`, `interfaceTrackingFvMesh` and the
`solidPointDisplacement` function object look up. Its presence is also how
`fsiInterface` tells which side of the interface is the solid.

## Case Setup

### Required dictionaries

```
constant/<solid>/mechanicalProperties   – rho, E, nu, planeStress, optional stabilisation
constant/<solid>/dynamicMeshDict        – staticFvMesh (the solid mesh is never moved)
system/<solid>/fvSolution               – D solver, relaxation, PICARD loop for D
system/<solid>/fvSchemes                – d2dt2, ddt, grad and the named laplacian schemes
```

### `mechanicalProperties`

```
planeStress     no;                         // required; yes => plane-stress lambda

rho     rho [1 -3 0 0 0 0 0]  1000;         // required
nu      nu  [0 0 0 0 0 0 0]   0.4;          // required
E       E   [1 -1 -2 0 0 0 0] 5.6e6;        // required

// Optional. If absent: { type RhieChow; scaleFactor 0.1; }
stabilisation
{
    type                 RhieChow;  // none | RhieChow | JamesonSchmidtTurkel | Laplacian
    scaleFactor          0.1;       // required unless type none or useCellZones yes
    useCellZones         false;     // default false; then <cellZoneName>ScaleFactor per zone
    interfaceScaleFactor 0.01;      // default 0.01; used on material-interface faces
}
```

### `multiRegionProperties` entry

```
regions
(
    ( fluid (pimpleFluid) )
    ( solid (solidStVK) )
);

DNA
{
    D
    {
        maxCoupleIter 100;

        residualControl
        {
            maxIntRes          1e-7;
            outputIntResField  yes;
        }
    }
}
```

`solidStVK` can be stacked with other region types on the same mesh, e.g.
`( solid (solidStVK conductTemperature) )` in `heatedBeamInCrossFlow`.

### Key `fvSolution` settings

```
solvers
{
    "D|DD"
    {
        solver          PCG;
        preconditioner  FDIC;
        tolerance       1e-09;
        relTol          0.1;
    }
}

relaxationFactors
{
    fields
    {
        "DD|D"  0.9;        // applied to D in postSolve()
    }
}

PICARD
{
    D
    {
        maxCorr         1000;   // Picard (outer) corrections per DNA iteration
        tolerance       1e-09;
        relTol          1e-09;
        infoFrequency   100;
    }
}
```

### Key `fvSchemes` settings

The momentum equation and the `RhieChow` stabilisation use the named
Laplacian scheme `laplacian(DD,D)`:

```
d2dt2Schemes        { default backward; }
ddtSchemes          { default backward; }
gradSchemes         { default leastSquares; }
divSchemes          { default Gauss linear; }

laplacianSchemes
{
    default            none;
    laplacian(DD,D)    Gauss linear corrected;
}

interpolationSchemes { default linear; }    // also used for interpolate(impK)
```

With `default none`, the `JamesonSchmidtTurkel` stabilisation also needs
`JSTinner` and `JSTouter` Laplacian entries. The `Laplacian` stabilisation
needs an entry for its default scheme name.

### Boundary conditions for `D`

| Patch role | BC type |
|------------|---------|
| FSI interface (solid side) | `regionCoupledTraction` (with `interfaceType (fsiInterface)`) |
| Clamped | `fixedValue` |
| Symmetry plane | `solidSymmetry` |
| Prescribed traction / pressure (uncoupled) | `solidTraction` |

```
interface
{
    type                regionCoupledTraction;
    neighbourFieldName  sigma;
    interfaceType       (fsiInterface);
    secondOrder         false;
    nonOrthCorr         false;
    value               uniform (0 0 0);
}
```

`regionCoupledTraction` and `solidTraction` look up `impK`, `rImpK`, `gradD`,
`sigma` and `Finv` from the solid region. All of these are provided by
`solidStVK`.

## Coupling

| Interface type | Exchanged fields | Coupling modes |
|----------------|-----------------|----------------|
| `fsiInterface` | Fluid traction `sigma·n` → solid (`regionCoupledTraction` on `D`); solid `pointD` → fluid mesh motion (`interfaceTrackingFvMesh`) | partitioned (DNA) |
| `heatTransferInterface` | `T`, heat flux (through a stacked `conductTemperature`, not `solidStVK` itself) | partitioned (DNA) |

Notes on `fsiInterface` with this region type:

- The interface mapping is built with the solid mesh temporarily moved to its
  deformed configuration (`points + pointD`). The mesh is then moved back, so
  the solid stays total Lagrangian.
- The fluid traction is mapped to the solid patch and its sign is flipped.
  The deformed normal is computed as $\mathbf{n} = \mathbf{F}^{-T}\mathbf{n}_0 / |\cdot|$.
- The interface residual for DNA convergence is the point-position mismatch
  between the deformed solid patch and the fluid patch.

## Known Limitations

- **Requires DNA coupling of `D`** — `solveRegion()` is a no-op and
  `selfCoupledFieldNames()` returns an empty list. If `D` is not listed in a
  `partitionedCoupledPatches` `coupledFields` entry (and in the `DNA` dict),
  the displacement equation is never solved. Single-region (stand-alone
  solid) runs are therefore not supported.
- **Partitioned only** — there is no monolithic FSI path (monolithic
  coupling is disabled in `multiRegionSystem`).
- **St. Venant-Kirchhoff only** — linear elastic in Green strain. It is not
  suited to large compressive strains. No other material laws, thermal
  expansion, or body forces (gravity) are implemented.
- **Constant, uniform material properties** — `rho`, `E`, `nu` are single
  `dimensionedScalar`s read from `mechanicalProperties`.
- **No time-step control** — `getMinDeltaT()` returns `GREAT`.
- **Static solid mesh** — `meshMotionCorrector()` is a no-op. The
  `staticFvMesh` reference configuration is used throughout, and deformation
  is represented only through `D`/`pointD`.
- **Unused code** — `solidStVKDdtCoeffs.H` (ddt-scheme coefficients) and
  the `threeK` member are present but not used anywhere.

## Tutorials

| Tutorial path | Description |
|---------------|-------------|
| `tutorials/fluidStructureInteraction/HronTurekFsi3` | Hron-Turek FSI3 benchmark (pimpleFluid + solidStVK) |
| `tutorials/fluidStructureInteraction/beamInCrossFlow` | Elastic beam in cross-flow (pimpleFluid + solidStVK) |
| `tutorials/fluidStructureInteraction/heatedBeamInCrossFlow` | Heated elastic beam in cross-flow; FSI plus conjugate heat transfer (pimpleFluid + transportTemperature / solidStVK + conductTemperature) |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
