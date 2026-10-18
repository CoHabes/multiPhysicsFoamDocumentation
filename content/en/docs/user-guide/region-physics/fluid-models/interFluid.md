---
mathjax: true
title: "interFluid"
description: "Transient, incompressible, isothermal two-phase flow of immiscible fluids using algebraic VOF interface capturing (MULES) and the PIMPLE algorithm."
tags: ["region-type", "fluid", "incompressible", "multiphase", "VOF", "moving-mesh", "fsi"]
weight: 10
---

## Overview

`interFluid` solves the flow of two incompressible, isothermal, immiscible
fluids (for example water and air) on one mesh. It captures the interface with
the algebraic **Volume-of-Fluid (VOF)** method: one phase-fraction field
$\alpha_1$ is transported with a MULES-limited, interface-compressed equation,
and a single mixture momentum equation is solved with the **PIMPLE**
algorithm. The momentum equation is written in terms of the pressure minus its
hydrostatic part, $p_{rgh}$. The region type is a port of OpenFOAM's
`interFoam` solver into the `regionType` framework, built on the same pattern
as `pimpleFluid`. It supports laminar, RANS and LES turbulence modelling,
surface tension (CSF), MRF, `fvOptions` and mesh motion through
`dynamicFvMesh`. It can run standalone, or as the fluid side of a
**fluid-structure interaction** (FSI) case. For FSI it provides the Cauchy
stress tensor `sigma` to `fsiInterface`, it supports mesh motion driven by the
solid's displacement, and it is laminar only. It cannot take part in any other
kind of inter-region coupling yet.

## Governing Equations

### Continuity

$$
\nabla \cdot \mathbf{U} = 0
$$

### Phase-fraction transport

$$
\frac{\partial \alpha_1}{\partial t} + \nabla \cdot \left(\alpha_1 \mathbf{U}\right) + \nabla \cdot \left(\alpha_1 (1 - \alpha_1)\,\mathbf{U}_r\right) = 0,
\qquad \alpha_2 = 1 - \alpha_1
$$

The third term is an artificial **interface-compression** term. It is non-zero
only inside the transition region ($0 < \alpha_1 < 1$) and works against the
numerical smearing of the interface. The compression velocity is never
assembled as a vector field. The code uses its face flux directly:

$$
\phi_r = c_\alpha \left| \frac{\phi}{|\mathbf{S}_f|} \right| \hat{\mathbf{n}}_f \cdot \mathbf{S}_f
$$

$\phi_r$ is set to zero on all non-coupled boundary faces, so compression
never moves phase across a domain boundary.

### Mixture properties

$$
\rho = \alpha_1 \rho_1 + \alpha_2 \rho_2,
\qquad
\mu = \alpha_1 \rho_1 \nu_1 + \alpha_2 \rho_2 \nu_2,
\qquad
\nu = \mu / \rho
$$

### Momentum

$$
\frac{\partial (\rho \mathbf{U})}{\partial t} + \nabla \cdot \left(\rho \mathbf{U}\mathbf{U}\right) - \nabla \cdot \mathbf{R}
= -\nabla p_{rgh} - (\mathbf{g} \cdot \mathbf{h})\,\nabla \rho + \mathbf{F}_\sigma + \mathbf{S}_U
$$

In the discrete momentum equation the convective term uses the mass flux
`rhoPhi`. This flux is rebuilt from the bounded MULES phase flux (see below),
so momentum and phase transport stay consistent.

### Surface tension (CSF)

$$
\mathbf{F}_\sigma = \sigma \kappa \nabla \alpha_1,
\qquad
\kappa = -\nabla \cdot \hat{\mathbf{n}},
\qquad
\hat{\mathbf{n}}_f = \frac{(\nabla \alpha_1)_f}{|(\nabla \alpha_1)_f| + \delta_N}
$$

At walls with an `alphaContactAngle`-type boundary condition, the interface
normal is rotated so that the interface meets the wall at the prescribed
contact angle.

### Gravity decomposition and pressure recovery

$$
p_{rgh} = p - \rho\, \mathbf{g} \cdot \mathbf{h},
\qquad
gh = \mathbf{g} \cdot \mathbf{x}_C - gh_\text{ref},
\qquad
gh_f = \mathbf{g} \cdot \mathbf{x}_{C_f} - gh_\text{ref}
$$

The solver solves for $p_{rgh}$. The static pressure $p = p_{rgh} + \rho\,gh$
is recovered after each pressure correction.

### Stress tensor (for FSI coupling)

After each `pressureCorrector()` call (once per outer corrector) the Cauchy
stress tensor is updated:

$$
\boldsymbol{\sigma} = -p\,\mathbf{I} + \mu\left[\nabla\mathbf{U} + (\nabla\mathbf{U})^T - \tfrac{2}{3}(\nabla\cdot\mathbf{U})\,\mathbf{I}\right]
$$

The bracket is OpenFOAM's `devTwoSymm(fvc::grad(U))`. Since
$\nabla\cdot\mathbf{U} = 0$, its trace term vanishes up to discretisation
error. The pressure is $p$,
which is already dimensional, so unlike `pimpleFluid` no $\rho$ factor is
needed. $\mu$ is the local **mixture** dynamic viscosity, so the viscous
traction switches between water and air values across the interface. The
field is `sigma`, and `fsiInterface` uses it to compute the traction on the
solid side (via `regionCoupledTraction`).

This is the stress formula of `pimpleFluid` with only the molecular viscosity
kept:
$\rho\left(-p_\text{kin}\mathbf{I} - \mathbf{R}^\text{dev}\right) = -p\,\mathbf{I} + \mu_\text{eff}\,\texttt{devTwoSymm}(\nabla\mathbf{U})$
with $\mu_\text{eff} \to \mu$. The turbulent (Reynolds-stress) part is
**not** included, because the VOF turbulence wrapper
(`incompressibleInterPhaseTransportModel`) has no public `devRhoReff()`. The
expression is exact for laminar flow.

```note
The field `sigma` (Cauchy stress, Pa) is not the same thing as the entry
`sigma` in `transportProperties` (surface tension coefficient, N/m). They
only share a name.
```

| Symbol | Description | Unit |
|--------|-------------|------|
| $\alpha_1$, $\alpha_2$ | Phase fraction of phase 1 / phase 2 | – |
| $\mathbf{U}$ | Mixture velocity (shared by both phases) | m/s |
| $\phi$ | Volumetric face flux | m³/s |
| $p$ | Static pressure | Pa |
| $p_{rgh}$ | Pressure minus hydrostatic contribution | Pa |
| $\rho_1$, $\rho_2$ / $\rho$ | Phase densities / mixture density | kg/m³ |
| $\nu_1$, $\nu_2$ / $\mu$ | Phase kinematic viscosities / mixture dynamic viscosity | m²/s / Pa·s |
| $\mathbf{R}$ | Effective stress tensor (viscous + turbulent) | Pa |
| $\mathbf{U}_r$, $\phi_r$ | Compression velocity / compression flux | m/s / m³/s |
| $c_\alpha$ | Interface-compression coefficient (`cAlpha`) | – |
| $\sigma$ | Surface tension coefficient | N/m |
| $\kappa$ | Interface curvature | 1/m |
| $\hat{\mathbf{n}}$ | Interface unit normal | – |
| $\mathbf{g}$ | Gravitational acceleration | m/s² |
| $h_\text{ref}$ | Reference height (`constant/hRef`, optional) | m |
| $\mathbf{F}_\sigma$ | Surface tension force | N/m³ |
| $\mathbf{S}_U$ | Momentum source (from `fvOptions`) | N/m³ |
| $\boldsymbol{\sigma}$ | Cauchy stress tensor (FSI) | Pa |

## Numerical Algorithm

The framework drives `interFluid` through the `regionType` PIMPLE phases.
For every outer corrector (`nOuterCorrectors`), `regionTypeList::solvePIMPLE()`
calls `prePredictor()`, then `momentumPredictor()`, then
`pressureCorrector()`. This is the same order as the loop body of stock
`interFoam`.

1. **Mesh motion** (`prePredictor()`, dynamic meshes only). Runs
   `mesh().controlledUpdate()` once per time step, or on every outer corrector
   if `moveMeshOuterCorrectors yes`. Then it recomputes `gh`/`ghf`, updates
   MRF and, if `correctPhi` is on, restores a divergence-free flux on the moved
   mesh. It also corrects `Uf`.
2. **Phase-fraction transport** (`prePredictor()`). Re-reads the alpha
   controls, then solves the $\alpha_1$ equation with MULES. There are two
   options:
   - **explicit MULES** (`MULESCorr no`), which applies a limited flux in one
     step, or
   - **semi-implicit MULES** (`MULESCorr yes`), which first solves an implicit,
     upwind-bounded predictor and then applies `nAlphaCorr` limited
     high-order correction passes.

   With `nAlphaSubCycles > 1` this step is sub-cycled and the mass flux
   `rhoPhi` is averaged over time across the sub-cycles. Afterwards
   $\rho = \alpha_1\rho_1 + \alpha_2\rho_2$ is updated, and
   `mixture.correct()` recomputes the normals, curvature and surface-tension
   force.
3. **Momentum predictor** (`momentumPredictor()`). Builds the momentum matrix
   with `rho`/`rhoPhi` and caches it. If `momentumPredictor yes`, it solves the
   matrix against the explicit surface-tension, buoyancy and $p_{rgh}$-gradient
   forcing, reconstructed from face fluxes.
4. **Pressure correction, PISO loop** (`pressureCorrector()`). Runs
   `nCorrectors` times. It builds $r_{AU} = 1/A$, $\mathbf{H}/A$ and a
   Rhie–Chow-consistent `phiHbyA` (including `ddtCorr`, which is aware of mesh
   motion). It adds the surface-tension and gravity flux:

$$
\phi_g = \left(\mathbf{F}_{\sigma,f} - gh_f\,\nabla^\perp_f \rho\right) r_{AU,f}\,|\mathbf{S}_f|
$$

   and solves the pressure equation:

$$
\nabla \cdot \left( r_{AU,f}\,\nabla p_{rgh} \right) = \nabla \cdot \phi_{HbyA}
$$

   The non-orthogonal corrector sub-loop is controlled by
   `nNonOrthogonalCorrectors`. On the final non-orthogonal iteration the
   solver corrects the flux and the velocity:

$$
\phi = \phi_{HbyA} - \left(r_{AU,f}\nabla p_{rgh}\right)_f,
\qquad
\mathbf{U} = \mathbf{H}/A + r_{AU}\,\mathcal{R}\!\left(\frac{\phi_g - \left(r_{AU,f}\nabla p_{rgh}\right)_f}{r_{AU,f}}\right)
$$

   where $\mathcal{R}$ is the face-to-cell reconstruction (`fvc::reconstruct`).
   It then reports continuity errors, corrects `Uf`, makes the flux relative to
   the mesh motion and recovers $p = p_{rgh} + \rho\,gh$. For closed domains it
   also applies the reference-pressure shift.
5. **Turbulence correction** (end of `pressureCorrector()`). Calls
   `turbulence->correct()` when `pimple.turbCorr()` is true.
6. **Stress update** (end of `pressureCorrector()`). Recomputes `sigma` from
   the new $p$ and $\mathbf{U}$ for FSI coupling. This runs on every outer
   corrector, including in standalone runs.

The scheme is fully segregated. $\alpha_1$, $\mathbf{U}$/$p_{rgh}$ and the
turbulence fields are each solved with the other fields frozen, and the outer
correctors provide the non-linear coupling between them.

### Time-step control

`getMinDeltaT()` computes the flow Courant number over all cells and the
interface Courant number $\mathrm{Co}_\alpha$ over cells with
$0.01 < \alpha_1 < 0.99$. It returns:

$$
\Delta t_\text{stab} = \Delta t \cdot \min\!\left(
\frac{\mathrm{Co}_\text{max}}{\mathrm{Co} + \epsilon},\;
\frac{\mathrm{Co}_{\alpha,\text{max}}}{\mathrm{Co}_\alpha + \epsilon}
\right)
$$

The solver's `setDeltaT.H` then limits how fast the step can grow (at most a
factor of 1.2 per step, damped) and caps it at `maxDeltaT`.

## Sub-Models

| Sub-model | Selection mechanism | Notes |
|-----------|---------------------|-------|
| Two-phase mixture | `constant/transportProperties` | `immiscibleIncompressibleTwoPhaseMixture`; `phases (<phase1> <phase2>)` plus a per-phase `transportModel`/`nu`/`rho` sub-dict |
| Surface tension | `constant/transportProperties` (`sigma`) | Run-time selectable `surfaceTensionModel`; a constant value is enough |
| Turbulence (laminar / RANS / LES) | `constant/turbulenceProperties` | `incompressibleInterPhaseTransportModel`, evaluated on mixture properties |
| Gravity / reference height | `constant/g`, `constant/hRef` (optional) | `hRef` defaults to 0 |
| Moving Reference Frame (MRF) | `constant/MRFProperties` | `IOMRFZoneList` |
| `fvOptions` | `system/fvOptions` | Explicit/implicit momentum sources, porosity |
| Mesh motion | `constant/dynamicMeshDict` | Any stock `dynamicFvMesh` + `motionSolver` (e.g. `dynamicMotionSolverFvMesh`) |

## Required Fields

| Field | Type | Description | Unit |
|-------|------|-------------|------|
| `alpha.<phase1>` | `volScalarField` | Phase fraction of the first phase listed in `phases` (e.g. `alpha.water`) | – |
| `U` | `volVectorField` | Velocity | m/s |
| `p_rgh` | `volScalarField` | Pressure minus hydrostatic contribution | Pa |
| Turbulence fields | – | As required by the selected model (e.g. `k`, `epsilon`, `omega`, `nut`) | – |

Additionally created at runtime:

| Field | Type | Description | Condition |
|-------|------|-------------|-----------|
| `phi` | `surfaceScalarField` | Volumetric face flux (written) | always |
| `p` | `volScalarField` | Static pressure $p_{rgh} + \rho\,gh$ (written) | always |
| `alphaPhi0.<phase1>` | `surfaceScalarField` | MULES phase flux (written; read back on restart) | always |
| `rho`, `rhoPhi` | vol / surface scalar | Mixture density and mass flux (not written) | always |
| `gh`, `ghf` | vol / surface scalar | $\mathbf{g}\cdot\mathbf{h}$ at cell and face centres (not written) | always |
| `alphaPhiUn` | `surfaceScalarField` | Unlimited compressed flux, registered for the `scalarTransport` function object (not written) | always |
| `sigma` | `volSymmTensorField` | Cauchy stress tensor $-p\mathbf{I} + \mu\,\texttt{devTwoSymm}(\nabla\mathbf{U})$ (written) | always (FSI) |
| `Uf` | `surfaceVectorField` | Face velocity | moving mesh only |
| `uniform/cumulativeContErr` | `uniformDimensionedScalarField` | Cumulative continuity error | always |

## Case Setup

### Required dictionaries

```
constant/transportProperties    – phases, per-phase nu/rho, sigma
constant/turbulenceProperties   – turbulence model selection
constant/g                      – gravitational acceleration
system/fvSolution               – alpha solver controls, PIMPLE sub-dict, linear solvers
system/fvSchemes                – discretisation schemes
```

Optional: `constant/hRef`, `constant/MRFProperties`, `system/fvOptions`,
`constant/dynamicMeshDict`.

### `transportProperties`

```
phases          (water air);

water
{
    transportModel  Newtonian;
    nu              1e-06;
    rho             1000;
}

air
{
    transportModel  Newtonian;
    nu              1.48e-05;
    rho             1;
}

sigma           0.07;
```

### `multiRegionProperties` entry

```
regions
(
    ( fluid (interFluid) )
);
```

For FSI, pair it with a solid region:

```
regions
(
    ( fluid (interFluid) )
    ( solid (solidStVK) )
);
```

### `controlDict` (singleRegionMode)

When running as a standalone single-region case via `singleRegionMode`:

```
application         multiPhysicsFoam;

singleRegionMode    yes;
regionType          interFluid;
```

### Key `fvSolution` settings

The alpha controls are read from the solver entry that matches
`alpha.<phase1>` on every time step. `nAlphaCorr` and `nAlphaSubCycles` are
**mandatory**.

```
solvers
{
    "alpha.water.*"
    {
        nAlphaCorr          2;      // MULES correction passes (mandatory)
        nAlphaSubCycles     1;      // alpha sub-cycles per time step (mandatory)
        cAlpha              1;      // interface compression (0 = off)

        MULESCorr           yes;    // semi-implicit MULES predictor-corrector
        nLimiterIter        5;      // MULES limiter iterations
        alphaApplyPrevCorr  no;     // reuse previous correction flux
        icAlpha             0;      // isotropic compression coefficient
        scAlpha             0;      // shear compression coefficient

        solver              smoothSolver;   // used for the MULESCorr predictor
        smoother            symGaussSeidel;
        tolerance           1e-8;
        relTol              0;
    }

    "pcorr.*"           { solver PCG; preconditioner DIC; tolerance 1e-5; relTol 0; }
    p_rgh               { solver PCG; preconditioner DIC; tolerance 1e-7; relTol 0.05; }
    p_rghFinal          { $p_rgh; relTol 0; }
    U                   { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-6; relTol 0; }
}

PIMPLE
{
    momentumPredictor        no;
    nOuterCorrectors         1;
    nCorrectors              3;
    nNonOrthogonalCorrectors 0;

    // Moving mesh
    correctPhi               yes;   // defaults to yes if mesh is dynamic
    checkMeshCourantNo       no;
    moveMeshOuterCorrectors  no;

    // Pressure reference (closed domains)
    p_rghRefCell             0;
    pRefValue                0;
}
```

### Key `fvSchemes` entries

```
divSchemes
{
    div(rhoPhi,U)    Gauss linearUpwind grad(U);
    div(phi,alpha)   Gauss vanLeer;     // bounded TVD scheme for advection
    div(phirb,alpha) Gauss linear;      // compression term
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
```

`ddt(alpha)` accepts `Euler` or `CrankNicolson` (with an off-centering
coefficient). `localEuler` (LTS) is not supported.

### Adjustable time stepping (`controlDict`)

```
adjustTimeStep  yes;
maxCo           1;          // maximum flow Courant number
maxAlphaCo      1;          // maximum interface Courant number
maxDeltaT       1;          // upper bound on time-step size
```

## Coupling

`interFluid` participates in the following multi-region interfaces:

| Interface type | Exchanged fields | Coupling modes |
|----------------|-----------------|----------------|
| `fsiInterface` | Displacement `D`, traction `sigma·n` | partitioned (DNA) |

FSI works the same way as for `pimpleFluid`, and nothing FSI-specific is
needed in the fluid's own boundary conditions:

- **Fluid → solid (traction).** The solid's `regionCoupledTraction` boundary
  condition looks up the fluid's registered `sigma` field by name on the
  neighbour mesh and takes `sigma·n` on the interface patch. `sigma` only has
  to exist under that name, so it does not need to be listed in
  `coupledFields`. Because $p$ includes the hydrostatic part ($p_{rgh} +
  \rho\,gh$) and $\mu$ is the mixture viscosity, the traction is correct
  whether the wall is wetted by water or by air.
- **Solid → fluid (displacement).** The fluid mesh follows the solid through
  a `dynamicMeshDict` (a `dynamicMotionSolverFvMesh` plus a motion solver),
  set up as in the `pimpleFluid` FSI tutorials (e.g. `HronTurekFsi3`). The DNA
  loop calls `meshMotionCorrector()` for the `D` field. The fluid interface patch uses
  standard boundary conditions, e.g. `movingWallVelocity` for `U` and
  `fixedFluxPressure` for `p_rgh`.

`setCoupledEqns` is still a no-op, and the region registers no
self-coupled fields of its own. Other interface types
(`heatTransferInterface`, `capillaryInterface`, `massTransferInterface`)
are not supported, since the region has no temperature field and no
interface-coupled boundary conditions.

## Known Limitations

- **FSI is the only coupling.** Only `fsiInterface` (partitioned, DNA) is
  supported. There is no monolithic coupling and no thermal or
  mass-transfer coupling (see [Coupling](#coupling)).
- **Laminar FSI traction.** `sigma` contains only the molecular viscous
  stress. With a RAS/LES model the fluid flow itself is turbulent, but the
  traction passed to the solid misses the Reynolds-stress part. That part is
  usually small compared with the pressure at a wall, but it is not zero.
- **No FSI tutorial yet.** Coupling with `solidStVK` has been shown to run,
  but there is no complete, validated VOF-FSI tutorial in the repository yet.
- **Two incompressible phases only.** Both densities are constant. There is
  no compressibility, no more than two phases, and no phase change (the
  `Su`/`Sp` alpha sources are zero). For phase change across a sharp
  interface, `interTrackFluid` with `massTransferInterface` is intended
  (currently disabled in `src/multiRegionSystem/Make/files`).
- **Isothermal.** No energy equation is solved.
- **No local time stepping (LTS).** The `localEuler` branch in the alpha
  equation is kept for fidelity with the stock code, but it is dead code.
- **Simplified `correctPhi` after mesh motion.** The flux correction uses a
  constant `rAUf = 1` instead of stock `interFoam`'s time-persisted
  `1/UEqn.A()`. The PISO loop still enforces continuity on the new mesh; only
  the initial flux guess is affected.
- **Mesh motion not yet validated.** The code path follows `pimpleFluid`, but
  no finished tutorial uses a `dynamicMeshDict` yet.
- **`sigma` is always computed and written.** This adds a `grad(U)`
  evaluation per outer corrector and one more field per write time, even in
  standalone runs where nothing reads it.
- **No wave generation.** `waveModels` is linked but not used by this region.
- **Diffuse interface.** As with any algebraic VOF/CSF method, the interface
  is smeared over a few cells and parasitic currents can appear near the
  interface in surface-tension-dominated flows.

## Tutorials

| Tutorial path | Description |
|---------------|-------------|
| `tutorials/singleRegion/interFluid/damBreak` | 2D laminar collapse of a water column over an obstacle, adapted from `interFoam/laminar/damBreak` (static mesh, single-region, adjustable time step) |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
- Rusche, H., *Computational Fluid Dynamics of Dispersed Two-Phase Flows at High Phase Fractions*, PhD thesis, Imperial College London, 2002.
- Brackbill, J. U., Kothe, D. B., Zemach, C., *A continuum method for modeling surface tension*, Journal of Computational Physics 100(2), 335–354, 1992.
  DOI: [10.1016/0021-9991(92)90240-Y](https://doi.org/10.1016/0021-9991(92)90240-Y)
- Deshpande, S. S., Anumolu, L., Trujillo, M. F., *Evaluating the performance of the two-phase flow solver interFoam*, Computational Science & Discovery 5, 014016, 2012.
  DOI: [10.1088/1749-4699/5/1/014016](https://doi.org/10.1088/1749-4699/5/1/014016)
