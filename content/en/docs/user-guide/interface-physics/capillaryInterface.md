---
mathjax: true
title: "capillaryInterface"
description: "Partitioned p-U coupling of two immiscible incompressible fluids across a sharp, mesh-tracked interface with constant surface tension."
tags: ["interface-type", "multiphase", "surface-tension", "interface-tracking", "moving-mesh"]
weight: 10
---

## Overview

`capillaryInterface` couples two incompressible fluid regions that share a
sharp, deforming interface, e.g. a gas bubble (one region) in a liquid (the
other region). It uses ALE interface tracking: each phase is its own mesh, and
the mesh boundary follows the interface. The interface type holds the interface
physics: a finite-area mesh (`faMesh`) on side A, a constant surface tension
$\sigma$, the surface velocity $\mathbf{u}_s$ and the curvature $\kappa$ taken
from the `faMesh`. The jump conditions are applied by four dedicated partitioned
p-U boundary conditions, which abort if they are not used together with a
`capillaryInterface`. Coupling runs through the special DNA field `pUPimple`.
No monolithic variant exists.

## Interface Conditions

Region $A$ is the first entry of `coupledPatchPair` and carries the `faMesh`.
Region $B$ is the second entry. $\mathbf{n}$ is the unit normal of the side-A
patch (outward from $A$), and $\mathbf{n}_B = -\mathbf{n}$. $\mathbf{u}_s$ is
the surface velocity on the `faMesh` (see [Implementation](#implementation)).
The conditions below are written the way the BCs evaluate them, after
interpolating the neighbour side onto the local patch.

### Jump conditions (Dirichlet sides)

**Velocity, set on side B** (`regionCoupledVelocityValue`):

$$
\mathbf{U}_B = \mathbf{U}_A - \mathbf{n}_B\,(\mathbf{n}_B\cdot\mathbf{u}_s) + \mathbf{n}_B\,\frac{\phi_{\text{mesh},B}}{|\mathbf{S}_f|}
\quad \text{on } \Gamma
$$

The tangential velocity is continuous. The normal component of the side-A
velocity is replaced by the side-B mesh normal velocity, corrected by the
normal part of $\mathbf{u}_s$ (which is the side-A mesh normal velocity, see
`updateUs()`).

**Pressure, set on side A** (`regionCoupledPressureValue`, field `p`):

$$
p_A = p_B + 2(\mu_B - \mu_A)\,\nabla_s\cdot\mathbf{u}_s - \sigma\kappa
+ (\rho_B - \rho_A)\,(\mathbf{x}_f - \mathbf{x}_\text{ref})\cdot\mathbf{g}
\quad \text{on } \Gamma
$$

The viscous term is the normal-stress jump, using
$\partial u_n/\partial n = -\nabla_s\cdot\mathbf{u}_s$. The gravity term
implies that `p` is the pressure without the hydrostatic part
$\rho\,\mathbf{g}\cdot(\mathbf{x}-\mathbf{x}_\text{ref})$ (inferred from the
formula). $\mathbf{x}_\text{ref}$ is `PISO.pRefPoint` of region A.
If the field is named `pKin`, the same jump is applied in kinematic form:
$p_{\text{kin},A} = \big(\rho_B\,p_{\text{kin},B} + \text{jump terms}\big)/\rho_A$.

### Flux conditions (Neumann sides)

**Velocity, set on side A** (`regionCoupledVelocityFlux`, diffusivity $\mu_A$):

$$
\mu_A \frac{\partial \mathbf{U}_A}{\partial n} =
-(\mathbf{I} - \mathbf{n}\mathbf{n})\cdot\big(\mu_B \partial_{n_B}\mathbf{U}_B\big)
+ \mathbf{f}_{\sigma,t}
- \mu_A\,\mathbf{n}\,(\nabla_s\cdot\mathbf{u}_s)
+ (\mu_B - \mu_A)\,(\nabla_s\mathbf{u}_s)\cdot\mathbf{n}
$$

with the tangential surface-tension force

$$
\mathbf{f}_{\sigma,t} = \sigma\,\frac{1}{S_f}\sum_{e} \mathbf{L}_e\,c_e \;-\; \sigma\kappa\,\mathbf{n}
$$

($\mathbf{L}_e$: edge length vectors, $c_e$: `edgeLengthCorrection`). With a
constant $\sigma$ this is only the discretisation difference between the
edge-integrated force and $\sigma\kappa\mathbf{n}$. The tangential shear stress
is continuous. The normal component comes from continuity, not from side B.

**Pressure, set on side B** (`regionCoupledPressureFlux`):

$$
k_B\,\frac{\partial p_B}{\partial n_B} = -\,\mathbf{n}_B\cdot\frac{\partial \mathbf{U}_B}{\partial t}
$$

The neighbour flux is added in `fluxJump()` and subtracted again in the generic
flux BC, so it cancels. The pressure gradient on side B depends only on the
local normal acceleration. $k_B$ is the BC's diffusivity (see
[Known Limitations](#known-limitations): with `diffusivity "v"` it actually
evaluates to $\rho_B$, not $1/\rho_B$; it is $1$ for `pKin`).

| Symbol | Description | Unit |
|--------|-------------|------|
| $\mathbf{U}$ | Velocity | m/s |
| $p$ | Pressure (static, without hydrostatic part) | Pa |
| $\mathbf{u}_s$ | Surface velocity on the `faMesh` (`<patchA>Us`) | m/s |
| $\sigma$ | Surface tension (`sigma`, constant) | N/m |
| $\kappa$ | Face curvature, `faMesh::faceCurvatures()` | 1/m |
| $\mu_A,\mu_B$ | Dynamic viscosities from each region's `transportProperties` | Pa·s |
| $\rho_A,\rho_B$ | Densities from each region's `transportProperties` | kg/m³ |
| $\mathbf{g}$ | Gravity (`constant/<regionA>/g`) | m/s² |
| $\phi_\text{mesh}$ | Mesh flux `fvc::meshPhi` on the patch | m³/s |

## Implementation

1. **Construction**: reads `capillaryInterfaceCoeffs/sigma` (required) and
   creates the area field `sigma<regionA><patchA>` (uniform $\sigma$,
   `zeroGradient`, not written). The `faMesh` of region A is looked up by the
   name `<patchA>FaMesh`. If none is registered, a new `faMesh(meshA)` is built
   from region A's finite-area definition and given that name. The
   `interfaceTrackingFvMesh` motion solver gives its own `faMesh` the same name
   so it can be shared (whether the lookup finds it in the current OpenFOAM
   versions has not been verified).
2. **`update()`**: called by **every** partitioned p-U BC in `updateCoeffs()`.
   This is the base-class `regionInterfaceType::update()`:
   - after a topology change (not pure motion), clear `Us`/`Phis`, re-fetch
     the `faMesh` and rebuild the global patches;
   - `updateK()`: on the edges of the `curvatureCorrectedSurfacePatches`,
     replace $\kappa$ by an average of neighbouring faces (wedge axis/fixed
     contact lines);
   - `correct()`: `updateUs()` sets
     $\mathbf{u}_s = \mathbf{U}_A - (\mathbf{U}_A\cdot\mathbf{n})\mathbf{n} + w_A\,\mathbf{n}$
     with $w_A = \phi_{\text{mesh},A}/|\mathbf{S}_f|$, then fixes the `Us`
     boundary values. `updatePhis()` sets the edge flux
     $\phi_s = (\mathbf{u}_s)_e\cdot\mathbf{L}_e$. $\sigma$ is not updated
     (`TODO` for a surfactant equation of state).
3. **DNA loop for `pUPimple`** (`multiRegionSystem::solve()`): every coupling
   iteration runs `solvePIMPLE()` on all regions (the BCs above are evaluated
   inside it), then `meshMotionCorrector()` (interface-tracking mesh motion),
   then `interfaces_->update()`. After the loop, `postSolvePIMPLE()` runs. The
   standalone `solvePIMPLE()` call is skipped when `pUPimple` is coupled.
   Convergence is checked from the jump/flux residuals of the fields listed
   under `residualControl` (`U`, `p`), which must exist in both regions.
4. **`interfaces_->update()`** (`regionInterfaceTypeList`; once in `preSolve()`
   and after every DNA iteration): rebuilds the interpolator if the mesh
   changes (subject to `interpolatorUpdateFrequency`), then calls `info()`.
   `info()` logs min/max/average curvature, the maximum capillary Courant
   number, the total surface-tension force and the total viscous + pressure
   force on side A.
5. **`getMinDeltaT()`**: capillary time-step limit (Brackbill-type), applied
   only if region A has a registered `volScalarField rho`:

   $$
   \Delta t_\sigma = \texttt{maxCapillaryCo}\,\sqrt{\frac{\rho_{A,\min}\,\Delta x_{\min}^3}{2\pi\,\sigma}}
   $$

   where $\Delta x_{\min}$ is the minimum `faMesh` `lPN` (distance between
   neighbouring face centres) and `maxCapillaryCo` is read from
   `system/controlDict` (default `1.0`). The solver uses it only with
   `adjustTimeStep yes` (`setDeltaT.H`). Only $\rho_A$ enters, not
   $\rho_A + \rho_B$.

### Coefficients and dictionary keywords

| Keyword | Location | Default | Description |
|---------|----------|---------|-------------|
| `sigma` | `capillaryInterfaceCoeffs` | **required** | Surface tension, `dimensionedScalar` |
| `coupled` | interface entry | `true` | Switch coupling on/off |
| `couplingStartTime` | interface entry | `0` | Enable coupling after this time |
| `curvatureCorrectedSurfacePatches` | top level of `regionInterfaceProperties` | **required** | `faPatch` names for curvature correction (may be `0()`) |
| `interpolatorUpdateFrequency` | top level | `1` | Interpolator rebuild interval in time steps (`0` = never rebuild) |
| `interfaceTransferMethod` | top level | `AMI` | `AMI` or `directMap`; requires a `<method>Coeffs` sub-dict |
| `maxCapillaryCo` | `system/controlDict` | `1.0` | Capillary Courant limit for `getMinDeltaT()` |

There is no keyword to pick a curvature method: $\kappa$ is always
`faMesh::faceCurvatures()`. The standalone utility `calcCurvature` writes the
same curvature (`curvature`) plus the normal component of the edge-integrated
force (`curvatureNew`) as volume fields for post-processing. It looks for a
patch named `freeSurface` or `interface`.

## Coupled Fields

`coupledFields (pUPimple);` couples `U` and `p` together:

| Field | Required boundary conditions (side A / side B) | Description |
|-------|-----------------------------------------------|-------------|
| `U` | `regionCoupledVelocityFlux` / `regionCoupledVelocityValue` | Velocity: shear-stress balance on A, continuity on B |
| `p` (or `pKin`) | `regionCoupledPressureValue` / `regionCoupledPressureFlux` | Pressure: normal-stress (Laplace) jump on A, acceleration-based gradient on B |

The side assignment is fixed by the code: the BCs on side A read $\mathbf{u}_s$,
$\kappa$ and the `faMesh` normals directly, and `regionCoupledVelocityValue`
interpolates $\mathbf{u}_s$ from its neighbour. The `faMesh` side must therefore
carry `U` flux / `p` value. All four BCs need `interfaceType (capillaryInterface);`.

Fields and objects required:

| Object | Region | Used by |
|--------|--------|---------|
| `U`, `p`, `phi` | A and B | BCs, `updateUs()`, `info()` |
| `transportProperties` with `rho`, `mu` | A and B | BCs, `info()`; also `nu` if a region solves `pKin` (checked: $\mu = \rho\nu$) |
| `g` (`uniformDimensionedVectorField`) | A | `regionCoupledPressureValue` |
| `PISO { pRefPoint ...; }` in `system/<regionA>/fvSolution` | A | `regionCoupledPressureValue` |
| `volScalarField rho` | A (optional) | enables `getMinDeltaT()`; also used for `meshPhi` |
| finite-area mesh on the interface patch | A | everything |

Created at run time (not written): `sigma<regionA><patchA>` (`areaScalarField`),
`<patchA>Us` (`areaVectorField`), `<patchA>Phis` (`edgeScalarField`).

## Compatible Region Types

| Region type (side A) | Region type (side B) | Notes |
|----------------------|----------------------|-------|
| `interTrackFluid` | `interTrackFluid` | Used by all tutorials, but **disabled** in `src/multiRegionSystem/Make/files` |
| `icoFluid` | `icoFluid` | Follows the same conventions (`p` in Pa, `rho` field, `PISO.pRefPoint`, `closedVolume`, `interfaceShadow`, `space`). Inferred, untested. No `g` is created. |
| `pimpleFluid` | `pimpleFluid` | BCs contain a `pKin` branch. `regionCoupledPressureValue` still needs a `PISO { pRefPoint; }` entry, which `pimpleFluid` itself does not read. Untested. |

The region needs a moving mesh whose interface patch follows the interface
(`interfaceTrackingFvMesh` from `libinterfaceTrackingFvMesh`, see
`src/dynamicFvMesh`). Otherwise the interface stays fixed. See
[icoFluid]({{< relref "icoFluid" >}}) and [pimpleFluid]({{< relref "pimpleFluid" >}}).

## Case Setup

### `regionInterfaceProperties` entry

From `tutorials/multiphaseFlow/2dRisingBubble/constant/regionInterfaceProperties`:

```
partitionedCoupledPatches
(
    fluidAfluidB
    {
        interfaceType (capillaryInterface);

        coupledPatchPair
        (
            (fluidA interface)
            (fluidB interfaceShadow)
        );

        coupledFields
        (
            pUPimple
        );

        capillaryInterfaceCoeffs
        {
            sigma       sigma [ 1 0 -2 0 0 0 0 ] 0.07232;
        }
    }
);

monolithicCoupledPatches
();

curvatureCorrectedSurfacePatches 0();

interpolatorUpdateFrequency 0;

interfaceTransferMethod directMap;

directMapCoeffs
{}
```

### `multiRegionProperties` DNA entry

```
DNA
{
    pUPimple
    {
        maxCoupleIter 10;

        residualControl
        {
            U
            {
                maxJumpRes 1e-5;
                outputJumpResField yes;
                maxFluxRes 1e-4;
                outputFluxResField yes;
            }

            p
            {
                maxJumpRes 1e-5;
                outputJumpResField yes;
                maxFluxRes 1e-5;
                outputFluxResField yes;
            }
        }
    }
}
```

### Boundary conditions

Side A (`0/fluidA/orig/partitioned/U` and `p`, patch `interface`, carries the `faMesh`):

```
interface
{
    type                    regionCoupledVelocityFlux;
    diffusivity             "mu";
    interfaceType           (capillaryInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform (0 0 0);
}

interface
{
    type                    regionCoupledPressureValue;
    diffusivity             "v";
    interfaceType           (capillaryInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform 0;
}
```

Side B (`0/fluidB/orig/partitioned/U` and `p`, patch `interfaceShadow`):

```
interfaceShadow
{
    type                    regionCoupledVelocityValue;
    diffusivity             "mu";
    interfaceType           (capillaryInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform (0 0 0);
}

interfaceShadow
{
    type                    regionCoupledPressureFlux;
    diffusivity             "v";
    interfaceType           (capillaryInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform 0;
}
```

### Other required files (from `2dRisingBubble`)

```
constant/fluidA/faMesh/faMeshDefinition  – faMesh on patch "interface" (Allrun: makeFaMesh -region fluidA)
constant/fluid{A,B}/transportProperties  – rho, v, mu
constant/fluidA/g                         – gravity
constant/fluidA/dynamicMeshDict           – dynamicFvMesh interfaceTrackingFvMesh; movingSurfacePatches (interface ...)
system/fluidA/fvSolution                  – PISO { pRefPoint; pRefValue; }  hasSpacePatch yes;
system/fluidB/fvSolution                  – closedVolume yes;
system/faSchemes, system/fluidA/faSchemes, faSolution
system/controlDict                        – libs ("libinterfaceTrackingFvMesh.so" ...)
```

## Known Limitations

- **No runnable tutorial**: all three tutorials select `interTrackFluid`,
  which is commented out in `src/multiRegionSystem/Make/files` (together with
  `pUCoupledIcoFluid`). They would fail when the region type is selected at run time (inferred).
- **Constant surface tension**: `sigma_` is created once from `sigma` and never
  updated (`correct()` has only a `TODO` for an equation of state). There is no
  Marangoni effect or surfactant transport.
- **No curvature model choice**: $\kappa$ = `faMesh::faceCurvatures()`. The only
  correction is the averaging on `curvatureCorrectedSurfacePatches`, and that
  keyword is mandatory even when it is empty.
- **`diffusivity "v"` gives $\rho$, not $1/\rho$** (code bug in
  `diffCoeff()` of both `genericRegionCoupledJump` and `genericRegionCoupledFlux`):
  `"v"` is renamed to `"rho"` before the `== "v"` test, so the reciprocal is
  never taken. For `regionCoupledPressureFlux` this makes
  $\partial p_B/\partial n = -\mathbf{n}_B\cdot\partial_t\mathbf{U}_B/\rho_B$
  instead of $-\rho_B\,\mathbf{n}_B\cdot\partial_t\mathbf{U}_B$. The effect on
  the tutorials was inferred, not tested.
- **Pressure flux ignores the neighbour**: the neighbour flux cancels in
  `regionCoupledPressureFlux`, so side B gets a purely local Neumann condition
  (convective and viscous normal terms are commented out).
- **Fixed side roles**: the `faMesh` must be on the first patch of
  `coupledPatchPair`, and that side must carry the `U` flux / `p` value BCs (see
  above). The hard-coded field names are `U`, `p`, `phi`, `rho`, `g`,
  `transportProperties`.
- **`update()` runs once per BC evaluation**: every coupled BC calls
  `rgInterface().update()`, so $\mathbf{u}_s$ and the curvature correction are
  recomputed up to four times per corrector (cost only).
- **Capillary time step**: uses only $\rho_A$ and $\min(\Delta x)$ of the
  `faMesh`. It is skipped when region A has no `rho` volume field. The
  `info()` capillary Courant number uses the `transportProperties` value instead.
- **`correctUsBoundaryConditions()`**: for `faPatch`es next to a `slip` or
  `symmetry` `U` patch it assigns the full internal `Us` field to the patch
  value (and computes an unused normal `N`). This looks like a bug, but none of
  the tutorials reach it, so it is untested.
- **`Phis` and `normal/tangentialSurfaceTensionForce()`** are computed or
  available but not used by any BC (`regionCoupledVelocityFlux` recomputes
  the forces itself).
- **Serial only (inferred)**: `2dRisingBubble/Allrun -P` prints "Running this
  case in parallel is currently not supported".
- **No monolithic coupling**: `monolithicCoupledPatches` must stay empty for
  this interface.

## Tutorials

No runnable tutorial exists: all of the cases below use `interTrackFluid` in
both regions, and that region type is disabled. They remain the reference
setups (`./Allrun -p -S`).

| Tutorial | Description |
|----------|-------------|
| `tutorials/multiphaseFlow/2dRisingBubble` | 2D gas bubble rising in liquid, $\sigma = 0.07232$ N/m, `directMap` transfer |
| `tutorials/multiphaseFlow/3dOscillatingBubble` | 3D oscillating bubble, $\sigma = 0.0727$ N/m |
| `tutorials/multiphaseFlow/3dDuineveldRisingBubble` | 3D rising bubble (Duineveld experiment), $\sigma = 0.0727$ N/m |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
