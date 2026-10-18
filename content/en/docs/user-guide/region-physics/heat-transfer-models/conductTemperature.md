---
mathjax: true
title: "conductTemperature"
description: "Transient heat conduction in a solid region, solved only through partitioned (DNA) coupling of the temperature field T."
tags: ["region-type", "solid", "heat-transfer", "conduction", "static-mesh"]
weight: 10
---

## Overview

`conductTemperature` models transient heat conduction in a solid. It solves a
single energy equation for the temperature `T`, using density `rho`, specific
heat `cv` and thermal conductivity `kappa`. There is no velocity field, no
advection and no pressure-velocity coupling. It is meant to be the solid side
of a conjugate heat transfer problem, coupled to a fluid region (usually
`pimpleFluid transportTemperature`) through a `heatTransferInterface`. It can
also be stacked with a solid mechanics region type (e.g. `solidStVK`) for
thermo-FSI.

## Governing Equations

### Energy (heat conduction)

$$
\frac{\partial \left( \rho \, c_v \, T \right)}{\partial t}
= \nabla \cdot \left( \kappa \, \nabla T \right)
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $T$ | Temperature | K |
| $\rho$ | Density | kg/m³ |
| $c_v$ | Specific heat capacity | J/(kg K) = m²/(s² K) |
| $\kappa$ | Thermal conductivity | W/(m K) = kg m/(s³ K) |

The equation is assembled as

```cpp
fvm::ddt(rho_*cv_, T()) == fvm::laplacian(kappa_(), T(), "laplacian(kappa,T)")
```

so `fvSchemes` must provide a `laplacian(kappa,T)` entry (or a `default`
laplacian scheme). The equation has no source term.

## Numerical Algorithm

`conductTemperature` does not run its own solve loop. `usesPIMPLE()` returns
`false`. `solveRegion()`, `prePredictor()`, `momentumPredictor()`,
`pressureCorrector()`, `meshMotionCorrector()` and `postSolve()` are all
no-ops. The temperature is solved only when `multiRegionSystem` asks for it:

1. The region does **not** override `selfCoupledFieldNames()`, so the base
   default (an empty list) applies. `T` is therefore solved **only** when it
   appears as a partitioned-coupled field, i.e. when the case has a
   `DNA { T { ... } }` entry in `multiRegionProperties` and `T` is listed in
   `coupledFields` of an interface in `regionInterfaceProperties`.
2. Inside each DNA coupling iteration,
   `multiRegionSystem::assembleAndSolveEqns<fvMatrix, scalar>("T")` calls
   `setCoupledEqns("T")`. This assembles the equation above and stores it
   under the key `T<regionName>MeshconductTemperatureTypeEqn`. The equation is
   then relaxed (`relaxationFactors` in `fvSolution`) and solved. The region
   Picard corrector (`maxCorr` / `relTol` / `tolerance` per field in
   `fvSolution`) can re-assemble and re-solve it within one DNA iteration.
   `setCoupledEqns` ignores its `fieldName` argument and always assembles the
   `T` equation.
3. `correct()` only calls `kappa.correctBoundaryConditions()`.
4. `getMinDeltaT()` returns `GREAT`, so this region never limits the global
   time step.

## Sub-Models

None. The region does not read `fvOptions`, a thermophysical model or
turbulence. All material properties come from `transportProperties` or from
optional field files.

## Required Fields

| Field | Type | Read mode | Description | Unit |
|-------|------|-----------|-------------|------|
| `T` | `volScalarField` | `MUST_READ`, `AUTO_WRITE` | Temperature | K |
| `kappa` | `volScalarField` | `READ_IF_PRESENT`, `AUTO_WRITE` | Thermal conductivity. If there is no file, a uniform field is created from `transportProperties` with `calculated` patches. | W/(m K) |
| `rho` | `volScalarField` | `READ_IF_PRESENT`, `NO_WRITE` | Density. If there is no file, it is uniform from `transportProperties`. | kg/m³ |
| `cv` | `volScalarField` | `READ_IF_PRESENT`, `NO_WRITE` | Specific heat. If there is no file, it is uniform from `transportProperties`. | m²/(s² K) |

`T` and `kappa` are created through `lookupOrRead`. If a sibling region type
on the same mesh has already registered a field with that name, that field is
reused. `rho` and `cv` are always constructed by this region.

## Case Setup

### Required dictionaries

```
constant/<region>/transportProperties   – rho, cv, kappa (all required, dimensioned)
system/<region>/fvSolution               – T solver, relaxationFactors
system/<region>/fvSchemes                – ddt and laplacian(kappa,T) schemes
0/<region>/T                             – temperature with coupled BC on the interface patch
```

### `multiRegionProperties` entry

```
regions
(
    ( fluid (pimpleFluid transportTemperature) )
    ( solid (conductTemperature) )
);

DNA
{
    T
    {
        maxCoupleIter 20;

        residualControl
        {
            maxJumpRes 1e-07;
            outputJumpResField no;
            maxFluxRes 1e-07;
            outputFluxResField no;
        }
    }
}
```

### `transportProperties`

All three entries are looked up as `dimensionedScalar` and are required, even
when a field file overrides them:

```
rho             rho   [1 -3 0 0 0 0 0]  1.0;
cv              cv    [0 2 -2 -1 0 0 0] 100.0;
kappa           kappa [1 1 -3 -1 0 0 0] 100.0;
```

### Key `fvSolution` settings

```
solvers
{
    T
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-15;
        relTol          0;
        minIter         2;
        maxIter         250;
    }
}

relaxationFactors
{
    T               1.0;
}
```

### Interface boundary condition on `T`

On the solid side of a partitioned `heatTransferInterface` (flux side; the
fluid side then uses `regionCoupledTemperatureJump`):

```
top
{
    type            regionCoupledHeatFlux;
    interfaceType   (heatTransferInterface);
    secondOrder     false;
    nonOrthCorr     false;
    value           $internalField;
}
```

The coupled patch fields look up the diffusivity field `kappa` by default (set
the patch keyword `diffusivity` to use another name). This matches the field
that `conductTemperature` registers.

## Coupling

| Interface type | Exchanged fields | Coupling mode |
|----------------|-----------------|---------------|
| `heatTransferInterface` | `T` (temperature jump / heat flux via `regionCoupledTemperatureJump` / `regionCoupledHeatFlux`, using `kappa`) | partitioned (DNA) |

`heatTransferInterface` only provides the `jumpCoefficient` (from
`heatTransferInterfaceCoeffs`). The patch fields do the actual transfer.
Monolithic coupling is commented out in `multiRegionSystem::solve()`, so
only partitioned coupling is available at the moment.

## Known Limitations

- **Coupled use only.** `selfCoupledFieldNames()` is not overridden, so `T`
  is never solved unless it is DNA-coupled. A standalone (single-region)
  `conductTemperature` case does not advance `T`. Use `laplacianFluid` for
  standalone conduction.
- **No source terms.** `fvOptions` are not read or applied.
- **Static mesh only.** `meshMotionCorrector()` is a no-op and the equation
  contains no mesh-flux terms.
- **No time-step control.** `getMinDeltaT()` returns `GREAT`; `deltaT` must be
  chosen by the user or limited by other regions.
- **Constant `rho` and `cv` in time.** They are read once and never updated
  (`NO_WRITE`). Only `kappa` gets a boundary refresh in `correct()`, and its
  values are never recomputed (for example as a function of `T`).
- **No non-orthogonal corrector loop.** The coupled solve path has no
  `correctNonOrthogonal()` loop, so `nNonOrthogonalCorrectors` in
  `fvSolution` is not used by this region.
- **No output/post-processing.** `postSolve()` is a no-op. No derived fields
  (e.g. heat flux) are written.

## Tutorials

| Tutorial path | Description |
|---------------|-------------|
| `tutorials/conjugateHeatTransfer/flowOverHeatedPlate` | Flow over a heated plate (`pimpleFluid transportTemperature` + `conductTemperature`) |
| `tutorials/conjugateHeatTransfer/ShellAndTubeHeatExchanger` | Shell-and-tube heat exchanger, two fluid regions + `conductTemperature` solid |
| `tutorials/fluidStructureInteraction/heatedBeamInCrossFlow` | Heated beam in cross-flow, thermo-FSI (`solidStVK conductTemperature` solid) |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
