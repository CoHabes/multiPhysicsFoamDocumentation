---
mathjax: true
title: "conductPotentialTemperature"
description: "Electrically and thermally conducting solid (e.g. battery current-collector tab) solving electric potential conduction and heat conduction with Joule heating."
tags: ["region-type", "electrochemistry", "battery", "thermal", "conduction"]
weight: 10
---

## Overview

`conductPotentialTemperature` models a solid, electrically conducting region
such as the copper or aluminium current-collector tab of a lithium-ion battery
cell. It solves a charge-conservation (Laplace) equation for the electric
potential and a transient heat-conduction equation for the temperature, with
ohmic (Joule) heating from the potential gradient as the heat source. It has no
electrochemistry of its own and is intended to be coupled to an
[`NTGK`]({{< relref "NTGK" >}}) battery region, which provides the current through the shared
interface.

The electrode polarity of the region is set by the `electrode` keyword in
`transportProperties`. It decides which potential field is used: `phiPos` for
`positive` and `phiNeg` for `negative`. These field names match the potential
fields of `NTGK`, so the same field name can be coupled across the
battery/tab interface. A single region type therefore covers both tabs.

## Governing Equations

### Electric potential

$$
\nabla \cdot \left( \sigma \, \nabla \phi \right) = 0
$$

where $\phi$ is `phiPos` or `phiNeg`, depending on `electrode`. The equation has
no transient term, so the potential is quasi-steady within each time step.

### Energy

$$
\rho c_p \frac{\partial T}{\partial t} - \nabla \cdot \left( k \, \nabla T \right) = S_T
$$

### Joule heating

$$
S_T = \dot{Q}_\text{ohm} = \sigma \, \nabla\phi \cdot \nabla\phi
$$

$S_T$ is evaluated explicitly with `fvc::grad` and stored in the field `ST`.

| Symbol | Field / keyword | Description | Unit |
|--------|-----------------|-------------|------|
| $\phi$ | `phiPos` / `phiNeg` | Electric potential | V |
| $T$ | `T` | Temperature | K |
| $\sigma$ | `sigma` | Electric conductivity | S/m |
| $\rho$ | `rho` | Density | kg/m³ |
| $c_p$ | `cp` | Specific heat capacity | J/(kg K) |
| $k$ | `k` | Thermal conductivity | W/(m K) |
| $S_T$ | `ST` | Volumetric heat source (Joule heating) | W/m³ |

All material properties are uniform `dimensionedScalar`s.

## Numerical Algorithm

`conductPotentialTemperature` has no region-local solve. `usesPIMPLE()` returns
`false`, and `solveRegion()`, `prePredictor()`, `momentumPredictor()`,
`pressureCorrector()` and `meshMotionCorrector()` are all no-ops. Each time step
goes through the coupled-matrix path of `multiRegionSystem`:

1. `setCoupledEqns(fieldName)` assembles one equation and registers it in
   `fvScalarMatrices` under the key
   `<field><region>MeshconductPotentialTemperatureTypeEqn`:
   - if `fieldName` is the region's potential field (`phiPos` or `phiNeg`):
     `fvm::laplacian(sigma, phi)`, using the laplacian scheme name
     `laplacian(sigma,fai)`;
   - if `fieldName` is `T`: `fvm::ddt(rho*cp, T) - fvm::laplacian(k, T)`
     (scheme name `laplacian(k,T)`), with `ST` as an explicit source;
   - any other field name, including the potential field of the opposite
     polarity, is ignored.
2. The equations are solved by the DNA (partitioned) coupling loop for each
   field listed in an interface's `coupledFields` (iteration controls come from
   the `DNA` block in `multiRegionProperties`).
3. `postSolve()` calls `correct()`, which recomputes `ST` from the updated
   potential. The new `ST` is first used in the `T` equation of the
   **next** time step.

`getMinDeltaT()` returns `GREAT`, so this region never limits the time step.

`selfCoupledFieldNames()` is not overridden. It returns an empty list, so the
region's equations are assembled and solved only when their field is listed
in the `coupledFields` of a partitioned interface in `regionInterfaceProperties`.

## Sub-Models

None. The region does not read or apply `fvOptions`, and there are no
turbulence, transport-model or MRF sub-models.

## Required Fields

| Field | Type | Description | Unit |
|-------|------|-------------|------|
| `phiPos` or `phiNeg` | `volScalarField` | Electric potential (`phiPos` if `electrode positive;`, `phiNeg` if `electrode negative;`) | V |
| `T` | `volScalarField` | Temperature | K |

Optional / created at runtime:

| Field | Type | Description | Unit |
|-------|------|-------------|------|
| `ST` | `volScalarField` | Joule heat source. Read if present, otherwise initialised to 0; written at write time | W/m³ |

## Case Setup

### Required dictionaries

```
constant/<region>/transportProperties   – electrode, sigma, rho, cp, k
system/<region>/fvSolution               – solvers for phiPos/phiNeg and T
system/<region>/fvSchemes                – ddt, grad (Joule heating), laplacian schemes
```

### `transportProperties`

Every keyword is required; there are no defaults. Any `electrode` value other
than `positive` or `negative` is a fatal error.

```
electrode       negative;                               // positive -> phiPos, negative -> phiNeg

rho             rho   [1 -3 0 0 0 0 0]   8933;          // kg/m3   (copper)
cp              cp    [0 2 -2 -1 0 0 0]  385;           // J/kg/K
k               k     [1 1 -3 -1 0 0 0]  398;           // W/m/K
sigma           sigma [-1 -3 3 0 0 2 0]  59.6e6;        // S/m
```

### `multiRegionProperties` entry

```
regions
(
    ( batteryRegion (NTGK) )
    ( negTabRegion  (conductPotentialTemperature) )
    ( posTabRegion  (conductPotentialTemperature) )
);

DNA
{
    phiPos { maxCoupleIter 50; residualControl { maxJumpRes 1e-05; maxFluxRes 1e-05; } }
    phiNeg { maxCoupleIter 25; residualControl { maxJumpRes 1e-05; maxFluxRes 1e-05; } }
    T      { maxCoupleIter 10; residualControl { maxJumpRes 1e-05; maxFluxRes 1e-05; } }
}
```

### Key `fvSolution` settings

```
solvers
{
    phiNeg      // or phiPos, depending on electrode
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-8;
        relTol          0;
        minIter         2;
        maxIter         250;
    }

    T
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-8;
        relTol          0;
        minIter         2;
        maxIter         250;
    }
}
```

### Interface boundary conditions

On the tab side of the battery/tab interface, the tutorial uses
`regionCoupledScalarFlux` (or `regionCoupledScalarJump`) with `diffusivity`
set to the name of the matching `transportProperties` entry: `sigma` for the
potential and `k` for `T`.

```
toCellFromNeg
{
    type            regionCoupledScalarFlux;
    interfaceType   (defaultInterface);
    diffusivity     "sigma";        // "k" for T
    secondOrder     false;
    nonOrthCorr     false;
    value           $internalField;
}
```

The external terminal is set by an ordinary boundary condition on the potential:
`fixedValue 0` on the negative tab and `fixedGradient` (applied current) on the
positive tab.

## Coupling

| Interface type | Partner region type | Exchanged fields | Coupling mode |
|----------------|--------------------|------------------|---------------|
| `defaultInterface` | `NTGK` | `phiPos` (positive tab) or `phiNeg` (negative tab), `T` | partitioned (DNA) |

No interface type contains code specific to this region type. Coupling goes
through the generic `regionCoupledScalarJump` / `regionCoupledScalarFlux`
patch fields, which read `sigma` / `k` from the region's `transportProperties`.
An interface's `coupledFields` may list both `phiPos` and `phiNeg`, because the
region ignores the potential field of the opposite polarity.

## Known Limitations

- **Static mesh only**: `meshMotionCorrector()` is a no-op.
- **No `fvOptions`**: additional sources cannot be added to the potential or
  temperature equation.
- **Uniform, constant properties**: `sigma`, `rho`, `cp` and `k` are
  `dimensionedScalar`s; temperature-dependent or spatially varying properties
  are not supported.
- **Quasi-steady potential**: the potential equation has no transient term. The
  member `C_` (a specific capacitance meant for a transient term) is declared
  and initialised to 1 but never used.
- **Explicit, lagged Joule heating**: `ST` is updated only in `postSolve()`, so
  the temperature equation uses the source from the previous time step. `ST` is
  not computed at construction, so the first time step uses the initial `ST`
  (0 unless an `ST` field is supplied).
- **Not self-coupled**: with an empty `selfCoupledFieldNames()`, the region
  cannot be run standalone (e.g. in singleRegionMode). Its fields are only
  solved through DNA coupling.
- **Monolithic coupling not supported**: the only tutorial runs partitioned
  coupling only (its `Allrun -m` reports monolithic coupling as unsupported).

## Tutorials

| Tutorial path | Description |
|---------------|-------------|
| `tutorials/electrochemistry/NTGKBatteryWithTabs` | Prismatic NTGK battery cell with copper negative tab and aluminium positive tab (`conductPotentialTemperature`), coupled via `defaultInterface` (`./Allrun -p -S`) |

## References

- Kwon, K. H. et al., *A two-dimensional modeling of a lithium-polymer battery*,
  Journal of Power Sources 163(1), 151–157 (2006).
  DOI: [10.1016/j.jpowsour.2005.12.051](https://doi.org/10.1016/j.jpowsour.2005.12.051)

- Kim, G.-H., Pesaran, A., Spotnitz, R., *A three-dimensional thermal abuse model
  for lithium-ion cells*, Journal of Power Sources 170(2), 476–489 (2007).
  DOI: [10.1016/j.jpowsour.2007.04.018](https://doi.org/10.1016/j.jpowsour.2007.04.018)

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
