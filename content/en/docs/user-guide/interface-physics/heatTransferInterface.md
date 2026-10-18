---
mathjax: true
title: "heatTransferInterface"
description: "Conjugate heat transfer interface: temperature continuity (with an optional constant jump coefficient) and heat flux continuity, coupled partitioned through the T boundary conditions."
tags: ["interface-type", "heat-transfer", "conjugate-heat-transfer"]
weight: 10
---

## Overview

`heatTransferInterface` couples the temperature field `T` of two regions that
share an interface, usually a fluid (`pimpleFluid transportTemperature`) and a
solid (`conductTemperature`). It is used for conjugate heat transfer (CHT)
and, stacked with `fsiInterface` on the same patch pair, for thermo-FSI. The
interface class itself holds only one parameter, the constant
`jumpCoefficient`. All transfer work happens in the two partitioned boundary
conditions `regionCoupledTemperatureJump` and `regionCoupledHeatFlux`, which
are solved inside the DNA coupling loop of `multiRegionSystem`.

## Interface Conditions

Let side $A$ carry the `regionCoupledTemperatureJump` BC (Dirichlet side) and
side $B$ carry the `regionCoupledHeatFlux` BC (Neumann side). These sides are
set by the BC types in `0/<region>/T`. They do not depend on the order in
`coupledPatchPair`.

### Jump condition

$$
T_A = c \, T_B \quad \text{on } \Gamma
$$

In the code, the jump BC sets its value to
$T_B + \texttt{valueJump}$ with $\texttt{valueJump} = T_B\,(c - 1)$, where
$T_B$ is the neighbour patch value interpolated onto side $A$. With the
default tutorial value $c = 1$ this is plain temperature continuity.

### Flux condition

$$
\kappa_A \nabla T_A \cdot \mathbf{n}_A + \kappa_B \nabla T_B \cdot \mathbf{n}_B = 0 \quad \text{on } \Gamma
$$

The flux BC sets its gradient to

$$
\left.\frac{\partial T}{\partial n}\right|_B = -\frac{\left(\kappa_A \, \partial_n T_A\right)_{A \to B}}{\kappa_B}
$$

where $\left(\kappa_A \, \partial_n T_A\right)_{A \to B}$ is the flux of the jump side
(`flux()` = `snGrad() * kappa`), interpolated onto side $B$. The flux
jump of `regionCoupledHeatFlux` is always zero (`fluxJump()` returns `*this * 0`).

| Symbol | Description | Unit |
|--------|-------------|------|
| $T$ | Temperature (coupled field) | K |
| $c$ | `jumpCoefficient` (dimensionless, uniform, constant) | – |
| $\kappa$ | Diffusivity looked up by the BC (default field name `kappa`) | W/(m K) (see note below) |
| $\mathbf{n}_A, \mathbf{n}_B$ | Outward patch normals of side $A$ / $B$ | – |
| $\Gamma$ | Interface (coupled patch pair) | – |

The flux balance only means energy conservation when $\kappa$ has the same
units on both sides. The code does not check this; see Known Limitations.

## Implementation

`heatTransferInterface` derives from `regionInterfaceType` and is registered in
its `IOdictionary` run-time selection table. At run time:

1. **Construction** (`regionInterfaceType` base): the object is named
   `<regionA><PatchA><RegionB><PatchB>HeatTransferInterface`, for example
   `fluidInterfaceSolidTopHeatTransferInterface`. The base class builds global
   patches for both sides and creates the interface-to-interface mapping
   immediately (`interfaceTransferMethod`, default `AMI`; see Case Setup).
   `heatTransferInterface` then reads
   `heatTransferInterfaceCoeffs/jumpCoefficient`.
2. **`correct()`**: no-op. It is called from `regionInterfaceType::update()`,
   which each coupled BC calls in `updateCoeffs()`. Because
   `heatTransferInterface` creates no `faMesh`, `update()` does no curvature work.
3. **`getMinDeltaT()`**: returns `GREAT`. The interface never limits the time
   step.
4. **`info()`**: empty.
5. **Mapping update**: when the interface reports `changing()`,
   `regionInterfaceTypeList::update()` calls
   `updateInterpolatorAndGlobalPatches()`. The mapping and global patches are
   rebuilt every `interpolatorUpdateFrequency` time steps (`0` disables this).
6. **Boundary conditions** (in each DNA iteration, when `T` is assembled and
   solved by `multiRegionSystem`):
   - `regionCoupledTemperatureJump::updateCoeffs()` looks up the neighbour `T`,
     interpolates its patch values, adds `valueJump()`, applies the
     acceleration model (`accType`, `relax`) and fixes the value.
   - `regionCoupledHeatFlux::updateCoeffs()` asks the neighbour
     jump BC for `flux()`, interpolates it, negates it, divides it by the local
     `kappa`, applies the acceleration model and fixes the gradient.
   - `regionCoupledTemperatureJump` looks up the interface found through its
     own `interfaceType` entry as a `heatTransferInterface` (to read
     `jumpCoefficient`) and stops with a `FatalError` if it is another type.
     `regionCoupledHeatFlux` defines the same check but never calls it.
7. **Residuals**: the DNA loop evaluates the jump residual
   $|T_A - T_B - \texttt{valueJump}|$ and the flux residual
   $|q_B + q_{A\to B} - 0|$ through the generic jump/flux BCs. These are checked
   against `maxJumpRes` / `maxFluxRes`.

### Diffusivity lookup

`regionCoupledTemperatureJump` and `regionCoupledHeatFlux` override
`diffusivityName()`. It returns the patch keyword `diffusivity` if that is
set, otherwise `kappa`. `diffCoeff()` then:

1. uses the patch values of the `volScalarField` with that name, if one is
   registered on the region mesh;
2. otherwise reads a uniform `dimensionedScalar` with that name from the
   region's `transportProperties` (the dictionary must be registered on the mesh).

## Coupled Fields

| Field | Required boundary conditions (side A / side B) | Description |
|-------|-----------------------------------------------|-------------|
| `T` | `regionCoupledTemperatureJump` / `regionCoupledHeatFlux` | Temperature. Value continuity (×`jumpCoefficient`) on the jump side, heat-flux continuity on the flux side. |

Fields needed but not created by the interface:

| Field | Where it must exist | Source |
|-------|---------------------|--------|
| `T` | both regions, same name | region type (e.g. `transportTemperature`, `conductTemperature`) |
| `kappa` (or the name given by `diffusivity`) | both regions, as `volScalarField` or as a `transportProperties` entry | region type or `transportProperties` |

The interface creates no fields. The BC keyword `neighbourFieldName` is read
and written by `interfaceToInterfaceCoupleManager`, but the jump/flux BCs always
look up the neighbour field by the local field name.

## Compatible Region Types

The region types below assemble a `T` equation in `setCoupledEqns` and/or
register a diffusivity the BCs can find. Rows marked *tutorial* are run in
the shipped tutorials. The other rows follow from the code and are untested.

| Region type (side A) | Region type (side B) | Notes |
|----------------------|----------------------|-------|
| `pimpleFluid transportTemperature` | `conductTemperature` | *Tutorial* (`flowOverHeatedPlate`, `ShellAndTubeHeatExchanger`). Both register `kappa` (W/(m K)). Fluid takes the jump BC, solid the flux BC. |
| `pimpleFluid transportTemperature` | `solidStVK conductTemperature` | *Tutorial* (`heatedBeamInCrossFlow`), with `fsiInterface` on the same patch pair. `T`/`kappa` come from `conductTemperature`, not `solidStVK`. |
| `rhoPimpleFluid` | `conductTemperature` | Untested. `rhoPimpleFluid` registers `kappa` from `turbulence->kappaEff()` and assembles a `T` equation in `setCoupledEqns`. |
| `laplacianFluid` | any of the above | Untested. `laplacianFluid` registers `DT` (m²/s), not `kappa`. Set `diffusivity DT;` on its patch, or provide `kappa` in its `transportProperties`. Mixing `DT` with a `kappa` in W/(m K) breaks flux consistency. |
| `conductPotentialTemperature` | any of the above | Untested. The conductivity is a `dimensionedScalar` `k` in `transportProperties` (not a field). Set `diffusivity k;` so the fallback lookup finds it. |

See `docs/regions/conductTemperature.md`, `docs/regions/transportTemperature.md`,
`docs/regions/pimpleFluid.md`, `docs/regions/rhoPimpleFluid.md`,
`docs/regions/solidStVK.md`, `docs/regions/laplacianFluid.md` and
`docs/regions/conductPotentialTemperature.md`.

## Case Setup

### `regionInterfaceProperties` entry

From `tutorials/conjugateHeatTransfer/flowOverHeatedPlate/constant/regionInterfacePropertiesPartitioned`:

```
partitionedCoupledPatches
(
    fluidsolid
    {
        interfaceType (heatTransferInterface);
        coupledPatchPair
        (
            (fluid interface)
            (solid top)
        );
        coupledFields
        (
            T
        );
        heatTransferInterfaceCoeffs
        {
            // Jump coefficient
            jumpCoefficient               jumpCoefficient [ 0 0 0 0 0 0 0 ] 1;
        }
    }
);

monolithicCoupledPatches
(
);

curvatureCorrectedSurfacePatches 0();

interpolatorUpdateFrequency 1;

interfaceTransferMethod directMap;

directMapCoeffs
{}
```

Thermo-FSI stacks both interface types on one patch pair
(`tutorials/fluidStructureInteraction/heatedBeamInCrossFlow/constant/regionInterfaceProperties`):

```
fluidsolid
{
    interfaceType (fsiInterface heatTransferInterface);
    coupledPatchPair
    (
        (fluid interfaceShadow)
        (solid interface)
    );
    coupledFields
    (
        D
        T
    );
    coupled             yes;
    fsiInterfaceCoeffs
    {}
    heatTransferInterfaceCoeffs
    {
        // Jump coefficient
        jumpCoefficient               jumpCoefficient [ 0 0 0 0 0 0 0 ] 1;
    }
}
```

#### Keywords

Per-interface entry (inside `partitionedCoupledPatches`):

| Keyword | Required | Default | Read by | Description |
|---------|----------|---------|---------|-------------|
| `interfaceType` | yes | – | `regionInterfaceTypeList` | List of interface types; one object is created per entry. |
| `coupledPatchPair` | yes | – | `regionInterfaceTypeList` | `((regionA patchA) (regionB patchB))`; must have exactly two entries. |
| `coupledFields` | yes | – | `regionInterfaceTypeList` | Fields coupled across this interface (`T`). |
| `heatTransferInterfaceCoeffs` | yes | – | `heatTransferInterface` | Sub-dictionary; `subDict` fails if it is missing. |
| `heatTransferInterfaceCoeffs/jumpCoefficient` | yes | – (no default in code) | `heatTransferInterface` | Dimensionless `dimensionedScalar` $c$, with $T_A = c\,T_B$. |
| `coupled` | no | `true` | `regionInterfaceType` | Turns the coupling on/off. When off, the coupled BCs skip `updateCoeffs()` and the residuals are zero. |
| `couplingStartTime` | no | `0` | `regionInterfaceType` | If `> 0` and `coupled` is off, coupling is turned on once `time >= couplingStartTime`. |

Top-level entries (shared by all interfaces in the file):

| Keyword | Required | Default | Description |
|---------|----------|---------|-------------|
| `interfaceTransferMethod` | no | `AMI` | Patch-to-patch mapping: `AMI` or `directMap` (`RBF` is not compiled; it needs Eigen). |
| `<interfaceTransferMethod>Coeffs` | yes | – | Must exist (`subDict`), even when empty. Neither `AMI` nor `directMap` reads any entries from it. |
| `interpolatorUpdateFrequency` | no | `1` | Rebuild period (in time steps) of the mapping for changing interfaces; `0` = never. |
| `curvatureCorrectedSurfacePatches` | yes | – | Read by the base class for every interface type; use `0()` for heat transfer. |

### `multiRegionProperties` DNA entry

From `tutorials/conjugateHeatTransfer/flowOverHeatedPlate/constant/multiRegionProperties`:

```
regions
(
    (fluid (pimpleFluid transportTemperature) )
    (solid (conductTemperature) )
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

### Boundary conditions

Jump side, fluid
(`tutorials/conjugateHeatTransfer/flowOverHeatedPlate/0/fluid/orig/partitioned/T`):

```
interface
{
    type                    regionCoupledTemperatureJump;
    interfaceType           (heatTransferInterface);

    accType                 aitken;
    relax                   0.75;

    secondOrder             false;
    nonOrthCorr             false;
    value                   $internalField;
}
```

Flux side, solid
(`tutorials/conjugateHeatTransfer/flowOverHeatedPlate/0/solid/orig/partitioned/T`):

```
top
{
    type                    regionCoupledHeatFlux;
    interfaceType           (heatTransferInterface);

    secondOrder             false;
    nonOrthCorr             false;
    value                   $internalField;
}
```

Patch keywords common to both BCs:

| Keyword | Required | Default | Description |
|---------|----------|---------|-------------|
| `interfaceType` | yes | – | Interface type list. Only the **first** entry is used to find the interface (matched by type plus local region/patch), so `heatTransferInterface` must come first. |
| `diffusivity` | no | `kappa` | Name of the diffusivity field / `transportProperties` entry. |
| `accType` | no | `fixed` (with a warning) | Acceleration model: `fixed`, `aitken`, `IQNILS`. |
| `relax` | no | `1.0` | Initial relaxation factor of the acceleration model. |
| `nonOrthCorr` | no | `false` | Non-orthogonal correction of the coupled `snGrad`. |
| `secondOrder` | no | `false` | Second-order correction of the coupled `snGrad`. |
| `neighbourFieldName` | no | local field name | Read and written, but not used by these two BCs (see above). |
| `value` | – | – | Initial patch value. |

The tutorials put the jump BC on the fluid side and the flux BC on the solid
side. The opposite assignment is not forbidden by the code, but the flux BC
`refCast`s the neighbour patch field to `genericRegionCoupledJumpFvPatchField`,
so the two sides must use one jump BC and one flux BC.

## Known Limitations

- **Partitioned only.** The monolithic solve in `multiRegionSystem::solve()`
  is commented out ("need to be ported"). The `*Monolithic` variants of
  `flowOverHeatedPlate` and `ShellAndTubeHeatExchanger` (BCs
  `monolithicTemperature`, `monolithicThermalDiffusivitySlave`) are compiled
  but do not drive a monolithic solve at present.
- **Constant, uniform jump coefficient.** `jumpCoefficient` is a single
  `dimensionedScalar`, read once. It cannot vary in space or time. The
  generic jump BC's `jumpCoefficient` patch keyword (field lookup) is not used
  by `regionCoupledTemperatureJump`, which overrides `valueJump()`.
- **No interface physics.** `correct()` and `info()` are empty, and
  `getMinDeltaT()` returns `GREAT`. There is no contact resistance, no
  radiation and no interfacial heat source (`fluxJump()` is always zero).
- **No unit/consistency check on the diffusivity.** The flux balance uses
  whatever field the `diffusivity` name resolves to on each side. For
  example, `laplacianFluid` provides `DT` in m²/s, while `conductTemperature`
  provides `kappa` in W/(m K). Mixing them silently gives a non-conservative
  coupling (inferred from the code).
- **Missing diffusivity is a fatal error.** If neither a `kappa` field nor a
  `kappa` entry in the region's `transportProperties` exists, the lookup in
  `diffCoeff()` stops the run.
- **Same patch names on both sides (inferred, untested).**
  `interfaceToInterfaceCoupleManager::nbrMesh()/nbrPatch()` and
  `interpolateFromNbrField()` decide the side by comparing **patch names only**
  with `patchA`. If both coupled patches have the same name in their regions
  (e.g. `interface`/`interface`), side B would resolve itself as its own
  neighbour. All shipped tutorials use different names.
- **`interfaceTransferMethod` is global.** It and its `Coeffs` sub-dictionary
  are read from the top level of `regionInterfaceProperties`, so all
  interfaces in a case share one mapping method.
- **`neighbourFieldName` has no effect** on `regionCoupledTemperatureJump` /
  `regionCoupledHeatFlux`. The neighbour field is always looked up by the
  local field name.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/conjugateHeatTransfer/flowOverHeatedPlate` | Laminar flow over a heated plate: `fluid (pimpleFluid transportTemperature)` / `solid (conductTemperature)`, `directMap` transfer. Partitioned (`-p`) and monolithic (`-m`) dictionaries are provided. |
| `tutorials/conjugateHeatTransfer/ShellAndTubeHeatExchanger` | Two fluid regions (`inner`, `outer`) coupled to one `solid` through two `heatTransferInterface`s, `AMI` transfer, `nonOrthCorr true`. Partitioned and monolithic dictionaries are provided. |
| `tutorials/fluidStructureInteraction/heatedBeamInCrossFlow` | Thermo-FSI: `heatTransferInterface` stacked with `fsiInterface` on one patch pair (`solidStVK conductTemperature` solid), `directMap` transfer. |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
