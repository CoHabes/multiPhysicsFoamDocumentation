---
mathjax: true
title: "defaultInterface"
description: "Generic pass-through interface without interface physics: couples a scalar field across two region patches through the generic regionCoupledScalarJump / regionCoupledScalarFlux boundary conditions."
tags: ["interface-type", "generic", "partitioned-coupling"]
weight: 10
---

## Overview

`defaultInterface` is the generic interface type. It adds no interface physics
of its own: `correct()` is a no-op, `getMinDeltaT()` returns `GREAT` and
`info()` prints nothing. All of its behaviour comes from the base class
`regionInterfaceType`, which every interface type inherits: it pairs two
patches of two regions, builds the patch-to-patch mapping and provides the
residuals for the DNA (partitioned) coupling loop. The coupling conditions
themselves are set by the boundary conditions on the two patches.

With the generic `regionCoupledScalarJump` / `regionCoupledScalarFlux` boundary
conditions, `defaultInterface` couples a scalar quantity that is continuous
(or continuous up to a constant partition coefficient) and whose diffusive flux
is conserved across the interface. Examples are temperature and electric
potential across a perfectly conducting contact, as between an
[`NTGK`]({{< relref "NTGK" >}}) battery region and its
[`conductPotentialTemperature`]({{< relref "conductPotentialTemperature" >}})
current-collector tabs.

This page also documents the base-class behaviour (interface naming,
dictionary keywords, interface-to-interface mapping, generic coupled boundary
conditions). It applies to every interface type unless the derived type
overrides it.

## Interface Conditions

The conditions are set by the boundary-condition pair, not by
`defaultInterface`. One side of the interface carries a jump (Dirichlet) BC,
side $J$, and the other side carries a flux (Neumann) BC, side $F$. Which
region gets which BC does not depend on the patch order in
`coupledPatchPair`. The tutorial uses different assignments for different
fields.

### Jump condition

`regionCoupledScalarJump` sets the patch value on side $J$ from the neighbour
value, interpolated onto the side-$J$ faces:

$$
\phi_J = c \, \phi_F \quad \text{on } \Gamma
$$

Internally this is $\phi_J = \phi_F + \Delta\phi$ with
`valueJump()` $\Delta\phi = (c - 1)\,\phi_F$. If the `jumpCoefficient` keyword
is not given, $c = 1$, which gives continuity.

### Flux condition

`regionCoupledScalarFlux` sets the normal gradient on side $F$ so that the
diffusive flux matches the flux computed on side $J$:

$$
\Gamma_F \nabla \phi_F \cdot \mathbf{n}_F + \Gamma_J \nabla \phi_J \cdot \mathbf{n}_J = 0 \quad \text{on } \Gamma
$$

In code this is `gradient() = (-flux_J + fluxJump())/k_F`, where
`flux_J = k_J * snGrad(phi_J)` is interpolated onto the side-$F$ faces.
`regionCoupledScalarFlux::fluxJump()` returns zero.

| Symbol | Keyword / field | Description | Unit |
|--------|-----------------|-------------|------|
| $\phi$ | field listed in `coupledFields` (e.g. `T`, `phiPos`, `phiNeg`) | Coupled scalar field (same name in both regions) | field-dependent |
| $\Gamma_J$, $\Gamma_F$ | `diffusivity` (BC keyword) | Diffusivity on each side | field-dependent |
| $c$ | `jumpCoefficient` (jump-BC keyword) | Value-jump (partition) coefficient, default 1 | — |
| $\mathbf{n}_J$, $\mathbf{n}_F$ | — | Outward patch normals of each side | — |

### Interface residuals

The DNA loop (`src/numerics/dnaSolutionControl`) measures convergence with
residuals computed by the BCs and collected by
`regionInterfaceType::normJumpRes` / `normFluxRes`:

$$
r_\text{jump} = \frac{\lVert \phi_J - \phi_F - \Delta\phi \rVert_2}
{\max\left(\min\left(\lVert \phi_J \rVert_2,\ \lVert \phi_F + \Delta\phi \rVert_2\right),\ \text{SMALL}\right)}
$$

$$
r_\text{flux} = \frac{\lVert \Gamma_F \partial_n \phi_F + q_J \rVert_2}
{\max\left(\min\left(\lVert \Gamma_F \partial_n \phi_F \rVert_2,\ \lVert q_J \rVert_2\right),\ \text{SMALL}\right)}
$$

with $q_J = \Gamma_J \partial_n \phi_J$ interpolated onto side $F$, and with
the neighbour quantity always interpolated onto the evaluating side. These are
compared with `maxJumpRes` / `maxFluxRes` in `multiRegionProperties`. If the
interface is not coupled (see `coupled` below), both residuals are 0.

## Implementation

`defaultInterface` only implements the pure virtuals of `regionInterfaceType`
(`correct()`, `getMinDeltaT()`, `info()`) as no-ops. Everything below is
base-class behaviour.

1. **Construction** (`regionInterfaceTypeList::reset`): for each entry of
   `partitionedCoupledPatches` (and `monolithicCoupledPatches`) in
   `constant/regionInterfaceProperties`, the two patches of `coupledPatchPair`
   are looked up. One interface object is created **per entry of the
   `interfaceType` list**, all on the same patch pair. The first pair entry
   becomes side A and the second becomes side B.
2. **Naming**: the interface name is
   `<regionA><PatchA><RegionB><PatchB><InterfaceType>`, where every part after
   the first is capitalised. For example, the tutorial gives
   `batteryRegionToNegTabNegTabRegionToCellFromNegDefaultInterface`. This name
   appears in the log and in the DNA residual output.
3. **interfaceKey**: the interface is also an `interfaceKey`, an unordered pair
   of words. Equality is symmetric (`(a,b) == (b,a)`) and the hash is the sum
   of the two word hashes. The interface object itself is keyed by the two
   **patch** names. The coupled-field tables of `regionInterfaceTypeList`
   are keyed by `(<regionA><patchA>, <regionB><patchB>)`.
4. **Global patches**: on construction, a `globalPolyPatch` is built for each
   side. It is a copy of the patch reconstructed in full on every processor,
   which is what makes the mapping work in parallel. The mapping object is
   created at the same time (`interfaceToInterface()`), so that a `directMap`
   can read its maps on restart.
5. **Field transfer**: the BCs call `interpolateFacesFromA/B` (and the point
   variants). These convert the patch field to the global patch, transfer it
   with the selected `interfaceToInterfaceMapping`, and scatter it back to the
   local patch. Only `scalar` and `vector` fields have specialisations. Other
   types hit a `notImplemented` fallback.
6. **Update per BC evaluation**: every coupled BC calls
   `regionInterfaceType::update()` in its `updateCoeffs()`. This function
   rebuilds the global patches after a topology change (mesh changing but not
   moving), updates the finite-area curvature if an `faMesh` exists, then calls
   `correct()`, which is a no-op for `defaultInterface`.
7. **Update per time step**: `multiRegionSystem::preSolve()` calls
   `regionInterfaceTypeList::update()`. For interfaces whose mesh A is changing,
   this rebuilds the mapping and global patches at the frequency set by
   `interpolatorUpdateFrequency` (and always after a topology change), then
   calls `info()`.
8. **Time step**: `getMinDeltaT()` returns `GREAT`, so `defaultInterface` never
   limits the time step.

### Interface-to-interface mapping

The mapping is selected **globally** for all interfaces by the top-level
keyword `interfaceTransferMethod` in `regionInterfaceProperties`. It is not set
per interface entry. The default is `AMI`. The matching `<method>Coeffs`
sub-dictionary is read with `subDict`, so it **must exist even when it is
empty**. Neither `AMI` nor `directMap` reads any keyword from it.

| Method | Class (`src/numerics/interfaceToInterfaceMappings`) | Behaviour |
|--------|------------------------------------------------------|-----------|
| `AMI` (default) | `amiInterfaceToInterfaceMapping` | Wraps `amiZoneInterpolation` (face-area weighted, `requireMatch` on, `lowWeightCorrection -1`) on the global patches; non-conformal patches allowed. Prints the max face-centre mapping error (`interface-to-interface face error`) on construction. |
| `directMap` | `directMapInterfaceToInterfaceMapping` | One-to-one face/point map for **conformal** patches; fatal error if face or point counts differ. Faces/points are matched by centre/position within `0.001 * min edge length`. Maps (`zoneAToZoneBFaceMap_<patchA>_to_<patchB>`, …) are written to and re-read from the time directory. |
| `RBF` | `rbfInterfaceToInterfaceMapping` | Radial basis function interpolation. **Disabled**: commented out in `src/multiRegionSystem/Make/files` (needs Eigen), so not selectable. |

### Generic coupled boundary conditions

`genericRegionCoupledJump<Type>` (derived from `fixedValue`) and
`genericRegionCoupledFlux<Type>` (derived from `fixedGradient`) are the base
classes of all partitioned coupling BCs. Both are registered for `scalar` and
`vector`, but their `valueJump()` / `fluxJump()` are `notImplemented`, so they
cannot be used directly in a case. `regionCoupledScalarJump` /
`regionCoupledScalarFlux` are the concrete scalar versions and the only
coupled BCs that do not need a specific interface type. All other
`regionCoupled*` BCs (`regionCoupledTemperatureJump`,
`regionCoupledVelocityValue`, `regionCoupledTraction`, …) cast the interface to
`heatTransferInterface`, `capillaryInterface`, `fsiInterface`, etc.

Behaviour shared by both generic BCs:

- The owning interface is found from the BC's `interfaceType` list (only the
  **first** entry is used) and the local region/patch name. It is a fatal
  error if no interface or more than one interface of that type contains the
  patch. The neighbour region and patch are taken from the matched interface.
- The neighbour field is looked up under the **same name** as the local field.
- If the interface is not `coupled`, `updateCoeffs()` returns early and the
  BC keeps its last value or gradient.
- After the new value (jump side) or gradient (flux side) is set, it is
  relaxed by the BC's acceleration model (`src/numerics/acceleration`).
- `diffCoeff()` takes the diffusivity named by `diffusivity` from a
  `volScalarField` of that name in the local region if one exists. Otherwise it
  reads a `dimensionedScalar` of that name from the region's
  `transportProperties`. For a field named `pKin` the diffusivity is forced
  to 1.
- Snapping the gradient (`snGrad`, jump side) and evaluating the patch value
  (flux side) use the cell gradient (`grad(<field>)` if registered, otherwise
  `fvc::grad`) for the optional non-orthogonal (`nonOrthCorr`) and
  second-order (`secondOrder`) corrections.

| Keyword | Side | Default | Description |
|---------|------|---------|-------------|
| `type` | both | — | `regionCoupledScalarJump` or `regionCoupledScalarFlux` |
| `interfaceType` | both | **required** | List of interface type names, e.g. `(defaultInterface)`; only the first is used |
| `diffusivity` | both | `word::null` | Name of the diffusivity (`volScalarField` or `transportProperties` entry). Effectively required: the flux side always needs it, and the jump side needs it to supply its flux to the flux side |
| `jumpCoefficient` | jump | `word::null` (→ $c = 1$) | Name of a `volScalarField` or `transportProperties` entry giving $c$ |
| `neighbourFieldName` | both | local field name | Read and written, but not used by the generic jump/flux evaluation, which always uses the local field name |
| `accType` | both | `fixed` | Acceleration model: `fixed`, `aitken`, `IQNILS` |
| `relax` | both | `1.0` | (Initial) relaxation factor |
| `couplingReuse` | both | `0` | `IQNILS` only: number of previous time steps reused |
| `nonOrthCorr` | both | `false` | Non-orthogonal correction |
| `secondOrder` | both | `false` | Second-order patch gradient / value |
| `value` | both | read if present | Initial patch value |
| `gradient` | flux | read if present | Initial patch gradient |

## Coupled Fields

| Field | Required boundary conditions (side A / side B) | Description |
|-------|-----------------------------------------------|-------------|
| any `volScalarField` present in both regions under the same name | `regionCoupledScalarJump` / `regionCoupledScalarFlux` (either order) | Generic continuity + flux conservation, optional partition coefficient |
| `T` (tutorial) | `regionCoupledScalarJump` (battery) / `regionCoupledScalarFlux` (tab) | Temperature, `diffusivity "k"` |
| `phiNeg` (tutorial) | `regionCoupledScalarJump` (battery) / `regionCoupledScalarFlux` (negative tab) | Negative-electrode potential, `diffusivity "sigmaNeg"` / `"sigma"` |
| `phiPos` (tutorial) | `regionCoupledScalarFlux` (battery) / `regionCoupledScalarJump` (positive tab) | Positive-electrode potential, `diffusivity "sigmaPos"` / `"sigma"` |

The flux side casts the neighbour patch field to
`genericRegionCoupledJump`, so the BCs must be one jump and one flux. Two flux
BCs fail with a bad cast. A field listed in `coupledFields` that is missing in
one of the regions is skipped by the DNA residual check (`couplesField()`).
That is why the tutorial can list `phiPos` and `phiNeg` on both tab
interfaces, although each tab region holds only one of them.

## Compatible Region Types

| Region type (side A) | Region type (side B) | Notes |
|----------------------|----------------------|-------|
| `NTGK` | `conductPotentialTemperature` | Tested in tutorial: `phiPos`/`phiNeg` and `T` |
| any region that assembles the field in `setCoupledEqns()` | any such region | The field must exist under the same name on both sides and the `diffusivity` must be resolvable in each region. Not covered by a tutorial |

[`laplacianFluid`]({{< relref "laplacianFluid" >}}) lists `defaultInterface` for
`T` in its coupling table, but no tutorial uses that combination (untested).

## Case Setup

### `regionInterfaceProperties` entry

From `tutorials/electrochemistry/NTGKBatteryWithTabs/constant/regionInterfaceProperties`:

```
partitionedCoupledPatches
(
    cellToNegTab
    {
        interfaceType (defaultInterface);

        coupledPatchPair
        (
            (batteryRegion toNegTab)
            (negTabRegion  toCellFromNeg)
        );

        coupledFields
        (
            phiPos
            phiNeg
            T
        );

        defaultInterfaceCoeffs
        {}
    }

    cellToPosTab
    {
        // same, with (batteryRegion toPosTab) (posTabRegion toCellFromPos)
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

Interface-entry keywords:

| Keyword | Default | Description |
|---------|---------|-------------|
| `interfaceType` | **required** | List of interface types; one interface object per entry is built on the same patch pair |
| `coupledPatchPair` | **required** | `((regionA patchA) (regionB patchB))`; exactly two entries, otherwise a fatal error |
| `coupledFields` | **required** | Field names coupled through this patch pair; their union across all entries defines the partitioned-coupled fields |
| `coupled` | `true` | Switch the interface coupling on or off |
| `couplingStartTime` | `0` | If `coupled false` and this value is > 0, coupling is switched on once `time >= couplingStartTime` |
| `<interfaceType>Coeffs` | — | Not read by `defaultInterface`. Kept for consistency with other interface types |

Top-level keywords (shared by all interfaces):

| Keyword | Default | Description |
|---------|---------|-------------|
| `partitionedCoupledPatches` | **required** | List of partitioned interface entries |
| `monolithicCoupledPatches` | **required** | List of monolithic interface entries (may be empty: `()`) |
| `curvatureCorrectedSurfacePatches` | **required** | `faPatch` names for curvature correction at wedge axes (used only when an `faMesh` exists); `0()` otherwise |
| `interpolatorUpdateFrequency` | `1` | Rebuild mapping/global patches every N time steps on changing meshes; `0` disables periodic rebuild |
| `interfaceTransferMethod` | `AMI` | `AMI` or `directMap` (see above) |
| `<method>Coeffs` | **required** | Sub-dictionary for the chosen method; may be empty |

`regionInterfaceProperties` is not read in `singleRegionMode`.

### `multiRegionProperties` DNA entry

From the same tutorial:

```
regions
(
    (batteryRegion (NTGK) )
    (negTabRegion  (conductPotentialTemperature) )
    (posTabRegion  (conductPotentialTemperature) )
);

DNA
{
    phiPos
    {
        maxCoupleIter 50;

        residualControl
        {
            maxJumpRes 1e-05;
            outputJumpResField no;
            maxFluxRes 1e-05;
            outputFluxResField no;
        }
    }

    phiNeg { maxCoupleIter 25; residualControl { maxJumpRes 1e-05; maxFluxRes 1e-05; } }
    T      { maxCoupleIter 10; residualControl { maxJumpRes 1e-05; maxFluxRes 1e-05; } }
}
```

`maxCoupleIter` defaults to 1. Either both `maxJumpRes` and `maxFluxRes`, or a
single `maxIntRes`, may be given. `defaultInterface` does not override
`normIntRes()`, which returns 0, so use `maxJumpRes` / `maxFluxRes`.

### Boundary conditions

Jump side, `0/batteryRegion/orig/partitioned/T`:

```
toNegTab
{
    type                    regionCoupledScalarJump;
    interfaceType           (defaultInterface);
    diffusivity             "k";
    accType                 aitken;
    relax                   0.75;
    secondOrder             false;
    nonOrthCorr             false;
    value                   $internalField;
}
```

Flux side, `0/negTabRegion/orig/partitioned/T`:

```
toCellFromNeg
{
    type                    regionCoupledScalarFlux;
    interfaceType           (defaultInterface);
    diffusivity             "k";
    secondOrder             false;
    nonOrthCorr             false;
    value                   $internalField;
}
```

Flux side with an initial gradient, `0/batteryRegion/orig/partitioned/phiPos`:

```
toPosTab
{
    type                    regionCoupledScalarFlux;
    interfaceType           (defaultInterface);
    diffusivity             "sigmaPos";
    secondOrder             false;
    nonOrthCorr             false;
    value                   $internalField;
    gradient                uniform 0.1266;
}
```

On a patch where a field is not coupled (e.g. `phiNeg` on the battery's
`toPosTab`), use an ordinary BC such as `zeroGradient`.

## Known Limitations

- **No interface physics**: `correct()` and `info()` are empty and the
  `defaultInterfaceCoeffs` dictionary is never read. Any jump or source at the
  interface must come from the BCs. With the generic BCs, this is limited to a
  constant or field-valued partition coefficient (`jumpCoefficient`) and a zero
  flux jump.
- **Scalar fields only in practice**: field transfer is specialised only for
  `scalar` and `vector`. The only concrete BCs usable with `defaultInterface`
  are `regionCoupledScalarJump` / `regionCoupledScalarFlux`. The generic
  `vector` BCs have `notImplemented` jumps, and all vector `regionCoupled*` BCs
  require `capillaryInterface`, `fsiInterface` or a mass-transfer interface.
  Coupling a vector field via `defaultInterface` would need a new BC.
- **Same field name on both sides**: the generic BCs always look up the
  neighbour field under the local field name. `neighbourFieldName` is ignored
  there.
- **Global mapping choice**: `interfaceTransferMethod` and its `Coeffs`
  sub-dictionary are top-level keywords. All interfaces in a case share one
  mapping method. `directMap` requires conformal patches. `RBF` is not
  compiled.
- **Only the first BC `interfaceType` is used**: if several interface types
  are listed on one patch, the BC binds to the first. A patch may belong to at
  most one interface of a given type.
- **Neighbour detection by patch name only**
  (`interfaceToInterfaceCoupleManager::nbrPatch()/nbrMesh()/interpolateFromNbrField()`):
  the local patch is compared with `patchA` by **name only**, without the
  region. If both sides of an interface use the same patch name, side B would
  identify itself as A and couple to itself (inferred from the code,
  untested). Give the two patches different names.
- **`diffusivity "v"` handling is ineffective** (inferred from the code):
  `diffCoeff()` renames `"v"` to `"rho"` before it checks for `"v"` to take
  the reciprocal, so the reciprocal is never applied.
- **Monolithic coupling not available**: `monolithicCoupledPatches` is parsed,
  but the monolithic solve in `multiRegionSystem::solve()` is commented out.
- **Residual control requires generic BCs**: if neither patch field derives
  from `genericRegionCoupledJump` / `genericRegionCoupledFlux`, a warning is
  printed and the residual is set to `GREAT`. The DNA loop then always runs to
  `maxCoupleIter` (inferred).

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/electrochemistry/NTGKBatteryWithTabs` | `NTGK` battery cell coupled to a negative and a positive `conductPotentialTemperature` tab via two `defaultInterface` entries (`phiPos`, `phiNeg`, `T`; `directMap`; `./Allrun -p -S`) |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
