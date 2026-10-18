---
mathjax: true
title: "massTransferInterface"
description: "Phase-change (evaporation/condensation) mass transfer across a sharp liquid-vapour interface, driven by the heat-flux imbalance at saturation temperature."
tags: ["interface-type", "phase-change", "mass-transfer", "multiphase"]
weight: 10
---

## Overview

`massTransferInterface` models phase change (evaporation/condensation) at a
sharp interface between two fluid regions, e.g. a liquid and its vapour. It is
**not** a species/concentration interface: there is no Henry or partition
coefficient and no species field. The interface is assumed to sit at a fixed
saturation temperature $T_{sat}$. The interfacial mass flux $\dot m$ follows
from the energy balance: the net conductive heat that reaches the interface
from both sides at $T = T_{sat}$, divided by the latent heat $h_{lv}$.

The interface type itself only computes $\dot m$ and the matching heat-flux
jump. These quantities enter the solution through dedicated partitioned
boundary conditions (`regionCoupledPhaseChangeHeatFlux` for `T`,
`regionCoupledMassTransfer*` for `U` and `p`). The `U`/`p` conditions also need
a `capillaryInterface` on the same patch pair (surface tension, curvature,
interface velocity). The two interface types are therefore used together.

> **Status:** No tutorial uses this interface type, and the interface-tracking
> fluid region it was written for (`interTrackFluid`) is disabled in
> `src/multiRegionSystem/Make/files`. In the current code `correct()` is never
> called, so $\dot m$ stays zero (see [Known Limitations](#known-limitations)).
> Treat this page as a description of the code, not of a validated feature.

## Interface Conditions

Below, $A$ and $B$ are the two regions of the `coupledPatchPair` and
$\mathbf{n}_A$, $\mathbf{n}_B$ the outward patch normals. In the BC formulas
"own" is the side the BC sits on and "nbr" is the other side. $\mathbf{n}_\varepsilon$
is the interface normal $\mathbf{n}$ on the `meshA` side of the
`capillaryInterface` and $-\mathbf{n}$ on the other side.

### Mass flux from the energy balance (`massTransferInterface::saturatedFlux`)

Both temperature fields are copied, the interface patch values are set to
$T_{sat}$ (`fixedValue`), and the normal gradients are evaluated:

$$
q_{sat} = k_A \left.\frac{\partial T_A}{\partial n_A}\right|_{T_\Gamma = T_{sat}}
        + k_B \left.\frac{\partial T_B}{\partial n_B}\right|_{T_\Gamma = T_{sat}}
$$

$$
\dot m = -\frac{q_{sat}}{h_{lv}} \quad \text{(set to } 0 \text{ at time index 1)}
$$

$-q_{sat}$ is the net heat conducted into the interface from both sides. It is
consumed as latent heat, so $\dot m > 0$ when the interface gains heat
(evaporation).

### Temperature flux condition (`regionCoupledPhaseChangeHeatFlux`)

$$
k_{own} \frac{\partial T_{own}}{\partial n_{own}}
= -\,k_{nbr} \frac{\partial T_{nbr}}{\partial n_{nbr}} + q_{sat}
$$

The other side of `T` must carry a jump-type BC. That BC provides the value
continuity $T_{own} = T_{nbr}$; this interface type adds no temperature jump.

### Velocity jump (`regionCoupledMassTransferVelocityValue`)

$$
\mathbf{u}_{own} = \mathbf{u}_{nbr}
+ \left(\frac{1}{\rho_{own}} - \frac{1}{\rho_{nbr}}\right)\dot m\, \mathbf{n}_\varepsilon
$$

### Pressure jump (`regionCoupledMassTransferPressureValue`)

$$
p_{own} = p_{nbr}
+ 2\,(\mu_{nbr} - \mu_{own})\,\nabla_s\!\cdot\mathbf{u}_s
- \sigma \kappa
+ (\rho_{nbr} - \rho_{own})\,\mathbf{g}\cdot(\mathbf{x}_f - \mathbf{x}_{ref})
- \left(\frac{1}{\rho_{own}} - \frac{1}{\rho_{nbr}}\right)\dot m^2
$$

If the pressure field is called `pKin` (kinematic pressure), the same jump is
built from $\rho_{nbr}\,p_{kin,nbr}$ and the result is divided by
$\rho_{own}$.

### Momentum flux jump (`regionCoupledMassTransferVelocityFlux`)

The diffusive flux $\mu_{own}\,\partial\mathbf{u}_{own}/\partial n$ is set to the
negated normal component of the neighbour flux plus

$$
\mathbf{f}_{\sigma,t}
- \mu_{own}\,\mathbf{n}\,\nabla_s\!\cdot\mathbf{u}_s
+ (\mu_{nbr} - \mu_{own})\,(\nabla_s \mathbf{u}_s)\cdot\mathbf{n}
+ \dot m \left[\left(\frac{\dot m}{\rho_{own}}\mathbf{n}_\varepsilon - \mathbf{u}_{own}\right)
            - \left(\frac{\dot m}{\rho_{nbr}}\mathbf{n}_\varepsilon - \mathbf{u}_{nbr}\right)\right]
$$

where $\mathbf{f}_{\sigma,t}$ is the tangential part of the surface-tension
force from the `capillaryInterface` finite-area mesh.

| Symbol | Description | Unit |
|--------|-------------|------|
| $T$ | Temperature (field `T` in both regions) | K |
| $T_{sat}$ | Saturation temperature (`TSat`) | K |
| $h_{lv}$ | Latent heat of vaporisation (`hlv`) | J/kg |
| $k_A$, $k_B$ | Thermal conductivity, entry `k` in each region's `transportProperties` | W/(m K) |
| $q_{sat}$ | Heat-flux imbalance at saturation (`<patchA>fluxMTJump`) | W/m² |
| $\dot m$ | Interfacial mass flux (`<patchA>mDotS`) | kg/(m² s) |
| $\rho$, $\mu$ | Density and dynamic viscosity from `transportProperties` | kg/m³, Pa s |
| $\sigma$, $\kappa$, $\mathbf{u}_s$ | Surface tension, curvature, interface velocity from `capillaryInterface` | N/m, 1/m, m/s |
| $\mathbf{g}$, $\mathbf{x}_{ref}$ | Gravity (`g`) and `PISO.pRefPoint` | m/s², m |

## Implementation

1. **Construction.** Reads `massTransferInterfaceCoeffs` from the interface
   entry. Creates on the interface finite-area mesh (`aMesh()`, built from
   region A if no `<patchA>FaMesh` is registered) the area fields
   `hlvInterface` (uniform $h_{lv}$) and `mDotInterface` (zero), both
   `NO_WRITE`.
2. **`correct()`** calls `updateMDotS()`:
   - `saturatedFlux()` looks up `T` in both regions and the scalar `k` from
     each region's `transportProperties` dictionary. It then computes
     $q_{sat}$ as above on copies of `T` (the solved fields are not changed).
   - The result is stored in the area field `<patchA>fluxMTJump`. The area
     field `<patchA>mDotS` is set to $-q_{sat}/h_{lv}$, or to zero when
     `timeIndex() == 1`. Both fields are created on demand and are
     `NO_WRITE`.
   - `mDotInterface` is set to `mDotS`. A diagnostic line (time, integrated
     $\dot m$, region volumes) is printed to the log.
3. **Boundary conditions** read the results: `regionCoupledPhaseChangeHeatFlux`
   uses `fluxMTJump()` as the flux jump of `T`. The
   `regionCoupledMassTransfer*` conditions interpolate `mDotS()` to their
   patch.
4. **`getMinDeltaT()`** returns `GREAT` (no time-step restriction).
   **`info()`** is empty.

## Coupled Fields

| Field | Required boundary conditions (side A / side B) | Description |
|-------|-----------------------------------------------|-------------|
| `T` | `regionCoupledPhaseChangeHeatFlux` with `interfaceType (massTransferInterface ...)` on one side / a `genericRegionCoupledJump`-derived BC on the other (*inferred*: `regionCoupledScalarJump`; `regionCoupledTemperatureJump` needs a `heatTransferInterface` listed first) | Temperature, flux jump $q_{sat}$ |
| `U` (via `pUPimple`) | `regionCoupledMassTransferVelocityValue` / `regionCoupledMassTransferVelocityFlux` | Velocity with mass-flux jump, `interfaceType (capillaryInterface massTransferInterface)` |
| `p` (via `pUPimple`) | `regionCoupledMassTransferPressureValue` / `regionCoupledPressureFlux` (*inferred*: there is no mass-transfer pressure-flux BC) | Pressure with recoil term $\dot m^2$ |

The BC's `interfaceType` list order matters. The first entry is the default
interface of the BC. `regionCoupledPhaseChangeHeatFlux` aborts unless its
first entry is `massTransferInterface`. `regionCoupledMassTransferPressureValue`
and `regionCoupledMassTransferVelocityFlux` abort unless their first entry is
`capillaryInterface`. Each BC finds the second interface type by name on the
same patch.

Fields and dictionaries needed in **both** regions:

| Name | Type | Used by |
|------|------|---------|
| `T` | `volScalarField` | `saturatedFlux()` |
| `transportProperties` | registered `IOdictionary` with `k` (and `rho`, `mu` for the `U`/`p` BCs, `nu` if `pKin` is used) | interface, BCs |
| `U`, `p` (or `pKin`), `phi` | vol/surface fields | `regionCoupledMassTransfer*` BCs |
| `g` | `uniformDimensionedVectorField` | `regionCoupledMassTransferPressureValue` |
| `PISO.pRefPoint` | entry in `fvSolution` | `regionCoupledMassTransferPressureValue` |

## Compatible Region Types

No compiled region type combination has been shown to work. The table lists
what the code requires (*inferred, untested*).

| Region type (side A) | Region type (side B) | Notes |
|----------------------|----------------------|-------|
| `interTrackFluid` | `interTrackFluid` | Intended pairing (moving-mesh interface tracking, as in the `capillaryInterface` tutorials). **Disabled** in `src/multiRegionSystem/Make/files`. |
| `icoFluid` + `transportTemperature` | `icoFluid` + `transportTemperature` | `icoFluid` provides `U`, `p`, `phi`, `transportProperties` (`rho`, `mu`) and documents these BCs (see `icoFluid`). `transportTemperature` provides `T`. Both read `kappa`, so `k` must be added to `transportProperties` as well. `g` must be provided separately. |
| any + `conductTemperature` | — | Only the `T` part can work (no `U`/`p`); `k` must be added next to `kappa`. |

The interface needs `k` as a single uniform `dimensionedScalar` per region. A
non-uniform `kappa` field (see `transportTemperature`, `conductTemperature`) is
not used.

## Case Setup

No tutorial exists. The snippets below are assembled from the keywords the
code reads and from the structure of
`tutorials/multiphaseFlow/2dRisingBubble` (which uses `capillaryInterface`
only). They have not been run.

### `regionInterfaceProperties` entry

```
partitionedCoupledPatches
(
    fluidAfluidB
    {
        interfaceType (capillaryInterface massTransferInterface);

        coupledPatchPair
        (
            (fluidA interface)
            (fluidB interfaceShadow)
        );

        coupledFields
        (
            pUPimple
            T
        );

        capillaryInterfaceCoeffs
        {
            sigma       sigma [ 1 0 -2 0 0 0 0 ] 0.07232;
        }

        massTransferInterfaceCoeffs
        {
            hlv         hlv  [ 0 2 -2 0 0 0 0 ] 2.257e6;   // required, J/kg
            TSat        TSat [ 0 0 0 1 0 0 0 ] 373.15;    // required, K
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

| Keyword | Where | Default | Description |
|---------|-------|---------|-------------|
| `hlv` | `massTransferInterfaceCoeffs` | — (required) | Latent heat $h_{lv}$ |
| `TSat` | `massTransferInterfaceCoeffs` | — (required) | Saturation temperature $T_{sat}$ |
| `coupled` | interface entry | `true` | Base class switch |
| `couplingStartTime` | interface entry | `0` | Base class |
| `interfaceTransferMethod` | top level | `AMI` | Patch-to-patch mapping used by the BCs. See also the conformity note in Known Limitations. |
| `curvatureCorrectedSurfacePatches` | top level | — (required) | Base class, needed by every interface type |

All interface types listed in `interfaceType` receive the same entry, so they
share `coupledPatchPair` and `coupledFields`.

### `multiRegionProperties` DNA entry

```
DNA
{
    pUPimple
    {
        maxCoupleIter 10;

        residualControl
        {
            U { maxJumpRes 1e-5; maxFluxRes 1e-4; }
            p { maxJumpRes 1e-5; maxFluxRes 1e-5; }
        }
    }

    T
    {
        maxCoupleIter   10;
        residualControl { maxJumpRes 1e-5; maxFluxRes 1e-5; }
    }
}
```

### Boundary conditions

Based on the `capillaryInterface` BCs in
`tutorials/multiphaseFlow/2dRisingBubble/0/fluidA/orig/partitioned/U` and
`.../fluidB/orig/partitioned/p`, with the mass-transfer variants substituted:

```
// 0/fluidA/U
interface
{
    type                    regionCoupledMassTransferVelocityFlux;
    diffusivity             "mu";
    interfaceType           (capillaryInterface massTransferInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform (0 0 0);
}

// 0/fluidB/U
interfaceShadow
{
    type                    regionCoupledMassTransferVelocityValue;
    diffusivity             "mu";
    interfaceType           (capillaryInterface massTransferInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform (0 0 0);
}

// 0/fluidA/p
interface
{
    type                    regionCoupledMassTransferPressureValue;
    diffusivity             "v";
    interfaceType           (capillaryInterface massTransferInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform 0;
}

// 0/fluidB/T
interfaceShadow
{
    type                    regionCoupledPhaseChangeHeatFlux;
    interfaceType           (massTransferInterface);
    secondOrder             false;
    nonOrthCorr             false;
    value                   uniform 373.15;
}
```

`regionCoupledPhaseChangeHeatFlux` fixes its diffusivity name to `k`. It uses a
`volScalarField k` if the region has one, otherwise `k` from
`transportProperties`.

## Known Limitations

- **`correct()` is never called (bug).** `regionInterfaceTypeList::update()`
  only refreshes the interpolators and calls `info()`. It does not call
  `regionInterfaceType::update()`, which is the only caller of `correct()`.
  Nothing else calls `updateMDotS()`. In the current code $\dot m$ and
  $q_{sat}$ therefore stay at their initial value of zero. *(Static code
  reading, not confirmed by a run.)* This affects the `correct()` of every
  interface type, not only this one.
- **No runnable tutorial.** No tutorial uses `massTransferInterface`. The
  `interTrackFluid` region type it was written for is commented out in
  `src/multiRegionSystem/Make/files`.
- **Constant properties.** $h_{lv}$, $T_{sat}$ and `k` are single uniform
  scalars. $T_{sat}$ does not depend on pressure.
- **Conformal interface assumed.** `saturatedFlux()` adds patch-A and
  patch-B face values index by index without mapping. Results are written into
  an area field on the finite-area mesh of region A. This only makes sense if
  both patches and the `aMesh` have the same faces in the same order
  (*inferred*). The BCs on side A also call `interpolateFromNbrField` on
  `mDotS`, which is already defined on A (*inferred to be harmless only for
  identical conformal patches*).
- **First time step forced to $\dot m = 0$.** The check is
  `timeIndex() == 1`, so the first step after a restart is presumably not
  affected, because the time index continues from `uniform/time` (*inferred*).
- **`regionCoupledMassTransferVelocityValue::updatePhi()` is never called.**
  The mass-flux correction of the boundary `phi` it implements is dead code.
- **No temperature constraint at the interface.** $T_{sat}$ is used only to
  compute $q_{sat}$. `T` on the patches is not fixed to $T_{sat}$, and the
  saturation check (`checkSaturation`) is computed and then discarded.
- **Name clash risk with `k`** (*inferred*). The heat-flux BC prefers a
  `volScalarField k`. In a turbulent region that is the turbulent kinetic
  energy, not the conductivity.
- **Diagnostics are always on.** `updateMDotS()` and the velocity BC write an
  `Info` line on every call. The output fields are `NO_WRITE`, so $\dot m$ is
  not written to disk.
- **Partitioned only.** Monolithic coupling is not available.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| — | No tutorial uses `massTransferInterface`. `tutorials/multiphaseFlow/2dRisingBubble` shows the companion `capillaryInterface` setup (also on the disabled `interTrackFluid`). |

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
