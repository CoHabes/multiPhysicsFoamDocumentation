---
title: "3D Duineveld Rising Bubble"
linkTitle: "3D Duineveld Rising Bubble"
weight: 2
mathjax: true
description: >
  3D rising bubble validation against Duineveld (1995) experimental data,
  using the ALE moving-mesh interface tracking method.
---

**Repository path:** `tutorials/multiphaseFlow/3dDuineveldRisingBubble`

**Reference:** Alkafri et al. (2024), *multiRegionFoam: A Unified Multiphysics Framework*, Sect. 7.3

## Overview

A single air bubble rises in still pure water (FluidA), following the
experimental setup of Duineveld (1995).  The moving-mesh ALE interface
tracking method (Hirt et al., originally implemented in OpenFOAM by Tuković
and Jasak) is used to track the deforming bubble surface.

The bubble starts as a sphere of radius $r_b$ at rest and accelerates to its
terminal rise velocity.  The computational domain consists of two regions:

- **FluidB** (bubble / air) — polyhedral mesh
- **FluidA** (surrounding water) — prismatic mesh with polyhedral base

The outer domain is a sphere of radius $20 r_b$.  The bubble mesh is bounded
by the `interfaceShadow` patch on the bubble side and the `interface` patch on
the liquid side where the coupled boundary conditions are applied.

**Bubble radii simulated:** $r_b = 0.5,\, 0.6,\, 0.7,\, 0.8,\, 0.9\,\text{mm}$

Time step: $\Delta t = 10^{-5}\,\text{s}$

## Computational domain

{{< figure src="/images/tutorials/multiphase-flow/3d-duineveld-rising-bubble/fig13-domain.png"
    caption="Fig. 13 — Sketch of the computational domain. FluidB (bubble) is enclosed in FluidA (water) with radius 20 r_b. The interfaceShadow patch on the bubble side coincides with the interface patch on the liquid side." >}}

## Mesh

{{< figure src="/images/tutorials/multiphase-flow/3d-duineveld-rising-bubble/fig14-mesh.png"
    caption="Fig. 14 — Mesh for the rising bubble simulation. Polyhedral cells for the bubble (FluidB) and prismatic cells with a polyhedral base for the water (FluidA)." >}}

## Boundary conditions

**Table: Boundary conditions for the rising bubble simulation**

| Boundary | Pressure | Velocity |
|---|---|---|
| **FluidA** | | |
| interface | `regionCoupledPressureValue` | `regionCoupledVelocityFlux` |
| space | zeroGradient | inletOutlet |
| **FluidB** | | |
| interfaceShadow | `regionCoupledPressureFlux` | `regionCoupledVelocityValue` |

## Physical properties

**Table: Physical properties of air (FluidB) and water (FluidA)**

| Property | Unit | FluidB (air) | FluidA (water) |
|---|---|---|---|
| Density | kg/m³ | 1.205 | 998.3 |
| Dynamic viscosity | kg/ms | 1.82 × 10⁻⁵ | 10⁻³ |
| Surface tension coefficient | N/m | 0.0727 | – |

## Numerical schemes

**Table: Numerical schemes**

| Scheme | Setting |
|---|---|
| `ddtScheme` | `backward` |
| `gradScheme` | `Gauss linear` |
| `divScheme div(phi,U)` | `Gauss GammaVDC 0.5` |
| `divScheme div(phi,T)` | `Gauss linearUpwind Gauss linear` |
| `laplacianScheme` | `Gauss linear corrected` |
| `interpolationScheme` | `linear` |
| `snGradScheme` | `corrected` |
| **Finite area schemes** | |
| `gradScheme` | `Gauss linear` |
| `divScheme` | `Gauss linear` |
| `interpolationScheme` | `linear` |

## Results summary

{{< figure src="/images/tutorials/multiphase-flow/3d-duineveld-rising-bubble/fig15-rise-velocity.png"
    caption="Fig. 15 — Rise velocity over time for a bubble with r_b = 0.5 mm converging to the experimental terminal velocity." >}}

{{< figure src="/images/tutorials/multiphase-flow/3d-duineveld-rising-bubble/fig16-terminal-velocity.png"
    caption="Fig. 16 — Terminal rise velocity for varying equivalent bubble radii. Comparison of simulation results with Duineveld (1995) experiments and Tuković & Jasak simulations." >}}

{{< figure src="/images/tutorials/multiphase-flow/3d-duineveld-rising-bubble/fig17-velocity-field.png"
    caption="Fig. 17 — Velocity field visualisation for the rising bubble with r_b = 0.9 mm at terminal state." >}}

The terminal rise velocity for $r_b = 0.5\,\text{mm}$ converges to the
experimental value from Duineveld.  Across all radii the simulated terminal
velocities show a high level of agreement with both the Duineveld experimental
data and the earlier OpenFOAM simulations by Tuković and Jasak, confirming that
the ALE moving-mesh implementation handles significant mesh deformation
robustly.

## Running the case

```bash
cd $FOAM_RUN/../multiPhysicsFoam/tutorials/multiphaseFlow/3dDuineveldRisingBubble
./Allrun
```
