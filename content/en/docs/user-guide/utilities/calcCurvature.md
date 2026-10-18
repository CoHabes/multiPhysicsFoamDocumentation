---
mathjax: true
title: "calcCurvature"
description: "Standalone pre-/post-processing application that evaluates the finite-area curvature of a free-surface or interface patch and writes it as volume fields."
tags: ["utility", "pre-processing", "post-processing", "finite-area", "curvature", "multiphase"]
weight: 10
---

## Overview

`calcCurvature` is a standalone application that builds the finite-area mesh
(`faMesh`) of a case and evaluates the mean curvature of that surface for a
single time. It writes two volume fields, `curvature` and `curvatureNew`. Both
are zero everywhere except on the boundary patch named `freeSurface` or
`interface`. The curvature it computes is the same `faMesh::faceCurvatures()`
that `capillaryInterface` uses for the surface-tension pressure jump. Use it to
check the curvature quality of a surface mesh (e.g. a spherical bubble, where
$\kappa$ should be uniform) before running a capillary case. It is a single-mesh
tool: it does not read `multiRegionProperties` and knows nothing about regions
beyond the standard `-region` option.

## Method

1. **Time and mesh.** The time is selected with the standard time options (see
   [Usage](#usage)). The run time is set to the **first** selected time, and only
   that time is processed. The mesh is created as a `dynamicFvMesh` (region
   from `-region`, default `region0`). If `constant/[<region>/]dynamicMeshDict`
   is absent, a `staticFvMesh` is used. The mesh is never moved.
2. **Finite-area mesh.** `faMesh(mesh)` is constructed for the default area
   name, i.e. it reads the `faMesh` written by `makeFaMesh`
   (`constant/[<region>/]faMesh`).
3. **Patch detection.** The boundary is scanned for a patch named
   `freeSurface` or `interface`. If both exist, the one with the higher patch
   index wins. If neither exists, the application aborts with a bare
   `FatalError` (no message).
4. **Curvature** (`curvature`). The `faMesh` face curvature is copied into the
   detected patch:

   $$
   \mathbf{k}_N = \frac{1}{S_f}\sum_{e \in f} \mathbf{L}_e\, c_e,
   \qquad
   \kappa = \operatorname{sign}\!\left(\mathbf{k}_N\cdot\mathbf{n}_f\right)\,|\mathbf{k}_N|
   $$

   with $S_f$ the face area, $\mathbf{L}_e$ the edge length vectors (`Le()`),
   $c_e$ the edge length correction (`edgeLengthCorrection()`,
   $c_e = \cos\!\big(\tfrac{1}{2}\arcsin|\mathbf{n}_a\times\mathbf{n}_b|\big)$ from
   the two edge point normals) and $\mathbf{n}_f$ the `faMesh` face area normal.
   $\kappa$ is therefore signed relative to the `faMesh` normal.
5. **Alternative curvature** (`curvatureNew`). The same edge sum
   $\mathbf{k}_N$ is computed again. A unit "normal" is built from it,
   $\hat{\mathbf{n}} = \mathbf{k}_N/|\mathbf{k}_N|$, falling back to the
   `fvPatch` face normal where $|S_f\,\mathbf{k}_N| \le$ `SMALL`. The written
   value is

   $$
   \kappa_\text{new} = \hat{\mathbf{n}}\cdot\mathbf{k}_N
   $$

   Because $\hat{\mathbf{n}}$ is parallel to $\mathbf{k}_N$, this is
   $|\mathbf{k}_N| = |\kappa|$ wherever the fallback is not active, i.e. the
   **unsigned** curvature.
6. **Report.** Min/max and arithmetic mean (unweighted by face area) of both
   patch fields are printed to the log.

## Usage

Command line (standalone application):

```bash
calcCurvature [-region <name>] [-time <ranges> | -latestTime] [-constant] [-noZero]
```

| Option | Description |
|--------|-------------|
| `-region <name>` | Mesh region to operate on. Default: `region0` (case root). |
| `-time <ranges>` | Time ranges to select (e.g. `0.1`, `:10,20 40:70`). Only the **first** selected time is processed. |
| `-latestTime` | Select the latest time. On its own, this processes the latest time; combined with other selections only the earliest one is processed. |
| `-constant` | Include `constant/` in the time list. |
| `-noZero` | Exclude `0/` from the time list. |

Without time options all time directories (including `0/`) are selected, so
the **earliest** time is processed. If no time directory matches, `constant/`
is used with a warning.

Typical sequence for a single-mesh case:

```bash
blockMesh
makeFaMesh
calcCurvature
```

## Input / Output

**Read**

| File | Description |
|------|-------------|
| `constant/[<region>/]polyMesh` | Volume mesh; must contain a patch named `freeSurface` or `interface` |
| `constant/[<region>/]faMesh` | Finite-area mesh from `makeFaMesh` on that patch |
| `constant/[<region>/]dynamicMeshDict` | Optional; selects the `dynamicFvMesh` type (static if absent) |
| `system/[<region>/]faSchemes`, `faSolution` | Required by the `faMesh` |
| `system/[<region>/]fvSchemes`, `fvSolution` | Required by the `fvMesh` |

**Written** to `<time>/[<region>/]` of the processed time:

| Field | Type | Unit | Content |
|-------|------|------|---------|
| `curvature` | `volScalarField` | 1/m | 0 in all cells and on all other patches; signed $\kappa$ on the free-surface/interface patch |
| `curvatureNew` | `volScalarField` | 1/m | 0 in all cells and on all other patches; $\kappa_\text{new}$ (unsigned, see [Method](#method)) on that patch |

All patches are of type `calculated`. No other fields are read or written.

## Known Limitations

- **Single time only.** The selected time list is not looped over. Only the
  earliest selected time is processed, even if `-time` selects several.
- **Fixed patch names.** The patch must be called exactly `freeSurface` or
  `interface`. There is no option to pick another name (e.g. `interfaceShadow`).
  A missing patch aborts with an empty error message.
- **faMesh must cover exactly that patch.** The `faMesh` face values are
  assigned directly to the patch, so the `faMesh` must be built on that patch
  only, with the same face count and ordering. A `faMesh` on a different or
  additional patch causes a size mismatch.
- **`curvatureNew` is unsigned.** Its "normal" is derived from the edge sum
  itself rather than from the surface normal, so it carries no sign. It is not
  the signed normal force density that its comments suggest.
- **Mean of face values.** The printed average is an unweighted face average,
  not an area-weighted one.
- **Libraries for the `dynamicMeshDict`.** The application links only
  `finiteVolume`, `finiteArea`, `dynamicFvMesh`, `dynamicMesh` and
  `meshTools`. A `dynamicMeshDict` that selects a mesh type or motion solver
  from another library (e.g. `interfaceTrackingFvMesh`, or a motion solver
  from `fvMotionSolvers`) needs that library loaded through the dictionary's
  library entries. Otherwise the selection fails.
- **No parallel-specific handling.** Min/max/average are evaluated per
  processor without reduction.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/multiphaseFlow/3dDuineveldRisingBubble/baseMeshes/meshA` | Base bubble mesh generation (`blockMesh`, `makeFaMesh` on patch `interface`). The `calcCurvature` call is present but **commented out** in `Allrun`. |
| `tutorials/multiphaseFlow/3dDuineveldRisingBubble/hexBaseMeshes/meshA` | Same for the hex base mesh; the call is likewise commented out. |
