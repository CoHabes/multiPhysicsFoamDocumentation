---
mathjax: true
title: "solidPointDisplacement"
description: "Function object that writes the displacement history of the solid mesh point closest to a given location, e.g. the tip of an elastic beam in FSI."
tags: ["utility", "function-object", "fsi", "solid-mechanics", "monitoring"]
weight: 10
---

## Overview

`solidPointDisplacement` is a function object (library `libmultiRegionFOs`)
that monitors the displacement of one point of a solid region. At start-up it
picks the mesh point closest to a user-given location. At every execution it
writes the time, the three components and the magnitude of the point
displacement field `pointD` at that point to a `.dat` file. `pointD` is
provided by the `solidStVK` region type. The typical use is the tip
displacement of an elastic structure in FSI benchmarks such as Hron-Turek FSI3.

## Method

1. **Region.** The `region` keyword names the solid mesh. If it is omitted,
   `solid` is used when a mesh of that name exists, otherwise `region0`.
2. **Point selection** (constructor, once). The vector `point` is read.
   On each processor the mesh point with the smallest distance
   $|\mathbf{x}_i - \mathbf{x}_\text{probe}|$ is found by a brute-force search
   over all points of the **reference (undeformed)** mesh. Only the processor(s)
   holding the global minimum distance keep their point. If several
   processors tie (e.g. a point on a processor boundary), the one with the
   highest processor number is kept. Each processor that keeps a point prints
   its distance to the given location. If no processor ends up with a point,
   the run aborts.

   Note: this search always uses the mesh `solid` if it exists, otherwise
   `region0`, regardless of the `region` keyword (see
   [Known Limitations](#known-limitations)).
3. **Sampling** (`execute()`, once per time step by default). In mesh
   `region`, the `pointVectorField` `pointD` is looked up. Its value at the
   chosen point is taken without interpolation. All other processors
   contribute zero, and the values are summed over processors. The master
   writes one line

   $$
   t \quad D_x \quad D_y \quad D_z \quad |\mathbf{D}|
   $$

   If `pointD` does not exist in that mesh, a message
   `pointVectorField pointD not found` is printed and nothing is written for
   that step.

`write()` does nothing and `read()` ignores the dictionary, so changes to the
entry during the run (`runTimeModifiable`) have no effect.

## Usage

`system/controlDict` entry (function object):

```
functions
{
    displacement
    {
        type            solidPointDisplacement;
        libs            (multiRegionFOs);
        region          solid;
        point           (0.6 0.2 0.0075);
    }
}
```

| Keyword | Description | Required | Default |
|---------|-------------|----------|---------|
| `type` | `solidPointDisplacement` | yes | — |
| `libs` | `(multiRegionFOs)` | yes | — |
| `point` | Location (m); the closest solid mesh point is monitored | yes | — |
| `region` | Name of the solid mesh that holds `pointD` | no | `solid` if that mesh exists, else `region0` |

The standard function object controls that OpenFOAM handles outside the
object itself (e.g. `executeControl`, `executeInterval`, `timeStart`,
`timeEnd`) can be added. The object itself reads no other keywords.

## Input / Output

**Read**

| Field | Type | Mesh | Description |
|-------|------|------|-------------|
| `pointD` | `pointVectorField` | `region` | Point displacement (m); registered by `solidStVK` |

**Written** (master processor only):

```
<case>/postProcessing/<startTime>/solidPointDisplacement_<name>.dat
```

`<name>` is the function object entry name (e.g. `displacement`) and
`<startTime>` is the case start time. In parallel runs the file is still
written to the undecomposed case directory (`processor*/../postProcessing`).
The file starts with the header `# Time Dx Dy Dz magD`, followed by one
space-separated line per execution. The file is opened anew at start-up, so a
restart from the same start time overwrites the earlier file.

## Known Limitations

- **`region` is only partly honoured.** The closest-point search always uses
  mesh `solid` (or `region0` if there is no `solid` mesh). The sampling uses
  `region`. In practice the solid region must be called `solid`. With another
  name the search falls back to `region0`, which does not exist in a
  multi-region case, and the run aborts.
- **Nearest mesh point, no interpolation.** The monitored location is the
  closest existing mesh point, not the exact `point`. The distance is printed
  at start-up, so check it in the log.
- **Requires `pointD`.** Only region types that register a `pointVectorField`
  called `pointD` (currently `solidStVK`) can be monitored. Other fields cannot
  be selected.
- **Reference configuration.** The point is selected on the solid mesh
  points at start-up. This is correct for the total-Lagrangian `solidStVK`,
  whose mesh never moves.
- **Brute-force search.** The point search is $O(N_\text{points})$ per
  processor. This is negligible since it runs only once.
- **Distributed cases.** With distributed data (no shared case directory),
  writing to `processor*/../postProcessing` can fail.
- **Header comment.** The class description in the header refers to
  "temperature" and least-squares vol-to-point interpolation. Both are wrong:
  the object reads `pointD` directly.

## Tutorials

| Tutorial | Description |
|----------|-------------|
| `tutorials/fluidStructureInteraction/HronTurekFsi3` | Beam tip displacement at `(0.6 0.2 0.0075)` (entry `displacement` in `system/functions`, included from `controlDict`); plotted by `deflection.gnuplot` from `postProcessing/0/solidPointDisplacement_displacement.dat` |
| `tutorials/fluidStructureInteraction/beamInCrossFlow` | Ships a `deflection.gnuplot` that reads the same file, but its `controlDict` defines **no** `solidPointDisplacement` entry, so the file is not produced |
