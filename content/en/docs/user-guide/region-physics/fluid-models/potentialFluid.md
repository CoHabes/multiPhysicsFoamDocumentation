---
mathjax: true
title: "potentialFluid"
description: "Inviscid, irrotational, incompressible potential flow; solves a Laplace equation for the velocity potential to generate a divergence-free initial flow field."
tags: ["region-type", "fluid", "incompressible", "inviscid", "initialisation"]
weight: 10
---

## Overview

`potentialFluid` solves for **potential flow**: an inviscid, irrotational,
incompressible idealisation in which the velocity field is the gradient of a
scalar velocity potential $\Phi$, reducing the flow problem to a single Laplace
equation. It is the multiPhysicsFoam port of OpenFOAM's `potentialFoam`.

Its primary use is **initialisation**: it produces a smooth, discretely
divergence-free velocity and flux field consistent with the boundary conditions,
from which a full Navier-Stokes region such as `pimpleFluid` can start without a
violent initial transient. It is also a cheap mesh-quality diagnostic. Because
the model contains no viscosity, it resolves no boundary layers, separation,
wakes, or drag, and is not intended as a final-answer solver.

The region re-solves on **every time step**. For the usual one-shot initialiser
workflow, limit the run to a single step via `controlDict` (`endTime 1`,
`deltaT 1`), as both tutorials do.

## Governing Equations

### Continuity and irrotationality

$$
\nabla \cdot \mathbf{U} = 0
\qquad\qquad
\nabla \times \mathbf{U} = 0
$$

### Velocity potential

Irrotationality permits $\mathbf{U} = \nabla \Phi$, which with continuity gives
Laplace's equation:

$$
\nabla^2 \Phi = 0
$$

Discretely this is solved in projection form. From the flux field $\phi^{*}$
implied by the boundary conditions,

$$
\nabla^2 \Phi = \nabla \cdot \phi^{*}
\qquad\text{then}\qquad
\phi = \phi^{*} - \nabla \Phi \cdot \mathbf{S}_f
$$

so that $\nabla \cdot \phi = 0$ to solver tolerance on the actual mesh. The cell
velocity is recovered by `fvc::reconstruct(phi)`.

### Euler pressure (optional)

Enabled by `writep`. Taking the divergence of the steady Euler equation, with
the source filtered along the flow direction by
$\mathbf{F} = \hat{\mathbf{U}} \otimes \hat{\mathbf{U}}$ to suppress spurious
cross-stream pressure variation:

$$
\nabla^2 p_\text{kin} + \nabla \cdot \left( \mathbf{F} \cdot \nabla \cdot (\mathbf{U}\mathbf{U}) \right) = 0
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $\Phi$ | Velocity potential | m²/s |
| $\mathbf{U}$ | Velocity | m/s |
| $\phi$ | Volumetric face flux | m³/s |
| $p_\text{kin} = p/\rho$ | Kinematic pressure | m²/s² |
| $\mathbf{S}_f$ | Face area vector | m² |
| $\mathbf{F}$ | Flow-direction filter tensor | – |

## Numerical Algorithm

There is no PIMPLE loop and no time integration: `usesPIMPLE()` returns `false`
and the whole algorithm runs in `solveRegion()`, which the orchestrator calls
before any PIMPLE stage. Per run:

1. **Make the flux relative** to any MRF zones.
2. **Balance the global flux** (`adjustPhi`). Applied only when $p_\text{kin}$
   needs a reference, i.e. no fixed-value pressure boundary exists; total inflow
   is then forced to equal total outflow.
3. **Non-orthogonal corrector loop** — assemble and solve
   $\nabla^2\Phi = \nabla\cdot\phi$ with `setReference` applied. On the final
   iteration only, correct the flux: `phi -= PhiEqn.flux()`.
4. **Restore the absolute flux.**
5. **Report the continuity error** $\overline{|\nabla \cdot \phi|}$, volume-weighted.
6. **Reconstruct** `U` from `phi` and correct its boundary conditions.
7. **Report the interpolated velocity error** — the drift between
   `fvc::flux(U)` and `phi`; a mesh-quality diagnostic.
8. **Solve the Euler pressure** in its own non-orthogonal loop, if `writep`.
9. **Update** the dimensional pressure $p = p_\text{kin}\rho$.

A reference cell is required because an all-Neumann Laplacian is singular
($\Phi$ and $\Phi + C$ are both solutions). `setReference` derives its keywords
from the field name, so the case dictionary uses `PhiRefCell`/`PhiRefValue` and
`pKinRefCell`/`pKinRefValue`.

Steps 1–9 execute on every time step; the solve is not latched. Since the
solution is fully determined by the boundary conditions, additional steps
recompute the same answer unless a sibling region has perturbed `U` or `phi`
in the meantime — so restrict the run to one step unless you want that
re-establishment.

## Sub-Models

| Sub-model | Selection mechanism | Notes |
|-----------|---------------------|-------|
| Moving Reference Frame (MRF) | `constant/MRFProperties` | `IOMRFZoneList`; multiple zones supported |
| Density | `constant/transportProperties` | `rho`, used only to form the dimensional `p` |
| Turbulence | — | Not supported; the model is inviscid |
| `fvOptions` | — | Not supported |

## Required Fields

| Field | Type | Description | Unit |
|-------|------|-------------|------|
| `U` | `volVectorField` | Velocity. Only its boundary conditions are used — the internal field is zeroed at construction. | m/s |
| `pKin` | `volScalarField` | Kinematic pressure $p/\rho$. Optional; BCs inferred from `U` if absent. | m²/s² |

Additionally created at runtime:

| Field | Type | Description | Condition |
|-------|------|-------------|-----------|
| `phi` | `surfaceScalarField` | Face-normal flux, from `fvc::flux(U)` | always |
| `Phi` | `volScalarField` | Velocity potential; BCs inferred from `pKin` if absent | always (written only if `writePhi`) |
| `p` | `volScalarField` | Dimensional pressure $p_\text{kin}\rho$, for sibling regions | always |

Where `pKin` and `Phi` are not supplied, their boundary types are inferred by
inverting Dirichlet and Neumann at each step: `pKin` is `zeroGradient` wherever
`U` fixes value and `fixedValue` elsewhere; `Phi` inverts that again. Supplying
`0/U` alone is therefore sufficient.

## Case Setup

### Required dictionaries

```
constant/transportProperties   – density (rho); required, read MUST_READ
system/fvSolution              – Phi/pKin linear solvers, potentialFlow sub-dict
system/fvSchemes               – discretisation schemes
```

Unlike stock `potentialFoam`, `constant/transportProperties` is mandatory:

```
rho             1000;
```

### `multiRegionProperties` entry

```
regions
(
    ( fluid (potentialFluid) )
);
```

To initialise a field that another region then solves, list `potentialFluid`
first so that it owns and constructs the shared fields:

```
regions
(
    ( fluid (potentialFluid pimpleFluid) )
);
```

### `controlDict` (singleRegionMode)

```
singleRegionMode    yes;
regionType          potentialFluid;

endTime             1;      // one step == one potential flow solve
deltaT              1;
adjustTimeStep      no;     // getMinDeltaT() returns GREAT
```

### Key `fvSolution` settings

```
solvers
{
    Phi
    {
        solver          GAMG;
        smoother        DIC;
        tolerance       1e-06;
        relTol          0.01;
    }

    pKin
    {
        $Phi;
    }
}

potentialFlow
{
    nNonOrthogonalCorrectors 3;

    writePhi        true;
    writep          true;

    // Reference cells (closed domains)
    PhiRefCell      0;
    PhiRefValue     0;
    pKinRefCell     0;
    pKinRefValue    0;
}
```

### `potentialFlow` dictionary entries

`potentialFoam`'s command-line flags become optional dictionary entries, all
read as `Foam::Switch` and defaulting to `false`:

| Entry | Upstream flag | Effect |
|-------|---------------|--------|
| `initialiseUBCs` | `-initialiseUBCs` | Call `U.correctBoundaryConditions()` and rebuild `phi` before solving |
| `writePhi` | `-writePhi` | AUTO_WRITE the velocity potential `Phi` |
| `writep` | `-writep` | **Calculate** the approximate Euler pressure; gates the calculation only |

Being `Switch`-valued, each accepts any of `true`/`false`, `yes`/`no`,
`on`/`off` or `none`/`any`; the tutorials use `true`/`false`.

Dropped: `-pName` (the field is always `pKin`), `-withFunctionObjects`
(function objects run automatically), and `-writephi` (`phi` is always written).

## Coupling

`potentialFluid` implements no interface-specific behaviour: `setCoupledEqns()`
is a no-op, it declares no `selfCoupledFieldNames()`, and it takes no part in
DNA partitioned coupling.

| Interface type | Exchanged fields | Coupling modes |
|----------------|-----------------|----------------|
| *(none)* | – | – |

Coupling is implicit, through the object registry: the fields it constructs
(`U`, `phi`, `pKin`, `p`) are visible by name to sibling regions on the same
mesh, which is what enables the initialiser role.

## Known Limitations

- **Re-solves every time step** — the solve is not latched, but the solution is
  fully determined by the boundary conditions, so a long `endTime` repeats the
  same computation. Limit the run to a single step in `controlDict`.
- **`writep false` writes a zero `pKin` and `p`** — the Euler pressure is never
  calculated, but both fields are AUTO_WRITE regardless.

## Tutorials

| Tutorial path | Description |
|---------------|-------------|
| `tutorials/singleRegion/potentialFluid/cylinder` | Flow around a cylinder; includes a `coded` function object computing the analytic solution and a relative-error field |
| `tutorials/singleRegion/potentialFluid/pitzDaily` | Backward-facing step; smoke test, and a demonstration of the model's limits (the real flow separates) |

Both reproduce the corresponding `$FOAM_TUTORIALS/basic/potentialFoam/` case
field-for-field, with `0/p` renamed to `0/pKin`, the command-line flags moved
into the `potentialFlow` dict, and `constant/transportProperties` added.

## References

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics Framework for Multi-Region Coupled Continuum-Physical Problems*, Engineering with Computers, 2023.
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
- Batchelor, G.K., *An Introduction to Fluid Dynamics*, CUP, 1967 — Ch. 6, irrotational flow theory.
- Upstream source: `$FOAM_APP/solvers/basic/potentialFoam/`
