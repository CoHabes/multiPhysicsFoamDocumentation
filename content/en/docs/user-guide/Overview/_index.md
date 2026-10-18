---
title: "Overview"
linkTitle: "Overview"
weight: 1
description: >
  How multiPhysicsFoam works: the multiPhysicsFoam executable, regions,
  interfaces and the coupling of transport equations between them.
---

multiPhysicsFoam is a library and solver for OpenFOAM (.com, v2512 and v2606)
that simulates **multi-region, coupled multiphysics problems**. Instead of
writing a new solver for every combination of physics, a problem is assembled
from building blocks:

- **Regions** – separate meshes, each carrying one or more physics models
  (*region types*), e.g. a fluid solved with `pimpleFluid` plus heat transport
  with `transportTemperature`, or a solid solved with `conductTemperature`.
- **Interfaces** – pairs of patches through which two regions are coupled, each
  with an *interface type* that defines the interface physics, e.g.
  `heatTransferInterface` or `fsiInterface`.
- **Coupled boundary conditions** – the patch fields on the interface patches
  that exchange values and fluxes between the regions.

All of this is configured in dictionaries at run time. The same executable,
`multiPhysicsFoam`, runs every case.

## The `multiPhysicsFoam` executable

`multiPhysicsFoam` itself only contains the time loop. On start-up it creates a
`multiRegionSystem`, which reads the case configuration and builds all regions
and interfaces. Then, in every time step, it

1. adjusts the time step to the smallest stable time step reported by all
   regions (adaptive time stepping, limited by `maxDeltaT`),
2. calls `preSolve()` and `solve()` of the multi-region system,
3. writes the results.

All physics lives in the `libmultiRegionSystem` library, so new region,
interface or boundary condition types extend the framework without touching
the solver.

### What happens in one time step

`multiRegionSystem::solve()` advances the coupled problem in the following
order:

1. **Region-local physics** – every region solves the equations that are not
   coupled to other regions (e.g. a VOF phase-fraction equation or the
   abuse reactions of a battery model).
2. **Pressure–velocity coupling** – regions that use the PIMPLE algorithm
   (fluid regions) run their momentum predictor and pressure corrector,
   including mesh motion where needed.
3. **Self-coupled fields** – fields that each region solves on its own, without
   exchange with a neighbour.
4. **Partitioned coupling** – for every field that is coupled across an
   interface (e.g. `T` in conjugate heat transfer), the equations of all
   regions sharing that field are solved repeatedly in a
   **Dirichlet–Neumann (DNA) loop**. One side of each interface receives a
   value (Dirichlet) and the other a flux (Neumann) from its partner. The loop
   repeats until the jump and flux residuals at the interfaces fall below the
   specified tolerances or `maxCoupleIter` is reached. Convergence can be
   accelerated with fixed relaxation, Aitken or IQN-ILS.

   Two coupled fields have their own loops: `pUPimple` for sharp-interface
   multiphase flows and `D` for fluid–structure interaction. They also run the
   PIMPLE algorithm, mesh-motion correction and the interface update inside
   each coupling iteration.

Monolithic coupling (all regions in one matrix) is planned; currently,
partitioned coupling is the supported mode.

## Setting up a case

A multi-region case looks like a standard OpenFOAM case with one sub-directory
per region:

```
case/
├── 0/<region>/                         initial and boundary conditions per region
├── constant/
│   ├── multiRegionProperties           regions, their types, DNA settings
│   ├── regionInterfaceProperties       interfaces and coupled fields
│   └── <region>/                       mesh and physical properties per region
└── system/
    ├── controlDict                     time control (shared by all regions)
    └── <region>/                       fvSchemes, fvSolution per region
```

### `multiRegionProperties` – regions

The `regions` list names every region and the region types stacked on it. The
`DNA` dictionary configures the partitioned coupling loop per coupled field:

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
            maxFluxRes 1e-07;
        }
    }
}
```

### `regionInterfaceProperties` – interfaces

Each entry of `partitionedCoupledPatches` defines one interface: its interface
type, the two patches it connects, the fields coupled across it and the
interface-specific coefficients:

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
            jumpCoefficient jumpCoefficient [ 0 0 0 0 0 0 0 ] 1;
        }
    }
);
```

Values are transferred between non-conformal interface patches by
interpolation (e.g. AMI or direct mapping).

### Coupled boundary conditions

On the interface patches, the coupled field uses a matching pair of
**coupled boundary conditions**: a value (jump) condition on one side and a
flux condition on the other. For the example above, `0/fluid/T` and
`0/solid/T` contain:

```
// 0/fluid/T
interface
{
    type            regionCoupledTemperatureJump;
    interfaceType   (heatTransferInterface);
    accType         aitken;
    relax           0.75;
    value           $internalField;
}
```

```
// 0/solid/T
top
{
    type            regionCoupledHeatFlux;
    interfaceType   (heatTransferInterface);
    value           $internalField;
}
```

### Running the case

Each region is meshed separately, with its mesh dictionaries in
`system/<region>/`, and the solver is started as usual:

```bash
blockMesh -region fluid
blockMesh -region solid
multiPhysicsFoam                              # serial

decomposePar -region fluid                    # parallel: decompose every region
decomposePar -region solid
mpirun -np 4 multiPhysicsFoam -parallel
```

The multi-region tutorials contain an `Allrun` script that does all of this,
e.g. `./Allrun -p -S` for partitioned coupling in serial (`-P` for parallel).

## Single-region mode

A single region type can also be run on a plain OpenFOAM case, without
`multiRegionProperties`, `regionInterfaceProperties` or region
sub-directories. Two entries in `system/controlDict` are enough:

```
singleRegionMode    yes;
regionType          pimpleFluid;
```

The case then behaves like the corresponding stock OpenFOAM solver (here
`pimpleFoam`), run through the multiPhysicsFoam infrastructure. This is useful
to verify a region type against the stock solver, or as a first step before
setting up a coupled case. The `tutorials/singleRegion/<regionType>/` cases
use this mode.

## Where to go next

- [Region Physics]({{< relref "region-physics" >}}) – available region types
  and their settings
- [Interface Physics]({{< relref "interface-physics" >}}) – available interface
  types
- [Boundary Conditions]({{< relref "boundary-conditions" >}}) – coupled and
  custom boundary conditions
- [Utilities]({{< relref "utilities" >}}) – pre- and post-processing tools
- [Tutorials]({{< relref "tutorials" >}}) – complete example cases
