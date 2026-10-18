---
title: "Single Region"
linkTitle: "Single Region"
weight: 5
description: >
  Single-region tutorial cases demonstrating individual physics solvers within multiPhysicsFoam.
---

Single-region cases use `multiPhysicsFoam` with a single physics region and no
interface coupling. They are useful as a first check that the solver infrastructure
works correctly and as a starting point before setting up multi-region cases.

## Cases

- [Pitz Daily](pimple-fluid/pitz-daily/) — turbulent flow over a backward-facing step
- [Mixer Vessel AMI2D](pimple-fluid/mixer-vessel-ami2d/) — rotating mixer with an arbitrary mesh interface (AMI)
- [Dam Break](inter-fluid/dam-break/) — two-phase collapse of a water column (`interFluid`)
- [Flange](laplacian-fluid/flange/) — transient heat conduction in a flange (`laplacianFluid`)
- [Cylinder](potential-fluid/cylinder/) — potential flow past a cylinder (`potentialFluid`)
- [Pitz Daily (potential flow)](potential-fluid/pitz-daily/) — potential flow over a backward-facing step (`potentialFluid`)
- [Oblique Shock](rho-pimple-fluid/oblique-shock/) — Mach 2.9 oblique shock (`rhoPimpleFluid`)
