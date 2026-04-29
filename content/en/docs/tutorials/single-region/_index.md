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
