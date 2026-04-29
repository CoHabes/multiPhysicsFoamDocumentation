---
title: "Pitz Daily"
linkTitle: "Pitz Daily"
weight: 1
description: >
  Turbulent flow over a backward-facing step — single-region pimpleFluid case.
---

**Repository path:** `tutorials/singleRegion/pimpleFluid/pitzDaily`

## Overview

The Pitz–Daily case simulates turbulent flow over a backward-facing step. It is
one of the standard OpenFOAM validation cases and is included here to verify
that the `pimpleFluid` region type within `multiPhysicsFoam` reproduces the
expected results.

## Running the case

```bash
cd $FOAM_RUN/../multiPhysicsFoam/tutorials/singleRegion/pimpleFluid/pitzDaily
./Allrun
```
