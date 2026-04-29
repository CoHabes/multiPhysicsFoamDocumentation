---
title: "Beam in Cross Flow"
linkTitle: "Beam in Cross Flow"
weight: 1
description: >
  Flexible beam deflected by a laminar cross-flow — FSI validation case.
---

**Repository path:** `tutorials/fluidStructureInteraction/beamInCrossFlow`

## Overview

A slender elastic beam is immersed in a laminar cross-flow channel. The fluid
exerts a drag and lift force on the beam which deflects and oscillates. This
case validates the partitioned FSI coupling with ALE moving meshes.

## Running the case

```bash
cd $FOAM_RUN/../multiPhysicsFoam/tutorials/fluidStructureInteraction/beamInCrossFlow
./Allrun
```
