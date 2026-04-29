---
title: "Quickstart"
linkTitle: "Quickstart"
weight: 2
description: >
  Download, compile, and verify multiPhysicsFoam.
---

## Quickstart

Source a [supported version of OpenFOAM](../installation-from-source/#supported-versions-of-openfoam),
then download, build and test multiPhysicsFoam-v0.1:

```bash
git clone --branch v0.1 https://github.com/multiPhysicsFoam/multiPhysicsFoam.git
cd multiPhysicsFoam && ./Allwmake -j && cd tutorials && ./Alltest
```

The `-j` flag instructs `Allwmake` to use all available cores.

Similarly, the latest  the `development` branch can be downloaded and compiled
 with

```bash
git clone --branch development https://github.com/multiPhysicsFoam/multiPhysicsFoam.git
cd multiPhysicsFoam && ./Allwmake -j && cd tutorials && ./Alltest
```