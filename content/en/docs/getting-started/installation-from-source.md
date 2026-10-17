---
title: "Installing from Source"
linkTitle: "Installing from Source"
weight: 3
description: >
  Download, compile, and verify multiPhysicsFoam.
---

## Detailed Installation Guide

---

### Supported Versions of OpenFOAM

multiPhysicsFoam requires a working version of OpenFOAM. Currently,
the following OpenFOAM versions are supported:

| multiPhysicsFoam version | OpenFOAM version                 |
| ------------------------ | -------------------------------- |
| multiPhysicsFoam-v0.1    | OpenFOAM-v2512                   |

multiPhysicsFoam is developed using the OpenFOAM.com variant of OpenFOAM
(i.e., OpenFOAM-v2512); consequently, this is the recommended variant for using
 multiPhysicsFoam.

multiPhysicsFoam builds against an unmodified OpenFOAM installation, e.g. the
precompiled packages from openfoam.com. No files of OpenFOAM itself need to be
changed or recompiled.

---

### Dependencies

> These dependencies are optional. You can skip them if you want to get up and
 running quickly.

Beyond a working version of OpenFOAM, multiPhysicsFoam does not have
any **mandatory** dependencies; however, several **optional** dependencies are
required to use the complete set of functionalities:

| Dependency | Functionality Provided                                           |
| ---------- | ---------------------------------------------------------------- |
| cfmesh     | Some tutorials use cfmesh for creating the meshes.               |
| gnuplot    | Some tutorials use Gnuplot to generate graphs after running the  |
|            | solver.                                                          |

#### cfmesh

A small number of tutorials require the `cartesianMesh` utility from cfmesh. If
the `cartesianMesh` executable is not found within the `$PATH` then the `Allrun`
script within these tutorials will exit.

```tip
- **OpenFOAM.com (OpenCFD/ESI version)**: cfmesh is provided as an optional
   submodule: see [https://develop.openfoam.com/Community/integration-cfmesh](https://develop.openfoam.com/Community/integration-cfmesh).
```

#### gnuplot

gnuplot can be installed on Ubuntu with:

```bash
> sudo apt-get install gnuplot
```

Or, on macOS with

```bash
> brew install gnuplot
```

---

### Downloading the multiPhysicsFoam Source Code

The multiPhysicsFoam directory can be downloaded to any reasonable location on your
computer; we suggest placing it in `$FOAM_RUN/..`.

#### Git repository: v0.1

`multiPhysicsFoam-v0.1` can be downloaded using

```bash
> git clone --branch v0.1 https://github.com/multiPhysicsFoam/multiPhysicsFoam.git
```

#### Git repository: latest development branch

The latest nightly build development branch can be downloaded with

```bash
> git clone --branch development git@github.com:multiPhysicsFoam/multiPhysicsFoam.git
```

---

### Building multiPhysicsFoam

Before building multiPhysicsFoam, a compatible version of OpenFOAM
should be sourced: see the [table above](#supported-versions-of-openfoam). `Allwmake`
stops with a message if no or an unsupported OpenFOAM version is sourced. To
build multiPhysicsFoam, enter the multiPhysicsFoam directory and execute the included
Allwmake script, e.g.

```bash
> cd multiPhysicsFoam
> ./Allwmake
```

Depending on your hardware, you can expect this build to last about 5-10 minutes.
To compile in parallel, set the number of cores beforehand, e.g.
`export WM_NCOMPPROCS=8`.

If multiPhysicsFoam is built successfully, you will be presented with the following
message:

```bash
There were no build errors: enjoy multiPhysicsFoam

To test the installation, run one of the tutorials, e.g.:
    cd tutorials/conjugateHeatTransfer/flowOverHeatedPlate
    ./Allrun -p -S
```

If the build encounters errors, you will receive the following message:

```bash
** BUILD ERROR **"
There were build errors in the following logs:
<LIST OF COMPILATION ERRORS>

Please examine these logs for additional details
```

Each build step writes its own log file, `log.Allwmake.<name>`, in the
multiPhysicsFoam directory; the message lists the ones containing errors.
Alternatively, if you believe you have encountered a bug, please create a new
issue at
[https://bitbucket.org/hmarschall/multiphysicsfoam/issues?status=new&status=open&status=submitted&is_spam=!spam](https://bitbucket.org/hmarschall/multiphysicsfoam/issues?status=new&status=open&status=submitted&is_spam=!spam).

---

## What Next?

Please see the [tutorial guide](../tutorials/README.md).
