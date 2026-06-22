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
| multiPhysicsFoam-v0.1    | OpenFOAM-v2506                   |

multiPhysicsFoam is developed using the OpenFOAM.com variant of OpenFOAM
(i.e., OpenFOAM-v2506); consequently, this is the recommended variant for using
 multiPhysicsFoam.

### Optional fixes for the OpenFOAM installation

The multiPhysicsFoam `Allwmake` script will ask you to fix files in the main
 OpenFOAM/foam installation. These checks will alwys be performed.
 The fixes will not break you OpenFOAM installation nor any other applications compiled
 against you OpenFOAM version. You will get clear information on what to do when running `Allwmake` for the first time
 e.g.

```bash
> Check filesToReplace

******** PLEASE FIX THIS ***********
You should replace the file '/usr/lib/openfoam/openfoam2512/src/OpenFOAM/fields/GeometricFields/GeometricField/GeometricField.H' with 'filesToReplace/GeometricField.H'

You can do it by running the following commands:
    cp filesToReplace/GeometricField.H /usr/lib/openfoam/openfoam2512/src/OpenFOAM/fields/GeometricFields/GeometricField/
    wmake libso /usr/lib/openfoam/openfoam2512/src/OpenFOAM
If you have not compiled OpenFOAM from scratch you might have to recompile it as a whole (sorry for that)

These actions will not break your working foam-extend installation.
It will just fix some minor bugs in the listed file which need to 
be fixed in order to make multiPhysicsFoam work propperly.
************************************

```

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
should be sourced: see the [table above](#supported-versions-of-openfoamfoam). To
build multiPhysicsFoam, enter the multiPhysicsFoam directory and execute the included
Allwmake script, e.g.

```bash
> cd multiPhysicsFoam
> ./Allwmake 2>&1 | tee log.Allwmake
```

Depending on your hardware, you can expect this build to last about 5-10 minutes.

If multiPhysicsFoam is built successfully, you will be presented with the following
message:

```bash
There were no build errors: enjoy multiPhysicsFoam

To test the installation test any of the provided tutorials, eg:
    > cd tutorials/conjugateHeatTransfer/flowOverHeatedPlate 
    > ./Allrun -p -P
```

If the build encounters errors, you will receive the following message:

```bash
** BUILD ERROR **"
There were build errors in the following logs:
<LIST OF COMPILATION ERRORS>

Please examine these logs for additional details
```

You can examine the source of the errors in the `log.Allwmake` file within the
multiPhysicsFoam parent directory.
Alternatively, if you believe you have encountered a bug, please create a new
issue at
[https://bitbucket.org/hmarschall/multiphysicsfoam/issues?status=new&status=open&status=submitted&is_spam=!spam](https://bitbucket.org/hmarschall/multiphysicsfoam/issues?status=new&status=open&status=submitted&is_spam=!spam).

---

## What Next?

Please see the [tutorial guide](../tutorials/README.md).
