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
|                          | OpenFOAM-v2512                   |

multiPhysicsFoam is primarily developed using the OpenFOAM.com variant of OpenFOAM
(i.e., OpenFOAM-v2512); consequently, this is the recommended variant for using
 multiPhysicsFoam.

### Optional fixes for the OpenFOAM installation

The multiPhysicsFoam `Allwmake` script can optionally ask you to fix files in the main
 OpenFOAM/foam installation; by default, these checks are **not** performed. To
 enable these checks, set the `S4F_NO_FILE_FIXES` environmental variable to `0`
 before running the Allwmake script when [building multiPhysicsFoam](#building-multiPhysicsFoam),
 e.g.

```bash
> export S4F_NO_FILE_FIXES=0 && ./Allwmake -j
```

To make these fixes, follow the instructions from the `Allwmake` script when
[building multiPhysicsFoam](#building-multiPhysicsFoam).

See the upstream
[Optional Fixes page](https://github.com/multiPhysicsFoam/multiPhysicsFoam/blob/master/optionalFixes/README.md)
for further details
on these optional changes.

---

### Dependencies

> These dependencies are optional. You can skip them if you want to get up and
 running quickly.

Beyond a working version of OpenFOAM or foam-extend, multiPhysicsFoam does not have
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
- **foam-extend**: cfmesh is included in foam-extend.
- **OpenFOAM.com (OpenCFD/ESI version)**: cfmesh is provided as an optional
   submodule: see [https://develop.openfoam.com/Community/integration-cfmesh](https://develop.openfoam.com/Community/integration-cfmesh).
- **OpenFOAM.org (Foundation version)**: the free version of cfmesh is currently
   not compatible with OpenFOAM.org.
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

Before building multiPhysicsFoam, a compatible version of OpenFOAM or foam-extend
should be sourced: see the [table above](#supported-versions-of-openfoamfoam). To
build multiPhysicsFoam, enter the multiPhysicsFoam directory and execute the included
Allwmake script, e.g.

```bash
> cd multiPhysicsFoam
> ./Allwmake -j 2>&1 | tee log.Allwmake
```

Depending on your hardware, you can expect this build to last about 5 minutes.

If multiPhysicsFoam is built successfully, you will be presented with the following
message:

```bash
There were no build errors: enjoy multiPhysicsFoam!
To test the installation, run:
    > cd tutorials && ./Alltest
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
[https://github.com/multiPhysicsFoam/multiPhysicsFoam/issues](https://github.com/multiPhysicsFoam/multiPhysicsFoam/issues).

---

### Testing the Installation

As instructed, after a successful build, you can test the tutorials using the
following commands, executed from the multiPhysicsFoam parent directory.

```bash
> cd tutorials && ./Alltest
```

If the tests pass, you will receive the message:

```bash
All tests passed: enjoy multiPhysicsFoam.
```

This means your multiPhysicsFoam installation is working as expected.

If any of the tests fail, you will receive the message:

```bash
The multiPhysicsFoam solver failed in the following cases:
<LIST OF FAILING TUTORIALS>
```

or if the errors do not come from the multiPhysicsFoam calls but elsewhere

```bash
The following commands failed:
<LIST OF FAILING COMMANDS AND TUTORIALS>
```

You can expect these _smoke-test_ tests to last a few minutes.
These _smoke tests_ just check that all cases can be solved for one time-step
 or iteration without an error; however, the values of the results are not
 checked.

---

## What Next?

Please see the [tutorial guide](../tutorials/README.md).
