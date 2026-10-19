---
title: "Contributing"
linkTitle: "Contributing"
weight: 5
description: >
  How to contribute to multiPhysicsFoam.
---

Thank you for considering contributing to **multiPhysicsFoam**! Whether you
fix a bug, add a new region or interface type, improve a tutorial or the
documentation, every contribution helps.

The source code is hosted on Bitbucket:
[bitbucket.org/hmarschall/multiphysicsfoam](https://bitbucket.org/hmarschall/multiphysicsfoam).

## Reporting Bugs and Requesting Features

```warning
**There is currently no issue tracker.** Atlassian has deprecated the built-in
Bitbucket issue tracker and is sunsetting Bitbucket Issues and Wikis completely
on August 20, 2026. The former multiPhysicsFoam issue tracker on Bitbucket is
therefore no longer available, and no replacement has been set up yet.

Until a new platform is announced on this page, bug reports and feature
requests cannot be filed publicly. If you already have a fix or an
implementation, please submit it directly as a pull request (see below):
pull requests are not affected by the sunset.
```

Once a new issue tracker is available, please include the following in a bug
report:

- A concise and descriptive title.
- Steps to reproduce the issue, and the behaviour you expected.
- Relevant details, such as logs, error messages or screenshots, and a
  description of your environment:

  ```plaintext
  - OS version: Ubuntu 24.04 LTS
  - OpenFOAM version: v2512
  - multiPhysicsFoam version/commit: <git commit hash>
  ```

- Ideally, a minimal case that reproduces the problem.

Feature requests should be titled **Feature request: \<feature title\>** and
describe the proposed feature and its benefit for multiPhysicsFoam.

## Submitting Code Changes

### Step 1: Fork and Clone

1. Fork the repository on Bitbucket.
2. Clone your fork and create a branch for your change:

   ```bash
   git clone https://bitbucket.org/<your-username>/multiphysicsfoam.git
   cd multiphysicsfoam
   git checkout -b feat/my-new-feature
   ```

### Step 2: Make Your Changes

- Follow the
  [OpenFOAM coding style guide](https://develop.openfoam.com/Development/openfoam/-/wikis/coding/style/style)
  and the structure and naming of the existing code.
- The code must compile against an **unmodified** OpenFOAM installation of
  every supported version (currently OpenFOAM v2512 and v2606). Do not require
  patches to OpenFOAM itself; missing core functionality belongs in
  multiPhysicsFoam.
- Build everything with `./Allwmake` and check that it finishes without errors.
- New region types, interface types and boundary conditions are compiled into
  `libmultiRegionSystem`: add the `.C` file to `src/multiRegionSystem/Make/files`.
  See [Adding Region Physics Types]({{< relref "adding-region-physics-types" >}}),
  [Adding Interface Physics Types]({{< relref "adding-interface-physics-types" >}})
  and [Adding Coupled BCs]({{< relref "adding-coupled-boundary-conditions" >}}).
- There is no unit-test suite: the **tutorials are the tests**. Run the
  tutorials affected by your change, and add or update a tutorial case if you
  add new functionality.
- Document your code with clear comments.

### Step 3: Document Your Changes

Every region type, interface type, boundary condition and utility has a page
in the `docs/` directory of the code repository. These pages are used in
this website, so please add or update them along with your code:

| Component | Documentation file | Template |
|-----------|--------------------|----------|
| Region type | `docs/regions/<regionType>.md` | `docs/regions/_template.md` |
| Interface type | `docs/interfaces/<interfaceType>.md` | `docs/interfaces/_template.md` |
| Boundary condition | `docs/boundaryConditions/{coupled,uncoupled}/<bcType>.md` | `docs/boundaryConditions/_template.md` |
| Utility | `docs/utilities/<name>.md` | `docs/utilities/_template.md` |

Keep the front matter of these files Hugo-compatible (`title`, `description`,
`tags`, `weight`).

### Step 4: Commit Your Changes

Write clear, descriptive commit messages following the
[Conventional Commits](https://www.conventionalcommits.org/) format
`<type>[(<scope>)]: <description>`, where `<type>` is one of `feat`, `fix`,
`docs`, `port`, `style`, `refactor`, `perf`, `test`, `build`, `ci`, `chore` or
`revert`:

```bash
git add src/regions/myRegion/myRegion.C
git commit -m "feat(regions): add myRegion region type"
```

### Step 5: Push and Open a Pull Request

Push your branch to your fork:

```bash
git push -u origin feat/my-new-feature
```

Then open a pull request against the `main` branch of the multiPhysicsFoam
repository on Bitbucket. Include:

- **A clear title** in the same format as the commit messages, e.g.
  `fix(tutorials): port HronTurekFsi3 to OpenFOAM.com`.
- **A detailed description** of your changes and the motivation behind them.
- **Instructions or example cases** to test your contribution, e.g. which
  tutorials you ran.

## Contributing to the Documentation

Documentation of individual components belongs in the code repository (see
[Step 3](#step-3-document-your-changes)). General pages of this website, such
as the getting-started guide or this page, live in the website repository on
GitHub:
[github.com/CoHabes/multiPhysicsFoamDocumentation](https://github.com/CoHabes/multiPhysicsFoamDocumentation).
Pull requests there are welcome as well.

Thank you for contributing to multiPhysicsFoam!
