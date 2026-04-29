# multiPhysicsFoam Website — Project Context

## Project Goal

This is the official documentation website for **multiPhysicsFoam**, a modular multiphysics framework for interface and volumetric coupled simulations based on OpenFOAM.

The website serves as the primary documentation hub for:
- Users running multiphysics simulations (installation, solver guides, boundary conditions, utilities)
- Developers extending the framework (architecture, adding solvers, contributing)
- Tutorial cases

## Tech Stack

- **Static site generator**: Hugo `0.124.1` (extended)
- **Theme**: [google/docsy](https://github.com/google/docsy) `v0.6.0` (git submodule in `themes/docsy`)
- **Node**: `18` (Bootstrap, FontAwesome via npm)
- **Deployment**: GitHub Pages via `.github/workflows/deploy.yml`
- **Theme mode**: Dark only (no theme switcher)

## Key Files

| File | Purpose |
|------|---------|
| `config.toml` | Site config: title, nav menu, GitHub links, params |
| `content/en/_index.html` | Homepage (hero section + feature cards) |
| `content/en/docs/` | All documentation content |
| `assets/scss/main.scss` | SCSS entry point (dark theme only) |
| `assets/scss/_variables_project.scss` | Color variables (`$mpf-*`), font, container widths |
| `assets/scss/_styles_project.scss` | Custom component styles |
| `assets/scss/dark/` | Dark theme Bootstrap and Docsy overrides |
| `layouts/partials/navbar.html` | Custom navbar (no theme toggler) |
| `.github/workflows/deploy.yml` | GitHub Pages deployment workflow |

## Content Structure

```
content/en/docs/
├── getting-started/       # Installation, prerequisites, quick start
├── user-guide/
│   ├── solvers/           # Available multiphysics solvers
│   ├── boundary-conditions/
│   └── utilities/
├── developer-guide/
│   ├── architecture/      # Modular framework design
│   ├── adding-solvers/    # How to implement new physics
│   └── contributing/
├── reference/             # API / configuration reference
├── examples/              # Tutorial cases
└── about/                 # Project info, license, community
```

## Important Conventions

- **Dark theme only**: `assets/scss/main.scss` imports `dark/bootstrap/main`, `dark/variables`, and `dark/styles` at the top level — no `[data-theme]` scoping
- **Color variables**: All color vars are prefixed `$mpf-` (e.g. `$mpf-blue: #4279f4`), defined in `_variables_project.scss` and `dark/bootstrap/_dark-variables.scss`
- **No version selector**: `config.toml` has no `[[params.versions]]` entries; the navbar version dropdown is hidden
- **No theme toggler**: The `theme-toggler` partial is not called in `navbar.html`
- **Docsy overrides**: Project-level overrides in `layouts/partials/` and `assets/scss/` take precedence over the theme — do not edit files inside `themes/docsy/`

## Owner

Constantin Habes
