# multiPhysicsFoam Website

Welcome to the GitHub repository for the multiPhysicsFoam documentation website.

multiPhysicsFoam is a modular multiphysics framework for interface and volumetric coupled simulations based on OpenFOAM.

The website is built with [Hugo](https://gohugo.io/) and the [google/docsy](https://github.com/google/docsy) theme, and deployed to [GitHub Pages](https://pages.github.com/).

## Local Development

### Prerequisites

- Hugo `0.124.1` (extended version)
- Node.js `18`

### Install Hugo

Download the extended Hugo binary for your platform from the [Hugo releases page](https://github.com/gohugoio/hugo/releases/tag/v0.124.1), or use a package manager:

```bash
# macOS / Linux via Homebrew (pin to the correct version)
HOMEBREW_COMMIT="9d025105a8be086b2eeb3b1b2697974f848dbaac" # 0.124.1
curl -fL -o "hugo.rb" "https://raw.githubusercontent.com/Homebrew/homebrew-core/${HOMEBREW_COMMIT}/Formula/h/hugo.rb"
brew install ./hugo.rb
brew pin hugo
```

### Install Node Dependencies

```bash
npm install
```

### Run the Local Server

1. Clone the repository and initialise the Docsy submodule:

   ```bash
   git clone git@github.com:<your-github-username>/multiPhysicsFoam-website.git
   cd multiPhysicsFoam-website/
   git submodule update --init --recursive
   ```

2. Install node packages:

   ```bash
   npm install
   ```

3. Start Hugo:

   ```bash
   hugo server
   ```

4. Open [http://localhost:1313/](http://localhost:1313/) in your browser.

## Contributing

1. Fork this repository and create a feature branch.
2. Make your changes.
3. Open a pull request against `master`.

See [`content/en/docs/developer-guide/contributing/`](content/en/docs/developer-guide/contributing/) for detailed contribution guidelines.

## Site Structure

| Path | Purpose |
|------|---------|
| `config.toml` | Hugo site configuration (title, nav, params) |
| `content/en/` | All documentation content (Markdown) |
| `layouts/partials/` | Hugo partial overrides (navbar, head, footer, …) |
| `assets/scss/` | Custom SCSS (variables, styles, dark theme) |
| `static/` | Static assets served at the root |
| `themes/docsy` | Docsy theme (git submodule — do not edit directly) |

### Navigation

- **Top bar** — defined by `[[menu.main]]` entries in `config.toml`
- **Left sidebar** — driven by the directory structure under `content/en/docs/`; page order is controlled by the `weight` front-matter field

### Customising Styles

- Override theme variables in [`assets/scss/_variables_project.scss`](assets/scss/_variables_project.scss)
- Add custom styles in [`assets/scss/_styles_project.scss`](assets/scss/_styles_project.scss)
- The site uses a dark-only theme; the active SCSS entry point is [`assets/scss/main.scss`](assets/scss/main.scss)

### Docsy Theme

The Docsy theme is pinned as a git submodule in `themes/docsy`. **Do not edit files inside that directory** — they belong to the upstream [google/docsy](https://github.com/google/docsy) repository.

To update the pinned Docsy version:

```bash
git -C themes/docsy fetch --tags
git -C themes/docsy checkout tags/v0.6.0
# WARNING: updating Docsy may require updating overrides in
#          layouts/partials/, assets/scss/, and the [[module.mounts]]
#          in config.toml
```

## Deployment

The site is deployed automatically to GitHub Pages on every push to `master` via the workflow defined in [`.github/workflows/deploy.yml`](.github/workflows/deploy.yml).

To enable GitHub Pages for this repository, go to **Settings → Pages → Source** and select **GitHub Actions**.

## License

This documentation website is licensed under the [GNU General Public License v3.0](LICENSE).

multiPhysicsFoam is built on [OpenFOAM](https://openfoam.org/), which is distributed under the GNU General Public License v3.0.
