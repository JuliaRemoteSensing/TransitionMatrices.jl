# CI Workflows

This repository uses GitHub Actions for package tests, documentation, and Julia
package maintenance.

## Main Pipeline

`CI.yml` runs on pushes to `main`, tags, pull requests, and manual dispatches.
It is split into:

- `detect-changes`: classifies PR changes so docs-only work can skip the test
  matrix.
- `test`: runs the main package test suite on Julia 1.10, current stable Julia,
  macOS, and Windows. Each matrix entry uploads coverage to Codecov.
- `docs`: runs doctests, builds Documenter HTML output, and uploads
  `docs/build/` as an artifact without deploying.
- `docs-deploy`: deploys documentation only for pushes to `main` and tags.

The local `just ci` command mirrors the core local checks: package tests,
doctests, and documentation build.

## Maintenance

- `CompatHelper.yml` opens compatibility-bound update pull requests.
- `TagBot.yml` creates release tags after Julia registry publication.
- `dependabot.yml` keeps GitHub Actions dependencies current.
