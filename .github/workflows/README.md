# CI Workflows

This repository uses GitHub Actions for package tests, documentation, and Julia
package maintenance.

## Main Pipeline

`CI.yml` runs on pushes to `main`, tags, pull requests, and manual dispatches.
It is split into:

- `detect-changes`: classifies PR changes so docs-only work can skip the test
  matrix.
- `test`: runs the main package test suite on Julia 1.10, current stable Julia,
  nightly, macOS, and Windows. Each matrix entry uploads coverage to Codecov.
  The nightly entry disables compiled modules during the test step to keep
  runtime compatibility coverage while upstream precompilation issues settle.
- `docs`: runs doctests, builds Documenter HTML output, and uploads
  `docs/build/` as an artifact without deploying.
- `docs-deploy`: deploys documentation only for pushes to `main` and tags.

The local `just ci` command mirrors the core local checks: package tests,
doctests, and documentation build.

## Benchmarks

`Benchmark.yml` runs the benchmark suite (`benchmark/benchmarks.jl`, the global
`SUITE`) via [AirspeedVelocity.jl](https://github.com/MilesCranmer/AirspeedVelocity.jl)
on pull requests that touch `src/`, `benchmark/`, or `Project.toml`. It
benchmarks the PR head commit against `main` — both revisions are run
back-to-back on the same runner (not against a stored baseline) — and writes a
`ratio` comparison table to the workflow run's **Summary** tab. It uses the
`pull_request` event with a read-only token and `job-summary: "true"`, so fork
PR code never runs with a write token and nothing is posted as a comment.
Shared CI runners are noisy, so treat the results as a guard for large
regressions rather than sub-10% changes. Run the same suite locally with
`just bench`.

## Maintenance

- `CompatHelper.yml` opens compatibility-bound update pull requests.
- `TagBot.yml` creates release tags after Julia registry publication.
- `dependabot.yml` keeps GitHub Actions dependencies current.
