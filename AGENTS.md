# Repository Guidelines

## Project Structure & Module Organization

PsyR is an R package for confidence intervals on ANOVA contrasts.

- `R/` contains package functions; keep one primary function or closely related implementation per file.
- `tests/testthat/` contains automated tests, with `tests/testthat.R` as the package test entry point.
- `man/` contains generated `.Rd` documentation and figures. Edit roxygen comments in `R/`, not generated files.
- `data/` contains package datasets (`.rda`), while `vignettes/` contains analysis tutorials and supporting references.
- `README.Rmd` is the source for `README.md`; update the R Markdown source when changing the README.

## Build, Test, and Development Commands

Run these from the repository root in R or RStudio:

```r
devtools::load_all()       # Load the package during development
devtools::test()           # Run the testthat suite
devtools::check()          # Run R package checks, examples, and documentation checks
devtools::document()       # Regenerate NAMESPACE and man/ from roxygen comments
devtools::build_vignettes()# Build the tutorial vignettes
```

Install dependencies listed in `DESCRIPTION` before testing. For a release-style build, use `devtools::build()` and install the resulting source package.

## Coding Style & Naming Conventions

Use standard R formatting: four-space indentation, snake_case names, and explicit namespaces for imported functions when clarity matters (for example, `stats::qt`). Use `<-` for assignment, keep functions focused, and validate user-facing arguments early. Add or update roxygen2 comments for exported functions, parameters, return values, and examples. No repository-wide formatter or linter is configured, so keep changes consistent with surrounding files.

## Testing Guidelines

Tests use testthat 3 (`Config/testthat/edition: 3`). Name files `test-<function>.R` and describe cases with `test_that("...", { ... })`. Add regression tests for numerical results, invalid inputs, and edge cases; use stable numeric tolerances where floating-point calculations are involved. Run `devtools::test()` locally and confirm `devtools::check()` is clean before submitting.

## Commit & Pull Request Guidelines

Existing commits use short, imperative-or-descriptive summaries such as `minor formatting tweaks` and `updated some padjusts`. Keep commits focused and explain the behavioral change. Pull requests should include a clear summary, testing performed, relevant issue or discussion link, and updated documentation, tests, README, or vignettes when behavior or public APIs change. Include rendered output or screenshots when changing vignette presentation.

## Data and Statistical Correctness

Preserve compatibility with `afex` and `emmeans` model and contrast objects. When changing confidence-interval or adjustment calculations, document the statistical rationale and compare representative results against the expected PsyR/Psy output where applicable.
