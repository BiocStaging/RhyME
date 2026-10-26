<table>
  <tr>
    <td width="170" align="center">
      <img src="man/figures/logo.png" width="150" alt="RhyME logo" />
    </td>
    <td>
      <h1>RhyME: Circadian rhythm detection for methylation array data</h1>
    </td>
  </tr>
</table>

## Introduction

`RhyME` is an R package for detecting circadian rhythmicity in
DNA methylation beta values. The package fits per-CpG beta-regression cosinor
models with sine and cosine terms for cyclic time effects, using penalized
quasi-likelihood (PQL).

The package supports fixed-effect designs and repeated-measures designs with a
subject-level random intercept. It also provides empirical Bayes shrinkage of
the beta-distribution precision parameter across CpG sites to stabilize
genome-wide inference.

## Installation

After RhyME becomes available through Bioconductor, install it with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RhyME")
```

## Usage

Input data should be a `SummarizedExperiment` with methylation beta values in an
assay and sample-level collection time in `colData()`.

```r
library(RhyME)
library(SummarizedExperiment)

example <- readRDS(system.file("extdata", "example_data.rds", package = "RhyME"))
se <- example$se

fit <- fitCosinor(
    se,
    time_col = "time",
    subject_col = "subject",
    assay_name = "beta",
    period = 24,
    prior = "auto"
)
```

If `subject_col` is supplied, `fitCosinor()` uses a mixed-effects model by
default. If `subject_col = NULL`, it uses a fixed-effects model. Optional
`covariates` are included in both the rhythmic and null models as nuisance
adjustment terms.

Rhythmicity is tested with a joint null-model score test of the sine and
cosine effects. A separate full model estimates rhythm parameters.
The default `prior = "auto"` estimates an EB prior and needs at least 1,000
usable CpGs; use `prior = "none"` for smaller datasets.

To adjust for covariates, add their columns to `colData(se)` and specify
their names, for example:

```r
fit_adjusted <- fitCosinor(
    se,
    time_col = "time",
    subject_col = "subject",
    covariates = c("sex", "age")
)
```

Results, including amplitude, phase, MESOR, p-values, and convergence, are
stored in `rowData(fit)` with the `RhyME_` prefix, such as `RhyME_pvalue` and
`RhyME_phase`. Existing annotations, including unprefixed columns such as
`pvalue`, are preserved. Refitting replaces only RhyME's result columns.
View the top CpGs by sorting the table:

```r
result <- as.data.frame(rowData(fit))
result <- result[order(result$RhyME_pvalue, na.last = NA), , drop = FALSE]
head(result, 10)
```

Fitted curves can be visualized with:

```r
plotCosinorFit(fit)
```

## Parallel Fitting

By default, `fitCosinor()` uses `BiocParallel::bpparam()`, respecting the
default backend registered with `BiocParallel::register()`. For example, to
register serial fitting for subsequent calls:

```r
BiocParallel::register(BiocParallel::SerialParam())
fit <- fitCosinor(se, time_col = "time", subject_col = "subject")
```

Pass `BPPARAM` explicitly to override the registered backend for one call.
To use eight forked workers on Linux or macOS:

```r
fit <- fitCosinor(
    se,
    time_col = "time",
    subject_col = "subject",
    BPPARAM = BiocParallel::MulticoreParam(workers = 8)
)
```

For controlled serial fitting:

```r
fit <- fitCosinor(
    se,
    time_col = "time",
    subject_col = "subject",
    BPPARAM = BiocParallel::SerialParam()
)
```

For socket-based parallel fitting, including Windows:

```r
fit <- fitCosinor(
    se,
    time_col = "time",
    subject_col = "subject",
    BPPARAM = BiocParallel::SnowParam(workers = 4, type = "SOCK")
)
```

## Main Functions

- `fitCosinor()` fits PQL beta-regression cosinor models, tests rhythmicity,
  and returns parameter estimates and rhythm summaries in `rowData()`.
- `plotCosinorFit()` plots raw beta values and fitted cosinor curves.

## Documentation

- [Package vignette](vignettes/RhyME.Rmd)
