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
models with sine and cosine terms for cyclic time effects.

The package supports fixed-effect designs and repeated-measures designs with a
subject-level random intercept. It also provides empirical Bayes shrinkage of
the beta-distribution precision parameter across CpG sites to stabilize
genome-wide inference.

The automatic prior starts from a no-prior pilot fit, summarizes log precision
robustly, and uses five common-precision parametric-bootstrap datasets by
default to correct finite-sample center bias and fitting-noise variance.

The main output is a `SummarizedExperiment` with fitted per-CpG results appended
to `rowData()`, including coefficient estimates, response-scale amplitude,
phase, MESOR, p-values, adjusted p-values, and convergence
status.

## Installation

Install the development version from GitHub:

```r
remotes::install_github("njuelsdke/RhyME")
```

## Usage

Input data should be a `SummarizedExperiment` with methylation beta values in an
assay and sample-level collection time in `colData()`.

```r
library(RhyME)
library(SummarizedExperiment)

fit <- fitCosinor(
    se,
    time_col = "time",
    subject_col = "subject",
    covariates = c("sex", "age")
)
```

If `subject_col` is supplied, `fitCosinor()` uses a mixed-effects model by
default. If `subject_col = NULL`, it uses a fixed-effects model. Optional
`covariates` are included in both the rhythmic and null models as nuisance
adjustment terms; the reported p-value tests only the sine/cosine effect.

By default, rhythmicity is tested with a covariance-aware joint Wald test of
the sine and cosine coefficients. The automatic prior uses 10,000 pilot CpGs
and five bootstrap replicates:

```r
set.seed(1)
fit <- fitCosinor(
    se,
    time_col = "time",
    subject_col = "subject",
    prior = "auto",
    test = "Wald"
)
```

Use `test = "LRT"` to request a likelihood-ratio test instead. Selecting
`prior = "auto"` always applies the bootstrap correction; use `prior = "none"`
to fit without a prior.
The default bootstrap-corrected prior requires at least 1,000 finite pilot
CpG fits.

Top rhythmic CpGs can be extracted with:

```r
topCpGs(fit, n = 10)
```

Fitted curves can be visualized with:

```r
plotCosinorFit(fit, cpgs = topCpGs(fit, n = 3)$cpg)
```

## Parallel Fitting

By default, `fitCosinor()` uses `BiocParallel::MulticoreParam()` for parallel
fitting on Linux systems. For controlled serial fitting:

```r
fit <- fitCosinor(
    se,
    time_col = "time",
    subject_col = "subject",
    BPPARAM = BiocParallel::SerialParam()
)
```

For socket-based parallel fitting:

```r
fit <- fitCosinor(
    se,
    time_col = "time",
    subject_col = "subject",
    BPPARAM = BiocParallel::SnowParam(workers = 4, type = "SOCK")
)
```

## Main Functions

- `fitCosinor()` fits circadian beta-regression models across CpG sites.
- `topCpGs()` returns top CpGs ranked by statistical significance.
- `plotCosinorFit()` plots raw beta values and fitted cosinor curves.
