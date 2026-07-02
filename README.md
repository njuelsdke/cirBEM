<table>
  <tr>
    <td width="170" align="center">
      <img src="man/figures/logo.png" width="150" alt="cirBEM logo" />
    </td>
    <td>
      <h1>CirBEM: Circadian rhythm detection for methylation array data</h1>
    </td>
  </tr>
</table>

## Introduction

`cirBEM` is an R package for detecting circadian rhythmicity in
DNA methylation beta values. The package fits per-CpG beta-regression cosinor
models with sine and cosine terms for cyclic time effects.

The package supports fixed-effect designs and repeated-measures designs with a
subject-level random intercept. It also provides empirical Bayes shrinkage of
the beta-distribution precision parameter across CpG sites to stabilize
genome-wide inference.

The main output is a `SummarizedExperiment` with fitted per-CpG results appended
to `rowData()`, including coefficient estimates, response-scale amplitude,
phase, MESOR, test statistics, p-values, adjusted p-values, and convergence
status.

## Installation

Install the package from source:

```r
install.packages("remotes")
remotes::install_github("njuelsdke/cirBEM")
```


## Usage

Input data should be a `SummarizedExperiment` with methylation beta values in an
assay and sample-level collection time in `colData()`.

```r
library(cirBEM)
library(SummarizedExperiment)

fit <- fitCosinor(
    se,
    time_col = "time",
    subject_col = "subject"
)
```

If `subject_col` is supplied, `fitCosinor()` uses a mixed-effects model by
default. If `subject_col = NULL`, it uses a fixed-effects model. The model type
can also be set explicitly:

```r
fit <- fitCosinor(
    se,
    time_col = "time",
    subject_col = "subject",
    mode = "mixed"
)
```

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
