
<img src="man/figures/logo.png" align="right" width="140" />

# cirBEM

Beta-regression cosinor modeling for detecting circadian rhythmicity in DNA
methylation beta values.

`cirBEM` fits per-CpG beta regression models with sine and cosine circadian
terms. It supports fixed-effect and repeated-measures designs, optional
empirical Bayes shrinkage of the precision parameter, and returns the fitted
results in a `SummarizedExperiment`.

## Installation

Install the package from source:

```r
install.packages("cirBEM", repos = NULL, type = "source")
```

## Basic Usage

Input data should be a `SummarizedExperiment` with methylation beta values in
an assay and sample-level time information in `colData`.

```r
library(cirBEM)
library(SummarizedExperiment)

fit <- fitCosinor(
    se,
    time_col = "time",
    subject_col = "subject"
)
```

If `subject_col` is supplied, `fitCosinor()` uses a mixed model with a random
subject intercept. If `subject_col = NULL`, it uses a fixed-effect model.

Top rhythmic CpGs can be extracted with:

```r
topCpGs(fit, n = 10)
```

Fitted curves can be visualized with:

```r
plotCosinorFit(fit, cpgs = topCpGs(fit, n = 3)$cpg)
```

## Output

Results are stored in:

```r
rowData(fit)
metadata(fit)$cirBEM
```

The per-CpG result table includes coefficient estimates, response-scale
amplitude, phase, MESOR, test statistics, p-values, adjusted p-values, and
convergence status. Coefficient standard errors are included when
`test_method = "Wald"`.

## Main Functions

- `fitCosinor()` fits circadian beta-regression models across CpG sites.
- `topCpGs()` returns the top CpGs ranked by significance.
- `plotCosinorFit()` plots raw beta values and fitted cosinor curves.
