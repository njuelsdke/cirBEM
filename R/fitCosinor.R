#' Fit CirBEM across CpG sites
#'
#' Fits a CirBEM to each CpG site in a methylation array
#' dataset. Supports both fixed-effects (cross-sectional) and mixed-effects
#' (repeated measures with random intercepts) designs, and optionally
#' applies empirical Bayes shrinkage to the precision parameter (phi).
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment] object. The
#'   specified assay should contain methylation beta values in (0, 1).
#' @param assay_name Character or NULL, name of the assay to use. The default
#'   NULL uses the first assay in `se`.
#' @param time_col Character, column name in `colData(se)` containing
#'   collection times in hours (default `"time"`).
#' @param subject_col Character or NULL. Column name in `colData(se)` for
#'   subject identifiers. If non-NULL, random intercepts per subject are
#'   included (mixed-effects mode). If NULL, fixed-effects mode is used.
#' @param mode Character or NULL. Set to `"mixed"` to fit random subject
#'   intercepts or `"fixed"` to fit a fixed-effects model. The default `NULL`
#'   uses `"mixed"` when `subject_col` is supplied and `"fixed"` otherwise.
#' @param period Numeric, the oscillation period in hours (default 24 for
#'   circadian rhythms). Can be set to other values for ultradian or
#'   infradian rhythms.
#' @param eb_shrink Logical, whether to apply empirical Bayes shrinkage to
#'   the precision parameter (default TRUE). This performs a two-pass fit:
#'   a pilot fit to estimate the prior, then a re-fit with the prior.
#' @param eb_trim Numeric, fraction to trim from each tail when estimating
#'   the EB prior (default 0, i.e. no trimming; the prior uses the full pilot
#'   log(phi) distribution). Increase (e.g. 0.02) to make the prior robust to
#'   pilot outliers. Set `eb_shrink = FALSE` to disable the prior entirely.
#' @param eb_pilot_n Integer, maximum number of CpGs to use for estimating the
#'   EB prior (default 10000). If the dataset has fewer CpGs, all CpGs are used.
#'   Random sampling uses R's current RNG state; call [set.seed()] before
#'   [fitCosinor()] for reproducible pilot selection.
#' @param min_cpgs Integer, minimum number of valid pilot CpGs required to
#'   estimate the EB prior (default 100). Lower this for small example or pilot
#'   datasets; increase it for large genome-wide analyses.
#' @param squeeze_bounds Numeric vector of length 2. Beta values are clamped
#'   to this open interval to avoid boundary issues (default `c(1e-4, 1-1e-4)`).
#' @param test_method Character, rhythmicity test to report. `"LRT"` (default)
#'   uses the likelihood ratio test computed during model fitting; `"Wald"`
#'   uses a joint Wald test that `b_sin = b_cos = 0`.
#' @param p.adjust.method Character, p-value adjustment method passed to
#'   [stats::p.adjust()] (default `"BH"`).
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object for execution.
#'   Defaults to [BiocParallel::MulticoreParam()] for parallel fitting on
#'   Unix-like systems. Pass [BiocParallel::SerialParam()] for serial fitting or
#'   [BiocParallel::SnowParam()] for socket-based parallelism.
#' @param verbose Logical, print progress messages (default TRUE).
#' @param progressbar Logical, show BiocParallel progress bars when supported
#'   by `BPPARAM` (default follows `verbose`).
#'
#' @return The input `SummarizedExperiment` with per-CpG results appended to
#'   `rowData(se)` and design info, EB prior, test settings, and fit summary
#'   stored in `metadata(se)$cirBEM`.
#'
#' @export
#'
#' @examples
#' ex <- readRDS(system.file("extdata", "example_data.rds",
#'                          package = "cirBEM"))
#' se <- ex$se[seq_len(5), ]
#'
#' fit <- fitCosinor(se, time_col = "time", subject_col = "subject",
#'                   eb_shrink = FALSE, verbose = FALSE)
#' head(SummarizedExperiment::rowData(fit))
fitCosinor <- function(se,
                       assay_name = NULL,
                       time_col = "time",
                       subject_col = NULL,
                       mode = NULL,
                       period = 24,
                       eb_shrink = TRUE,
                       eb_trim = 0,
                       eb_pilot_n = 10000L,
                       min_cpgs = 100L,
                       squeeze_bounds = c(1e-4, 1 - 1e-4),
                       test_method = c("LRT", "Wald"),
                       p.adjust.method = "BH",
                       BPPARAM = BiocParallel::MulticoreParam(),
                       verbose = TRUE,
                       progressbar = verbose) {

    test_method <- match.arg(test_method)
    compute_se <- test_method == "Wald"
    BPPARAM <- .setBPProgressBar(BPPARAM, progressbar)
    if (eb_shrink) {
        if (!is.numeric(eb_pilot_n) || length(eb_pilot_n) != 1L ||
            !is.finite(eb_pilot_n) || eb_pilot_n < 1) {
            stop("'eb_pilot_n' must be a positive integer.", call. = FALSE)
        }
        eb_pilot_n <- as.integer(eb_pilot_n)
        if (!is.numeric(min_cpgs) || length(min_cpgs) != 1L ||
            !is.finite(min_cpgs) || min_cpgs < 1) {
            stop("'min_cpgs' must be a positive integer.", call. = FALSE)
        }
        min_cpgs <- as.integer(min_cpgs)
    }

    # Resolve model type from the supplied design unless explicitly requested.
    mode <- .resolveFitMode(mode, subject_col)
    subject_col_model <- if (mode == "mixed") subject_col else NULL

    # Validate input
    assay_name <- .validateSE(se, assay_name, time_col, subject_col_model)

    # Extract data
    beta_mat <- .getAssay(se, assay_name)
    cd <- as.data.frame(colData(se))
    cpg_ids <- rownames(beta_mat)
    n_cpg <- length(cpg_ids)

    if (verbose) {
        message(sprintf("cirBEM: %d CpGs, %d samples, mode=%s, period=%g",
                        n_cpg, ncol(beta_mat), mode, period))
    }

    # Build per-sample data frame with cosinor terms
    cos_terms <- .makeCosinorTerms(cd[[time_col]], period = period)
    sample_df <- data.frame(
        sin_term = cos_terms$sin_term,
        cos_term = cos_terms$cos_term
    )
    if (mode == "mixed") {
        sample_df$subject <- cd[[subject_col_model]]
    }

    eb_prior <- list()
    pilot_log_phi <- rep(NA_real_, n_cpg)

    if (eb_shrink) {
        # ---- Pass 1: Pilot fit on a random subset (no prior) ----
        n_pilot <- min(eb_pilot_n, n_cpg)
        if (n_pilot < min_cpgs) {
            stop(sprintf(
                paste0("EB prior needs at least %d pilot CpGs, but only %d ",
                       "are available. Lower 'min_cpgs' or set ",
                       "'eb_shrink = FALSE'."),
                min_cpgs, n_pilot
            ), call. = FALSE)
        }
        pilot_idx <- if (n_pilot < n_cpg) {
            sample(seq_len(n_cpg), n_pilot)
        } else {
            seq_len(n_cpg)
        }

        if (verbose) {
            message(sprintf("Pass 1: Pilot fit on %d CpGs for EB prior...",
                            n_pilot))
        }

        pilot_results <- BiocParallel::bplapply(
            pilot_idx,
            function(i) {
                .fitOneCpG(
                    cpg_id = cpg_ids[i],
                    y = as.numeric(beta_mat[i, ]),
                    df = sample_df,
                    mode = mode,
                    prior_phi = NULL,
                    eps = squeeze_bounds[1],
                    need_lrt = FALSE,
                    compute_se = FALSE
                )
            },
            BPPARAM = BPPARAM
        )

        pilot_log_phi_subset <- vapply(pilot_results, function(r) r$log_phi,
                                       numeric(1))
        pilot_log_phi[pilot_idx] <- pilot_log_phi_subset

        if (verbose) message("Estimating EB prior on log(phi)...")
        eb_prior <- ebShrinkPhi(pilot_log_phi_subset, trim = eb_trim,
                                min_success = min_cpgs)
        eb_prior$n_pilot <- n_pilot
        eb_prior$min_cpgs <- min_cpgs
        eb_prior$pilot_fraction <- n_pilot / n_cpg
        if (verbose) {
            message(sprintf(
                "  EB prior: mean=%.3f, sd=%.3f (n_pilot=%d, n_used=%d)",
                eb_prior$mean, eb_prior$sd, n_pilot, eb_prior$n_used))
        }

        if (verbose) message("Pass 2: Fitting with EB prior...")

        prior_spec <- list(mean = eb_prior$mean, sd = eb_prior$sd)

        final_results <- BiocParallel::bplapply(
            seq_len(n_cpg),
            function(i) {
                .fitOneCpG(
                    cpg_id = cpg_ids[i],
                    y = as.numeric(beta_mat[i, ]),
                    df = sample_df,
                    mode = mode,
                    prior_phi = prior_spec,
                    eps = squeeze_bounds[1],
                    need_lrt = test_method == "LRT",
                    compute_se = compute_se
                )
            },
            BPPARAM = BPPARAM
        )
    } else {
        if (verbose) message("Fitting...")

        final_results <- BiocParallel::bplapply(
            seq_len(n_cpg),
            function(i) {
                .fitOneCpG(
                    cpg_id = cpg_ids[i],
                    y = as.numeric(beta_mat[i, ]),
                    df = sample_df,
                    mode = mode,
                    prior_phi = NULL,
                    eps = squeeze_bounds[1],
                    need_lrt = test_method == "LRT",
                    compute_se = compute_se
                )
            },
            BPPARAM = BPPARAM
        )
    }

    # ---- Assemble results into the SE ----
    if (verbose) message("Assembling results...")

    se <- .assembleResults(se, final_results, eb_prior,
                           eb_shrink, mode, period, assay_name, time_col,
                           subject_col_model, ncol(beta_mat),
                           if (mode == "mixed") {
                               length(unique(cd[[subject_col_model]]))
                           } else {
                               NA_integer_
                           },
                           test_method, p.adjust.method)

    if (verbose) .messageFitSummary(se)
    se
}

#' Assemble per-CpG results into the SummarizedExperiment
#' @noRd
.assembleResults <- function(se, results, eb_prior, eb_shrink,
                             mode, period, assay_name, time_col, subject_col,
                             n_samples, n_subjects, test_method,
                             p.adjust.method) {

    # Extract vectors from list of lists
    extract_field <- function(field, type = "numeric") {
        if (type == "numeric") {
            vapply(results, function(r) {
                v <- r[[field]]
                if (is.null(v)) NA_real_ else as.numeric(v)
            }, numeric(1))
        } else if (type == "logical") {
            vapply(results, function(r) {
                v <- r[[field]]
                if (is.null(v)) NA else as.logical(v)
            }, logical(1))
        } else {
            vapply(results, function(r) {
                v <- r[[field]]
                if (is.null(v)) NA_character_ else as.character(v)
            }, character(1))
        }
    }

    cpg_names <- extract_field("cpg", "character")

    # Coefficients
    b0     <- extract_field("b0")
    b_sin  <- extract_field("b_sin")
    b_cos  <- extract_field("b_cos")

    tau_id <- extract_field("tau_id")

    # Circadian parameters (computed from coefficients; response-scale
    # summaries are marginal over random intercepts when tau_id is available).
    circ <- .computeCircadianParams(b0, b_sin, b_cos, period, tau_id = tau_id)

    if (test_method == "LRT") {
        statistic <- extract_field("lrt_chisq")
        df <- extract_field("lrt_df")
        pvalue <- extract_field("lrt_p")
    } else {
        se_sin <- extract_field("se_sin")
        se_cos <- extract_field("se_cos")
        statistic <- (b_sin / se_sin)^2 + (b_cos / se_cos)^2
        df <- rep(2, length(statistic))
        pvalue <- stats::pchisq(statistic, df = 2, lower.tail = FALSE)
    }

    # Build the per-CpG results DataFrame. Hessian-based coefficient SEs are
    # only included when they were actually computed for Wald testing.
    res_df <- DataFrame(
        b0     = b0,
        b_sin  = b_sin,
        b_cos  = b_cos,
        amplitude_response  = circ$amplitude_response,
        phase               = circ$phase,
        mesor_response      = circ$mesor_response,
        statistic = statistic,
        df        = df,
        pvalue    = pvalue,
        padj      = stats::p.adjust(pvalue, method = p.adjust.method),
        converged = extract_field("converged", "logical"),
        row.names = cpg_names
    )
    if (test_method == "Wald") {
        res_df <- cbind(
            res_df[, c("b0", "b_sin", "b_cos"), drop = FALSE],
            DataFrame(
                se_b0  = extract_field("se_b0"),
                se_sin = extract_field("se_sin"),
                se_cos = extract_field("se_cos"),
                row.names = cpg_names
            ),
            res_df[, setdiff(names(res_df), c("b0", "b_sin", "b_cos")),
                   drop = FALSE]
        )
    }

    # Store per-CpG results as ordinary rowData columns.
    rd <- rowData(se)
    keep_cols <- setdiff(names(rd), c("cirBEM", names(res_df)))
    rowData(se) <- cbind(rd[, keep_cols, drop = FALSE], res_df)

    n_total <- nrow(res_df)
    n_converged <- sum(res_df$converged %in% TRUE)
    fit_summary <- list(
        total = n_total,
        converged = n_converged,
        failed = n_total - n_converged
    )

    # Store design info and EB prior in metadata
    metadata(se)$cirBEM <- list(
        design = list(
            period      = period,
            mode        = mode,
            assay_name  = assay_name,
            time_col    = time_col,
            subject_col = subject_col,
            n_samples   = n_samples,
            n_subjects  = n_subjects
        ),
        ebPrior = eb_prior,
        test = list(
            method = test_method,
            p.adjust.method = p.adjust.method
        ),
        fitSummary = fit_summary
    )

    se
}

#' Enable BiocParallel progress bars when supported
#' @noRd
.setBPProgressBar <- function(BPPARAM, progressbar) {
    if (!isTRUE(progressbar)) {
        return(BPPARAM)
    }
    tryCatch({
        BiocParallel::bpprogressbar(BPPARAM) <- TRUE
        BPPARAM
    }, error = function(e) BPPARAM)
}

#' Print a compact fit summary
#' @noRd
.messageFitSummary <- function(se) {
    fs <- metadata(se)$cirBEM$fitSummary
    message("cirBEM fit summary:")
    message(sprintf("  converged: %d / %d", fs$converged, fs$total))
    message(sprintf("  failed: %d", fs$failed))
    invisible(NULL)
}
