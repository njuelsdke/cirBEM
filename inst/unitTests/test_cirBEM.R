make_test_se <- function(n_cpg = 50, n_subjects = 5, n_timepoints = 6,
                         n_rhythmic = 10, seed = 42) {
    set.seed(seed)

    times <- seq(0, 20, length.out = n_timepoints)
    n_samples <- n_subjects * n_timepoints
    omega <- 2 * pi / 24

    meta <- data.frame(
        time = rep(times, n_subjects),
        subject = rep(paste0("S", seq_len(n_subjects)), each = n_timepoints),
        stringsAsFactors = FALSE
    )
    meta$sample_id <- paste0(meta$subject, "_T", meta$time)

    sin_t <- sin(omega * meta$time)
    cos_t <- cos(omega * meta$time)

    beta_mat <- matrix(NA_real_, nrow = n_cpg, ncol = n_samples)
    rownames(beta_mat) <- paste0("cg", sprintf("%04d", seq_len(n_cpg)))
    colnames(beta_mat) <- meta$sample_id

    for (i in seq_len(n_cpg)) {
        b0 <- rnorm(1, mean = 0, sd = 0.5)
        if (i <= n_rhythmic) {
            b_sin <- rnorm(1, mean = 0.3, sd = 0.1)
            b_cos <- rnorm(1, mean = 0.3, sd = 0.1)
        } else {
            b_sin <- 0
            b_cos <- 0
        }
        tau <- 0.2
        phi <- exp(rnorm(1, mean = 3, sd = 0.3))

        u <- rnorm(n_subjects, 0, tau)
        eta <- b0 + b_sin * sin_t + b_cos * cos_t +
            rep(u, each = n_timepoints)
        mu <- plogis(eta)

        y <- rbeta(n_samples, mu * phi, (1 - mu) * phi)
        beta_mat[i, ] <- pmin(pmax(y, 1e-4), 1 - 1e-4)
    }

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(beta = beta_mat),
        colData = S4Vectors::DataFrame(meta)
    )
    list(se = se)
}

make_test_se_fixed <- function(n_cpg = 30, n_samples = 20,
                               n_rhythmic = 5, seed = 123) {
    set.seed(seed)

    times <- runif(n_samples, 0, 24)
    omega <- 2 * pi / 24

    meta <- data.frame(
        time = times,
        sample_id = paste0("S", seq_len(n_samples)),
        stringsAsFactors = FALSE
    )

    sin_t <- sin(omega * times)
    cos_t <- cos(omega * times)

    beta_mat <- matrix(NA_real_, nrow = n_cpg, ncol = n_samples)
    rownames(beta_mat) <- paste0("cg", sprintf("%04d", seq_len(n_cpg)))
    colnames(beta_mat) <- meta$sample_id

    for (i in seq_len(n_cpg)) {
        b0 <- rnorm(1, mean = 0, sd = 0.5)
        if (i <= n_rhythmic) {
            b_sin <- rnorm(1, mean = 0.4, sd = 0.1)
            b_cos <- rnorm(1, mean = 0.4, sd = 0.1)
        } else {
            b_sin <- 0
            b_cos <- 0
        }
        phi <- exp(rnorm(1, mean = 3, sd = 0.3))

        eta <- b0 + b_sin * sin_t + b_cos * cos_t
        mu <- plogis(eta)

        y <- rbeta(n_samples, mu * phi, (1 - mu) * phi)
        beta_mat[i, ] <- pmin(pmax(y, 1e-4), 1 - 1e-4)
    }

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(beta = beta_mat),
        colData = S4Vectors::DataFrame(meta)
    )
    list(se = se)
}

test_cirBEM_mixed_workflow <- function() {
    dat <- make_test_se(n_cpg = 60, n_subjects = 5, n_timepoints = 6,
                        n_rhythmic = 10, seed = 42)

    set.seed(1)
    se <- fitCosinor(
        dat$se,
        time_col = "time",
        subject_col = "subject",
        mode = "mixed",
        eb_shrink = TRUE,
        eb_pilot_n = 50,
        min_cpgs = 30,
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    )

    res <- SummarizedExperiment::rowData(se)
    meta <- S4Vectors::metadata(se)$cirBEM

    checkTrue(inherits(se, "SummarizedExperiment"))
    checkEquals(nrow(res), 60)
    checkEquals(meta$design$mode, "mixed")
    checkEquals(meta$ebPrior$n_pilot, 50)
    checkEquals(meta$ebPrior$min_cpgs, 30)
    checkEquals(meta$fitSummary$total, 60)
    checkEquals(meta$fitSummary$converged, sum(res$converged %in% TRUE))
    checkEquals(meta$fitSummary$failed, 60 - sum(res$converged %in% TRUE))
    checkTrue(!"relative_amplitude" %in% names(res))
    hidden_cols <- c("phi", "log_phi", "log_phi_pilot", "tau_id", "error",
                     "amplitude_link", "mesor_link",
                     "se_b0", "se_sin", "se_cos")
    checkTrue(!any(hidden_cols %in% names(res)))

    top <- topCpGs(se, n = 5)
    checkTrue(inherits(top, "DataFrame"))
    checkTrue(nrow(top) <= 5)
    checkEquals(names(top), c("cpg", "amplitude_response", "phase",
                              "mesor_response", "pvalue", "padj"))
    if (nrow(top) > 1) {
        checkTrue(all(diff(top$pvalue) >= 0, na.rm = TRUE))
    }

    p <- plotCosinorFit(se, cpgs = top$cpg[1], n_grid = 20)
    checkTrue(inherits(p, "ggplot"))
}

test_cirBEM_fixed_wald <- function() {
    dat <- make_test_se_fixed(n_cpg = 25, n_samples = 20,
                              n_rhythmic = 5, seed = 123)

    se <- fitCosinor(
        dat$se,
        time_col = "time",
        subject_col = NULL,
        mode = "fixed",
        eb_shrink = FALSE,
        test_method = "Wald",
        p.adjust.method = "bonferroni",
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    )

    res <- SummarizedExperiment::rowData(se)
    meta <- S4Vectors::metadata(se)$cirBEM

    checkEquals(meta$design$mode, "fixed")
    checkEquals(meta$test$method, "Wald")
    checkEquals(meta$test$p.adjust.method, "bonferroni")
    checkTrue(all(c("se_b0", "se_sin", "se_cos") %in% names(res)))
    checkTrue(all(res$df == 2, na.rm = TRUE))
    checkTrue(all(res$pvalue >= 0 & res$pvalue <= 1, na.rm = TRUE))
}

test_cirBEM_default_assay <- function() {
    dat <- make_test_se_fixed(n_cpg = 20, n_samples = 18,
                              n_rhythmic = 4, seed = 99)
    SummarizedExperiment::assayNames(dat$se) <- "methylation"

    se <- fitCosinor(
        dat$se,
        time_col = "time",
        eb_shrink = FALSE,
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    )

    meta <- S4Vectors::metadata(se)$cirBEM
    checkEquals(meta$design$assay_name, "methylation")
    checkEquals(meta$design$mode, "fixed")
}

test_cirBEM_mode_resolution <- function() {
    dat <- make_test_se(n_cpg = 12, n_subjects = 4, n_timepoints = 4,
                        n_rhythmic = 3, seed = 7)

    fixed <- fitCosinor(
        dat$se,
        time_col = "time",
        subject_col = "subject",
        mode = "fixed",
        eb_shrink = FALSE,
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    )
    checkEquals(S4Vectors::metadata(fixed)$cirBEM$design$mode, "fixed")
    checkTrue(is.null(S4Vectors::metadata(fixed)$cirBEM$design$subject_col))

    checkException(fitCosinor(
        dat$se,
        time_col = "time",
        subject_col = NULL,
        mode = "mixed",
        eb_shrink = FALSE,
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    ))
}

test_cirBEM_internal_summaries <- function() {
    terms <- cirBEM:::.makeCosinorTerms(c(0, 6, 12), period = 24)
    checkEquals(terms$sin_term, c(0, 1, 0), tolerance = 1e-10)
    checkEquals(terms$cos_term, c(1, 0, -1), tolerance = 1e-10)

    params <- cirBEM:::.computeCircadianParams(
        b0 = 0, b_sin = 0.3, b_cos = 0.4, period = 24
    )
    checkEquals(names(params), c("amplitude_link", "amplitude_response",
                                 "phase", "mesor_link", "mesor_response"))
    checkEquals(params$amplitude_link, 0.5, tolerance = 1e-10)

    marginal <- cirBEM:::.computeCircadianParams(
        b0 = -1, b_sin = 0.5, b_cos = 0.3, period = 24, tau_id = 1
    )
    checkTrue(!isTRUE(all.equal(params$mesor_response,
                                marginal$mesor_response)))

    set.seed(42)
    log_phi <- rnorm(500, mean = 3.0, sd = 0.5)
    prior <- cirBEM:::ebShrinkPhi(log_phi, trim = 0)
    checkEquals(prior$mean, 3.0, tolerance = 0.15)
    checkEquals(prior$sd, 0.5, tolerance = 0.15)

    checkException(cirBEM:::ebShrinkPhi(rep(NA_real_, 10)))
}
