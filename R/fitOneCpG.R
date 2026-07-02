#' Fit a single CpG with the compiled TMB beta-regression cosinor model
#'
#' Internal workhorse function called per CpG by [fitCosinor()]. Fits the full
#' (sin + cos) and null (intercept) models with the package's compiled TMB
#' template (`beta_glmm` for mixed / repeated-measures designs, `beta_fe` for
#' fixed / cross-sectional designs) and derives the likelihood-ratio test.
#' Replaces the former glmmTMB backend.
#'
#' @param cpg_id Character, the CpG identifier.
#' @param y Numeric vector of beta values for this CpG.
#' @param df A data.frame with columns: `sin_term`, `cos_term`, and optionally
#'   `subject`. Must be aligned with `y`.
#' @param mode Character, either `"fixed"` or `"mixed"`.
#' @param prior_phi Named list with `mean` and `sd` for the normal prior on
#'   log(phi), or NULL for no prior (pilot pass / EB disabled).
#' @param eps Numeric, squeeze bound for beta values.
#' @param need_lrt Logical, whether to fit the null model and compute the LRT.
#' @param compute_se Logical, whether to compute Hessian-based coefficient SEs.
#'
#' @return A named list with elements: `cpg`, `n`, `b0`, `b_sin`, `b_cos`,
#'   `se_b0`, `se_sin`, `se_cos`, `phi`, `log_phi`, `tau_id`, `lrt_chisq`,
#'   `lrt_df`, `lrt_p`, `converged`, `error`.
#' @noRd
.fitOneCpG <- function(cpg_id, y, df, mode = "fixed",
                       prior_phi = NULL, eps = 1e-4,
                       need_lrt = TRUE, compute_se = FALSE) {

    na_result <- list(
        cpg = cpg_id, n = NA_integer_,
        b0 = NA_real_, b_sin = NA_real_, b_cos = NA_real_,
        se_b0 = NA_real_, se_sin = NA_real_, se_cos = NA_real_,
        phi = NA_real_, log_phi = NA_real_, tau_id = NA_real_,
        lrt_chisq = NA_real_, lrt_df = NA_real_, lrt_p = NA_real_,
        converged = FALSE, error = NA_character_
    )

    # Prepare data
    df$beta_val <- y
    keep <- !is.na(df$beta_val) & is.finite(df$beta_val)
    df <- df[keep, , drop = FALSE]

    df$beta_val <- .squeezeValues(df$beta_val, c(eps, 1 - eps))
    yv <- df$beta_val

    # Dispersion prior (EB shrinkage). use_prior = 0 -> plain MLE.
    use_prior <- if (is.null(prior_phi)) 0L else 1L
    pm <- if (use_prior) prior_phi$mean else 0
    ps <- if (use_prior) prior_phi$sd   else 1

    # Design matrices: full (intercept + sin + cos) and null (intercept)
    X1 <- cbind(1, df$sin_term, df$cos_term)
    X0 <- matrix(1, nrow(df), 1)

    # Moment-based starting value for log(phi)
    mu0 <- mean(yv); v0 <- stats::var(yv)
    phi0 <- if (is.finite(v0) && v0 > 0) mu0 * (1 - mu0) / v0 - 1 else 10
    lphi_init <- if (use_prior) pm else log(max(phi0, 1))

    if (mode == "mixed") {
        id0 <- as.integer(factor(df$subject)) - 1L
        nid <- max(id0) + 1L
    }

    .fitTMB <- function(X) {
        if (mode == "mixed") {
            dat <- list(model = "beta_glmm", y = yv, X = X, id = id0,
                        prior_mean = pm, prior_sd = ps, use_prior = use_prior)
            par <- list(beta = c(stats::qlogis(mu0), rep(0, ncol(X) - 1)),
                        log_phi = lphi_init, log_tau = log(0.1),
                        u = rep(0, nid))
            obj <- TMB::MakeADFun(dat, par, random = "u",
                                  DLL = "cirBEM", silent = TRUE)
        } else {
            dat <- list(model = "beta_fe", y = yv, X = X,
                        prior_mean = pm, prior_sd = ps, use_prior = use_prior)
            par <- list(beta = c(stats::qlogis(mu0), rep(0, ncol(X) - 1)),
                        log_phi = lphi_init)
            obj <- TMB::MakeADFun(dat, par, DLL = "cirBEM", silent = TRUE)
        }
        opt <- stats::nlminb(obj$par, obj$fn, obj$gr)
        list(obj = obj, opt = opt, ll = -opt$objective)
    }

    tryCatch({
        m1 <- .fitTMB(X1)
        m0 <- NULL

        # Likelihood-ratio test (full vs null), df = 2. The EB pilot pass
        # and Wald-only fits skip this second optimization.
        if (isTRUE(need_lrt)) {
            m0 <- .fitTMB(X0)
            lrt_chisq <- 2 * (m1$ll - m0$ll)
            lrt_p <- stats::pchisq(lrt_chisq, df = 2, lower.tail = FALSE)
        } else {
            lrt_chisq <- NA_real_
            lrt_p <- NA_real_
        }

        # Fixed-effect coefficients (X1 columns: intercept, sin, cos)
        bidx <- which(names(m1$opt$par) == "beta")
        bhat <- m1$opt$par[bidx]

        # Standard errors via inverse Hessian of the fixed effects (Wald).
        se <- if (isTRUE(compute_se)) {
            tryCatch({
                H <- stats::optimHess(m1$opt$par, m1$obj$fn, m1$obj$gr)
                sqrt(diag(solve(H)))[bidx]
            }, error = function(e) rep(NA_real_, length(bidx)))
        } else {
            rep(NA_real_, length(bidx))
        }

        # Dispersion and random-intercept SD
        log_phi <- unname(m1$opt$par[names(m1$opt$par) == "log_phi"])
        phi_hat <- exp(log_phi)
        tau_id <- NA_real_
        if (mode == "mixed") {
            tau_id <- exp(unname(m1$opt$par[names(m1$opt$par) == "log_tau"]))
        }

        converged <- m1$opt$convergence == 0 &&
            (!isTRUE(need_lrt) || m0$opt$convergence == 0)

        list(
            cpg = cpg_id, n = nrow(df),
            b0 = unname(bhat[1]), b_sin = unname(bhat[2]), b_cos = unname(bhat[3]),
            se_b0 = unname(se[1]), se_sin = unname(se[2]), se_cos = unname(se[3]),
            phi = phi_hat, log_phi = log_phi,
            tau_id = tau_id,
            lrt_chisq = lrt_chisq, lrt_df = 2, lrt_p = lrt_p,
            converged = converged, error = NA_character_
        )
    }, error = function(e) {
        na_result$n <- nrow(df)
        na_result$error <- conditionMessage(e)
        na_result
    })
}
