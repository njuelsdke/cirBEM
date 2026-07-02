#' Empirical Bayes shrinkage for beta regression precision parameter
#'
#' Estimates a Normal prior on log(phi) from a vector of per-CpG log(phi)
#' estimates. The prior is used in a second-pass TMB fit to shrink
#' CpG-specific precision estimates toward a global mean, analogous to
#' limma's moderated variance approach.
#'
#' @param log_phi_hat Numeric vector of per-CpG log(phi) estimates from the
#'   pilot (first-pass) fit.
#' @param trim Numeric, fraction of observations to trim from each tail
#'   before estimating the prior (default 0).
#' @param min_success Integer, minimum number of valid (finite) estimates
#'   required. Throws an error if fewer are available.
#'
#' @return A list with components:
#'   \describe{
#'     \item{mean}{Estimated prior mean of log(phi).}
#'     \item{sd}{Estimated prior standard deviation of log(phi).}
#'     \item{n_used}{Number of CpGs used for prior estimation (after trimming).}
#'     \item{trim}{The trim fraction that was applied.}
#'   }
#'
#' @noRd
#'
#' @examples
#' set.seed(42)
#' log_phi <- rnorm(500, mean = 3, sd = 0.5)
#' prior <- ebShrinkPhi(log_phi)
#' prior$mean  # close to 3
#' prior$sd    # close to 0.5
ebShrinkPhi <- function(log_phi_hat, trim = 0, min_success = 50) {

    # Filter to valid values
    valid <- log_phi_hat[is.finite(log_phi_hat)]
    n_ok <- length(valid)

    if (n_ok < min_success) {
        stop(sprintf(
            paste0("ebShrinkPhi: only %d valid log(phi) estimates ",
                   "(min_success = %d). Prior would be unstable."),
            n_ok, min_success
        ), call. = FALSE)
    }

    # Quantile trimming
    trim_used <- 0
    trimmed <- valid

    if (is.finite(trim) && trim > 0) {
        trim2 <- min(trim, 0.49)
        qs <- stats::quantile(valid, probs = c(trim2, 1 - trim2),
                              na.rm = TRUE, names = FALSE, type = 7)
        keep <- (valid >= qs[1]) & (valid <= qs[2])
        if (sum(keep) >= min_success) {
            trimmed <- valid[keep]
            trim_used <- trim2
        }
    }

    mu_0 <- mean(trimmed)
    sigma_0 <- stats::sd(trimmed)

    list(
        mean = mu_0,
        sd = sigma_0,
        n_used = length(trimmed),
        trim = trim_used
    )
}
