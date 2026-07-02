#' Build cosinor terms from time vector
#'
#' @param time Numeric vector of collection times (in hours).
#' @param period Numeric scalar, the period in hours (default 24).
#' @return A data.frame with columns `sin_term` and `cos_term`.
#' @noRd
.makeCosinorTerms <- function(time, period = 24) {
    omega <- 2 * pi / period
    data.frame(
        sin_term = sin(omega * time),
        cos_term = cos(omega * time)
    )
}
