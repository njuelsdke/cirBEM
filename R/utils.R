# Internal utility functions

.invlogit <- function(x) {
    1 / (1 + exp(-x))
}

#' Squeeze beta values away from 0 and 1
#'
#' @param x Numeric vector of beta values.
#' @param bounds Length-2 numeric vector with lower and upper bounds.
#' @return Numeric vector with values clamped to `(bounds[1], bounds[2])`.
#' @noRd
.squeezeValues <- function(x, bounds = c(1e-5, 1 - 1e-5)) {
    pmin(pmax(x, bounds[1]), bounds[2])
}

#' Validate a SummarizedExperiment for cirBEM
#'
#' @param se A SummarizedExperiment.
#' @param assay_name Character or NULL, name of the assay to use.
#' @param time_col Character, column name in colData for time.
#' @param subject_col Character or NULL, column name in colData for subject.
#' @return Invisible NULL. Throws errors on invalid input.
#' @noRd
.validateSE <- function(se, assay_name, time_col, subject_col) {
    if (!inherits(se, "SummarizedExperiment")) {
        stop("'se' must be a SummarizedExperiment object.", call. = FALSE)
    }
    assay_name <- .resolveAssayName(se, assay_name)
    cd <- as.data.frame(colData(se))
    if (!time_col %in% names(cd)) {
        stop(sprintf("Column '%s' not found in colData(se).", time_col),
             call. = FALSE)
    }
    if (!is.numeric(cd[[time_col]])) {
        stop(sprintf("Column '%s' in colData must be numeric (hours).",
                     time_col), call. = FALSE)
    }
    if (!is.null(subject_col)) {
        if (!subject_col %in% names(cd)) {
            stop(sprintf("Column '%s' not found in colData(se).",
                         subject_col), call. = FALSE)
        }
    }
    assay_name
}

#' Resolve fixed versus mixed model mode
#' @noRd
.resolveFitMode <- function(mode, subject_col) {
    if (is.null(mode)) {
        return(if (!is.null(subject_col)) "mixed" else "fixed")
    }
    if (!is.character(mode) || length(mode) != 1L ||
        is.na(mode) || !nzchar(mode)) {
        stop("'mode' must be NULL, 'mixed', or 'fixed'.", call. = FALSE)
    }
    mode <- match.arg(mode, c("mixed", "fixed"))
    if (mode == "mixed" && is.null(subject_col)) {
        stop("'subject_col' must be supplied when mode = 'mixed'.",
             call. = FALSE)
    }
    mode
}

#' Resolve a requested assay name
#' @noRd
.resolveAssayName <- function(se, assay_name) {
    an <- assayNames(se)
    if (is.null(assay_name)) {
        if (length(an) == 0L) {
            stop("'se' must contain at least one assay.", call. = FALSE)
        }
        if (!is.na(an[1]) && nzchar(an[1])) {
            return(an[1])
        }
        return(NULL)
    }
    if (!is.character(assay_name) || length(assay_name) != 1L ||
        is.na(assay_name) || !nzchar(assay_name)) {
        stop("'assay_name' must be a single assay name or NULL.",
             call. = FALSE)
    }
    if (!assay_name %in% an) {
        stop(sprintf("Assay '%s' not found in se. Available: %s",
                     assay_name, paste(an, collapse = ", ")),
             call. = FALSE)
    }
    assay_name
}

#' Extract an assay by name, or the default assay when name is NULL
#' @noRd
.getAssay <- function(se, assay_name) {
    if (is.null(assay_name)) {
        assay(se)
    } else {
        assay(se, assay_name)
    }
}
