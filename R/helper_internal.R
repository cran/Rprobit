#' Read Model Formula
#'
#' @description
#' `read_formula()` parses a model formula into the response variable, the
#' three supported covariate groups, and the ASC flag.
#'
#' `read_formula_lhs()` returns the response variable.
#'
#' `read_formula_rhs()` returns the parsed right hand side.
#'
#' `read_formula_asc()` returns whether ASCs are included.
#'
#' @param formula \[`formula`\]\cr
#' A symbolic description of the choice model, see details.
#'
#' @section Specifying the model formula:
#' The structure of `formula` is `choice ~ A | B | C`, i.e., a standard
#' \code{\link[stats]{formula}} object but with (up to) three parts on the
#' right-hand side, separated by `|`, where
#' \itemize{
#'   \item `choice` is the name of the discrete response variable,
#'   \item `A` are names of \strong{alternative-specific covariates} with
#'   \strong{a coefficient that is constant across alternatives},
#'   \item `B` are names of \strong{covariates that are constant across
#'   alternatives},
#'   \item and `C` are names of \strong{alternative-specific covariates}
#'   with \strong{alternative-specific coefficients}.
#' }
#'
#' The following rules apply:
#' \enumerate{
#'   \item By default, intercepts (referred to as alternative-specific
#'   constants, ASCs) are added to the model. They can be removed by adding
#'   `+ 0` in the second part, e.g., `choice ~ A | B + 0 | C`. To not include
#'   any covariates of the second type but to estimate ASCs, add `1` in the
#'   second part, e.g., `choice ~ A | 1 | C`. The expression
#'   `choice ~ A | 0 | C` is interpreted as no covariates of the second type and
#'   no ASCs.
#'   \item To not include covariates of any type, add `0` in the respective
#'   part, e.g., `choice ~ 0 | B | C`.
#'   \item Some parts of the formula can be omitted when there is no ambiguity.
#'   For example, `choice ~ A` is equivalent to `choice ~ A | 1 | 0`.
#'   \item Multiple covariates in one part are separated by a `+` sign, e.g.,
#'   `choice ~ A1 + A2`.
#' }
#'
#' @return
#' For `read_formula()`, a `list` with the following elements:
#'
#' - `choice_variable`: a `character(1)`, the response variable.
#' - `covariate_types`: a `list` of length three. Its elements contain the
#'   covariate names from the three right-hand-side parts.
#' - `ASC`: a `logical(1)` that indicates whether alternative-specific
#'   constants are included.
#'
#' For `read_formula_lhs()`, a `character(1)` with the response variable.
#'
#' For `read_formula_rhs()`, a `list` of length three. Its elements contain
#' the covariate names from the three right-hand-side parts.
#'
#' For `read_formula_asc()`, a `logical(1)` that indicates whether
#' alternative-specific constants are included.
#'
#' @keywords internal

read_formula <- function(formula) {
  choice_variable <- read_formula_lhs(formula)
  covariate_types <- read_formula_rhs(formula)
  ASC <- read_formula_asc(formula)
  list(
    "choice_variable" = choice_variable,
    "covariate_types" = covariate_types,
    "ASC" = ASC
  )
}

#' @rdname read_formula
#' @keywords internal

read_formula_lhs <- function(formula) {
  if (!inherits(formula, "formula")) {
    stop("Argument `formula` must be a formula.", call. = FALSE)
  }
  if (length(formula) != 3) {
    stop("Argument `formula` must be a two-sided formula.", call. = FALSE)
  }
  lhs <- all.vars(formula[[2]])
  if (length(lhs) != 1) {
    stop(
      "Argument `formula` must have exactly one left-hand-side variable.",
      call. = FALSE
    )
  }
  lhs
}

#' @rdname read_formula
#' @keywords internal

read_formula_rhs <- function(formula) {
  read_formula_lhs(formula)
  fs <- as.character(formula)
  covariate_types <- strsplit(fs[3], split = "|", fixed = TRUE)[[1]] |>
    strsplit(split = "+", fixed = TRUE) |>
    lapply(trimws)
  var_levels <- length(covariate_types)
  if (var_levels > 3) {
    stop(
      "Argument `formula` must have at most three covariate parts.",
      call. = FALSE
    )
  }
  if (var_levels == 1) {
    covariate_types[[2]] <- character(0)
    covariate_types[[3]] <- character(0)
  }
  if (var_levels == 2) {
    covariate_types[[3]] <- character(0)
  }
  for (i in seq_along(covariate_types)) {
    keep_covariates <- !(covariate_types[[i]] %in% c("", "0", "1"))
    covariate_types[[i]] <- covariate_types[[i]][keep_covariates]
  }
  covariate_types
}

#' @rdname read_formula
#' @keywords internal

read_formula_asc <- function(formula) {
  read_formula_rhs(formula)
  fs <- as.character(formula)
  rhs_parts <- strsplit(fs[3], split = "|", fixed = TRUE)[[1]]
  if (length(rhs_parts) == 1) {
    return(TRUE)
  }
  asc_part <- strsplit(rhs_parts[[2]], split = "+", fixed = TRUE)[[1]]
  asc_part <- trimws(asc_part)
  has_no_asc <- any(asc_part %in% "0")
  has_asc <- any(asc_part %in% "1")
  if (has_no_asc && has_asc) {
    stop(
      "Terms `0` and `1` cannot both occur in the second formula part.",
      call. = FALSE
    )
  }
  !has_no_asc || has_asc
}

#' Reconstruct Cholesky Root from Half-Vectorization
#'
#' @description
#' This function reconstructs a lower triangular Cholesky root from its
#' half-vectorization, e.g., for a covariance matrix \eqn{M}, it maps
#' \eqn{\mathrm{vech}(\mathrm{chol}(M))} to \eqn{\mathrm{chol}(M)}.
#'
#' @param vech_chol \[`numeric()`\]\cr
#' Half-vectorized Cholesky root in column-wise order.
#'
#' @return
#' A lower triangular Cholesky root.
#'
#' @keywords internal

vech_chol_to_chol <- function(vech_chol) {
  size <- round(-0.5 + sqrt(2 * length(vech_chol) + 0.25))
  if (size * (size + 1) / 2 != length(vech_chol)) {
    stop("Length of `vech_chol` is not triangular.", call. = FALSE)
  }
  if (size == 0) {
    return(matrix(0, 0, 0))
  }
  if (size == 1) {
    return(matrix(vech_chol, 1, 1))
  }
  chol <- matrix(0, size, size)
  chol[lower.tri(chol, diag = TRUE)] <- vech_chol
  chol
}

#' Check Covariance Matrix
#'
#' @description
#' This function checks whether the input is a valid covariance matrix.
#'
#' @param x \[any\]\cr
#' Object to check.
#'
#' @param dim \[`NULL` | `integer(1)`\]\cr
#' The matrix dimension.
#'
#' @param tolerance \[`numeric(1)`\]\cr
#' A non-negative tolerance value.
#'
#' @param warning \[`logical(1)`\]\cr
#' If `TRUE`, warn about the failed check.
#'
#' @return
#' Either `TRUE` or `FALSE`.
#'
#' @keywords internal

check_covariance_matrix <- function(
    x, dim = NULL, tolerance = sqrt(.Machine$double.eps), warning = FALSE
) {
  problem <- NULL
  if (!is.numeric(x)) {
    problem <- "The value must be numeric."
  } else if (!is.matrix(x)) {
    problem <- "The value must be a matrix."
  } else if (nrow(x) != ncol(x)) {
    problem <- "The matrix must be square."
  } else if (any(is.na(x))) {
    problem <- "The matrix must not contain missing values."
  } else if (any(!is.finite(x))) {
    problem <- "The matrix must contain only finite values."
  } else if (any(abs(x - t(x)) > tolerance)) {
    problem <- "The matrix must be symmetric."
  } else if (nrow(x) == 0) {
    if (!is.null(dim) && nrow(x) != dim) {
      problem <- paste0("The matrix must have dimension ", dim, ".")
    }
  } else if (any(eigen(x)$value < -tolerance)) {
    problem <- "The matrix must have only positive eigenvalues."
  } else if (!is.null(dim) && nrow(x) != dim) {
    problem <- paste0("The matrix must have dimension ", dim, ".")
  }
  if (!is.null(problem)) {
    if (isTRUE(warning)) {
      warning(problem, call. = FALSE)
    }
    return(FALSE)
  }
  return(TRUE)
}

#' Adjust Covariance Matrix
#'
#' @description
#' This helper function symmetrizes a covariance matrix estimate and replaces
#' eigenvalues below a positive tolerance to obtain a numerically stable
#' inverse. A matrix containing non-finite values is replaced by an identity
#' matrix.
#'
#' @param x \[`matrix`\]\cr
#' Numeric square matrix to adjust.
#'
#' @param label \[`character(1)`\]\cr
#' Matrix label used in warning messages.
#'
#' @param tolerance \[`numeric(1)`\]\cr
#' Positive lower bound for adjusted eigenvalues.
#'
#' @return
#' A `list` with the adjusted matrix `mat` and its inverse `inverse`.
#'
#' @keywords internal

adjust_covariance_matrix <- function(
    x, label = "covariance matrix", tolerance = 1e-6
) {
  checkmate::assert_matrix(x, mode = "numeric", any.missing = TRUE)
  if (nrow(x) != ncol(x)) {
    stop("Argument `x` must be a square matrix.", call. = FALSE)
  }
  checkmate::assert_string(label)
  checkmate::assert_number(tolerance, lower = 0, finite = TRUE)
  if (tolerance == 0) {
    stop("Argument `tolerance` must be positive.", call. = FALSE)
  }
  p <- nrow(x)
  if (p == 0L) {
    return(list(mat = x, inverse = x))
  }
  if (any(!is.finite(x))) {
    warning(
      sprintf(
        paste0(
          "The fitted %s contained non-finite values and was replaced ",
          "by an identity matrix."
        ),
        label
      ),
      call. = FALSE
    )
    x <- diag(p)
  }
  x <- (x + t(x)) / 2
  eig <- eigen(x, symmetric = TRUE)
  if (p == 1L && eig$values > 0) {
    return(list(mat = x, inverse = 1 / x))
  }
  if (any(eig$values < tolerance)) {
    warning(
      sprintf(
        "The fitted %s was not positive definite and was adjusted.",
        label
      ),
      call. = FALSE
    )
  }
  values <- pmax(eig$values, tolerance)
  adjusted <- eig$vectors %*% diag(values, nrow = p) %*% t(eig$vectors)
  inverse <- eig$vectors %*% diag(1 / values, nrow = p) %*% t(eig$vectors)
  list(mat = adjusted, inverse = inverse)
}

#' Complete Long-Format Choice Data
#'
#' @description
#' This helper function cycles through the decisions and adds rows to a
#' long-format data frame, if alternatives are not available. The choice is not
#' taken and characteristics are added as NA.
#'
#' @param df \[`data.frame`\]\cr
#' Choice data in long format.
#'
#' @param column_choice \[`character(1)`\]\cr
#' Name of the choice indicator column.
#'
#' @param column_decider_id \[`character(1)` | `NULL`\]\cr
#' Name of the decider identifier column.
#'
#' @param column_occasion_id \[`character(1)` | `NULL`\]\cr
#' Name of the occasion identifier column.
#'
#' If `NULL` and a decider identifier is supplied, a column named `occasion` is
#' added and numbered within each decider.
#'
#' @param alt_names \[`character()`\]\cr
#' Alternative labels to include in every completed choice set.
#'
#' @param column_alternative \[`character(1)`\]\cr
#' Name of the alternative label column.
#'
#' @return
#' The completed `data.frame`.
#'
#' @keywords internal

complete_long_data_frame <- function(
    df, column_choice, column_decider_id, column_occasion_id, alt_names,
    column_alternative
) {
  if (nrow(df) == 0L || is.null(column_decider_id)) {
    return(df)
  }
  if (is.null(column_occasion_id)) {
    column_occasion_id <- "occasion"
    if (!column_occasion_id %in% names(df)) {
      df[[column_occasion_id]] <- integer(nrow(df))
      ids <- df[[column_decider_id]]
      index_deciders <- split(seq_along(ids), match(ids, unique(ids)))
      for (index_decider in index_deciders) {
        n_rows <- length(index_decider)
        n_alt <- length(alt_names)
        if (n_rows <= n_alt) {
          occasions <- rep(1L, n_rows)
        } else if (n_rows %% n_alt == 0L) {
          n_occasions <- n_rows %/% n_alt
          occasions <- rep(seq_len(n_occasions), each = n_alt)
        } else {
          stop(
            paste(
              "Long data without an occasion identifier must contain",
              "consecutive complete choice sets within each decider."
            ),
            call. = FALSE
          )
        }
        df[[column_occasion_id]][index_decider] <- occasions
      }
    }
  }
  df_ext <- df[0, , drop = FALSE]
  template_row <- df[1, , drop = FALSE]
  # group rows by decider and occasion
  ids <- df[[column_decider_id]]
  index_deciders <- split(seq_along(ids), match(ids, unique(ids)))
  choice_sets <- unlist(lapply(index_deciders, function(index_decider) {
    occasions <- df[[column_occasion_id]][index_decider]
    split(index_decider, match(occasions, unique(occasions)))
  }), recursive = FALSE, use.names = FALSE)
  # pre-allocated list to save computation time
  rows <- vector("list", length(choice_sets) * length(alt_names))
  row_index <- 0L
  for (index_rows in choice_sets) {
    alts_occurring <- df[[column_alternative]][index_rows]
    jd <- df[[column_decider_id]][index_rows[1]]
    jo <- df[[column_occasion_id]][index_rows[1]]
    alt_match <- match(alts_occurring, alt_names)
    alt_count <- tabulate(alt_match, nbins = length(alt_names))
    if (any(alt_count > 1L)) {
      stop(
        "Each long-data choice set must contain each alternative once.",
        call. = FALSE
      )
    }
    ind_alt <- match(alt_names, alts_occurring)
    for (ja in seq_along(alt_names)) {
      row_index <- row_index + 1L
      if (!is.na(ind_alt[ja])) {
        index_alt <- index_rows[ind_alt[ja]]
        rows[[row_index]] <- df[index_alt, , drop = FALSE]
      } else {
        # add row with missing covariates for this alternative
        added_row <- template_row
        added_row[, names(df) != column_alternative] <- NA
        added_row[[column_decider_id]] <- jd
        added_row[[column_occasion_id]] <- jo
        added_row[[column_alternative]] <- alt_names[ja]
        added_row[[column_choice]] <- 0
        rows[[row_index]] <- added_row
      }
    }
  }
  df_ext <- do.call(rbind, rows)
  return(df_ext)
}
