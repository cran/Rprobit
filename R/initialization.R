# The following non-base functions are used in initialization.R:
#
# # from `checkmate`:
#
# - `assert_count()`
# - `assert_int()`
# - `assert_integerish()`
# - `assert_matrix()`
# - `assert_number()`
# - `assert_true()`
#
# # from `oeli`:
#
# - `assert_covariance_matrix()`
# - `assert_list_of_lists()`
# - `assert_numeric_vector()`
# - `assert_one_hot_matrix()`
# - `assert_probability_vector()`
# - `chunk_vector()`
# - `cov_to_chol()`
# - `delta()`
# - `insert_matrix_column()`
#
# # from `Rprobit`:
#
# - `SJ_grad()`
# - `SJ()`
#
# # from `RprobitB`:
#
# - `update_coefficient()`
# - `update_Omega()`
# - `update_Sigma()`
# - `update_U()`
#
# # from `stats`:
#
# - `cov()`
# - `kmeans()`
# - `lm.fit()`
# - `median()`
# - `nlm()`

#' Initialize nMMNP Model Parameters
#'
#' @description
#' This function computes initial values for all parameters of a normally mixed
#' multinomial probit (nMMNP) model.
#'
#' @param data \[`list`\]\cr
#' A list, each element representing choice occasions by one decider, containing
#' - `X`, a list of design matrices,
#' - `y`, an integer vector, representing the chosen alternatives.
#'
#' @param J \[`integer(1)`\]\cr
#' The number of choice alternatives.
#'
#' @param P \[`integer(3)`\]\cr
#' The numbers of type-1, type-2, and type-3 coefficients.
#'
#' @param ref_alt_index \[`integer(1)`\]\cr
#' Position of the reference alternative for type-2 coefficients.
#'
#' @param re_cols \[`integer`\]\cr
#' The design matrix column positions connected to Gaussian random effects.
#'
#' @param cluster_max \[`integer(1)`\]\cr
#' The maximum number of clusters for localization in type-2 initialization.
#'
#' @param cluster_min_weight \[`numeric(1)`\]\cr
#' The minimum relative cluster size for localization in type-2 initialization.
#'
#' @param cluster_value_tolerance \[`numeric(1)`\]\cr
#' The minimum absolute cluster value for localization in type-2 initialization.
#'
#' @param R \[`integer(1)`\]\cr
#' The number of Gibbs iterations for covariance initialization.
#'
#' @param B \[`integer(1)`\]\cr
#' The length of the burn-in period for covariance initialization.
#'
#' @param Q \[`integer(1)`\]\cr
#' The thinning factor for covariance initialization.
#'
#' @param N_max \[`integer(1)`\]\cr
#' The maximum number of deciders used for covariance initialization.
#'
#' @param check_inputs \[`logical(1)`\]\cr
#' If `TRUE`, inputs are checked.
#'
#' @return
#' A `numeric` vector containing
#'
#' - the coefficient initial values `b`, scaled to have Euclidean norm one,
#'   which fixes the common utility scale,
#' - followed by the Cholesky-root parameters of the Gaussian random-effect
#'   covariance matrix `Omega` if any `re_cols` are set,
#' - and the Cholesky-root parameters of the differenced `(J - 1) x (J - 1)`
#'   error-term covariance matrix `Sigma` relative to alternative 1 for level
#'   identification.
#'
#' The computation time in seconds is attached as attribute `sec`.

initialize <- function(
    data, J, P, ref_alt_index = 1L, re_cols = integer(0),
    cluster_max = 10, cluster_min_weight = 0.1, cluster_value_tolerance = 1e-3,
    R = 200, B = R %/% 2, Q = 1, N_max = 1000, check_inputs = TRUE
) {
  if (isTRUE(check_inputs)) {
    oeli::assert_list_of_lists(data)
    checkmate::assert_int(J, lower = 2)
    checkmate::assert_integerish(P, lower = 0, any.missing = FALSE, len = 3)
    checkmate::assert_int(ref_alt_index, lower = 1L, upper = J)
    checkmate::assert_int(cluster_max, lower = 1)
    checkmate::assert_number(
      cluster_min_weight, lower = 0.01, upper = 1, finite = TRUE
    )
    checkmate::assert_number(cluster_value_tolerance, lower = 0, finite = TRUE)
    checkmate::assert_count(R, positive = TRUE)
    checkmate::assert_int(B, lower = 0, upper = R - 1)
    checkmate::assert_int(Q, lower = 1, upper = R - B)
    checkmate::assert_int(N_max, lower = 1)
  }
  J <- as.integer(J)
  P <- as.integer(P)
  ref_alt_index <- as.integer(ref_alt_index)
  b_len <- P[1] + P[2] * (J - 1) + P[3] * J
  if (isTRUE(check_inputs)) {
    checkmate::assert_integerish(
      re_cols, any.missing = FALSE, unique = TRUE
    )
    if (length(re_cols) > 0L) {
      checkmate::assert_integerish(
        re_cols, lower = 1, upper = b_len, any.missing = FALSE, unique = TRUE
      )
    }
  }
  re_cols <- as.integer(re_cols)

  # put the reference row first and group blocks by regressor type
  X <- unlist(lapply(data, `[[`, "X"), recursive = FALSE)
  if (sum(P > 0L) == 1L) {
    block_types <- rep(which(P > 0L), P[P > 0L])
  } else {
    block_types <- regressor_types(
      X = X, J = J, P = P, ref_alt_index = ref_alt_index,
      check_inputs = check_inputs
    )
  }
  alt_order <- c(ref_alt_index, setdiff(seq_len(J), ref_alt_index))
  widths <- c(1L, J - 1L, J)
  block_widths <- widths[block_types]
  block_ends <- cumsum(block_widths)
  block_starts <- block_ends - block_widths + 1L
  blocks <- Map(seq.int, from = block_starts, to = block_ends)
  type_3 <- which(block_types == 3L)
  blocks[type_3] <- lapply(blocks[type_3], `[`, alt_order)
  block_order <- order(block_types, seq_along(block_types))
  col_order <- unlist(blocks[block_order], use.names = FALSE)
  data <- lapply(data, function(decider) {
    decider$X <- lapply(decider$X, function(X_i) {
      X_i[alt_order, col_order, drop = FALSE]
    })
    decider$y <- match(decider$y, alt_order)
    decider
  })
  re_cols <- match(re_cols, col_order)

  ### get initial values
  t0 <- Sys.time()
  init_b <- initialize_b(
    data = data, J = J, P = P, ref_alt_index = 1L,
    cluster_max = cluster_max, cluster_min_weight = cluster_min_weight,
    cluster_value_tolerance = cluster_value_tolerance,
    check_inputs = check_inputs
  )

  # undo the type-2 row order so covariance differences use alternative 1
  inverse_alt_order <- order(alt_order)
  data <- lapply(data, function(decider) {
    decider$X <- lapply(decider$X, function(X_i) {
      X_i[inverse_alt_order, , drop = FALSE]
    })
    decider$y <- alt_order[decider$y]
    decider
  })
  init_cov <- initialize_covariance(
    data = data, J = J, b_fixed = init_b, re_cols = re_cols,
    R = R, B = B, Q = Q, N_max = N_max,
    check_inputs = check_inputs
  )
  init_o <- init_cov$init_o
  init_l <- init_cov$init_l

  # restore coefficients to the original design order
  original_b <- numeric(length(init_b))
  original_b[col_order] <- init_b
  init_b <- original_b
  out <- c(init_b, init_o, init_l)
  structure(out, "sec" = difftime(Sys.time(), t0, units = "secs"))
}

# Orthogonalize regressors via the FWL theorem
# y: numeric matrix
# X: numeric matrix (nrow(y) x p)
# ind: indices of X columns to orthogonalize, returned in this order
# check_inputs: whether to check inputs
# Returns: list(Ry, RX), i.e., y-residuals and X[, ind]-residuals wrt X[, -ind]
orthogonalize_regressors <- function(y, X, ind, check_inputs = TRUE) {
  if (isTRUE(check_inputs)) {
    checkmate::assert_matrix(
      y, mode = "numeric", any.missing = TRUE, min.rows = 1, min.cols = 1
    )
    checkmate::assert_true(all(is.finite(y[!is.na(y)])))
    n <- nrow(y)
    checkmate::assert_matrix(
      X, mode = "numeric", any.missing = TRUE, nrows = n, min.cols = 2
    )
    checkmate::assert_true(all(is.finite(X[!is.na(X)])))
    p <- ncol(X)
    checkmate::assert_integerish(
      ind, lower = 1, upper = p, min.len = 1, max.len = p - 1, unique = TRUE
    )
  } else {
    n <- nrow(y)
    p <- ncol(X)
  }
  ind <- as.integer(ind)
  mind <- setdiff(seq_len(p), ind)
  keep <- stats::complete.cases(y, X)
  y <- y[keep, , drop = FALSE]
  X <- X[keep, , drop = FALSE]
  if (nrow(X) == 0L) {
    return(list(
      Ry = matrix(numeric(), nrow = 0L, ncol = ncol(y)),
      RX = matrix(numeric(), nrow = 0L, ncol = length(ind))
    ))
  }
  Y <- cbind(y, X[, ind, drop = FALSE])
  fit <- stats::lm.fit(x = X[, mind, drop = FALSE], y = Y)
  res <- fit$residuals
  Ry <- res[, seq_len(ncol(y)), drop = FALSE]
  RX <- res[, ncol(y) + seq_along(ind), drop = FALSE]
  list(Ry = Ry, RX = RX)
}

# Localize regressor values using k-means
# values: numeric vector
# cluster_max: maximum number of clusters
# cluster_min_weight: minimum relative cluster size (at least 0.01)
# check_inputs: whether to check inputs
# Returns: list(loc, ids, w), i.e., cluster centers, allocation ids, weights
localize_regressors <- function(
    values, cluster_max = 10, cluster_min_weight = 0.1, check_inputs = TRUE
) {
  if (isTRUE(check_inputs)) {
    oeli::assert_numeric_vector(
      values, finite = TRUE, any.missing = FALSE, min.len = 1
    )
    checkmate::assert_int(cluster_max, lower = 1)
    checkmate::assert_number(cluster_min_weight, lower = 0.01, upper = 1)
  }
  values <- as.numeric(values)
  n <- length(values)
  # cap centers by distinct values, minimum weight, and the sample size
  centers <- as.integer(min(
    cluster_max,
    length(unique(values)),
    floor(1 / cluster_min_weight),
    max(1L, n - 1L)
  ))
  x <- matrix(values, ncol = 1)
  repeat {
    if (centers == 1L) {
      loc <- mean(values)
      ids <- rep(1L, n)
      w <- 1
      break
    }
    cl <- NULL
    for (method in c("Hartigan-Wong", "MacQueen")) {
      cl_try <- tryCatch(
        suppressWarnings(stats::kmeans(
          x = x, centers = centers, iter.max = 100, nstart = 100,
          algorithm = method
        )),
        error = function(e) NULL
      )
      valid <- !is.null(cl_try) &&
        (is.null(cl_try$ifault) || isTRUE(cl_try$ifault == 0L)) &&
        length(cl_try$cluster) == n &&
        length(cl_try$centers) == centers &&
        all(is.finite(cl_try$centers)) &&
        all(cl_try$cluster %in% seq_len(centers))
      if (valid) {
        cl <- cl_try
        break
      }
    }
    if (is.null(cl)) {
      centers <- centers - 1L
      next
    }
    ids <- cl$cluster
    loc <- as.numeric(cl$centers)
    w <- tabulate(ids, nbins = centers) / n
    if (any(w < cluster_min_weight)) {
      centers <- centers - 1L
    } else {
      break
    }
  }
  # order according to descending weight
  o <- order(w, decreasing = TRUE)
  w <- w[o]
  loc <- loc[o]
  new_id <- integer(length(o))
  new_id[o] <- seq_along(o)
  ids <- new_id[ids]
  list(loc = loc, ids = ids, w = w)
}

# Compute choice probabilities from systematic utilities and covariance
# sys: length-J utility vector (J >= 2)
# cov: J x J covariance matrix
# check_inputs: whether to check inputs
# Returns: numeric vector of probabilities with J x J Jacobian attribute
sys_to_prob <- function(sys, cov, check_inputs = TRUE) {
  if (isTRUE(check_inputs)) {
    oeli::assert_numeric_vector(
      sys, finite = TRUE, any.missing = FALSE, min.len = 2
    )
    J <- length(sys)
    oeli::assert_covariance_matrix(cov, dim = J)
  } else {
    J <- length(sys)
  }
  sys <- as.numeric(sys)
  jac <- matrix(0, J, J)
  prob <- numeric(J)
  for (j in seq_len(J)) {
    D <- oeli::delta(ref = j, dim = J)
    sys_d <- as.numeric(D %*% sys)
    cov_d <- D %*% cov %*% t(D)
    scale <- 1 / sqrt(diag(cov_d))
    cov_t <- cov_d * tcrossprod(scale)
    sys_t <- sys_d * scale
    p <- exp(SJ(-sys_t, cov_t))
    prob[j] <- p
    gr <- as.numeric(SJ_grad(-sys_t, cov_t))
    jac[j, ] <- -as.numeric((gr[seq_len(J - 1)] * scale) %*% D) * p
  }
  structure(prob, jac = jac)
}

# Compute SSR between choice probabilities and observed frequencies
# sys_m1: length-(J - 1) utility vector excluding the first alternative
# cov: J x J covariance matrix
# freq: length-J probability vector of observed choice frequencies
# fix_sys1: fixed systematic utility for the first alternative
# check_inputs: whether to check inputs
# Returns: scalar SSR with gradient attribute with respect to sys_m1
ssr_prob_freq <- function(
    sys_m1, cov, freq, fix_sys1 = 0, check_inputs = TRUE
  ) {
  J <- length(sys_m1) + 1
  if (isTRUE(check_inputs)) {
    oeli::assert_numeric_vector(
      sys_m1, finite = TRUE, any.missing = FALSE, min.len = 1
    )
    oeli::assert_covariance_matrix(cov, dim = J)
    oeli::assert_probability_vector(freq, len = J)
    checkmate::assert_number(fix_sys1, na.ok = FALSE, finite = TRUE)
  }
  sys_m1 <- as.numeric(sys_m1)
  freq <- as.numeric(freq)
  fix_sys1 <- as.numeric(fix_sys1)
  prob <- sys_to_prob(
    sys = c(fix_sys1, sys_m1), cov = cov, check_inputs = FALSE
  )
  res <- freq - prob
  crit <- sum(res^2)
  jac <- attr(prob, "jac")[, -1, drop = FALSE]
  gradient <- drop(-2 * crossprod(jac, res))
  structure(crit, gradient = gradient)
}

# Check whether a design column is a type-1 regressor
# X: list of numeric design matrices with alternatives in rows
# i: index of the candidate column
# Returns: TRUE if column i varies across alternatives in any matrix
is_type_1 <- function(X, i) {
  any(vapply(
    X,
    function(X_i) {
      values <- X_i[, i]
      length(unique(values[!is.na(values)])) > 1L
    },
    logical(1)
  ))
}

# Check whether a design block is a type-2 regressor
# X: list of numeric design matrices with alternatives in rows
# J: number of choice alternatives (J >= 2)
# i: index of the first candidate column
# ref_alt_index: position of the reference alternative for type-2 coefficients
# Returns: TRUE if every block has type-2 structure and at least one matrix
# has a nonzero shared diagonal value
is_type_2 <- function(X, J, i, ref_alt_index) {
  nonzero_diag <- FALSE
  for (X_i in X) {
    B <- X_i[, i:(i + J - 2L), drop = FALSE]
    S <- B[-ref_alt_index, , drop = FALSE]
    d <- diag(S)
    ref <- B[ref_alt_index, ]
    off_diag <- S[row(S) != col(S)]
    ref <- ref[!is.na(ref)]
    off_diag <- off_diag[!is.na(off_diag)]
    d <- d[!is.na(d)]
    valid <- all(ref == 0) && all(off_diag == 0) && length(unique(d)) <= 1L
    if (!valid) return(FALSE)
    nonzero_diag <- nonzero_diag || any(d != 0)
  }
  nonzero_diag
}

# Check whether a design block is a type-3 regressor
# X: list of numeric design matrices with alternatives in rows
# J: number of choice alternatives (J >= 2)
# i: index of the first candidate column
# Returns: TRUE if every block is diagonal and every alternative has a
# nonzero diagonal entry in at least one matrix
is_type_3 <- function(X, J, i) {
  nonzero_diag <- logical(J)
  for (X_i in X) {
    C <- X_i[, i:(i + J - 1L), drop = FALSE]
    off_diag <- C[row(C) != col(C)]
    off_diag <- off_diag[!is.na(off_diag)]
    if (!all(off_diag == 0)) return(FALSE)
    d <- diag(C)
    nonzero_diag[!is.na(d)] <- nonzero_diag[!is.na(d)] | d[!is.na(d)] != 0
  }
  all(nonzero_diag)
}

# Detects regressor types from design matrices
# The matrices are scanned from left to right:
# - Single columns with alternative variation are type-1,
# - J - 1 reference-coded diagonal columns are type-2,
# - and J diagonal columns are type-3.
# - For J = 2, type-1 and type-2 columns are indistinguishable; ambiguous
#   columns are assigned to type-1 until P[1] type-1 regressors have been found.
# X: J x (P[1] + P[2] * (J - 1) + P[3] * J) design matrix, or a list thereof
# J: number of choice alternatives (J >= 2)
# P: integer vector of length 3 with regressor-type counts
# ref_alt_index: position of the reference alternative for type-2 coefficients
# check_inputs: whether to check inputs
# Returns: integer vector with entries 1, 2, or 3, one per detected regressor
regressor_types <- function(X, J, P, ref_alt_index, check_inputs = TRUE) {
  if (isTRUE(check_inputs)) {
    checkmate::assert_int(J, lower = 2)
    checkmate::assert_integerish(P, lower = 0, any.missing = FALSE, len = 3)
    checkmate::assert_int(ref_alt_index, lower = 1L, upper = J)
    P <- as.integer(P)
    J <- as.integer(J)
    ref_alt_index <- as.integer(ref_alt_index)
    pX <- P[1] + P[2] * (J - 1) + P[3] * J
    if (is.matrix(X)) X <- list(X)
    checkmate::assert_list(X, min.len = 1L)
    for (X_i in X) {
      checkmate::assert_matrix(
        X_i, mode = "numeric", any.missing = TRUE, nrows = J, ncols = pX
      )
      checkmate::assert_true(all(is.finite(X_i[!is.na(X_i)])))
    }
  } else {
    P <- as.integer(P)
    J <- as.integer(J)
    ref_alt_index <- as.integer(ref_alt_index)
    if (is.matrix(X)) {
      X <- list(X)
    }
  }
  # scan from left to right and advance by the detected block width
  nc <- ncol(X[[1L]])
  i <- 1L
  k <- 0L
  widths <- c(1L, J - 1L, J)
  out <- integer(sum(P))
  while (i <= nc) {
    found <- tabulate(out[seq_len(k)], nbins = 3L)
    remaining <- P - found
    # test structured blocks first
    if (i + J - 2L <= nc) {
      if (remaining[2L] > 0L && is_type_2(X, J, i, ref_alt_index)) {
        if (J == 2L && is_type_1(X, i)) {
          # use P to resolve type-1 and type-2 ambiguity for J = 2
          if (remaining[1L] > 0L) {
            k <- k + 1L
            out[k] <- 1L
            i <- i + 1L
            next
          }
        }
        k <- k + 1L
        out[k] <- 2L
        i <- i + J - 1L
        next
      }
    }
    if (i + J - 1L <= nc) {
      if (remaining[3L] > 0L && is_type_3(X, J, i)) {
        k <- k + 1L
        out[k] <- 3L
        i <- i + J
        next
      }
    }
    # finally, test for simple type-1 blocks
    if (remaining[1L] > 0L && is_type_1(X, i)) {
      k <- k + 1L
      out[k] <- 1L
      i <- i + 1L
      next
    }
    # assign unrecognized blocks from the remaining counts
    type <- which(remaining > 0L & widths <= nc - i + 1L)[1L]
    if (is.na(type)) {
      warning(
        "Conflict in detecting regressor types; guessing them.",
        call. = FALSE
      )
      out <- rep(seq_along(P), P)
      k <- length(out)
      break
    }
    last <- i + widths[type] - 1L
    zero <- all(vapply(
      X,
      function(X_i) {
        values <- X_i[, i:last, drop = FALSE]
        all(values[!is.na(values)] == 0)
      },
      logical(1)
    ))
    warning(
      sprintf(
        "The %s block at column %d was assigned to type-%d.",
        if (zero) "zero" else "unrecognized", i, type
      ),
      call. = FALSE
    )
    k <- k + 1L
    out[k] <- type
    i <- last + 1L
  }
  out <- if (k > 0L) out[seq_len(k)] else integer()
  detP <- tabulate(out, nbins = 3L)
  if (any(detP != P)) {
    warning(
      "Conflict in detecting regressor types; guessing them.",
      call. = FALSE
    )
    out <- rep(seq_along(P), P)
  }
  as.integer(out)
}

# Initialize type-1 coefficients via a linear probability model
# X: N x (J * P[1]) numeric matrix of centered type-1 regressors
# y: N x J one-hot encoded choice matrix
# check_inputs: whether to check inputs
# Returns: numeric vector of length P[1], rescaled to have Euclidean norm one,
# with computation time attribute
initialize_b_type_1 <- function(X, y, check_inputs = TRUE) {
  if (isTRUE(check_inputs)) {
    checkmate::assert_matrix(
      y, mode = "numeric", any.missing = TRUE, min.rows = 1L, min.cols = 2L
    )
    checkmate::assert_true(all(is.finite(y[!is.na(y)])))
    J <- ncol(y)
    N <- nrow(y)
    checkmate::assert_int(J, lower = 2)
    checkmate::assert_matrix(
      X, nrows = N, min.cols = 1, mode = "numeric", any.missing = TRUE
    )
    checkmate::assert_true(all(is.finite(X[!is.na(X)])))
    P1 <- ncol(X) / J
    checkmate::assert_int(P1)
  } else {
    J <- ncol(y)
    N <- nrow(y)
    P1 <- ncol(X) / J
  }
  keep <- stats::complete.cases(y, X)
  y <- y[keep, , drop = FALSE]
  X <- X[keep, , drop = FALSE]
  if (isTRUE(check_inputs) && nrow(y) > 0L) {
    oeli::assert_one_hot_matrix(y, ncols = J)
  }
  t0 <- Sys.time()
  if (nrow(y) == 0L) {
    warning(
      "Type-1 initialization failed; using unit fallback.",
      call. = FALSE
    )
    return(structure(
      numeric(P1), sec = difftime(Sys.time(), t0, units = "secs")
    ))
  }
  N <- nrow(y)
  X <- array(t(X), dim = c(J, P1, N))
  X <- matrix(aperm(X, c(1, 3, 2)), ncol = P1)
  y <- as.numeric(t(y))
  col_norm <- sqrt(colSums(X^2))
  active <- is.finite(col_norm) & col_norm > 1e-8 * max(col_norm, 0)
  est <- numeric(P1)
  if (any(active)) {
    est[active] <- stats::lm.fit(x = X[, active, drop = FALSE], y = y)$coefficients
  }
  est[!is.finite(est)] <- 0
  norm_est <- sqrt(sum(est^2))
  if (!is.finite(norm_est) || norm_est == 0) {
    warning(
      "Type-1 initialization failed; using unit fallback.",
      call. = FALSE
    )
    est <- numeric(P1)
  } else {
    est <- as.numeric(est / norm_est)
  }
  structure(est, sec = difftime(Sys.time(), t0, units = "secs"))
}

# Initialize one type-2 coefficient block via choice frequency matching
# X: N x (J - 1) numeric matrix of one type-2 regressor
# y: N x J one-hot encoded choice matrix
# cluster_max: maximum number of clusters
# cluster_min_weight: minimum relative cluster size (at least 0.01)
# cluster_value_tolerance: minimum absolute cluster value
# check_inputs: whether to check inputs
# Returns: numeric vector of length J - 1, with computation time attribute
initialize_b_type_2_single <- function(
    X, y, cluster_max = 10, cluster_min_weight = 0.1,
    cluster_value_tolerance = 1e-3, check_inputs = TRUE
) {
  if (isTRUE(check_inputs)) {
    checkmate::assert_matrix(
      y, mode = "numeric", any.missing = TRUE, min.rows = 1L, min.cols = 2L
    )
    checkmate::assert_true(all(is.finite(y[!is.na(y)])))
    J <- ncol(y)
    N <- nrow(y)
    checkmate::assert_int(J, lower = 2)
    checkmate::assert_matrix(
      X, nrows = N, ncols = J - 1L, mode = "numeric",
      any.missing = TRUE
    )
    checkmate::assert_true(all(is.finite(X[!is.na(X)])))
    checkmate::assert_int(cluster_max, lower = 1)
    checkmate::assert_number(
      cluster_min_weight, lower = 0.01, upper = 1, finite = TRUE
    )
    checkmate::assert_number(cluster_value_tolerance, lower = 0, finite = TRUE)
  } else {
    J <- ncol(y)
  }
  t0 <- Sys.time()
  keep <- stats::complete.cases(y, X)
  y <- y[keep, , drop = FALSE]
  X <- X[keep, , drop = FALSE]
  if (isTRUE(check_inputs) && nrow(y) > 0L) {
    oeli::assert_one_hot_matrix(y, ncols = J)
  }
  zero <- rep(0, J - 1L)
  if (nrow(y) == 0L) {
    warning(
      "Type-2 initialization failed; using unit fallback.",
      call. = FALSE
    )
    return(structure(zero, sec = difftime(Sys.time(), t0, units = "secs")))
  }
  # initialize binary type-2 block as type-1 block
  if (J == 2L) {
    X_type_1 <- cbind(0, X)
    est <- initialize_b_type_1(X = X_type_1, y = y, check_inputs = FALSE)
    return(structure(
      as.numeric(est), sec = difftime(Sys.time(), t0, units = "secs")
    ))
  }
  out <- localize_regressors(
    values = rowMeans(X), cluster_max = cluster_max,
    cluster_min_weight = cluster_min_weight, check_inputs = FALSE
  )
  loc <- out$loc
  ids <- out$ids
  w <- out$w
  cl <- which(abs(loc) > cluster_value_tolerance)
  w <- w[cl]
  cov <- diag(J)
  ests <- matrix(0, nrow = length(cl), ncol = J - 1L)
  valid <- logical(length(cl))
  cluster_failed <- FALSE
  for (k in seq_along(cl)) {
    i <- cl[k]
    ind <- which(ids == i)
    freq <- colSums(y[ind, , drop = FALSE]) / length(ind)
    fit <- tryCatch(
      stats::nlm(
        ssr_prob_freq, p = zero, cov = cov, freq = freq,
        fix_sys1 = 0, check_inputs = FALSE
      ),
      error = function(e) NULL
    )
    fit_ok <- !is.null(fit) &&
      is.numeric(fit$estimate) &&
      length(fit$estimate) == J - 1L &&
      all(is.finite(fit$estimate)) &&
      is.numeric(fit$minimum) &&
      length(fit$minimum) == 1L &&
      is.finite(fit$minimum) &&
      isTRUE(fit$code %in% 1:2)
    if (!fit_ok) {
      cluster_failed <- TRUE
      next
    }
    est <- fit$estimate
    est <- est / loc[i]
    norm_est <- sqrt(sum(est^2))
    if (!is.finite(norm_est) || norm_est == 0) {
      cluster_failed <- TRUE
      next
    }
    valid[k] <- TRUE
    ests[k, ] <- as.numeric(est / norm_est)
  }
  if (any(valid)) {
    w <- w[valid] / sum(w[valid])
    est <- as.numeric(w %*% ests[valid, , drop = FALSE])
    if (cluster_failed) {
      warning(
        "Type-2 initialization failed for some clusters; ignoring them.",
        call. = FALSE
      )
    }
  } else {
    est <- zero
  }
  norm_est <- sqrt(sum(est^2))
  if (!is.finite(norm_est) || norm_est == 0) {
    warning(
      "Type-2 initialization failed; using unit fallback.",
      call. = FALSE
    )
    est <- zero
  }
  structure(est, sec = difftime(Sys.time(), t0, units = "secs"))
}

# Initialize type-2 coefficients via choice frequency matching
# X: N x (P[2] * (J - 1)) numeric matrix of type-2 regressors
# y: N x J one-hot encoded choice matrix
# cluster_max: maximum number of clusters
# cluster_min_weight: minimum relative cluster size (at least 0.01)
# cluster_value_tolerance: minimum absolute cluster value
# check_inputs: whether to check inputs
# Returns: numeric vector of length P[2] * (J - 1), rescaled to have Euclidean
# norm one, with computation time attribute
initialize_b_type_2 <- function(
    X, y, cluster_max = 10, cluster_min_weight = 0.1,
    cluster_value_tolerance = 1e-3, check_inputs = TRUE
) {
  if (isTRUE(check_inputs)) {
    checkmate::assert_matrix(
      y, mode = "numeric", any.missing = TRUE, min.rows = 1L, min.cols = 2L
    )
    checkmate::assert_true(all(is.finite(y[!is.na(y)])))
    J <- ncol(y)
    N <- nrow(y)
    checkmate::assert_int(J, lower = 2)
    checkmate::assert_matrix(
      X, nrows = N, min.cols = 1, mode = "numeric", any.missing = TRUE
    )
    checkmate::assert_true(all(is.finite(X[!is.na(X)])))
    P2 <- ncol(X) / (J - 1)
    checkmate::assert_int(P2)
    checkmate::assert_int(cluster_max, lower = 1)
    checkmate::assert_number(
      cluster_min_weight, lower = 0.01, upper = 1, finite = TRUE
    )
    checkmate::assert_number(cluster_value_tolerance, lower = 0, finite = TRUE)
  } else {
    J <- ncol(y)
    N <- nrow(y)
    P2 <- ncol(X) / (J - 1)
  }
  keep <- stats::complete.cases(y, X)
  y <- y[keep, , drop = FALSE]
  X <- X[keep, , drop = FALSE]
  if (isTRUE(check_inputs) && nrow(y) > 0L) {
    oeli::assert_one_hot_matrix(y, ncols = J)
  }
  N <- nrow(y)
  t0 <- Sys.time()
  if (N == 0L) {
    warning(
      "Type-2 initialization failed; using unit fallback.",
      call. = FALSE
    )
    return(structure(
      numeric(P2 * (J - 1L)),
      sec = difftime(Sys.time(), t0, units = "secs")
    ))
  }
  # initialize all binary type-2 blocks jointly as type-1 blocks
  if (J == 2L) {
    X_type_1 <- matrix(0, nrow = N, ncol = 2L * P2)
    X_type_1[, 2L * seq_len(P2)] <- X
    est <- initialize_b_type_1(X = X_type_1, y = y, check_inputs = FALSE)
    return(structure(
      as.numeric(est), sec = difftime(Sys.time(), t0, units = "secs")
    ))
  }
  if (P2 == 1) {
    est <- initialize_b_type_2_single(
      X = X, y = y, cluster_max = cluster_max,
      cluster_min_weight = cluster_min_weight,
      cluster_value_tolerance = cluster_value_tolerance,
      check_inputs = FALSE
    )
  } else {
    idx <- oeli::chunk_vector(seq_len(P2 * (J - 1)), n = J - 1, type = 2)
    est <- numeric((J - 1) * P2)
    for (p in seq_len(P2)) {
      ind <- idx[[p]]
      R <- orthogonalize_regressors(
        y = y, X = X, ind = ind, check_inputs = check_inputs
      )
      est[ind] <- initialize_b_type_2_single(
        X = R$RX, y = y, cluster_max = cluster_max,
        cluster_min_weight = cluster_min_weight,
        cluster_value_tolerance = cluster_value_tolerance,
        check_inputs = FALSE
      )
    }
    Xb <- sweep(X, 2, est, "*")
    zero_pos <- (J - 1) * (seq_len(P2) - 1L)
    Xb0 <- oeli::insert_matrix_column(Xb, 0, zero_pos)
    norm <- initialize_b_type_1(
      Xb0, y, check_inputs = check_inputs
    )
    est <- as.vector(
      sweep(matrix(est, J - 1, P2), 2, norm, "*")
    )
  }
  norm_est <- sqrt(sum(est^2))
  if (!is.finite(norm_est) || norm_est == 0) {
    if (P2 > 1L) {
      warning(
        "Type-2 initialization failed; using unit fallback.",
        call. = FALSE
      )
    }
    est <- numeric(length(est))
  } else {
    est <- as.numeric(est / norm_est)
  }
  structure(est, sec = difftime(Sys.time(), t0, units = "secs"))
}

# Initialize one type-3 coefficient block via constant utility direction
# X: N x J numeric matrix of one type-3 regressor
# y: N x J numeric response matrix
# check_inputs: whether to check inputs
# Returns: numeric vector of length J with Euclidean norm one and computation
# time attribute
initialize_b_type_3_single <- function(X, y, check_inputs = TRUE) {
  J <- ncol(y)
  N <- nrow(y)
  if (isTRUE(check_inputs)) {
    checkmate::assert_matrix(
      y, mode = "numeric", any.missing = TRUE, min.rows = 1L, min.cols = 2L
    )
    checkmate::assert_matrix(
      X, nrows = N, ncols = J, mode = "numeric", any.missing = TRUE
    )
    checkmate::assert_true(all(is.finite(y[!is.na(y)])))
    checkmate::assert_true(all(is.finite(X[!is.na(X)])))
  }
  t0 <- Sys.time()
  keep <- stats::complete.cases(y, X)
  y <- y[keep, , drop = FALSE]
  X <- X[keep, , drop = FALSE]
  if (nrow(y) == 0L) {
    warning(
      "Type-3 initialization failed; using unit fallback.",
      call. = FALSE
    )
    return(structure(
      numeric(J), sec = difftime(Sys.time(), t0, units = "secs")
    ))
  }
  M <- stats::cov(y, X)
  decomp <- tryCatch(svd(M, nu = 0, nv = J), error = function(e) NULL)
  tol <- sqrt(.Machine$double.eps)
  finite_svd <- !is.null(decomp) && all(is.finite(c(decomp$d, decomp$v)))
  unique_direction <- finite_svd && decomp$d[1L] > 0 &&
    decomp$d[J - 1L] - decomp$d[J] > tol * decomp$d[1L]
  if (unique_direction) {
    v <- decomp$v[, J]
    v_tol <- tol * max(abs(v))
    k <- v / pmax(v^2, v_tol^2)
    k <- as.numeric(k / sqrt(sum(k^2)))
    clear_direction <- isTRUE(
      all(is.finite(k)) && abs(sum(diag(M) * k)) > tol * sum(abs(diag(M) * k))
    )
  } else {
    clear_direction <- FALSE
  }
  if (!clear_direction) {
    warning(
      "Type-3 initialization failed; using unit fallback.",
      call. = FALSE
    )
    est <- numeric(J)
  } else {
    # determine direction
    est <- sign(sum(diag(M) * k)) * k
  }
  structure(est, sec = difftime(Sys.time(), t0, units = "secs"))
}

# Initialize type-3 coefficients via constant utility direction
# X: N x (P[3] * J) numeric matrix of type-3 regressors
# y: N x J one-hot encoded choice matrix
# check_inputs: whether to check inputs
# Returns: numeric vector of length P[3] * J, rescaled to have Euclidean norm
# one, with computation time attribute
initialize_b_type_3 <- function(X, y, check_inputs = TRUE) {
  if (isTRUE(check_inputs)) {
    checkmate::assert_matrix(
      y, mode = "numeric", any.missing = TRUE, min.rows = 1L, min.cols = 2L
    )
    checkmate::assert_true(all(is.finite(y[!is.na(y)])))
    J <- ncol(y)
    N <- nrow(y)
    checkmate::assert_int(J, lower = 2)
    checkmate::assert_matrix(
      X, nrows = N, min.cols = 1, mode = "numeric", any.missing = TRUE
    )
    checkmate::assert_true(all(is.finite(X[!is.na(X)])))
    P3 <- ncol(X) / J
    checkmate::assert_int(P3)
  } else {
    J <- ncol(y)
    N <- nrow(y)
    P3 <- ncol(X) / J
  }
  keep <- stats::complete.cases(y, X)
  y <- y[keep, , drop = FALSE]
  X <- X[keep, , drop = FALSE]
  if (isTRUE(check_inputs) && nrow(y) > 0L) {
    oeli::assert_one_hot_matrix(y, ncols = J)
  }
  N <- nrow(y)
  t0 <- Sys.time()
  if (N == 0L) {
    warning(
      "Type-3 initialization failed; using unit fallback.",
      call. = FALSE
    )
    return(structure(
      numeric(P3 * J), sec = difftime(Sys.time(), t0, units = "secs")
    ))
  }
  if (P3 == 1) {
    est <- initialize_b_type_3_single(
      X = X, y = y, check_inputs = FALSE
    )
  } else {
    idx <- oeli::chunk_vector(seq_len(P3 * J), n = J, type = 2)
    est <- numeric(P3 * J)
    for (p in seq_len(P3)) {
      ind <- idx[[p]]
      R <- orthogonalize_regressors(
        y = y, X = X, ind = ind, check_inputs = check_inputs
      )
      est[ind] <- initialize_b_type_3_single(
        X = R$RX, y = R$Ry, check_inputs = FALSE
      )
    }
    active <- vapply(idx, function(ind) any(est[ind] != 0), logical(1))
    if (any(active)) {
      Xb <- sweep(X, 2, est, "*")
      cols <- rep(active, each = J)
      norm <- numeric(P3)
      norm[active] <- initialize_b_type_1(
        Xb[, cols, drop = FALSE], y, check_inputs = check_inputs
      )
      est <- as.vector(
        sweep(matrix(est, J, P3), 2, norm, "*")
      )
    }
  }
  norm_est <- sqrt(sum(est^2))
  if (!is.finite(norm_est) || norm_est == 0) {
    est <- numeric(length(est))
  } else {
    est <- as.numeric(est / norm_est)
  }
  structure(est, sec = difftime(Sys.time(), t0, units = "secs"))
}

# Initialize all coefficients
# data: list, each element is a list of
# - X, a list, where the t-th element is the covariate matrix at occasion t,
# - y, integer vector, t-th element is the chosen alternative at occasion t.
# J: number of choice alternatives (J >= 2)
# P: integer vector of length 3 with regressor-type counts
# ref_alt_index: position of the reference alternative for type-2 coefficients
# cluster_max: maximum number of clusters
# cluster_min_weight: minimum relative cluster size (at least 0.01)
# cluster_value_tolerance: minimum absolute cluster value
# check_inputs: whether to check inputs
# Returns: numeric vector of length P[1] + P[2] * (J - 1) + P[3] * J, rescaled
# to have Euclidean norm one, with computation time attribute
initialize_b <- function(
    data, J, P, ref_alt_index, cluster_max = 10, cluster_min_weight = 0.1,
    cluster_value_tolerance = 1e-3, check_inputs = TRUE
) {
  if (isTRUE(check_inputs)) {
    oeli::assert_list_of_lists(data)
    checkmate::assert_int(J, lower = 2)
    checkmate::assert_integerish(P, lower = 0, any.missing = FALSE, len = 3)
    checkmate::assert_int(ref_alt_index, lower = 1L, upper = J)
    checkmate::assert_int(cluster_max, lower = 1)
    checkmate::assert_number(
      cluster_min_weight, lower = 0.01, upper = 1, finite = TRUE
    )
    checkmate::assert_number(cluster_value_tolerance, lower = 0, finite = TRUE)
  }
  J <- as.integer(J)
  P <- as.integer(P)
  ref_alt_index <- as.integer(ref_alt_index)
  if (sum(P) == 0) {
    out <- structure(numeric(), "sec" = as.difftime(0, units = "secs"))
    return(out)
  }
  ### prepare data
  data_flat <- unlist(
    Map(function(d) Map(function(x, y) list(X = x, y = y), d$X, d$y), data),
    recursive = FALSE
  )
  NT <- length(data_flat)

  if (sum(P > 0L) == 1L) {
    regressor_order <- rep(which(P > 0L), P[P > 0L])
  } else {
    regressor_order <- regressor_types(
      X = lapply(data_flat, function(d) d$X),
      J = J, P = P, ref_alt_index = ref_alt_index,
      check_inputs = check_inputs
    )
  }
  if (isTRUE(check_inputs)) {
    stopifnot(
      "Detected regressor types must match `P`." =
        all(tabulate(regressor_order, nbins = 3) == P)
    )
  }
  regressor_position <- unlist(
    lapply(
      regressor_order, switch, "1" = 1, "2" = rep(2, J - 1), "3" = rep(3, J)
    ),
    use.names = FALSE
  )
  ind_t1 <- which(regressor_position == 1)
  ind_t2 <- which(regressor_position == 2)
  ind_t3 <- which(regressor_position == 3)
  reg <- P * c(J, J - 1, J)
  X <- matrix(NA_real_, nrow = NT, ncol = sum(reg))
  y <- matrix(NA_real_, nrow = NT, ncol = J)
  for (nt in seq_len(NT)) {
    X_nt <- data_flat[[nt]]$X
    if (length(ind_t1) > 0L) {
      X_nt[, ind_t1] <- sweep(
        X_nt[, ind_t1, drop = FALSE], 2, X_nt[1L, ind_t1], "-"
      )
    }
    val_t1 <- as.numeric(X_nt[, ind_t1, drop = FALSE])
    val_t2 <- colSums(X_nt[, ind_t2, drop = FALSE])
    val_t3 <- colSums(X_nt[, ind_t3, drop = FALSE])
    X[nt, ] <- c(val_t1, val_t2, val_t3)
    choice <- data_flat[[nt]]$y
    valid_choice <- length(choice) == 1L && !is.na(choice) &&
      is.finite(choice) && choice %in% seq_len(J)
    if (valid_choice) {
      y[nt, ] <- 0
      y[nt, choice] <- 1
    }
  }
  keep <- stats::complete.cases(y, X)
  X <- X[keep, , drop = FALSE]
  y <- y[keep, , drop = FALSE]
  ### build initial values
  t0 <- Sys.time()
  if (nrow(y) == 0L) {
    warning(
      "Coefficient initialization failed; using unit fallback.",
      call. = FALSE
    )
    out <- numeric(length(regressor_position))
    return(structure(
      out, sec = difftime(Sys.time(), t0, units = "secs")
    ))
  }
  NT <- nrow(y)
  n_types <- sum(P > 0)
  if (P[1] > 0) {
    if (sum(P) > P[1]) {
      ind <- seq_len(reg[1])
      ortho <- orthogonalize_regressors(
        y = y, X = X, ind = ind, check_inputs = check_inputs
      )
      RX <- ortho$RX
      Ry <- ortho$Ry
    } else {
      RX <- X
      Ry <- y
    }
    b_t1 <- initialize_b_type_1(
      X = RX, y = Ry,
      check_inputs = isTRUE(check_inputs) && n_types == 1L
    )
  } else {
    b_t1 <- numeric()
  }
  if (P[2] > 0) {
    if (sum(P) > P[2]) {
      ind <- reg[1] + seq_len(reg[2])
      ortho <- orthogonalize_regressors(
        y = y, X = X, ind = ind, check_inputs = check_inputs
      )
      RX <- ortho$RX
    } else {
      RX <- X
    }
    b_t2 <- initialize_b_type_2(
      X = RX, y = y, cluster_max = cluster_max,
      cluster_min_weight = cluster_min_weight,
      cluster_value_tolerance = cluster_value_tolerance,
      check_inputs = check_inputs
    )
  } else {
    b_t2 <- numeric()
  }
  if (P[3] > 0) {
    if (sum(P) > P[3]) {
      ind <- reg[1] + reg[2] + seq_len(reg[3])
      ortho <- orthogonalize_regressors(
        y = y, X = X, ind = ind, check_inputs = check_inputs
      )
      RX <- ortho$RX
      Ry <- ortho$Ry
    } else {
      RX <- X
      Ry <- y
    }
    b_t3 <- initialize_b_type_3(
      X = RX, y = Ry,
      check_inputs = isTRUE(check_inputs) && n_types == 1L
    )
  } else {
    b_t3 <- numeric()
  }
  est <- c(b_t1, b_t2, b_t3)
  ### find common utility scale
  if (n_types > 1L) {
    scale_t1 <- rep(b_t1, each = J)
    scale_t2 <- numeric()
    if (P[2] > 0L) {
      scale_t2 <- as.vector(rbind(0, matrix(b_t2, nrow = J - 1)))
    }
    scale_t3 <- b_t3
    scale <- c(scale_t1, scale_t2, scale_t3)
    active <- colSums(matrix(abs(scale), nrow = J)) > 0
    X_scale <- X
    if (P[2] > 0) {
      X_scale <- oeli::insert_matrix_column(
        X_scale, 0, reg[1] + (seq_len(P[2]) - 1L) * (J - 1L)
      )
    }
    Xb <- sweep(X_scale, 2, scale, "*")
    for (p in seq_len(sum(P))) {
      block_cols <- (p - 1L) * J + seq_len(J)
      block <- Xb[, block_cols, drop = FALSE]
      Xb[, block_cols] <- block - rowMeans(block)
    }
    norm <- numeric(sum(P))
    if (any(active)) {
      cols <- rep(active, each = J)
      norm[active] <- initialize_b_type_1(
        X = Xb[, cols, drop = FALSE], y = y - 1 / J, check_inputs = FALSE
      )
    }
    norm_t2 <- numeric()
    if (P[2] > 0L) {
      norm_t2 <- norm[P[1] + seq_len(P[2])]
      if (P[2] > 1L) {
        idx_t2 <- oeli::chunk_vector(
          reg[1] + seq_len(reg[2]), n = J - 1, type = 2
        )
        constant_t2 <- vapply(idx_t2, function(ind) {
          all(X[, ind, drop = FALSE] == matrix(
            X[1L, ind], nrow = NT, ncol = length(ind), byrow = TRUE
          ))
        }, logical(1))
        active_t2 <- active[P[1] + seq_len(P[2])]
        constant <- constant_t2 & active_t2
        varying <- !constant_t2 & active_t2
        if (any(constant) && any(varying)) {
          norm_t2[constant] <- stats::median(norm_t2[varying])
        }
      }
    }
    norm_factor <- c(
      norm[seq_len(P[1])],
      rep(norm_t2, each = J - 1),
      rep(norm[P[1] + P[2] + seq_len(P[3])], each = J)
    )
    est <- est * norm_factor
  }
  out <- numeric(length(est))
  out[order(regressor_position, seq_along(regressor_position))] <- est
  norm_out <- sqrt(sum(out^2))
  if (!is.finite(norm_out)) {
    warning(
      "Coefficient initialization failed; using unit fallback.",
      call. = FALSE
    )
    out <- numeric(length(out))
  } else if (norm_out > 0) {
    out <- as.numeric(out / norm_out)
  }
  structure(out, "sec" = difftime(Sys.time(), t0, units = "secs"))
}

# Initialize covariance via Gibbs sampling
# data: list, each element is a list of
# - X, a list, where the t-th element is the covariate matrix at occasion t,
# - y, integer vector, t-th element is the chosen alternative at occasion t.
# J: number of choice alternatives (J >= 2)
# b_fixed: fixed numeric coefficient vector
# re_cols: integer positions of coefficients with Gaussian random effects
# R: number of Gibbs iterations
# B: length of burn-in period
# Q: thinning factor
# N_max: maximum number of deciders used
# check_inputs: whether to check inputs
# Returns: list of (with computation time (sec) as attribute)
# - init_l, numeric vector of length J * (J - 1) / 2 (Cholesky root elements
#   for the differenced (J - 1) x (J - 1) error-term covariance matrix relative
#   to alternative 1 for level identification)
# - init_o, numeric vector of length length(re_cols) * (length(re_cols) + 1) / 2
#   (Cholesky root elements for the Gaussian random effect covariance matrix)
initialize_covariance <- function(
    data, J, b_fixed, re_cols = seq_along(b_fixed), R = 200,
    B = R %/% 2, Q = 1, N_max = 1000, check_inputs = TRUE
) {
  if (isTRUE(check_inputs)) {
    oeli::assert_list_of_lists(data)
    checkmate::assert_int(J, lower = 2)
    oeli::assert_numeric_vector(
      b_fixed, finite = TRUE, any.missing = FALSE
    )
    checkmate::assert_integerish(
      re_cols, any.missing = FALSE, unique = TRUE
    )
    if (length(re_cols) > 0L) {
      checkmate::assert_integerish(
        re_cols, lower = 1, upper = length(b_fixed), any.missing = FALSE,
        unique = TRUE
      )
    }
    checkmate::assert_count(R, positive = TRUE)
    checkmate::assert_int(B, lower = 0, upper = R - 1)
    checkmate::assert_int(Q, lower = 1, upper = R - B)
    checkmate::assert_int(N_max, lower = 1)
  }
  J <- as.integer(J)
  R <- as.integer(R)
  B <- as.integer(B)
  Q <- as.integer(Q)
  N_max <- as.integer(N_max)
  re_ind <- as.integer(re_cols)
  fe_ind <- setdiff(seq_along(b_fixed), re_ind)
  any_fe <- length(fe_ind) > 0L
  any_re <- length(re_ind) > 0L
  t0 <- Sys.time()
  # discard incomplete data
  data <- lapply(data, function(decider) {
    keep <- Map(function(X_i, y_i) {
      valid_y <- length(y_i) == 1L && !is.na(y_i) && is.finite(y_i) &&
        y_i %in% seq_len(J)
      valid_y && all(is.finite(X_i))
    }, decider$X, decider$y)
    keep <- unlist(keep, use.names = FALSE)
    decider$X <- decider$X[keep]
    decider$y <- decider$y[keep]
    decider
  })
  data <- Filter(function(decider) length(decider$y) > 0L, data)
  if (length(data) == 0L) {
    warning(
      paste(
        "Covariance initialization failed because no complete observations",
        "were available; using identity fallback."
      ),
      call. = FALSE
    )
    est_l <- oeli::cov_to_chol(diag(J - 1L))
    est_o <- if (any_re) {
      oeli::cov_to_chol(diag(length(re_ind)))
    } else {
      numeric()
    }
    return(structure(
      list("init_l" = est_l, "init_o" = est_o),
      "sec" = difftime(Sys.time(), t0, units = "secs")
    ))
  }
  ### difference covariates
  N <- min(length(data), N_max)
  data <- data[round(seq(1L, length(data), length.out = N))]
  Tp <- sapply(data, function(n) length(n$y))
  cs_Tp <- c(0, cumsum(Tp))
  NT <- sum(Tp)
  X <- unlist(lapply(data, `[[`, "X"), recursive = FALSE)
  Y <- unlist(lapply(data, `[[`, "y"))
  U <- lapply(Y, function(x) {out <- numeric(J); out[x] <- 1; return(out)})
  D <- oeli::delta(ref = 1, dim = J)
  X_diff <- lapply(X, function(x) D %*% x)
  U_diff <- lapply(U, function(u) D %*% u)
  if (length(b_fixed) == 0L) {
    sys <- lapply(X_diff, function(x) numeric(nrow(x)))
  } else {
    sys <- lapply(X_diff, function(x) x %*% b_fixed)
  }
  ### start Gibbs sampling
  Sigma <- diag(J - 1)
  draws_ind <- seq(B + 1, R, Q)
  draws_id <- 0
  Sigma_draws <- vector(mode = "list", length = length(draws_ind))
  Omega_draws <- vector(mode = "list", length = length(draws_ind))
  beta <- matrix(b_fixed[re_ind], ncol = N, nrow = length(re_ind))
  for (r in seq_len(R)) {
    Sigma_inv <- solve(Sigma)
    U_diff <- mapply(
      function(u, s, y) {
        U_diff_n <- RprobitB::update_U(
          U = u, y = y - 1, sys = s, Sigma_inv = Sigma_inv
        )
        if (any(is.na(U_diff_n) | is.infinite(U_diff_n))) {
          U_diff_n <- u
        }
        return(U_diff_n)
      },
      U_diff, sys, Y, SIMPLIFY = FALSE
    )
    Sigma <- RprobitB::update_Sigma(
      n_Sigma_0 = J + 1,
      V_Sigma_0 = diag(J - 1),
      N = NT,
      S = Reduce(`+`, mapply(
        function(u, s) (u - s) %*% t(u - s),
        U_diff, sys, SIMPLIFY = FALSE
      ))
    )
    if (any_re) {
      Sigma_inv <- solve(Sigma)
      Omega <- RprobitB::update_Omega(
        beta = beta,
        b = matrix(b_fixed[re_ind]),
        z = rep(1, N),
        m = N,
        n_Omega_0 = length(re_ind) + 2,
        V_Omega_0 = diag(length(re_ind))
      )
      Omega <- matrix(Omega, length(re_ind), length(re_ind))
      for (n in seq_len(N)) {
        ind_n <- (cs_Tp[n] + 1):(cs_Tp[n + 1])
        X_diff_n <- lapply(X_diff[ind_n], function(X) X[, re_ind, drop = FALSE])
        U_diff_n <- U_diff[ind_n]
        if (any_fe) {
          W_diff_n <- lapply(
            X_diff[ind_n], function(X) X[, fe_ind, drop = FALSE]
          )
          U_diff_n <- mapply(
            function(u, w) u - w %*% b_fixed[fe_ind], U_diff_n, W_diff_n,
            SIMPLIFY = FALSE
          )
        }
        XSigX_n <- Reduce(`+`, lapply(X_diff_n, function(X) {
          crossprod(X, Sigma_inv %*% X)
        }))
        XSigU_n <- Reduce(`+`, mapply(function(X, U) {
          crossprod(X, Sigma_inv %*% U)
        }, X_diff_n, U_diff_n, SIMPLIFY = FALSE))
        beta_n <- RprobitB::update_coefficient(
          mu_beta_0 = matrix(b_fixed[re_ind]),
          Sigma_beta_0_inv = solve(Omega),
          XSigX = XSigX_n,
          XSigU = XSigU_n
        )
        for (i in ind_n) {
          b_n <- b_fixed
          b_n[re_ind] <- beta_n
          sys[[i]] <- X_diff[[i]] %*% b_n
        }
        beta[, n] <- beta_n
      }
    }
    if (r %in% draws_ind) {
      draws_id <- draws_id + 1
      Sigma_draws[[draws_id]] <- Sigma
      if (isTRUE(any_re)) {
        Omega_draws[[draws_id]] <- Omega
      }
    }
  }
  Sigma <- Reduce(`+`, Sigma_draws) / length(draws_ind)
  est_l <- tryCatch(
    oeli::cov_to_chol(Sigma), error = function(e) NULL
  )
  len_l <- J * (J - 1L) / 2L
  if (!is.numeric(est_l) || length(est_l) != len_l || any(!is.finite(est_l))) {
    warning(
      "Sigma initialization failed; using identity fallback.",
      call. = FALSE
    )
    est_l <- oeli::cov_to_chol(diag(J - 1L))
  }
  if (isTRUE(any_re)) {
    Omega <- Reduce(`+`, Omega_draws) / length(draws_ind)
    est_o <- tryCatch(
      oeli::cov_to_chol(Omega), error = function(e) NULL
    )
    len_o <- length(re_ind) * (length(re_ind) + 1L) / 2L
    if (!is.numeric(est_o) ||length(est_o) != len_o || any(!is.finite(est_o))) {
      warning(
        "Omega initialization failed; using identity fallback.",
        call. = FALSE
      )
      est_o <- oeli::cov_to_chol(diag(length(re_ind)))
    }
  } else {
    est_o <- numeric()
  }
  structure(
    list("init_l" = est_l, "init_o" = est_o),
    "sec" = difftime(Sys.time(), t0, units = "secs")
  )
}
