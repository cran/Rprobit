# Compare DGP parameters with initial parameter vector
# theta_dgp: numeric vector, DGP parameters
# theta_init: numeric vector, initial parameters
# len: expected vector length
# rel_tol: maximum allowed relative Euclidean distance
# Returns: no return value, called for side effects
check_init <- function(theta_dgp, theta_init, len, rel_tol = 1) {
  checkmate::expect_count(len, positive = TRUE)
  checkmate::expect_number(rel_tol, lower = 0, finite = TRUE)
  checkmate::expect_numeric(
    theta_dgp, finite = TRUE, any.missing = FALSE, len = len
  )
  checkmate::expect_numeric(
    theta_init, finite = TRUE, any.missing = FALSE, len = len
  )
  sec <- attr(theta_init, "sec")
  checkmate::expect_class(sec, classes = "difftime", null.ok = FALSE)
  expect_identical(attr(sec, "unit"), "secs")
  rel_dist <- sqrt(sum((theta_dgp - theta_init)^2)) /
    max(sqrt(sum(theta_dgp^2)), .Machine$double.eps)
  expect_lte(rel_dist, rel_tol)
  rbind("theta_dgp" = theta_dgp, "theta_init" = theta_init) |> round(2)
}

test_that("regressors can be orthogonalized", {
  set.seed(1)
  N <- 100
  X <- matrix(stats::rnorm(N * 3), nrow = N)
  X[, 2] <- 0.8 * X[, 1] + stats::rnorm(N, sd = 0.2)
  beta <- stats::rnorm(3)
  y <- matrix(X %*% beta + stats::rnorm(N))
  c_orig <- stats::lm(y ~ X - 1) |> stats::coef()
  ortho <- orthogonalize_regressors(
    y = y, X = X, ind = 2:3, check_inputs = TRUE
  )
  c_m1 <- stats::lm(ortho$Ry ~ ortho$RX - 1) |> stats::coef()
  expect_equal(as.numeric(c_orig)[-1], as.numeric(c_m1))
  ortho <- orthogonalize_regressors(
    y = y, X = X, ind = 1, check_inputs = TRUE
  )
  c_m23 <- stats::lm(ortho$Ry ~ ortho$RX - 1) |> stats::coef()
  expect_equal(as.numeric(c_orig)[1], as.numeric(c_m23))
})

test_that("regressors can be localized", {
  set.seed(1)
  n <- 120
  cluster_max <- 10
  cluster_min_weight <- 0.1
  values <- c(
    stats::rnorm(80, mean = 0, sd = 0.1),
    stats::rnorm(40, mean = 3, sd = 0.1)
  )
  out <- localize_regressors(values, check_inputs = TRUE)
  checkmate::expect_list(out, len = 3)
  checkmate::expect_names(names(out), identical.to = c("loc", "ids", "w"))
  expect_length(out$ids, n)
  expect_true(all.equal(sum(out$w), 1))
  expect_lte(length(out$loc), cluster_max)
  expect_gte(min(out$w), cluster_min_weight)
  expect_true(all(out$ids %in% seq_along(out$loc)))
  expect_true(all(diff(out$w) <= 0))
  id_factor <- factor(out$ids, levels = seq_along(out$loc))
  expect_equal(out$loc, as.numeric(tapply(values, id_factor, mean)))
  expect_equal(out$w, as.numeric(table(id_factor) / n))
})

test_that("choice probabilities can be calculated", {
  J <- 3
  sys <- rep(0, J)
  cov <- diag(J)
  out <- sys_to_prob(sys = sys, cov = cov, check_inputs = TRUE)
  expect_equal(as.numeric(out), rep(1, J) / J)
})

test_that("SSR between probabilities and frequencies can be computed", {
  J <- 3
  sys_m1 <- rep(0, J - 1)
  fix_sys1 <- 0
  cov <- diag(J)
  freq <- rep(1, J) / 3
  out <- ssr_prob_freq(
    sys_m1 = sys_m1, cov = cov, freq = freq, fix_sys1 = fix_sys1,
    check_inputs = TRUE
  )
  checkmate::expect_number(out, lower = 0, finite = TRUE)
  gradient <- attr(out, "gradient", exact = TRUE)
  expect_true(oeli::test_numeric_vector(gradient, len = J - 1))
  expect_equal(as.numeric(out), 0)
  expect_equal(gradient, rep(0, J - 1), tolerance = 1e-12)
})

test_that("type-1 regressors can be identified", {
  J <- 3L
  zero <- matrix(0, nrow = J, ncol = 1L)
  varying <- matrix(c(0, 1, 2), ncol = 1L)
  expect_true(is_type_1(list(zero, varying), 1L))
  expect_false(is_type_1(list(zero), 1L))
})

test_that("type-2 regressors can be identified", {
  J <- 3L
  zero <- matrix(0, nrow = J, ncol = J - 1L)
  type_2 <- diag(J)[, -1L, drop = FALSE] * 2
  expect_true(is_type_2(list(zero, type_2), J, 1L, 1L))
  expect_false(is_type_2(list(zero), J, 1L, 1L))
  invalid <- type_2
  invalid[2L, 2L] <- 1
  expect_false(is_type_2(list(type_2, invalid), J, 1L, 1L))
  type_2 <- diag(J)[, -2L, drop = FALSE] * 2
  expect_true(is_type_2(list(type_2), J, 1L, 2L))
})

test_that("type-3 regressors can be identified", {
  J <- 3L
  first <- diag(c(1, 0, 0))
  second <- diag(c(0, 2, 3))
  expect_true(is_type_3(list(first, second), J, 1L))
  expect_false(is_type_3(list(first), J, 1L))
  invalid <- second
  invalid[1L, 2L] <- 1
  expect_false(is_type_3(list(first, invalid), J, 1L))
})

test_that("regressor types can be detected", {
  J <- 4L
  P <- c(2L, 1L, 1L)
  type_1a <- c(0, 3, 4, 5)
  type_1b <- c(0, -1, 2, -3)
  type_2 <- diag(J)[, -1, drop = FALSE] * 2
  type_3 <- diag(c(1, 1, 2, 3))
  X <- cbind(type_2, type_1a, type_3, type_1b)
  out <- regressor_types(X, J, P, ref_alt_index = 1L, check_inputs = TRUE)
  expect_identical(out, c(2L, 1L, 3L, 1L))

  # uses the reference alternative when identifying type-2 blocks
  J <- 3L
  P <- c(1L, 1L, 1L)
  type_1 <- c(1, 2, 4)
  type_3 <- diag(c(2, 5, 7))
  for (ref_alt_index in seq_len(J)) {
    type_2 <- diag(J)[, -ref_alt_index, drop = FALSE] * 3
    X <- cbind(type_2, type_1, type_3)
    out <- regressor_types(X, J, P, ref_alt_index = ref_alt_index)
    expect_identical(out, c(2L, 1L, 3L))
  }
})

test_that("initialize b type-1 works", {
  set.seed(1)
  Tp <- rep(3L, 300)
  J <- 3
  P <- c(8L, 0L, 0L)
  type_1_terms <- paste0("x", seq_len(P[1]), collapse = " + ")
  form <- stats::as.formula(paste("choice ~", type_1_terms, "| 0 | 0"))
  m <- mod$new(formula = form, alt = J)
  sim <- choice_data$new()
  sim$simulate(mod = m, N = length(Tp), Tp = Tp)
  data <- sim$build_design(m)
  data_flat <- data |>
    Map(f = function(d) Map(function(x, y) list(X = x, y = y), d$X, d$y)) |>
    unlist(recursive = FALSE)
  X <- do.call(rbind, lapply(data_flat, function(d) as.numeric(d$X)))
  y <- t(sapply(data_flat, function(d) {
    vec <- numeric(J)
    vec[d$y] <- 1
    vec
  }))
  theta_init <- initialize_b_type_1(X, y, check_inputs = TRUE)
  theta_dgp <- as.numeric(m$build_parameter(sim$theta_dgp)$b)
  theta_dgp <- theta_dgp / sqrt(sum(theta_dgp^2))
  check_init(theta_dgp = theta_dgp, theta_init = theta_init, len = P[1])
})

test_that("initialize b type-2 single works", {
  set.seed(1)
  Tp <- sample(5, size = 200, replace = TRUE)
  ASC <- FALSE
  P <- c(0L, 1L, 0L)
  J <- 3L
  type_2_rhs <- paste(
    c(if (ASC) "1" else "0", paste0("x", seq_len(P[2] - ASC))),
    collapse = " + "
  )
  form <- stats::as.formula(paste("choice ~ 0 |", type_2_rhs, "| 0"))
  m <- mod$new(formula = form, alt = J)
  sim <- choice_data$new()
  sim$simulate(mod = m, N = length(Tp), Tp = Tp)
  data <- sim$build_design(m)
  data_flat <- data |>
    Map(f = function(d) Map(function(x, y) list(X = x, y = y), d$X, d$y)) |>
    unlist(recursive = FALSE)
  X <- do.call(rbind, lapply(data_flat, function(d) {
    as.numeric(rep(1, J) %*% d$X)
  }))
  y <- t(sapply(data_flat, function(d) {
    vec <- numeric(J)
    vec[d$y] <- 1
    vec
  }))
  theta_init <- initialize_b_type_2_single(X, y, check_inputs = TRUE)
  theta_dgp <- as.numeric(m$build_parameter(sim$theta_dgp)$b)
  theta_dgp <- theta_dgp / sqrt(sum(theta_dgp^2))
  check_init(
    theta_dgp = theta_dgp, theta_init = theta_init, len = P[2] * (J - 1L)
  )
})

test_that("initialize b type-2 works", {
  set.seed(1)
  Tp <- sample(5, size = 200, replace = TRUE)
  ASC <- TRUE
  P <- c(0L, 4L, 0L)
  J <- 3L
  type_2_rhs <- paste(
    c(if (ASC) "1" else "0", paste0("x", seq_len(P[2] - ASC))),
    collapse = " + "
  )
  form <- stats::as.formula(paste("choice ~ 0 |", type_2_rhs, "| 0"))
  m <- mod$new(formula = form, alt = J)
  sim <- choice_data$new()
  sim$simulate(mod = m, N = length(Tp), Tp = Tp)
  data <- sim$build_design(m)
  data_flat <- data |>
    Map(f = function(d) Map(function(x, y) list(X = x, y = y), d$X, d$y)) |>
    unlist(recursive = FALSE)
  X <- do.call(rbind, lapply(data_flat, function(d) {
    as.numeric(rep(1, J) %*% d$X)
  }))
  y <- t(sapply(data_flat, function(d) {
    vec <- numeric(J)
    vec[d$y] <- 1
    vec
  }))
  theta_init <- initialize_b_type_2(X, y, check_inputs = TRUE)
  theta_dgp <- as.numeric(m$build_parameter(sim$theta_dgp)$b)
  theta_dgp <- theta_dgp / sqrt(sum(theta_dgp^2))
  check_init(
    theta_dgp = theta_dgp, theta_init = theta_init, len = P[2] * (J - 1)
  )
})

test_that("initialize b type-3 single works", {
  set.seed(1)
  Tp <- sample(2, size = 200, replace = TRUE)
  P <- c(0L, 0L, 1L)
  J <- 3L
  type_3_terms <- paste0("x", seq_len(P[3]), collapse = " + ")
  form <- stats::as.formula(paste("choice ~ 0 | 0 |", type_3_terms))
  m <- mod$new(formula = form, alt = J)
  sim <- choice_data$new()
  sim$simulate(mod = m, N = length(Tp), Tp = Tp)
  data <- sim$build_design(m)
  data_flat <- data |>
    Map(f = function(d) Map(function(x, y) list(X = x, y = y), d$X, d$y)) |>
    unlist(recursive = FALSE)
  X <- do.call(rbind, lapply(data_flat, function(d) {
    as.numeric(rep(1, J) %*% d$X)
  }))
  y <- t(sapply(data_flat, function(d) {
    vec <- numeric(J)
    vec[d$y] <- 1
    vec
  }))
  theta_init <- initialize_b_type_3_single(X, y, check_inputs = TRUE)
  theta_dgp <- as.numeric(m$build_parameter(sim$theta_dgp)$b)
  theta_dgp <- theta_dgp / sqrt(sum(theta_dgp^2))
  check_init(
    theta_dgp = theta_dgp, theta_init = theta_init, len = P[3] * J
  )
})

test_that("initialize b type-3 works", {
  set.seed(1)
  Tp <- sample(2, size = 300, replace = TRUE)
  P <- c(0L, 0L, 4L)
  J <- 3L
  type_3_terms <- paste0("x", seq_len(P[3]), collapse = " + ")
  form <- stats::as.formula(paste("choice ~ 0 | 0 |", type_3_terms))
  m <- mod$new(formula = form, alt = J)
  sim <- choice_data$new()
  sim$simulate(mod = m, N = length(Tp), Tp = Tp)
  data <- sim$build_design(m)
  data_flat <- data |>
    Map(f = function(d) Map(function(x, y) list(X = x, y = y), d$X, d$y)) |>
    unlist(recursive = FALSE)
  X <- do.call(rbind, lapply(data_flat, function(d) {
    as.numeric(rep(1, J) %*% d$X)
  }))
  y <- t(sapply(data_flat, function(d) {
    vec <- numeric(J)
    vec[d$y] <- 1
    vec
  }))
  theta_init <- initialize_b_type_3(X, y, check_inputs = TRUE)
  theta_dgp <- as.numeric(m$build_parameter(sim$theta_dgp)$b)
  theta_dgp <- theta_dgp / sqrt(sum(theta_dgp^2))
  check_init(theta_dgp = theta_dgp, theta_init = theta_init, len = P[3] * J)
})

test_that("initialize b works", {
  set.seed(1)
  Tp <- sample(2, size = 600, replace = TRUE)
  ASC <- TRUE
  P <- c(2L, 2L, 1L)
  J <- 3L
  type_1_terms <- paste0("x", seq_len(P[1]), collapse = " + ")
  type_2_rhs <- paste(
    c(if (ASC) "1" else "0", paste0("z", seq_len(P[2] - ASC))),
    collapse = " + "
  )
  type_3_terms <- paste0("w", seq_len(P[3]), collapse = " + ")
  formula_rhs <- paste(type_1_terms, type_2_rhs, type_3_terms, sep = " | ")
  form <- stats::as.formula(paste("choice ~", formula_rhs))
  m <- mod$new(formula = form, alt = J)
  sim <- choice_data$new()
  sim$simulate(mod = m, N = length(Tp), Tp = Tp)
  data <- sim$build_design(m)
  perm <- c(3L, 4L, 1L, 7L, 8L, 9L, 2L, 5L, 6L)
  theta_dgp <- as.numeric(m$build_parameter(sim$theta_dgp)$b)[perm]
  data <- lapply(data, function(d) {
    d$X <- lapply(d$X, function(X) X[, perm, drop = FALSE])
    d
  })
  theta_init <- initialize_b(
    data = data, J = J, P = P, ref_alt_index = 1L,
    check_inputs = TRUE
  )
  theta_dgp <- theta_dgp / sqrt(sum(theta_dgp^2))
  check_init(
    theta_dgp = theta_dgp, theta_init = theta_init, len = length(theta_dgp)
  )
})

test_that("initialize covariances works", {
  set.seed(1)
  Tp <- sample(2, size = 200, replace = TRUE)
  J <- 3L
  P <- c(5L, 0L, 0L)
  type_1_terms <- paste0("x", seq_len(P[1]), collapse = " + ")
  form <- stats::as.formula(paste("choice ~", type_1_terms, "| 0 | 0"))
  m <- mod$new(
    formula = form,
    re = c("x1", "x3"),
    corr_re = TRUE,
    alt = J
  )
  sim <- choice_data$new()
  sim$simulate(mod = m, N = length(Tp), Tp = Tp)
  data <- sim$build_design(m)
  param <- m$build_parameter(sim$theta_dgp)
  b <- as.numeric(param$b)
  cols <- colnames(data[[1]]$X[[1]])
  re_cols <- match(m$var_names$random, cols)
  theta_init <- initialize_covariance(
    data = data, J = J, b_fixed = b, re_cols = re_cols,
    R = 30L, B = 5L, Q = 5L, N_max = 80L, check_inputs = TRUE
  )
  theta_init_vec <- structure(
    c(theta_init$init_o, theta_init$init_l), sec = attr(theta_init, "sec")
  )
  Omega_root <- t(chol(param$Omega))
  Sigma_root <- t(chol(param$Sigma[-1, -1, drop = FALSE]))
  theta_dgp <- c(
    Omega_root[lower.tri(Omega_root, diag = TRUE)],
    Sigma_root[lower.tri(Sigma_root, diag = TRUE)]
  )
  check_init(
    theta_dgp = theta_dgp, theta_init = theta_init_vec,
    len = length(theta_dgp)
  )
})

test_that("initialize works", {
  set.seed(1)
  Tp <- sample(5, size = 100, replace = TRUE)
  J <- 3L
  P <- c(2L, 2L, 1L)
  type_1_terms <- paste0("x", seq_len(P[1]), collapse = " + ")
  type_2_terms <- paste0("z", seq_len(P[2] - 1L), collapse = " + ")
  type_3_terms <- paste0("w", seq_len(P[3]), collapse = " + ")
  formula_rhs <- paste(type_1_terms, type_2_terms, type_3_terms, sep = " | ")
  form <- stats::as.formula(paste("choice ~", formula_rhs))
  m <- mod$new(
    formula = form,
    re = c("x1", "ASC"),
    corr_re = TRUE,
    alt = J
  )
  sim <- choice_data$new()
  sim$simulate(mod = m, N = length(Tp), Tp = Tp)
  data <- sim$build_design(m)
  param <- m$build_parameter(sim$theta_dgp)
  cols <- colnames(data[[1]]$X[[1]])
  re_cols <- match(m$var_names$random, cols)
  theta_init <- initialize(
    data = data, J = J, P = P, ref_alt_index = 1L,
    re_cols = re_cols,
    R = 20L, B = 5L, Q = 5L, N_max = 80L, check_inputs = TRUE
  )
  b <- as.numeric(param$b)
  norm_b <- sqrt(sum(b^2))
  Omega_root <- t(chol(param$Omega))
  Sigma_root <- t(chol(param$Sigma[-1, -1, drop = FALSE]))
  theta_dgp <- c(
    b / norm_b,
    Omega_root[lower.tri(Omega_root, diag = TRUE)] / norm_b,
    Sigma_root[lower.tri(Sigma_root, diag = TRUE)] / norm_b
  )
  check_init(
    theta_dgp = theta_dgp, theta_init = theta_init,
    len = length(theta_dgp)
  )
})
